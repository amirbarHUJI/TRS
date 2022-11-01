import sys
import argparse
import yaml
import logging
import datetime
import numpy as np
import pandas as pd

from pathlib import Path
from collections import defaultdict
from scipy.stats import nbinom

from statsmodels.stats.multitest import multipletests

from ARES.MergedPeakCaller import  MergedPeakCaller
from ARES.Globals import Strand
from ARES.Lib import Lib

SIGNIFICANT_PEAKS_FILE = "significant_peaks.csv"
ALL_MEAN_PEAKS_FILE = "all_mean_peaks.csv"


def process_command_line(argv=None):
    """
    Return a 2-tuple: (settings object, args list).
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object, replace the description
    parser = argparse.ArgumentParser(
        description="Uses read starts to determine 3' termini of transcripts.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'scheme',
        help="Path for file containing info of the libs to process.")

    parser.add_argument(
        '-w',
        '--workdir',
        help='Working directory.',
        default="./runs")

    parser.add_argument(
        '-t',
        '--threshold',
        help='The signficance level used by the model.',
        type=float,
        default=0.01)

    def process_validator(value):
        ivalue = int(value)

        if ivalue <= 0:
            raise argparse.ArgumentTypeError(
                "Process count must be a positive integer."
            )

        return ivalue

    parser.add_argument(
        '-b',
        '--force_bam',
        action='store_true',
        help='Forces the program to reprocess the bam files.')

    parser.add_argument(
        '-c',
        '--min_count',
        default=10,
        type=int,
        help='The minimal coverage for a region to be considered.')

    parser.add_argument(
        '--min_height',
        default=None,
        type=float,
        help='The minimal ratio to consider as a peak.')

    parser.add_argument(
        '--window_margin',
        default=3,
        type=int,
        help='Defines the region which will be used to count the '
        'local number of reads starts. The margin is the number '
        'of nucleotides (upstream/downstream) that will be added '
        'to the region around the considered site.')

    parser.add_argument(
        '--merge_distance',
        default=0,
        type=int,
        help='The distance which below it, peaks will be merged together.')

    parser.add_argument(
        '--rel_height',
        default=0.75,
        type=float,
        help='The relative height of peaks to use in the '
        'scipy find_peaks function.')

    parser.add_argument(
        '--chr_list',
        default='',
        help='Conversion of chromosome files from the bam to new name '
        'in the following format: bam1:new1,bam2:new2,...')

    parser.add_argument(
        '-d',
        '--ds_distance',
        type=int,
        default=70,
        help='The distance (downstream) used to compute the read starts ratio.')

    parser.add_argument(
        '--signif_min_lib_count',
        type=int,
        default=2,
        help="The minimal number of libraries in which the 3' terminus should be found as significant to report it depending on the threshold.")

    parser.add_argument(
        '--insignif_min_ratio',
        type=float,
        default=0.5,
        help="The minimal mean ratio to accept if the 3' terminus wasn't significant in all repeats.")

    parser.add_argument(
        '-l',
        '--log_level',
        default="debug",
        help='The logging level to report to the log file.')

    settings = parser.parse_args(argv)

    return settings


def process_config(config_file):
    with open(config_file, "rt") as fl:
        config_dct = yaml.safe_load(fl)

    lib_list = [
        Lib(lib_dct["name"], lib_dct["group"], lib_dct["count_file"])
        for lib_dct in config_dct["libs"]
    ]

    group_dct = defaultdict(list)

    for lib in lib_list:
        group_dct[lib.Group].append(lib)

    return lib_list, group_dct


def configure_log(log_dir, log_level):
    time = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    log_file = log_dir.joinpath("log_%s.log" % time)
    log_format = "%(asctime)s [%(levelname)s] %(message)s"

    log_level = getattr(logging, log_level.upper(), None)

    if not isinstance(log_level, int):
        raise ValueError('Invalid log level: %s' % log_level)

    logging.basicConfig(
        filename=str(log_file),
        format=log_format,
        level=log_level)


def estimate_dispersion(lib_list, df, window_margin, ds_distance):
    mean_list = []
    var_list = []

    for i, lib in enumerate(lib_list):
        for idx, peak in df.iterrows():
            chromosome = peak.Chromosome
            strand = Strand(peak.Strand)
            coordinate = peak.Dominant

            if strand == Strand.FWD:
                ds_start = coordinate + window_margin + 1
                ds_end = coordinate + ds_distance

            else:
                ds_start = coordinate - ds_distance
                ds_end = coordinate - window_margin - 1

            ds_ends = lib.Counts.Position[chromosome, strand, ds_start, ds_end]

            mu = ds_ends.mean()
            var = ds_ends.var() * len(ds_ends) / (len(ds_ends) - 1)

            if var != 0 and mu != 0:
                mean_list.append(mu)
                var_list.append(var)

    mean_list = np.log(mean_list)
    var_list = np.log(var_list)

    return np.polyfit(mean_list, var_list, 2)


def get_peak_coordinate(chromosome, strand, start, end, counts, window_margin):
    window_end_count = np.zeros(end - start + 1)

    for i in range(-window_margin, window_margin + 1):
        window_end_count += counts.Position[
            chromosome,
            strand,
            start + i,
            end + i
        ]

    return start + np.argmax(window_end_count)


def test_peak(
    peak,
    counts,
    ds_distance,
    window_margin,
    min_count,
    coef
):
    chromosome = peak.Chromosome
    strand = Strand(peak.Strand)

    # Assign the coordinate in the peak as the highest ends position
    # coordinate = peak.Start + np.argmax(counts.Position[chromosome, strand, peak.Start, peak.End])
    coordinate = get_peak_coordinate(
        chromosome=chromosome,
        strand=strand,
        start=peak.Start,
        end=peak.End,
        counts=counts,
        window_margin=window_margin
    )

    if strand == Strand.FWD:
        ds_start = coordinate + window_margin + 1
        ds_end = coordinate + ds_distance

    else:
        ds_start = coordinate - ds_distance
        ds_end = coordinate - window_margin - 1

    ds_ends = counts.Position[chromosome, strand, ds_start, ds_end]
    window_ends = counts.Position[chromosome, strand, coordinate - window_margin, coordinate + window_margin]
    total_ends = ds_ends.sum() + window_ends.sum()

    # If the position is below coverage ratio return p-value 1
    if total_ends < min_count:
        return 1, 0

    position_ratio = counts.Ratio[chromosome, strand, peak.Start, peak.End].max()

    mu = ds_ends.mean()

    # No read starts downstream
    if mu == 0:
        p = 1
        n = 1
        var = None

    else:
        var = np.exp(np.polyval(coef, np.log(mu)))

        p = mu / var
        n = mu ** 2 / (var - mu)

    pval = nbinom.sf(
        k=window_ends.sum(),
        n=n * (2 * window_margin + 1),
        p=p
    )

    return pval, position_ratio


def main():
    settings = process_command_line()
    log_dir = Path(settings.workdir).joinpath("logs")
    data_dir = Path(settings.workdir).joinpath("data")
    results_dir = Path(settings.workdir).joinpath("results")

    log_dir.mkdir(exist_ok=True)
    data_dir.mkdir(exist_ok=True)
    results_dir.mkdir(exist_ok=True)

    configure_log(
        log_dir=log_dir,
        log_level=settings.log_level
    )

    logging.info("Running from: '%s'" % str(Path().cwd()))
    logging.info("Exectuing: '%s'" % " ".join(sys.argv))

    lib_list, group_dct = process_config(settings.scheme)

    if settings.chr_list:
        chr_conv_dct = dict([
            tuple(entry.split(":"))
            for entry in settings.chr_list.split(",")
        ])

    else:
        chr_conv_dct = {}

    for lib in lib_list:
        logging.info("Loading library '%s'" % lib.Name)
        lib.load(
            workdir=data_dir,
            force_bam=settings.force_bam,
            chr_conv_dct=chr_conv_dct,
            min_count=settings.min_count,
            window_width=settings.window_margin,
            ds_distance=settings.ds_distance
        )

    # If no minimal ratio threshold was defined, take the expected mean of the distribution
    if settings.min_height:
        min_height = settings.min_height

    else:
        min_height = (1 + 2 * settings.window_margin) / (1 + settings.window_margin + settings.ds_distance)

    mpc = MergedPeakCaller()

    # Process each group of libraries separately
    for group, group_lib_list in group_dct.items():
        logging.info("Processing group '%s'" % group)
        group_dir = results_dir.joinpath(group)
        group_dir.mkdir(exist_ok=True)

        df_group_peaks = mpc.call_peaks(
            lib_list=group_lib_list,
            window_width=settings.window_margin,
            merge_distance=settings.merge_distance,
            rel_height=settings.rel_height,
            min_height=settings.min_height,
        )

        pvals_list = []
        observed_r_list = []

        logging.info("Estimating parameters...")

        coef = estimate_dispersion(
            lib_list=group_lib_list,
            df=df_group_peaks,
            window_margin=settings.window_margin,
            ds_distance=settings.ds_distance
        )
        logging.info("Done")

        for idx, peak in df_group_peaks.iterrows():
            pvals, observed_ratios = zip(*[
                test_peak(
                    peak=peak,
                    counts=lib.Counts,
                    ds_distance=settings.ds_distance,
                    window_margin=settings.window_margin,
                    min_count=settings.min_count,
                    coef=coef
                )
                for lib in group_lib_list
            ])

            pvals_list.append(pvals)
            observed_r_list.append(observed_ratios)

        pvals_list = np.array(pvals_list)
        observed_r_list = np.array(observed_r_list)

        # Define passed criteria
        padj_list = np.array([
            multipletests(pvals=lst, alpha=settings.threshold, method='bonferroni')[1]
            for lst in pvals_list.T
        ]).T

        for pval_lst, lib in zip(padj_list.T, group_lib_list):
            df_group_peaks[lib.Name] = pval_lst

        idx_in_all = (padj_list <= settings.threshold).sum(axis=1) == len(group_lib_list)
        idx_at_least_two = (padj_list <= settings.threshold).sum(axis=1) >= settings.signif_min_lib_count
        idx_above_ratio_threshold = observed_r_list.mean(axis=1) >= settings.insignif_min_ratio
        passed_idx = idx_in_all | (idx_at_least_two & idx_above_ratio_threshold)

        logging.info("Total ends: %d" % passed_idx.sum())
        logging.info("Significant in all repeats: %d" % idx_in_all.sum())
        logging.info("Peaks above minimal ratio: %d" % idx_above_ratio_threshold.sum())

        significant_peaks_file = str(group_dir.joinpath(SIGNIFICANT_PEAKS_FILE))
        all_mean_peaks_file = str(group_dir.joinpath(ALL_MEAN_PEAKS_FILE))

        logging.info("Generating significant peaks output '%s'" % significant_peaks_file)

        df_passed = df_group_peaks.iloc[np.array([i for i, is_passd in enumerate(passed_idx) if is_passd])]
        df_passed = df_passed.sort_values(["Chromosome", "Dominant"], ascending=[True, True])
        df_passed.to_csv(significant_peaks_file, sep="\t", index=False)

        logging.info("Generating all mean peaks output '%s'" % all_mean_peaks_file)
        df_group_peaks.to_csv(all_mean_peaks_file, sep="\t", index=False)

    return 0


if __name__ == "__main__":
    exit(main())
