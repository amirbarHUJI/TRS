import logging
import itertools
import numpy as np
import pandas as pd

from multiprocessing import Pool
from intervaltree import IntervalTree

from TRS.LibPeakCaller import LibPeakCaller
from TRS.Globals import Strand


def run_process(args):
    return args[0].get_paired_significant_peaks(**args[1])


class MergedPeakCaller(LibPeakCaller):

    def __init__(self):
        pass

    def call_peaks(
        self,
        lib_list,
        window_width=3,
        merge_distance=0,
        rel_height=0.75,
        min_height=0.1,
    ):

        lib_signal_dct = {}
        lib_signal_peak_dct = {}

        for lib in lib_list:
            lib_signal_dct[lib.Name] = self.get_signal_dct(
                signal_counts=lib.Counts.Ratio)

            lib_signal_peak_dct[lib.Name] = self.get_peak_dct(
                signal_dct=lib_signal_dct[lib.Name],
                prominence=(None, None),
                width=(window_width, None),
                rel_height=rel_height,
                distance=10,
                min_height=min_height,
            )

        mean_signal_dct = self.get_mean_signal(
            lib_list=lib_list,
            lib_signal_dct=lib_signal_dct,
        )

        mean_signal_peak_dct = self.get_peak_dct(
            signal_dct=mean_signal_dct,
            prominence=(None, None),
            width=(window_width, None),
            rel_height=rel_height,
            distance=10,
            min_height=min_height,
        )

        mean_peak_list = [
            peak
            for sublist in mean_signal_peak_dct.values()
            for peak in sublist
        ]
        mean_peak_list.sort(key=lambda peak: peak.Signal, reverse=True)

        mean_peak_row_list = []

        for i, peak in enumerate(mean_peak_list, 1):
            dominant, signal = self.find_dominant_coord(
                peak.Chromosome,
                Strand(peak.Strand),
                peak.Start,
                peak.End,
                mean_signal_dct)
            mean_peak_row_list.append([
                peak.Chromosome,
                peak.Strand,
                peak.Start + 1,
                peak.End + 1,
                dominant,
                signal,
                i
            ])

        df_all_mean_peaks = pd.DataFrame(
            mean_peak_row_list,
            columns=[
                "Chromosome",
                "Strand",
                "Start",
                "End",
                "Dominant",
                "Signal",
                "Rank"
            ]
        )
        # df_all_mean_peaks = df_all_mean_peaks.sort_values(["Chromosome", "Dominant"], ascending=[True, True])
        df_all_mean_peaks = df_all_mean_peaks.sort_values(["Signal"], ascending=[False])

        if merge_distance > 0:
            df_significant_signals = self.merge_signals(
                df=df_all_mean_peaks,
                padding=merge_distance
            )

        return df_all_mean_peaks

    def get_mean_signal(
        self,
        lib_list,
        lib_signal_dct,
    ):
        mean_signal_dct = {}

        for chr_name in self.get_libs_chromosomes(lib_list):
            mean_signal_dct[chr_name] = {}

            for strand in (Strand.FWD, Strand.REV):
                signal_array = np.array([
                    lib_signal_dct[lib.Name][chr_name][strand]
                    for lib in lib_list
                ])

                mean_signal_dct[chr_name][strand] = \
                    np.mean(signal_array, axis=0)

        return mean_signal_dct

    def get_libs_chromosomes(self, lib_list):
        chr_dct = {}

        for lib in lib_list:
            for chr_name, chr_len in lib.Chromosomes:
                if chr_name not in chr_dct:
                    chr_dct[chr_name] = chr_len

                else:
                    if chr_dct[chr_name] != chr_len:
                        raise ValueError("Chromosomes lengths should match between libraries")

        return chr_dct

    def find_dominant_coord(
        self,
        chr_name,
        strand,
        start,
        end,
        mean_signal_dct
    ):
        x = mean_signal_dct[chr_name][strand][start:end + 1]
        max_index_array = np.where(x == x.max())[0]

        # Take the most downstream if equal ratios
        if strand == Strand.FWD:
            i = max_index_array[-1]

        else:
            i = max_index_array[0]

        return start + i, x[i]

    def merge_signals(self, df, padding=50):
        tree_dct = {}

        for idx, row in df.iterrows():
            chr_name = row["Chromosome"]
            strand = Strand(row["Strand"])
            start = row["Start"]
            end = row["End"]
            dominant = row["Dominant"]
            signal = row["Signal"]

            if chr_name not in tree_dct:
                tree_dct[chr_name] = {
                    Strand.FWD: IntervalTree(),
                    Strand.REV: IntervalTree()
                }

            tree_dct[chr_name][strand].addi(
                start - padding,
                end + padding + 1,
                (dominant, signal))

        row_list = []

        for chr_name in tree_dct:
            for strand in tree_dct[chr_name]:
                tree = tree_dct[chr_name][strand]
                tree.merge_overlaps(
                    strict=False,
                    data_initializer=[],
                    data_reducer=lambda x, y: x + [y])

                for interval in tree:
                    dominant_list, signal_list = zip(*sorted(interval.data))
                    x = np.array(signal_list)

                    max_index_array = np.where(x == x.max())[0]

                    # Take the most downstream if equal ratios
                    if strand == Strand.FWD:
                        i = max_index_array[-1]

                    else:
                        i = max_index_array[0]

                    dominant_string = ",".join([
                        str(val) for val in dominant_list
                    ])
                    signal_string = ",".join([str(val) for val in signal_list])

                    row_list.append([
                        chr_name,
                        strand.value,
                        interval.begin,
                        interval.end - 1,
                        dominant_list[i],
                        signal_list[i],
                        dominant_string,
                        signal_string])

        df = pd.DataFrame(
            row_list,
            columns=[
                "Chromosome",
                "Strand",
                "Start",
                "End",
                "Dominant",
                "Signal",
                "Dominants Merged",
                "Signals Merged"])

        df = df.sort_values(["Chromosome", "Dominant"], ascending=[True, True])

        return df
