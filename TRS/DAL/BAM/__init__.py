import gzip
import pickle
import numpy as np

from TRS.Globals import Strand
from TRS.DAL.BAM.BamReader import BamReader
from TRS.DAL.Track.ChrCounts import ChrCounts
from TRS.DAL.Track.Counts import Counts


def bam_to_coverage(
    path,
    is_reverse=False,
    skip_clipped=False,
    paired_as_single=False,
    ignore_read1=False,
    ignore_read2=False,
    chr_conv_dct={}
):
    coverage_dct = {}

    with BamReader(path=path,
                   is_reverse=is_reverse,
                   skip_clipped=skip_clipped,
                   paired_as_single=paired_as_single,
                   ignore_read1=ignore_read1,
                   ignore_read2=ignore_read2) as br:

        # Generate chromosome structure
        for chromosome, length in br.Chromosomes.items():
            coverage_dct[chromosome] = {Strand.FWD: np.zeros(length),
                                        Strand.REV: np.zeros(length)}

        # Generate coverage
        for chromosome, strand, start, end in br:
            coverage_dct[chromosome][strand][start:end] += 1

        coverage = Counts()

        for chromosome, length in br.Chromosomes.items():
            chr_counts = ChrCounts(
                chr_len=length,
                fwd_counts=coverage_dct[chromosome][Strand.FWD],
                rev_counts=coverage_dct[chromosome][Strand.REV])

            chrom_name = chromosome

            if chromosome in chr_conv_dct:
                chrom_name = chr_conv_dct[chromosome]

            coverage.add_chromosome(
                chromosome=chrom_name,
                chr_counts=chr_counts)

    return coverage


def bam_to_read_ends(
    path,
    is_reverse=False,
    skip_clipped=False,
    paired_as_single=False,
    ignore_read1=False,
    ignore_read2=False,
    chr_conv_dct={}
):
    coverage_dct = {}

    with BamReader(path=path,
                   is_reverse=is_reverse,
                   skip_clipped=skip_clipped,
                   paired_as_single=paired_as_single,
                   ignore_read1=ignore_read1,
                   ignore_read2=ignore_read2) as br:

        # Generate chromosome structure
        for chromosome, length in br.Chromosomes.items():
            coverage_dct[chromosome] = {Strand.FWD: np.zeros(length),
                                        Strand.REV: np.zeros(length)}

        # Generate coverage
        for chromosome, strand, start, end in br:
            if strand == Strand.FWD:
                coverage_dct[chromosome][strand][end - 1] += 1

            else:
                coverage_dct[chromosome][strand][start] += 1

        coverage = Counts()

        for chromosome, length in br.Chromosomes.items():
            chr_counts = ChrCounts(
                chr_len=length,
                fwd_counts=coverage_dct[chromosome][Strand.FWD],
                rev_counts=coverage_dct[chromosome][Strand.REV])

            chrom_name = chromosome

            if chromosome in chr_conv_dct:
                chrom_name = chr_conv_dct[chromosome]

            coverage.add_chromosome(
                chromosome=chrom_name,
                chr_counts=chr_counts)

    return coverage


def read_bam_coverage_and_pos(
    path,
    is_reverse=False,
    skip_clipped=False,
    paired_as_single=False,
    ignore_read1=False,
    ignore_read2=False,
    chr_conv_dct={}
):
    coverage_dct = {}
    pos_dct = {}

    with BamReader(path=path,
                   is_reverse=is_reverse,
                   skip_clipped=skip_clipped,
                   paired_as_single=paired_as_single,
                   ignore_read1=ignore_read1,
                   ignore_read2=ignore_read2) as br:

        # Generate chromosome structure
        for chromosome, length in br.Chromosomes.items():
            pos_dct[chromosome] = {
                Strand.FWD: np.zeros(length),
                Strand.REV: np.zeros(length)
            }

            coverage_dct[chromosome] = {
                Strand.FWD: np.zeros(length),
                Strand.REV: np.zeros(length)
            }

        # Generate per chromosome coverage and pos counts
        for chromosome, strand, start, end in br:
            if strand == Strand.FWD:
                pos_dct[chromosome][strand][end - 1] += 1

            else:
                pos_dct[chromosome][strand][start] += 1

            coverage_dct[chromosome][strand][start:end] += 1

        coverage = Counts()
        pos = Counts()

        for chromosome, length in br.Chromosomes.items():
            chr_pos_counts = ChrCounts(
                chr_len=length,
                fwd_counts=pos_dct[chromosome][Strand.FWD],
                rev_counts=pos_dct[chromosome][Strand.REV])

            chr_coverage_counts = ChrCounts(
                chr_len=length,
                fwd_counts=coverage_dct[chromosome][Strand.FWD],
                rev_counts=coverage_dct[chromosome][Strand.REV])

            chrom_name = chromosome

            if chromosome in chr_conv_dct:
                chrom_name = chr_conv_dct[chromosome]

            pos.add_chromosome(
                chromosome=chrom_name,
                chr_counts=chr_pos_counts)

            coverage.add_chromosome(
                chromosome=chrom_name,
                chr_counts=chr_coverage_counts)

    return pos, coverage


def save_processed_bam(path, content):
    with gzip.open(path, "wb") as fl:
        pickle.dump(content, fl)


def restore_processed_bam(path):
    with gzip.open(path, "rb") as fl:
        return pickle.load(fl)
