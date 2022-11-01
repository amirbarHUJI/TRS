import gzip
import pickle
import numpy as np

from pathlib import Path
from collections import namedtuple

from ARES.Globals import Strand
from ARES.DAL.BAM import bam_to_read_ends
from ARES.DAL.Track.Counts import Counts
from ARES.DAL.Track.ChrCounts import ChrCounts
from ARES.DAL.Track import save, load


CountsTuple = namedtuple("CountsTuple", ["Position", "Coverage", "Ratio"])


class Lib(object):
    def __init__(self, name, group, count_file):
        self.__name = name
        self.__group = group
        self.__count_file = count_file
        self.__counts = None

    def load(
        self,
        workdir,
        force_bam,
        min_count,
        window_width,
        ds_distance=70,
        chr_conv_dct={}
    ):
        bam_path = Path(self.__count_file["path"])

        counts_path = workdir.joinpath(bam_path.name + ".counts.npz")

        if force_bam or \
           not counts_path.exists() or \
           not counts_path.exists():
            pos_counts = bam_to_read_ends(
                chr_conv_dct=chr_conv_dct,
                **self.__count_file
            )

            save(counts=pos_counts, path=counts_path)

        else:
            pos_counts = load(path=counts_path)

        coverage_counts = self.compute_pseudo_coverage(
            pos_counts=pos_counts,
            ds_distance=ds_distance
        )

        ratio_counts = self.calculate_ratio(
            pos_counts=pos_counts,
            min_count=min_count,
            window_width=window_width,
            ds_distance=ds_distance
        )

        self.__counts = CountsTuple(
            pos_counts,
            coverage_counts,
            ratio_counts
        )

    def compute_pseudo_coverage(
        self,
        pos_counts,
        ds_distance
    ):
        coverage = Counts()
        coverage_dct = {}

        for chr_name, chr_len in pos_counts.Chromosomes:
            coverage_dct[chr_name] = {}

            for strand in (Strand.FWD, Strand.REV):
                ds_ends_list = np.zeros(chr_len)

                # Compute pseudo coverage efficiently
                # (all array at once instead of per position)
                start, end = (0, ds_distance) if strand == Strand.FWD else (-ds_distance, 0)
                for i in range(start, end + 1):
                    ds_ends_list += pos_counts[
                        chr_name,
                        strand,
                        1 + i,
                        chr_len + i
                    ]

                coverage_dct[chr_name][strand] = ds_ends_list

            coverage.add_chromosome(
                chromosome=chr_name,
                chr_counts=ChrCounts(
                    chr_len=chr_len,
                    fwd_counts=coverage_dct[chr_name][Strand.FWD],
                    rev_counts=coverage_dct[chr_name][Strand.REV])
            )

        return coverage

    def calculate_ratio(
        self,
        pos_counts,
        min_count,
        window_width,
        ds_distance
    ):
        ratio_dct = {}
        ratio = Counts()

        for chr_name, chr_len in pos_counts.Chromosomes:
            ratio_dct[chr_name] = {}

            for strand in (Strand.FWD, Strand.REV):
                pos_list = np.zeros(chr_len)
                ds_ends_list = np.zeros(chr_len)

                # Set ranges for downstream coverage
                if strand == Strand.FWD:
                    start, end = window_width + 1, ds_distance
                else:
                    start, end = -ds_distance, -window_width - 1

                for i in range(start, end + 1):
                    ds_ends_list += pos_counts[
                        chr_name,
                        strand,
                        1 + i,
                        chr_len + i
                    ]

                # Compute local explained ratio efficiently
                # (all array at once instead of per position)
                for i in range(-window_width, window_width + 1):
                    pos_list += pos_counts[chr_name,
                                           strand,
                                           1 + i,
                                           chr_len + i]

                denominator = ds_ends_list + pos_list
                pos_list[denominator < min_count] = 0.

                ratio_array = np.divide(
                    pos_list,
                    denominator,
                    out=np.zeros_like(pos_list),
                    where=(denominator != 0))

                ratio_dct[chr_name][strand] = ratio_array

            ratio.add_chromosome(
                chromosome=chr_name,
                chr_counts=ChrCounts(
                    chr_len=chr_len,
                    fwd_counts=ratio_dct[chr_name][Strand.FWD],
                    rev_counts=ratio_dct[chr_name][Strand.REV])
            )

        return ratio

    def __repr__(self):
        return "Lib(Name='%s', Group='%s')" % (self.Name, self.Group)

    @property
    def Name(self):
        return self.__name

    @property
    def Group(self):
        return self.__group

    @property
    def Counts(self):
        return self.__counts

    @property
    def Chromosomes(self):
        return self.__counts.Position.Chromosomes
