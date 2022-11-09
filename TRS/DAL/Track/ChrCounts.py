import numpy as np

from TRS.Globals import Strand


class ChrCounts(object):

    def __init__(self, chr_len, fwd_counts, rev_counts):

        self.__chr_len = chr_len

        if len(fwd_counts) != len(rev_counts):
            raise BaseException("Unmatched chromosome length on strands")

        if self.__chr_len != len(fwd_counts):
            raise BaseException("Chromosome len is different then declared")

        self.__chr = {Strand.FWD: fwd_counts,
                      Strand.REV: rev_counts}

    @property
    def ChromosomeLength(self):
        return self.__chr_len

    def get_region(self, strand, start, end):

        def handle_exceeding(coord):
            if coord < 0 or coord > self.__chr_len:
                return coord % self.__chr_len

            return coord

        # Zero based
        start = start - 1

        start = handle_exceeding(start)
        end = handle_exceeding(end)

        if start < end:
            return self.__chr[strand][start:end]

        else:
            return np.concatenate((self.__chr[strand][start:], self.__chr[strand][:end]))

    def __getitem__(self, key):
        strand, start, end = key
        return self.get_region(strand=strand, start=start, end=end)
