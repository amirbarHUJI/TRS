import numpy as np
from collections import defaultdict

from TRS.Globals import Strand
from TRS.DAL.Track.Counts import Counts
from TRS.DAL.Track.ChrCounts import ChrCounts


def save(counts, path):
    np_dct = {}

    for chromosome, length in counts.Chromosomes:
        for strand in (Strand.FWD, Strand.REV):
            data_id = "%s_%s" % (chromosome, strand.value)
            np_dct[data_id] = counts[chromosome, strand, 1, length]

    np.savez_compressed(str(path), **np_dct)


def load(path):
    try:
        np_dct = np.load(path)
        counts = Counts()
        chr_dct = defaultdict(dict)

        for identifier, strand_counts in np_dct.items():
            chromosome, strand = tuple(identifier.rsplit("_", 1))
            strand = Strand(strand)

            chr_dct[chromosome][strand] = strand_counts

        for chromosome in chr_dct:
            fwd_counts = chr_dct[chromosome][Strand.FWD]
            rev_counts = chr_dct[chromosome][Strand.REV]

            counts.add_chromosome(
                chromosome=chromosome,
                chr_counts=ChrCounts(
                    chr_len=len(fwd_counts),
                    fwd_counts=fwd_counts,
                    rev_counts=rev_counts
                )
        )

    # TODO: raise specific exception to handle above
    except KeyError:
        return None

    except OSError:
        return None

    return counts
