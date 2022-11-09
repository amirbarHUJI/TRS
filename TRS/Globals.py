from enum import Enum
from collections import namedtuple

Peak = namedtuple("Peak", ["Chromosome", "Strand", "Start", "End", "Signal", "Dominant"])


class Strand(Enum):
    FWD = "+"
    REV = "-"
    ANY = "*"

    @classmethod
    def flip(cls, strand):
        """Return the opposite strand

        Args:
            strand (Strand): The strand to flip

        Returns:
            Strand: The strand on the opposite side
        """
        if strand == Strand.FWD:
            return Strand.REV

        elif strand == Strand.REV:
            return Strand.FWD

        return Strand.ANY

    @classmethod
    def has_value(cls, value):
        return value in [e.value for e in Strand]

    @classmethod
    def values(cls):
        return [e.value for e in Strand]
