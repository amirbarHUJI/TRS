import pysam

from ARES.Globals import Strand


class BamReader(object):

    def __init__(
        self,
        path,
        is_reverse=False,
        skip_clipped=False,
        paired_as_single=False,
        ignore_read1=False,
        ignore_read2=False
    ):

        self.__bam_path = path
        self.__is_reverse = is_reverse
        self.__skip_clipped = skip_clipped
        self._paired_as_single = paired_as_single
        self._ignore_read1 = ignore_read1
        self._ignore_read2 = ignore_read2

        self.__bamfile = None

    @property
    def Chromosomes(self):
        return {
            chromosome: self.__bamfile.lengths[i]
            for i, chromosome in enumerate(self.__bamfile.references)
        }

    def _single_pos(self, read):
        strand = Strand.FWD

        if read.is_reverse != self.__is_reverse:
            strand = Strand.REV

        start = read.reference_start
        end = read.reference_end
        chrom = read.reference_name

        return chrom, strand, start, end

    def _paired_pos(self, read1, read2):
        strand = self._get_paired_read_strand(read1)

        start = min(read1.reference_start, read2.reference_start)
        end = max(read1.reference_end, read2.reference_end)

        chrom = read1.reference_name

        return chrom, strand, start, end

    def _get_paired_read_strand(self, read):
        strand = Strand.FWD

        # If the read is on the reverse set it to reverse
        if read.is_reverse:
            strand = Strand.flip(strand)

        # Flip the direction if its read 2
        # (usually read2 is the reverse complement)
        if read.is_read2:
            strand = Strand.flip(strand)

        # Reverse again if the sequencing is Livny
        if self.__is_reverse:
            strand = Strand.flip(strand)

        return strand

    def _paired_as_single_pos(self, read):
        strand = self._get_paired_read_strand(read)

        start = read.reference_start
        end = read.reference_end
        chrom = read.reference_name

        return chrom, strand, start, end

    def __enter__(self):
        self.__bamfile = pysam.AlignmentFile(self.__bam_path, "rb")

        return self

    def __exit__(self, type, value, traceback):
        self.__bamfile.close()

    def __iter__(self):

        paired_dct = {}

        for read in self.__bamfile:

            # Ignore non-primary alignments
            if read.is_secondary or read.is_supplementary:
                continue

            # discard all reads with soft clipped reads.
            # it can produce incorrect read start positions
            if read.cigarstring and ("S" in read.cigarstring) and \
               self.__skip_clipped:
                continue

            # Check if need to ignore certain reads
            if self._paired_as_single:
                if (self._ignore_read1 and read.is_read1) or \
                   (self._ignore_read2 and read.is_read2):
                    continue

            # Single end
            if not read.is_paired:
                # Invalid single read
                if read.is_unmapped:
                    continue

                yield self._single_pos(read)

            # Paired end
            else:
                # Invalid paired end read
                if read.is_unmapped or \
                   read.mate_is_unmapped or \
                   read.is_reverse == read.mate_is_reverse or \
                   not read.is_proper_pair:
                    continue

                if self._paired_as_single:
                    yield self._paired_as_single_pos(read)

                else:
                    other = paired_dct.get(read.query_name)

                    if other:
                        del paired_dct[read.query_name]
                        yield self._paired_pos(read, other)

                    else:
                        paired_dct[read.query_name] = read
