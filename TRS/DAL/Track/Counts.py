class Counts(object):

    def __init__(self):
        self.__counts = {}

    @property
    def Chromosomes(self):
        return [(name, self.__counts[name].ChromosomeLength) for name in self.__counts]

    def add_chromosome(self, chromosome, chr_counts):
        self.__counts[chromosome] = chr_counts

    def rename_chr(self, old, new):
        self.__counts[new] = self.__counts[old]
        del self.__counts[old]

    def get_region(self, chromosome, strand, start, end):
        return self.__counts[chromosome][strand, start, end]

    def __getitem__(self, key):
        chromosome, strand, start, end = key

        return self.get_region(chromosome=chromosome,
                               strand=strand,
                               start=start,
                               end=end)
