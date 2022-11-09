class Coverage(object):

    def __init__(self):
        self.__coverage = {}

    @property
    def Chromosomes(self):
        return [(name, self.__coverage[name].ChromosomeLength) for name in self.__coverage]

    def add_chromosome(self, chromosome, chr_counts):
        self.__coverage[chromosome] = chr_counts

    def rename_chr(self, old, new):
        self.__coverage[new] = self.__coverage[old]
        del self.__coverage[old]

    def get_region(self, chromosome, strand, start, end):
        return self.__coverage[chromosome][strand, start, end]

    def __getitem__(self, key):
        chromosome, strand, start, end = key

        return self.get_region(chromosome=chromosome,
                               strand=strand,
                               start=start,
                               end=end)
