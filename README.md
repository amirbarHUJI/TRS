# ReadEndAnalysis

<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#description">Description</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#scheme-file">Scheme file</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
    <li><a href="#references">References</a></li>
  </ol>
</details>

## Description
This algorithm is intended to determine 3' termini from RNATag-seq [[1]](#1) sequencing data. We observed that in RNATag-seq data, sequencing reads accumulate at the 3' termini of transcripts that enables their identification. To do so, we compute for each position in the genome the ..... The computed ratio is acomparable measure between libraries and it signifies the positions where drastic loss in expression occur. We expect positions subjected to random shearing of the RNA in the experimental procedure to resolve in less reproducible signal. To retrieve the reproducible signal, we follow Adams et al. 2021, [[2]](#2) method to analyse Term-seq [[3]](#3) data, which uses peak calling followed by the irreproducibility discovery rate (IDR) procedure [[4]](#4).
The pipeline is fully described in [].


## Getting Started

This section briefly describes how to setup the package and shows how test it works. 

### Prerequisites

This package requires the installation of the following packages:
- numpy
- pandas
- scipy
- statsmodels
- pysam
- pyaml
- intervaltree

### Installation
Either ``pip install ARES`` or ``conda install ARES``

### Files required for running the program
The package requires the following files to determine 3' termini:
- BAM files containing the mapped reads of the sequencing libraries.
- A scheme file describing how to read the BAM files, assigning them to groups, which enables processing sequencing libraries of multiple conditions in a single run (<a href="#scheme-file">see here for further detail</a>).

### Running an example file
COMPLETE

## Scheme File
COMPLETE

## Usage
```
usage: peakcaller.py [-h] [-w WORKDIR] [-t THRESHOLD] [-b] [-c MIN_COUNT] [--min_height MIN_HEIGHT] [--window_margin WINDOW_MARGIN] [--merge_distance MERGE_DISTANCE] [--rel_height REL_HEIGHT]
                     [--chr_list CHR_LIST] [-d DS_DISTANCE] [--signif_min_lib_count SIGNIF_MIN_LIB_COUNT] [--insignif_min_ratio INSIGNIF_MIN_RATIO] [-l LOG_LEVEL]
                     scheme

Uses read starts to determine 3' termini of transcripts.

positional arguments:
  scheme                Path for file containing info of the libs to process.

optional arguments:
  -h, --help            show this help message and exit
  -w WORKDIR, --workdir WORKDIR
                        Working directory. (default: ./runs)
  -t THRESHOLD, --threshold THRESHOLD
                        The signficance level used by the model. (default: 0.01)
  -b, --force_bam       Forces the program to reprocess the bam files. (default: False)
  -c MIN_COUNT, --min_count MIN_COUNT
                        The minimal coverage for a region to be considered. (default: 10)
  --min_height MIN_HEIGHT
                        The minimal ratio to consider as a peak. (default: None)
  --window_margin WINDOW_MARGIN
                        Defines the region which will be used to count the local number of reads starts. The margin is the number of nucleotides (upstream/downstream) that will be added to the region
                        around the considered site. (default: 3)
  --merge_distance MERGE_DISTANCE
                        The distance which below it, peaks will be merged together. (default: 0)
  --rel_height REL_HEIGHT
                        The relative height of peaks to use in the scipy find_peaks function. (default: 0.75)
  --chr_list CHR_LIST   Conversion of chromosome files from the bam to new name in the following format: bam1:new1,bam2:new2,... (default: )
  -d DS_DISTANCE, --ds_distance DS_DISTANCE
                        The distance (downstream) used to compute the read starts ratio. (default: 70)
  --signif_min_lib_count SIGNIF_MIN_LIB_COUNT
                        The minimal number of libraries in which the 3' terminus should be found as significant to report it depending on the threshold. (default: 2)
  --insignif_min_ratio INSIGNIF_MIN_RATIO
                        The minimal mean ratio to accept if the 3' terminus wasn't significant in all repeats. (default: 0.5)
  -l LOG_LEVEL, --log_level LOG_LEVEL
                        The logging level to report to the log file. (default: debug)
```

## Contact
Amir Bar - amir.bar@mail.huji.ac.il

## Acknowledgements
We are thankful to:
- Term-seq peak-caller [[2]](#2) - https://github.com/NICHD-BSPC/termseq-peaks

## License
TODO: add MIT license (and a copy of Adams et al)

## References
- <a id="1">[1]</a>  Shishkin, Alexander A., et al. "Simultaneous generation of many RNA-seq libraries in a single reaction." Nature methods 12.4 (2015): 323-325.
- <a id="2">[2]</a>  Adams, Philip P., et al. "Regulatory roles of Escherichia coli 5'UTR and ORF-internal RNAs detected by 3'end mapping." Elife 10 (2021): e62438.
- <a id="3">[3]</a>  Dar, Daniel, et al. "Term-seq reveals abundant ribo-regulation of antibiotics resistance in bacteria." Science 352.6282 (2016).

