# motif_scan.py

A Python script for scanning RNA-binding protein (RBP) motifs in a given set of
sequences. The PWMs are downloaded from the CISBP website:
http://cisbp-rna.ccbr.utoronto.ca/.

This program has been developed for scanning motifs under three modes:
    1. DNA/RNA motifs
    2. Contextual secondary structure motifs
    3. RNA motifs *and* secondary structure (RNA+structure)

# Prerequisites

## Python libraries
This program was written for Python 2.7 and uses the following Python libraries
(which needs to be installed if you haven't already):
 - [`pandas`](http://pandas.pydata.org) (v0.17 or higher): for handling the results using DataFrames
 - [`numpy`](http://www.numpy.org/) (v1.10 or higher): for numerical computations
 - [`biopython`](http://biopython.org) (v1.66 or higher): for parsing FASTA, PWMs, and performing motif scanning

Alternatively, both of the above pacakges can be installed via the
[Anaconda](https://www.continuum.io/why-anaconda) distribution.

# Installation

Download the latest source code and install using provided `setup.py` script:

```
git clone git@github.com:kcha/motif_scan.git
cd motif_scan
python setup.py install
```

Or for user-specific installation:
```
python setup.py install --user
```

# Usage

## Motif Scanning

For full documentation of options, please refer to the help message:

```
motif_scan -h
```

### Quick Usage
A minimal usage command for scanning RNA (default) motifs:

```
motif_scan -d pwm_dir sequences.fasta > hits.tab
```

Parallelization is implemented via Python's [`multiprocessing`](https://docs.python.org/2/library/multiprocessing.html) module:
```
motif_scan -c 8 -d pwm_dir sequences.fasta > hits.tab
```

To run a test sequence:

```
motif_scan -d pwm_dir -s AGTTCCGGTCCGGCAGAGATCGCG > hits.tab
```

For scanning DNA and secondary structure motifs, use the option `-t` to change
the mode to `DNA` or `SS`, respectively. For RNA+structure, see below.

### Scanning RNA+structure PFMs

For RNA+structure motif scanning, a new PFM must
be computed from the given RNA PFM and secondary structure PFM. This can
be done using the command `combine_pfms`:

```
combine_pfms rna_pfm.txt secondary_structure_pfm.txt > combined_pwm_dir/pfm.txt
```

Next, supply two FASTA sequences when calling `motif_scan`:

    1. RNA sequences
    2. Contextual secondary structure sequences

```
motif_scan -d combined_pwm_dir rna_sequences.fa secondary_structure_sequences.fa
```

# References

 - Ray D, Kazan H, Cook KB, Weirauch MT, Najafabadi HS, Li X, Gueroussov S, Albu
   M, Zheng H, Yang A, Na H, Irimia M, Matzat LH, Dale RK, Smith SA, Yarosh CA,
   Kelly SM, Nabet B, Mecenas D, Li W, Laishram RS, Qiao M, Lipshitz HD, Piano
   F, Corbett AH, Carstens RP, Frey BJ, Anderson RA, Lynch KW, Penalva LO, Lei
   EP, Fraser AG, Blencowe BJ, Morris QD, Hughes TR. [A compendium of RNA-binding
   motifs for decoding gene
   regulation](http://www.nature.com/nature/journal/v499/n7457/full/nature12311.html). Nature. 2013 Jul 11;499(7457):172-7.
   doi: 10.1038/nature12311. PubMed PMID: 23846655.
 - [Biopython tutorial on sequence motif analysis](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc213)
