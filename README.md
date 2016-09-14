# motif_scan.py

A Python script for scanning RNA-binding protein (RBP) motifs in a given set of
sequences. The PWMs are downloaded from the CISBP website:
http://cisbp-rna.ccbr.utoronto.ca/.

This program has been develop for scanning RNA motifs *and* secondary structure motifs. This is done by combining RNA and secondary struture PFMs to
generate a combined scoring for a given sequence.

# Prerequisites

## Python libraries
This program was written for Python 2.7 and uses the following Python libraries
(which needs to be installed if you haven't already):
 - [`pandas`](http://pandas.pydata.org) (v0.17.1 or higher): for handling the results using DataFrames
 - [`biopython`](http://biopython.org) (v1.66 or higher): for parsing FASTA, PWMs, and performing motif scanning

Alternatively, both of the above pacakges can be installed via the
[Anaconda](https://www.continuum.io/why-anaconda) distribution.

## PWMs

CISBP PWMs and RBP Info files can be downloaded from
http://cisbp-rna.ccbr.utoronto.ca/bulk.php.

# Installation

Download the latest source code and install using provided `setup.py` script:

```
git clone git@github.com:kcha/motif_scan.git

cd motif_scan

python setup.py install
# or python setup.py install --user
```

# Usage

For full documentation of options, please refer to the help message:

```
motif_scan -h
```

A minimal usage command:

```
motif_scan -d pwm_dir sequences.fasta > hits.tab
```

Parallelization is implemented via Python's [`multiprocessing`](https://docs.python.org/2/library/multiprocessing.html) module:
```
motif_scan -c 8 -d pwm_df sequences.fasta > hits.tab
```

To run a test sequence:

```
motif_scan -d pwm_dir -s AGTTCCGGTCCGGCAGAGATCGCG > hits.tab
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
