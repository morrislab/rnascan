# motif_scan

A Python script for scanning RNA-binding protein motifs in a given set of
sequences. The PWMs are downloaded from the CISBP website:
http://cisbp-rna.ccbr.utoronto.ca/

Has been tuned to scan for RNA motifs, but can be adjusted to scan for DNA
motifs as well. 

# Prerequisites

This program was written for Python 2.7 and uses the following Python libraries
(which needs to be install if you haven't already):
 - `pandas`: for handling the results using DataFrames
 - `Biopython`: for parsing FASTA, PWMs, and performing motif scanning

# Input Requirements

 1. Directory containing PWMs from CISBP (default directory: `db/pwms_all_motifs/`)
 1. RBP metadata (default: `db/RBP_Information.txt`)
 1. Sequence(s) in FASTA format

# Usage

```
# Results sent to standard output
python motif_scan.py sequences.fasta
```

The script is also parallelized via Python's `multiprocessing` module:
```
python motif_scan.py -c 8 sequences.fasta
```

To run a test sequence:

```
python motif_scan.py -s AGTTCCGGTCCGGCAGAGATCGCGGAGAGACGCAGAACGCAGCCCGCTCCT
```

