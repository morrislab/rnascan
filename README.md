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
 - `Biopython`: for parsing FASTAr, PWMs, and performing motif scanning

# Input Requirements

 1. Directory containing PWMs from CISBP (default directory: `db/pwms_all_motifs/`)
 1. RBP metadata (default: `db/RBP_Information.txt`)
 1. Sequence(s) in FASTA format

# Steps

 1. Load PWMs
 1. Convert to a PSSM (log odds) using pseudocount to avoid infinite values
 1. For each input sequence, scan using each PSSMs and get the log odds ratio
 1. Combine results with metadata to generate complete table of results

# Usage

```
python motif_scan.py sequences.fasta
```


