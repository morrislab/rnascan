# motif_scan

A Python script for scanning RNA-binding protein motifs in a given set of
sequences. The PWMs are downloaded from the CISBP website:
http://cisbp-rna.ccbr.utoronto.ca/

Has been tuned to scan for RNA motifs, but can be adjusted to scan for DNA
motifs as well. 

*This script is experimental and not fully tested.* 

# Prerequisites

## Python libraries
This program was written for Python 2.7 and uses the following Python libraries
(which needs to be installed if you haven't already):
 - [`pandas`](http://pandas.pydata.org): for handling the results using DataFrames
 - `[`biopython`](http://biopython.org): for parsing FASTA, PWMs, and performing motif scanning

 Alternatively, both of the above pacakges can be installed via the 
 [Anaconda](https://www.continuum.io/why-anaconda) distribution.

## PWMs

CISBP PWMs and RBP Info files can be downloaded from
http://cisbp-rna.ccbr.utoronto.ca/bulk.php and saved in the default directory:
`db`. Expects the following files:
     1. PWM files (default save location: `db/pwms`)
     1. RBP metadata (default: `db/RBP_Information_all_motifs.txt`)

# Usage

```
python motif_scan.py sequences.fasta > hits.tab
```
Parallelization is implemented via Python's [`multiprocessing`](https://docs.python.org/2/library/multiprocessing.html) module:
```
python motif_scan.py -c 8 sequences.fasta > hits.tab
```

To run a test sequence:

```
python motif_scan.py -s AGTTCCGGTCCGGCAGAGATCGCGGAGAGACGCAGAACGCAGCCCGCTCCT >
hits.tab
```

# Output

Results are sent to standard output by default.

The output is tab-delimited and the format resembles the output of the CISBP RNA
Scan tool:
```
RBP_ID	Motif_ID	DBID	RBP_Name	RBP_Status	Family_Name	RBDs	RBP_Species	Start	End	Sequence	Score
T37299_0.6	M352_0.6	ENSG00000161547	SRSF2	D	RRM	RRM_1	Homo_sapiens	12	20	GGCAGAGAU	7.421
T37328_0.6	M352_0.6	ENSG00000180771	ENSG00000180771	I	RRM	RRM_1	Homo_sapiens	12	20	GGCAGAGAU	7.421
...
```

The score is in log (base 2) odds. 

# References

 - Ray D, Kazan H, Cook KB, Weirauch MT, Najafabadi HS, Li X, Gueroussov S, Albu
   M, Zheng H, Yang A, Na H, Irimia M, Matzat LH, Dale RK, Smith SA, Yarosh CA,
   Kelly SM, Nabet B, Mecenas D, Li W, Laishram RS, Qiao M, Lipshitz HD, Piano
   F, Corbett AH, Carstens RP, Frey BJ, Anderson RA, Lynch KW, Penalva LO, Lei
   EP, Fraser AG, Blencowe BJ, Morris QD, Hughes TR. [A compendium of RNA-binding
   motifs for decoding gene
   regulation](http://www.nature.com/nature/journal/v499/n7457/full/nature12311.html). Nature. 2013 Jul 11;499(7457):172-7.
   doi: 10.1038/nature12311. PubMed PMID: 23846655.
