# motif_scan.py

A Python script for scanning RNA-binding protein motifs in a given set of
sequences. The PWMs are downloaded from the CISBP website:
http://cisbp-rna.ccbr.utoronto.ca/.

Has been tuned to scan for RNA motifs from human and mouse, but can be adjusted to scan for DNA
motifs as well. 

*This script is experimental and not fully tested.* 

# Prerequisites

## Python libraries
This program was written for Python 2.7 and uses the following Python libraries
(which needs to be installed if you haven't already):
 - [`pandas`](http://pandas.pydata.org): for handling the results using DataFrames
 - [`biopython`](http://biopython.org): for parsing FASTA, PWMs, and performing motif scanning

Alternatively, both of the above pacakges can be installed via the 
[Anaconda](https://www.continuum.io/why-anaconda) distribution.

## PWMs

CISBP PWMs and RBP Info files can be downloaded from
http://cisbp-rna.ccbr.utoronto.ca/bulk.php and saved in the default sub-directory:
`db`. 

Steps:

 1. In the `motif_scan` program folder, create a new sub-directory called `db`:

	```
	> mkdir db
	```
 2. Download PWMs and RBP Info files from http://cisbp-rna.ccbr.utoronto.ca/bulk.php (choose "Download Entire Datasets Archive" or alternatively choose a specific species).
 3. Unzip the downloaded file and save the contents inside the newly created `db` directory:
 
 	```
 	> unzip -d db entiredata_2016_01_30_11_56_pm.zip 
 	> ls -l db
	total 49296
	-rw-r--r--@   1 kevinha  staff    10M 10 Apr  2013 RBP_Information.txt
	-rw-r--r--@   1 kevinha  staff    14M 10 Apr  2013 RBP_Information_all_motifs.txt
	-rw-r--r--@   1 kevinha  staff   3.9K 30 Jan 23:56 README.txt
	drwxr-xr-x@ 359 kevinha  staff    12K 30 Jan 23:56 pwms
	```
	
 4. If you chose a specific species (note: currently only support human and mouse), the PWMs sub-folder might not be named the same way as above (e.g. it might be `pwms_all_motifs`). If so, you should rename it to the default (or use `-d` to point to this PWMs folder):
 	```
 	> cd db
 	> mv pwms_all_motifs pwms
 	```

By default, the script will look for the PWMs saved in `db/pwms` (specified by `-d`) as
 well as the file `RBP_Information_all_motifs.txt` (specified by `-r`).

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
python motif_scan.py -s AGTTCCGGTCCGGCAGAGATCGCGGAGAGACGCAGAACGCAGCCCGCTCCT > hits.tab
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

The score is in log (base 2) odds ratio. 

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
