#!/usr/bin/env python
# Jan 29, 2016
#
# Calculates motif scores for all PWMs in a given set of sequences in FASTA format
#
# This tool was motivated by the RNA RBP motif scanning tool from CISBP-RNA:
# http://cisbp-rna.ccbr.utoronto.ca/TFTools.php

import sys
import time
import glob
import os.path as op
from optparse import OptionParser, OptionGroup
from collections import defaultdict
#import multiprocessing
import pandas as pd
from Bio import motifs, SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def getoptions():
    usage = "usage: python %prog [options] sequences.fa"
    desc = "Scan sequence for potential RBP binding sites."
    parser = OptionParser(usage = usage, description = desc)
    parser.add_option('-d', type = 'string', dest = "pwm_dir",
    	default = op.dirname(sys.argv[0]) + "/db/pwms_all_motifs",
        help = "Directory of PWMs [%default]")
    parser.add_option('-p', '--pseudocount', type = "float", dest = "pseudocount",
    	default = 0.5,
    	help = "Pseudocount for normalizing PWM. [%default]")
    parser.add_option('-r', '--rbpinfo', type = 'string', dest = 'rbpinfo',
    	default = op.dirname(sys.argv[0]) + "/db/RBP_Information.txt",
    	help = "RBP info for adding meta data to results. [%default]")
    parser.add_option('-t', '--type', type = 'string', dest = 'seqtype',
    	default = "DNA", 
    	help = "Alphabet of input sequence (DNA or RNA). [%default]")
    # parser.add_option('-c', type = "int", default = 1,
    #     dest = "cores", metavar = "CORES",
    #     help = "Number of processing cores [%default]")
    (opts, args) = parser.parse_args()
    
    if len(args) < 1: 
        print >> sys.stderr, "Error: missing input FASTA file\n"
        parser.print_help()
        exit(-1)

    return (opts, args)

def load_motifs(db, *args):
	"""
	Load all motifs from given directory. Will look for *.txt files
	"""
	motifs = {}
	print >> sys.stderr, "Loading motifs ",
	tic = time.time()
	for file in glob.glob(db + "/*.txt")[1:3]:
		try:
			id = op.splitext(op.basename(file))[0]
			motifs[id] = pwm2pssm(file, args[0])
		except:
			continue
		print >> sys.stderr, "\b.",
		sys.stderr.flush()
	toc = time.time()
	print >> sys.stderr, "done in %0.2f seconds!" % (float(toc - tic))

	return motifs

def pwm2pssm(file, pseudocount):
	pwm = pd.read_table(file)
	# Assuming we are doing RNA motif scanning. Need to replace U with T
	# as Biopython's motif scanner only does DNA
	pwm.rename(columns = {'U':'T'}, inplace=True)
	pwm = pwm.drop("Pos", 1).to_dict(orient = 'list')
	pwm = motifs.Motif(alphabet = IUPAC.IUPACUnambiguousDNA(), counts = pwm)
	pwm = pwm.counts.normalize(pseudocount)

	# Can optionally add background, but for now assuming uniform probability
	pssm = pwm.log_odds()

	return(pssm)

def header():
	return ("Input Sequence", "RBP_ID", "Name", "Motif ID", "Gene ID", "Family",
		"Sequence", "From", "To", "Score")

def collect(x, db):
	"""
	Finilize results into a DataFrame for output
	"""

	# Get metadata
	columns = ["RBP_ID", "Motif_ID", "DBID", "RBP_Name", "Family_Name", "RBDs"]
	meta = pd.read_table(db).loc[:, columns]


	hits = pd.DataFrame(x, columns = ['Motif_ID', 'Start', 'End', 'Sequence', 'Score'])

	# Update rows with new data
	final = pd.merge(meta, hits).sort_values(['Start', 'Motif_ID'])

	# Return
	return final


def main():
	(opts, args) = getoptions()

	results = []

	# Load PWMs
	pssms = load_motifs(opts.pwm_dir, opts.pseudocount)

	# Print header
	# print "\t".join(header())

	# Scan in sequence
	print >> sys.stderr, "Scanning sequences ",
	tic = time.time()
	for seqrecord in SeqIO.parse(open(args[0]), "fasta"):

		seq = seqrecord.seq
		if opts.seqtype == "RNA":
			seq = seq.back_transcribe()
		seq.alphabet = IUPAC.IUPACUnambiguousDNA()

		for motif_id, pssm in pssms.iteritems():
			for position, score in pssm.search(seq, threshold = 3, both = False):

				end_position = position + len(pssm.consensus)
				values = [motif_id,
					position, end_position,
					str(seq[position:end_position].transcribe()), 
					score]
				results.append(values)

		print >> sys.stderr, "\b.",
		sys.stderr.flush()
	toc = time.time()
	print >> sys.stderr, "done in %0.2f seconds!" % (float(toc - tic))
	
	# Collect results
	print >> sys.stderr, "Getting metadata and finalizing... ",
	final = collect(results, opts.rbpinfo)
	print final.to_csv(sep="\t", index = False)

if __name__ == '__main__':
	main()
