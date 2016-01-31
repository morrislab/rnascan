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
import multiprocessing
import pandas as pd
from Bio import motifs, SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def getoptions():
    usage = "usage: python %prog [options] sequences.fa"
    desc = "Scan sequence for potential RBP binding sites."
    parser = OptionParser(usage = usage, description = desc)
    parser.add_option('-d', type = 'string', dest = "pwm_dir",
    	default = op.dirname(sys.argv[0]) + "/db/pwms",
        help = "Directory of PWMs [%default]")
    parser.add_option('-p', '--pseudocount', type = "float", dest = "pseudocount",
    	default = 0,
    	help = "Pseudocount for normalizing PWM. [%default]")
    parser.add_option('-r', '--rbpinfo', type = 'string', dest = 'rbpinfo',
    	default = op.dirname(sys.argv[0]) + "/db/RBP_Information_all_motifs.txt",
    	help = "RBP info for adding meta data to results. [%default]")
    parser.add_option('-t', '--type', type = 'string', dest = 'seqtype',
    	default = "DNA", 
    	help = "Alphabet of input sequence (DNA or RNA). [%default]")
    parser.add_option('-m', '--minscore', type = 'float', dest = 'minscore',
    	default = 6,
    	help = "Minimum score for motif hits. [%default]")
    parser.add_option('-s', '--seq', type = 'string', dest = 'testseq',
    	default = None,
    	help = "Supply a test sequence to scan. FASTA files will be ignored.")
    parser.add_option('-c', type = "int", default = 8,
        dest = "cores", metavar = "CORES",
        help = "Number of processing cores [%default]")
    (opts, args) = parser.parse_args()
    
    if opts.testseq is None and len(args) < 1: 
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
	for file in glob.glob(db + "/*.txt"):
		try:
			id = op.splitext(op.basename(file))[0]
			motifs[id] = pwm2pssm(file, *args)
		except:
			continue
		print >> sys.stderr, "\b.",
		sys.stderr.flush()
	toc = time.time()
	print >> sys.stderr, "done in %0.2f seconds!" % (float(toc - tic))

	return motifs

def pwm2pssm(file, pseudocount):
	"""
	Convert load PWM and covernt it to PSSM (take the log_odds)
	"""
	pwm = pd.read_table(file)
	# Assuming we are doing RNA motif scanning. Need to replace U with T
	# as Biopython's motif scanner only does DNA
	pwm.rename(columns = {'U':'T'}, inplace=True)
	pwm = pwm.drop("Pos", 1).to_dict(orient = 'list')
	pwm = motifs.Motif(alphabet = IUPAC.IUPACUnambiguousDNA(), counts = pwm)
	pwm = pwm.counts.normalize(pseudocount)

	# Can optionally add background, but for now assuming uniform probability
	pssm = pwm.log_odds()

	# Replace negative infinity values with very low number
	for letter, odds in pssm.iteritems():
		pssm[letter] = [-10**6 if x == -float("inf") else x for x in odds]

	return(pssm)

def collect(x, db):
	"""
	Finilize results into a DataFrame for output
	"""

	# Get metadata
	columns = ["RBP_ID", "Motif_ID", "DBID", "RBP_Name", "RBP_Status", "Family_Name", "RBDs", "RBP_Species"]
	meta = pd.read_table(db).loc[:, columns]
	meta = meta[meta['RBP_Species'].isin(['Homo_sapiens', 'Mus_musculus'])]
	meta['RBP_Name'] = meta['RBP_Name'].str.upper()
	
	# Create DataFrame from motif hits
	hits = pd.DataFrame(x, columns = ['Motif_ID', 'Start', 'End', 'Sequence', 'Score'])

	# Merge metadata with hits
	return pd.merge(meta, hits).sort_values(['Start', 'Motif_ID'])

def scan(pssm, seq, minscore, motif_id):
	results = []
	for position, score in pssm.search(seq, threshold = minscore, both = False):
		end_position = position + len(pssm.consensus)
		values = [motif_id,
			position + 1, end_position,
			str(seq[position:end_position].transcribe()), 
			round(score, 3)]
		results.append(values)
	return results

def scan_all(pssms, seq, opts):
	"""
	Scan seq for all motifs in pssms
	"""
	hits = []
	tasks = []

	p = multiprocessing.Pool(opts.cores)
	for motif_id, pssm in pssms.iteritems():
		tasks.append((pssm, seq, opts.minscore, motif_id,))
	results = [p.apply_async(scan, t) for t in tasks]

	for r in results:
		hits.extend(r.get())

	# Collect results
	print >> sys.stderr, "Getting metadata and finalizing... ",
	final = collect(hits, opts.rbpinfo)
	print >> sys.stderr, "done"
	return final

def main():
	(opts, args) = getoptions()

	# Load PWMs
	pssms = load_motifs(opts.pwm_dir, opts.pseudocount)

	if opts.testseq is not None:
		if opts.seqtype == 'RNA':
			seq = Seq(opts.testseq, IUPAC.IUPACUnambiguousRNA()).back_transcribe()
		else:
			seq = Seq(opts.testseq, IUPAC.IUPACUnambiguousDNA())
		final = scan_all(pssms, seq, opts)
		print final.to_csv(sep="\t", index = False)
	else:
		# Scan in sequence
		print >> sys.stderr, "Scanning sequences ",
		tic = time.time()
		for seqrecord in SeqIO.parse(open(args[0]), "fasta"):

			seq = seqrecord.seq
			if opts.seqtype == "RNA":
				seq = seq.back_transcribe()
			seq.alphabet = IUPAC.IUPACUnambiguousDNA()

			final = scan_all(pssms, seq, opts)
			print final.to_csv(sep="\t", index = False)

		toc = time.time()
		print >> sys.stderr, "done in %0.2f seconds!" % (float(toc - tic))
		
if __name__ == '__main__':
	main()
