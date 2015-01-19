# Copyright (C) 2014-2015 Kate Cook
#
# This file is part of rnascan.
# 
# rnascan is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# rnascan is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with rnascan.  If not, see <http://www.gnu.org/licenses/>.


from Bio.Seq import Seq
from Bio import SeqIO
import sys,csv
import argparse
import numpy as np
sys.path.append("..")
import rnascan
from rnascan.pfmutil import read_pfm
from rnascan.scanner import sequence_scan,structure_scan

parser = argparse.ArgumentParser(description='Run sequence scanning on RNA sequence and structure profile.')
parser.add_argument("--in_seq", type=str, help="name of file containing input RNA sequence to fold, in fasta format (required)", required=True)
parser.add_argument("--in_struct", type=str, help="name of file containing input RNA structure profile (required)", required=True)
parser.add_argument("-o","--outfile", type=str, help="filename for output (required)", required=True)
parser.add_argument("--pfm_seq", type=str, help="name of file containing sequence PFM (required)", required=True)
parser.add_argument("--pfm_struct", type=str, help="name of file containing structure PFM (required)", required=True)
parser.add_argument("--bg_seq_probs", type=str, help="name of file containing background sequence probabilities (optional, default is 0.25/all bases)", required=False)
parser.add_argument("--bg_struct_probs", type=str, help="name of file containing background structural probabilities (required)", required=True)
parser.add_argument("-p","--p_pfm", type=float, default=0.001, help="prior probability of pfm binding (optional, default: 0.001)", required=False)
parser.add_argument("-w","--weight_seq", type=float, default=0.5, help="weighting of sequence (optional, default: 0.5)", required=False)

args = parser.parse_args()

# input

pfm_seq = read_pfm(args.pfm_seq)
pfm_struct = read_pfm(args.pfm_struct)

seq_bg_probs = {}
if args.bg_seq_probs is None:
    seq_bg_probs = {'A':0.25, 'C':0.25, 'G':0.25, 'U':0.25, 'T':0.25}
else:
    bgprobreader = csv.reader(open(args.bg_seq_probs,'r'),delimiter='\t')
    for row in bgprobreader:
        seq_bg_probs[row[0]] = float(row[1])
if 'T' in seq_bg_probs:
    seq_bg_probs['U'] = seq_bg_probs['T']
elif 'U' in seq_bg_probs:
    seq_bg_probs['T'] = seq_bg_probs['U']
else:
    raise Exception("Background sequence file contains neither T or U, check that this isn't the structural context background")


bgprobreader = csv.reader(open(args.bg_struct_probs,'r'),delimiter='\t')
struct_bg_probs = {}
for row in bgprobreader:
    struct_bg_probs[row[0]] = float(row[1])

seq_weight = args.weight_seq
if seq_weight < 0 or seq_weight > 1:
    raise Exception("seq weight out of range: "+str(seq_weight))
struct_weight = 1 - seq_weight

p_pfm = args.p_pfm
if p_pfm < 0 or p_pfm > 1:
    raise Exception("PFM probability out of range: "+str(p_pfm))
p_bg = 1 - p_pfm

# file handles 

inhandle_seq = open(args.in_seq)
ofhandle = open(args.outfile,'w')

# run sequence scan

record = next(SeqIO.parse(inhandle_seq, "fasta"))
seqobj = record.seq
rna = seqobj.transcribe()
scores_seq = np.array(sequence_scan(pfm_seq,rna,seq_bg_probs,p_pfm,p_bg))

inhandle_seq.close()

# run structure scan

rna_structure = read_pfm(args.in_struct)

scores_struct = np.array(structure_scan(pfm_struct,rna_structure,struct_bg_probs,p_pfm,p_bg))

# combine scores using weighting

scores_seqstruct = scores_seq * seq_weight + scores_struct * struct_weight

ofhandle.write("\t".join(str(x) for x in scores_seqstruct)+"\n")

ofhandle.close()