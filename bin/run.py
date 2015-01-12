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
import numpy as np
sys.path.append("..")
import rnascan
from rnascan.pfmutil import read_pfm
from rnascan.scanner import sequence_scan,structure_scan

infile_seq = sys.argv[1]
infile_struct = sys.argv[2]
outfile = sys.argv[3]
pfm_file_seq = sys.argv[4]
pfm_file_struct = sys.argv[5]
struct_bg_probs_file = sys.argv[6]
p_pfm = float(sys.argv[7])
p_bg = float(sys.argv[8])
seq_weight = float(sys.argv[9])

# input

pfm_seq = read_pfm(pfm_file_seq)
pfm_struct = read_pfm(pfm_file_struct)

seq_bg_probs = {'A':0.25, 'C':0.25, 'G':0.25, 'U':0.25, 'T':0.25}

bgprobreader = csv.reader(open(struct_bg_probs_file,'r'),delimiter='\t')
bgprobreader.next()
struct_bg_probs = {}
for row in bgprobreader:
    struct_bg_probs[row[0]] = float(row[1])

if seq_weight < 0 or seq_weight > 1:
    raise Exception("seq weight out of range: "+str(seq_weight))

struct_weight = 1 - seq_weight

# file handles 

inhandle_seq = open(infile_seq)
ofhandle = open(outfile,'w')

# run sequence scan


record = next(SeqIO.parse(inhandle_seq, "fasta"))
seqobj = record.seq
rna = seqobj.transcribe()
scores_seq = np.array(sequence_scan(pfm_seq,rna,seq_bg_probs,p_pfm,p_bg))

inhandle_seq.close()

# run structure scan

rna_structure = read_pfm(infile_struct)

scores_struct = np.array(structure_scan(pfm_struct,rna_structure,struct_bg_probs,p_pfm,p_bg))

scores_seqstruct = scores_seq * seq_weight + scores_struct * struct_weight

ofhandle.write("\t".join(str(x) for x in scores_seqstruct)+"\n")

ofhandle.close()