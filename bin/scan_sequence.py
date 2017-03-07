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
sys.path.append("..")
import rnascan
from rnascan.pfmutil import read_pfm
from rnascan.scanner import sequence_scan,structure_scan

infile_seq = sys.argv[1]
outfile_seq = sys.argv[2]
pfm_file_seq = sys.argv[3]
p_pfm = float(sys.argv[4])
p_bg = float(sys.argv[5])

# run sequence scan

pfm_seq = read_pfm(pfm_file_seq)

seq_bg_probs = {'A':0.25, 'C':0.25, 'G':0.25, 'U':0.25, 'T':0.25}


ofhandle_seq = open(outfile_seq,'w')
inhandle_seq = open(infile_seq)

for record in SeqIO.parse(inhandle_seq, "fasta"):
    seqobj = record.seq
    rna = seqobj.transcribe()
    scores = sequence_scan(pfm_seq,rna,seq_bg_probs,p_pfm,p_bg)
    ofhandle_seq.write("\t".join(str(x) for x in scores)+"\n")

ofhandle_seq.close()
inhandle_seq.close()

