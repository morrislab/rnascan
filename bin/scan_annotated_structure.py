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

infile_struct = sys.argv[1]
outfile_struct = sys.argv[2]
pfm_file_struct = sys.argv[3]
struct_bg_probs_file = sys.argv[4]
p_pfm = float(sys.argv[5])
p_bg = float(sys.argv[6])

# run structure scan

pfm_struct = read_pfm(pfm_file_struct)

bgprobreader = csv.reader(open(struct_bg_probs_file,'r'),delimiter='\t')
bgprobreader.next()
struct_bg_probs = {}
for row in bgprobreader:
    struct_bg_probs[row[0]] = float(row[1])

ofhandle_struct = open(outfile_struct,'w')
inhandle_struct = open(infile_struct)

for record in SeqIO.parse(inhandle_struct, "fasta"):
    rna_struct = str(record.seq)
    scores = sequence_scan(pfm_struct,rna_struct,struct_bg_probs,p_pfm,p_bg)
    ofhandle_struct.write("\t".join(str(x) for x in scores)+"\n")

ofhandle_struct.close()
inhandle_struct.close()

