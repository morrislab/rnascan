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

import sys
import argparse
from Bio import SeqIO
sys.path.append("..")
import rnascan
from rnascan.average_structure import get_structure_probability_matrix_for_sequence
from rnascan.pfmutil import write_pfm

# old arguments
# infile = sys.argv[1]
# outfile = sys.argv[2]
# window_size = sys.argv[3]
# overlap_size = sys.argv[4]
parser = argparse.ArgumentParser(description='Calculate the average structural profile of an RNA sequence.')
parser.add_argument("infile",help="name of file containing input RNA sequence to fold, in fasta format")
parser.add_argument("outfile",help="filename for output")
parser.add_argument("-w","--window_size", type=int, default=100, help="size of folding window (in nt) default: 100")
parser.add_argument("-o","--overlap_size", type=int, default=95, help="overlap between folding windows (in nt) default: 95")


args = parser.parse_args()

#print "window size: "+str(args.window_size)
#print "overlap size: "+str(args.overlap_size)

for seq_record in SeqIO.parse(args.infile, "fasta"):
    seq_id = seq_record.id
    seq = seq_record.seq
    structure = get_structure_probability_matrix_for_sequence(seq_id,seq,args.window_size,args.overlap_size)
    write_pfm(structure,args.outfile)

