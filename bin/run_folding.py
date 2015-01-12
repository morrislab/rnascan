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



import os,sys,csv
sys.path.append("..")
import rnascan
from rnascan.average_structure import get_structure_probabilities_for_sequence
from rnascan.pfmutil import write_pfm
from Bio import SeqIO

infile = sys.argv[1]
outfile = sys.argv[2]

for seq_record in SeqIO.parse(infile, "fasta"):
    seq_id = seq_record.id
    seq = seq_record.seq
    structure = get_structure_probability_matrix_for_sequence(seq_id,seq)
    write_pfm(structure,outfile)

