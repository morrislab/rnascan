import os,sys,csv
sys.path.append("..")
import rnascan
from rnascan.average_structure import get_structure_probabilities_for_sequence
from rnascan.pfmutil import write_pfm
from Bio import SeqIO

infile = sys.argv[1]

for seq_record in SeqIO.parse(infile, "fasta"):
    seq_id = seq_record.id
    seq = seq_record.seq
    structure = get_structure_probability_matrix_for_sequence(seq_id,seq)
    write_pfm(structure,seq_id+"_structure.tmp")
