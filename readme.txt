============
Introduction
============


rnascan is a (mostly) python tool to scan RNA sequences and secondary structures with sequence and secondary structure PFMs. Secondary structure is represented as weights in different secondary structure contexts, similar to how a PFM represents weights of different nucleotides or amino acids. This allows representation and use of secondary structures in a way that is similar to how PFMs are used to scan nucleotide sequences, and also allows for some flexibility in the structure, as you might find in the boltzmann distribution of secondary structures.

The secondary structure alphabet is as follows:
B - bulge loop
E - external (unpaired) RNA
H - hairpin loop
L - left paired RNA (i.e., a '(' in dot-bracket format)
M - multiloop
R - right paired RNA (i.e., a ')' in dot-bracket format)
T - internal loop

============
Installation
============

To install rnascan:

1. Install BioPython if not already installed
2. Compile secondary structure parser:

cd lib/
g++ parse_secondary_structure.cpp -o parse_secondary_structure

=========
Execution
=========

rnascan consists of two tools:

1. a tool (run_folding.py which calls average_structure.py) to calculate the structural context profile of an RNA sequence by folding overlapping 100nt subsequences, then averaging across 

To run the secondary structure averaging:

cd bin/
python run_folding.py ../example/HIST2H3C_3p_end.fa HIST2H3C_3p_structure.txt

To run a test script which scans the HIST2H3C sequence with sequence and structural PFMs representing the binding specificity of SLBP:

(in bin/)
python run.py <infile_sequence> <infile_structure> <outfile_sequence_scan> <outfile_structure_scan> <pfm_file_seq> <pfm_file_struct> <structural_background> <prior_pfm> <prior_bg>

for example:

python run.py ../example/HIST2H3C_3p_end.fa ../example/HIST2H3C_3p_end_structure.txt test_output_seq.tmp test_output_struct.tmp ../example/SLBP_pfm_assembled_normalized_seq.txt ../example/SLBP_pfm_assembled_normalized_struct.txt ../example/3p_UTR_background_structural_context.txt 0.001 0.999

