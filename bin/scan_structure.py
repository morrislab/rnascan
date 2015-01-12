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
rna_structure = read_pfm(infile_struct)

scores = structure_scan(pfm_struct,rna_structure,struct_bg_probs,p_pfm,p_bg)
ofhandle_struct.write("\t".join(str(x) for x in scores)+"\n")

ofhandle_struct.close()