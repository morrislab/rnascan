from Bio.Seq import Seq
from Bio import SeqIO
import sys,csv
sys.path.append("..")
import rnascan
from rnascan.pfmutil import read_pfm
from rnascan.scanner import sequence_scan,structure_scan

infile_seq = sys.argv[1]
infile_struct = sys.argv[2]
outfile_seq = sys.argv[3]
outfile_struct = sys.argv[4]
pfm_file_seq = sys.argv[5]
pfm_file_struct = sys.argv[6]
struct_bg_probs_file = sys.argv[7]
p_pfm = float(sys.argv[8])
p_bg = float(sys.argv[9])

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