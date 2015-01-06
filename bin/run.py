from Bio.Seq import Seq
from Bio import SeqIO
import sys
sys.path.append("..")
import rnascan
from rnascan.pfmutil import read_pfm
from rnascan.scanner import sequence_scan

seq_file = sys.argv[1]
outfile = sys.argv[2]
pwm_file = sys.argv[3]
p_pfm = float(sys.argv[4])
p_bg = float(sys.argv[5])

bg_probs = {'A':0.25, 'C':0.25, 'G':0.25, 'U':0.25, 'T':0.25}

pwm = read_pfm(pwm_file)

ofhandle = open(outfile,'w')
inhandle = open(seq_file)

for record in SeqIO.parse(inhandle, "fasta"):
	seqobj = record.seq
	rna = seqobj.transcribe()
	scores = sequence_scan(pwm,rna,bg_probs,p_pfm,p_bg)
	ofhandle.write("\t".join(str(x) for x in scores)+"\n")


ofhandle.close()
