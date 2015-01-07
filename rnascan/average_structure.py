import os,sys,csv
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from rnascan.pfmutil import norm_pfm
import subprocess

def struct_pfm_from_aligned(sequences):
    #alphabet = ['B','E','H','I','L','M','R','T']
    alphabet = ['B','E','H','L','M','R','T']
    length = len(sequences[0])
    
    counts = {}
    for base in alphabet:
        counts[base] = [0] * length

    for seq in sequences:
        for index, char in enumerate(seq):
            if char != '-':
                counts[char][index] = counts[char][index] + 1
    
    return counts

def get_structure_probability_matrix_for_sequence(id,seq):
    frag_length = 100
    overlap = 95
    
    aligned_annotated_sequences = []
    
    for i in xrange(-frag_length/2,len(seq)-frag_length/2,frag_length-overlap):
        temphandle = tempfile.NamedTemporaryFile(delete=False,mode='r+b') # for centroid structure
        temphandle2 = tempfile.NamedTemporaryFile(delete=False,mode='r+b') # for output of structure parser
        
        # set up sequence fragment record & format using SeqIO
        realstart = i
        if i<0:
            realstart = 0
        
        subseq = seq[realstart:(i + frag_length)]
        #print i,i+frag_length,subseq
        frag_id = id+"_frag_"+str(i)
        input_record = SeqRecord(subseq,id=frag_id,description="")
        
        # call RNAfold and pipe in sequence fragment as fasta
        rnafold_args = ["RNAfold","-p","--noPS"]
        rnafold_proc = subprocess.Popen(rnafold_args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        rnafold_output = rnafold_proc.communicate(input_record.format("fasta"))[0]
        #print rnafold_output
        
        # process output to get the actual centroid structure & write it to a temp file
        os.system("rm -f *.ps") # remove ps files created by RNAfold because it is dumb
        centroid_struct = get_centroid_from_RNAfold_output(rnafold_output)
        temphandle.write(centroid_struct+'\n')
        temphandle.close()
        
        # translate the centroid structure into structural context alphabet
        parse_args = ['../lib/parse_secondary_structure', temphandle.name, temphandle2.name]
        print parse_args
        parse_structure_proc = subprocess.Popen(parse_args)
        parse_structure_proc.wait()
        
        annotated_struct = temphandle2.readline()
        annotated_struct.rstrip()
        print annotated_struct
        
        os.remove(temphandle.name) # remove centroid structure file
        os.remove(temphandle2.name) # remove annotated structure file
        
        # generate aligned fragment string with - for gaps
        start_gap = "-"*i
        end_gap = "-"*( len(seq) - (i + frag_length) )
        aligned = start_gap + annotated_struct.rstrip() + end_gap    
        aligned_annotated_sequences.append(aligned)

    
    # make count pfm & then normalize
    counts = struct_pfm_from_aligned(aligned_annotated_sequences)
    normalized = norm_pfm(counts)
    
    return normalized


def get_centroid_from_RNAfold_output(rnafold_output):
    lines = rnafold_output.split('\n')
    centroid_line = lines[2]
    structetc = centroid_line.split(' ')
    return structetc[0]