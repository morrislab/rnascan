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
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from rnascan.pfmutil import norm_pfm
import subprocess
import numpy as np


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

def get_structure_probability_matrix_for_sequence(id,seq,frag_length,overlap):
    aligned_annotated_sequences = []
    
    for i in xrange(-frag_length/2,len(seq)-frag_length/2,frag_length-overlap):
        temphandle = tempfile.NamedTemporaryFile(delete=False,mode='r+b') # for centroid structure
        temphandle2 = tempfile.NamedTemporaryFile(delete=False,mode='r+b') # for output of structure parser
        
        # set up sequence fragment record & format using SeqIO
        realstart = i
        if i<0:
            realstart = 0
        
        subseq = seq[realstart:(i + frag_length)]
        #print >> sys.stderr, subseq
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
        #print >> sys.stderr, centroid_struct
        temphandle.write(centroid_struct+'\n')
        temphandle.close()
        
        # translate the centroid structure into structural context alphabet
        parse_args = ['parse_secondary_structure', temphandle.name, temphandle2.name]
        #print parse_args
        parse_structure_proc = subprocess.Popen(parse_args)
        parse_structure_proc.wait()
        
        annotated_struct = temphandle2.readline()
        annotated_struct.rstrip()
        #print >> sys.stderr, annotated_struct
        
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

def get_structure_probability_matrix_from_probabilities(id,seq,frag_length):
    input_record = SeqRecord(seq,id=id,description="")
    
    alphabet = ['E','H','I','M']
    programs = {'E':'E_RNAplfold_nolunp', 'H':'H_RNAplfold_nolunp', 'I':'I_RNAplfold_nolunp', 'M':'M_RNAplfold_nolunp'}
    
    plfold_args = ["-L",str(frag_length),"-W ",str(frag_length),"-u","1"]
    
    probabilities = {}
    
    for alph,p in programs.iteritems():
        args = [p] + plfold_args
        plfold_proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=None)
        #plfold_output = plfold_proc.communicate(input_record.format("fasta"))[0]
        plfold_output = plfold_proc.communicate(str(seq))[0]
        data = np.fromstring(plfold_output, dtype=float, sep="\t")
        probabilities[alph] = data
    
    # calculated paired probability
    length = len(probabilities[alph[0]])
    sum = np.zeros(length)
    for a in alphabet:
        sum = sum + probabilities[a]
    paired = np.subtract(np.ones(length),sum)
    probabilities['P'] = paired
    
    for a in alphabet+['P']:
        probabilities[a] = probabilities[a].tolist()
    
    return probabilities



def get_centroid_from_RNAfold_output(rnafold_output):
    lines = rnafold_output.split('\n')
    
    #print >> sys.stderr, "::".join(lines)
    
    centroid_line = lines[4]
    structetc = centroid_line.split(' ')
    return structetc[0]


if __name__ == "__main__":
    seq = Seq('GUACUCGAAAAAAUGUCAUGGACCCCUUAAAAUUACUGAGGGGUUCAGAAAAUACCGUGCAAAAGACGAAAAAAGACGAAUUUCAUUUGAUUUAUAUUUUAUAAAUGACUGUUGCAUUAAACAAUAGACCAAUUAUUUCAAUUUAAUAUUCUUUGCAGGAAACUUUCACAAUGGAAUAACGCCACAUAUUCAUUGUAAAGAUGUUGCGUACUUCUCUUACUAAAGGGGCACGGCUAACUGGGACAAGAUUUGUUCAAACAAAGGCCCUUUCGAAGGCAACAUUGACAGAUCUGCCCGAAAGAUGGGAAAAUAUGCCAAACUUAGAACAGAAAGAGAUUGCAGAUAAUUUGACAGAACGUCAAAAGCUUCCAUGGAAAACUCUCAAUAACGAGGAAAUCAAAGCAGCUUGGUACAUAUCCUACGGCGAGUGGGGACCUAGAAGACCUGUACACGGAAAAGGCGAUGUUGCAUUUAUAACUAAAGGAGUAUUUUUAGGGUUAGGAAUCUCAUUUGGGCUCUUUGGUUUAGUGAGACUAUUAGCCAAUCCUGAAACUCCAAAGACUAUGAACAGGGAAUGGCAGUUGAAAUCAGACGAGUAUCUGAAGUCAAAAAAUGCCAAUCCUUGGGGAGGUUAUUCUCAAGUUCAAUCUAAAUAAGUAGACGAGGAAAAUAAAAUUGUUUCGUAUAUUCCGUGUUUGGGGUAUAAGUAGAUUGUUUUCAUAUAUACGCAUUUGGUCUUAGUUCAGUAGGUUGAUUACUUAGUUCCUUGUACCUUCUUCUGCAAAUAUCAUUCAUUGUUACUUCGAAGAAGAAAAAAAAUAAUCAUGGAAAAUUGGAAAAAAAAAAAGUCCAAUCU')
    id = 'seq_1'
    
    get_structure_probability_matrix_from_probabilities(id,seq,40)