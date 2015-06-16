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


import sys,csv

from math import log

IUPAC_to_pfm = {'A' : {'A':1., 'C':0., 'G':0., 'U':0.},
         'C' : {'A':0., 'C':1., 'G':0., 'U':0.},
         'G' : {'A':0., 'C':0., 'G':1., 'U':0.},
         'U' : {'A':0., 'C':0., 'G':0., 'U':1.},
         'R' : {'A':0.5, 'C':0., 'G':0.5, 'U':0.},
         'Y' : {'A':0., 'C':0.5, 'G':0., 'U':0.5},
         'S' : {'A':0., 'C':0.5, 'G':0.5, 'U':0.},
         'W' : {'A':0.5, 'C':0., 'G':0., 'U':0.5},
         'K' : {'A':0., 'C':0., 'G':0.5, 'U':0.5},
         'M' : {'A':0.5, 'C':0.5, 'G':0., 'U':0.},
         'B' : {'A':0., 'C':(1./3.), 'G':(1./3.), 'U':(1./3.)},
         'D' : {'A':(1./3.), 'C':0., 'G':(1./3.), 'U':(1./3.)},
         'H' : {'A':(1./3.), 'C':(1./3.), 'G':0., 'U':(1./3.)},
         'V' : {'A':(1./3.), 'C':(1./3.), 'G':(1./3.), 'U':0.},
         'N' : {'A':0.25, 'C':0.25, 'G':0.25, 'U':0.25} }

RNA_ALPHABET = ['A','C','G','U']

def read_pfm(pfmfile):
    pfm = {}
    with open(pfmfile) as f:
        reader = csv.reader(f, delimiter='\t')
        headerline = reader.next()
        alphabet = headerline[1:]
        for base in alphabet:
            pfm[base] = []
        for row in reader:
            for i in range(len(row)):
                if i==0:
                    continue
                pfm[alphabet[i-1]].append(float(row[i]))
    
    f.close()
    return(pfm)

def write_pfm(pfm,pfmoutfile):
    of = open(pfmoutfile, 'w')
        
    alphabet = sorted(pfm.keys())
    pfm_length = len(pfm[alphabet[0]])
    
    of.write('PO')
    for base in alphabet:
        of.write('\t' + base)
    of.write('\n')

    for pos in range(0,pfm_length):
        if pfm[alphabet[0]][pos] is None:
            break
        of.write(str(pos))
        for base in alphabet:
            of.write('\t' + str(pfm[base][pos]))
        of.write('\n')

    of.close()

def norm_pfm(pfm):
    alphabet = sorted(pfm.keys())
    pfm_length = len(pfm[alphabet[0]])

    pfm_normalized = {}
    for base in alphabet:
        pfm_normalized[base] = [None] * (pfm_length)

    for pos in range(0,pfm_length):
        sum = 0
        for base in alphabet:
            sum = sum + pfm[base][pos]
        for base in alphabet:
            normalized = pfm[base][pos] / float(sum)
            pfm_normalized[base][pos] = normalized
    return pfm_normalized


def is_normalized(pfm,epsilon = 1e-6):
    alphabet = sorted(pfm.keys())
    pfm_length = len(pfm[alphabet[0]])
    
    for pos in range(0,pfm_length):
        sum = 0
        for base in alphabet:
            sum = sum + pfm[base][pos]
        if abs(1-sum) > epsilon:
            return False
    return True

def pfm_from_IUPAC(iupac):
    pfm = {}
    alphabet = RNA_ALPHABET
    for base in alphabet:
        pfm[base] = [0.]*len(iupac)
    for i,char in enumerate(iupac):
        freqs = IUPAC_to_pfm[char]
        for base in alphabet:
            pfm[base][i] = freqs[base]
    return(pfm)

def pfm_from_string(string, alphabet):
    pfm = {}
    for base in alphabet:
        pfm[base] = [0.]*len(string)
    for i,base in enumerate(string):
        if base not in alphabet:
            raise Exception("char "+base+" not in alphabet "+str(alphabet))
        pfm[base][i] = 1.
    return(pfm)

def pfm_to_pwm(pfm, num_sites):
    # assume equal background model for all bases (ie .25/.25/.25/.25 for nucleotides)
    alphabet = sorted(pfm.keys())
    num_bases = len(alphabet)
    
    pfm_length = len(pfm[alphabet[0]])
    
    pwm = {}
    for base in alphabet:
        pwm[base] = [None] * (pfm_length)
    
    for pos in range(0,pfm_length):
        for base in alphabet:
            corrected_p = (float(pfm[base][pos])*num_sites + 1./num_bases) / (num_sites+1)
            pwm[base][pos] = log(corrected_p / (1./num_bases), 2)
    
    return pwm

def pwm_scan_fwd(pwm,seq):
    #score forward strand only
    alphabet = sorted(pwm.keys())
    pwm_length = len(pwm[alphabet[0]])
    
    pos = 0
    FwdScores = []
    
    while pos < (len(seq) - pwm_length + 1):
            score = 0
            subseq = seq[pos:(pos + pwm_length)]
            for subpos in range(0, len(subseq)):
                    base = subseq[subpos]
                    score = score + pwm[base][subpos]
            FwdScores.append(score)
            pos = pos + 1
    return FwdScores

if __name__ == "__main__":
    iupac = "AUGN"
    pfm = pfm_from_IUPAC(iupac)
    #pfm2 = pfm_from_string(iupac,RNA_ALPHABET)
    
    print pfm
    #print pfm2