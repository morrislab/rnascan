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

def multi_pfm_iter(filename):
    fh = open(filename)
    
    id = ''
    alphabet = []
    pfm = {}
    
    for isheader,group in groupby(fh, lambda line: line[0] == "#"):
        if isheader:
            headerlines = [x.rstrip()[1:] for x in group]
            id = headerlines[0]
            alphabet = str.split(headerlines[1],"\t")
        else:
            for base in alphabet:
                pfm[base] = []
            # read the pfm from the rest of the non-header
            rows = [x.rstrip() for x in group]
            for row in rows:
                row = str.split(row,'\t')
                for i in range(len(row)):
                    if i==0:
                        continue #skip the position column
                    pfm[alphabet[i-1]].append(float(row[i]))
            yield id,pfm

def write_multi_pfm(idlist,pfmlist,outfile):
    of = open(outfile,'w')
    for (id,pfm) in zip(idlist,pfmlist):
        alphabet = sorted(pfm.keys())
        pfm_length = len(pfm[alphabet[0]])
        of.write('#'+id+'\n')
        of.write('#PO')
        for base in alphabet:
            of.write('\t'+base)
        of.write('\n')
        for pos in range(0,pfm_length):
            if pfm[alphabet[0]][pos] is None:
                break
            of.write(str(pos))
            for base in alphabet:
                of.write('\t' + str(pfm[base][pos]))
            of.write('\n')
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
