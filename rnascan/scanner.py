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
import rnascan.pfmutil as pfmutil
import numpy as np

NUM_STRUCT = 7 # number of secondary structure contexts

def pfm_likelihood_seq(pfm,subseq):
    """ calculates likelihood that a subsequence was generated from a pfm """
    score = 1.
    for subpos in range(0, len(subseq)):
        base = subseq[subpos]
        score = score * pfm[base][subpos]
    return score

def bg_likelihood_seq(bg_probs,subseq):
    """ calculates likelihood that a subsequence was generated from a background distribution """
    score = 1.
    for subpos in range(0, len(subseq)):
        base = subseq[subpos]
        score = score * bg_probs[base]
    return score

def pfm_likelihood_struct(pfm,substructure):
    """ calculates likelihood that a substructure was generated from a pfm """
    pfm_length = pfm.shape[1]
    score = 1.
    for i in range(0, pfm_length):
        pfm_column = pfm[:,i]
        substructure_column = substructure[:,i]
        d = np.dot(pfm_column,substructure_column)
        score = score * d
    return score

def pfm_compare_seq(pfm,subseq,bg_probs,prob_pfm,prob_bg):
    """ calculates probability of a sequence hit given prior probabilities """
    llhood = pfm_likelihood_seq(pfm,subseq)
    bg_llhod = bg_likelihood_seq(bg_probs,subseq)
    probability = llhood * prob_pfm / (llhood*prob_pfm + bg_llhod*prob_bg)
    return probability

def pfm_compare_struct(pfm,substructure,bg_pfm,prob_pfm,prob_bg):
    """ calculates probability of a structure hit given prior probabilities """
    llhood = pfm_likelihood_struct(pfm,substructure)
    bg_llhood = pfm_likelihood_struct(bg_pfm,substructure)
    probability = llhood * prob_pfm / (llhood*prob_pfm + bg_llhood*prob_bg)
    return probability

def sequence_scan(pfm,seq,bg_probs,prob_pfm,prob_bg):
    """ Scan a sequence with a sequence pfm, returning a list of hit probabilities """
    #score forward strand only
    alphabet = sorted(pfm.keys())
    pfm_length = len(pfm[alphabet[0]])
    
    pos = 0
    scores = []
    
    while pos < (len(seq) - pfm_length + 1):
            subseq = seq[pos:(pos + pfm_length)]
            score = pfm_compare_seq(pfm,subseq,bg_probs,prob_pfm,prob_bg)
            scores.append(score)
            pos = pos + 1
    return scores

def structure_scan(pwm,structure,bg_model,prob_pfm,prob_bg):
    """ scan a structure with a structure pfm, returning a list of hit probabilities """
    #score forward strand only
    alphabet = sorted(pwm.keys())
    alphabet_struct = sorted(structure.keys())
    
    assert alphabet == alphabet_struct

    pwm_length = len(pwm[alphabet[0]])
    sequence_length = len(structure[alphabet[0]])
    
    pwm_np = convert_pfm_to_numpy(pwm)
    structure_np = convert_pfm_to_numpy(structure)
    
    pos = 0
    scores = []
    
    bg_probs = [bg_model[x] for x in sorted(bg_model.keys())]
    bg_pwm = np.reshape(np.repeat(bg_probs,pwm_length),[NUM_STRUCT,pwm_length])
    
    while pos < (sequence_length - pwm_length + 1):
        substructure = structure_np[:,pos:(pos+pwm_length)]
        p_hit = pfm_compare_struct(pwm_np,substructure,bg_pwm,prob_pfm,prob_bg)
        scores.append(p_hit)
        pos = pos + 1
    
    return scores

def convert_pfm_to_numpy(pfm):
    """ helper function to convert a (structure) pfm to a numpy array format """
    alphabet = sorted(pfm.keys())
    pfm_length = len(pfm[alphabet[0]])
    
    pfm_array = np.empty([len(alphabet),pfm_length])
    for i,base in enumerate(alphabet):
        for j in range(0, pfm_length):
            pfm_array[i,j] = pfm[base][j]
    
    return pfm_array


