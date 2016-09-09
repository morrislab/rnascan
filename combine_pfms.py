#!/usr/bin/env python
"""
Script for combining probabilities of RNA sequence and secondary structure PFMs

Outputs another n x m PFM, where n is the number of positions and
m = m_i * m_j where m_i is the number of letters in RNA sequence PFM and m_j is
the number of letters in secondary struture PFMs.
"""
from BioAddons.Alphabet import *
import sys
import pandas as pd
import numpy as np
import argparse


def getoptions():
    desc = "Multiply two PFMs element-wise"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('PFMfiles', metavar='PFMs', nargs=2,
                        help="Two PFM files")
    return parser.parse_args()

def load_pfm(pfmfile):
    """Load the PFM
    """
    df = pd.read_table(pfmfile)
    df = df.drop(df.columns[0], 1)
    return df


def combine_pfms(pfm1, pfm2):
    """Given two PFMs (e.g. sequence and secondary structure), combine their
    probabilities by taking the outer product of each element and position
    """
    assert pfm1.shape[0] == pfm2.shape[0]

    results = []
    # Go through each position
    for i in xrange(0, pfm1.shape[0]):
        p = np.outer(pfm1.iloc[i], pfm2.iloc[i])
        results.append(np.ravel(p))

    alphabet = [ContextualSequenceSecondaryStructure.convert(i,j) \
                for i in pfm1.columns for j in pfm2.columns]
    return pd.DataFrame(results, columns=alphabet)


def main():
    args = getoptions()

    # Load first PFM
    seqpfm = load_pfm(args.PFMfiles[0])

    # Load second PFM
    structpfm = load_pfm(args.PFMfiles[1])

    combinedpfm = combine_pfms(seqpfm, structpfm)

    # Print new PFM to STDOUT
    combinedpfm.to_csv(sys.stdout, sep="\t", index_label='PO')


if __name__ == '__main__':
    main()
