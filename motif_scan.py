#!/usr/bin/env python
#
# Calculates motif scores for all PWMs in a given set of sequences in FASTA
# format
#
# This tool was motivated by the RNA RBP motif scanning tool from CISBP-RNA:
# http://cisbp-rna.ccbr.utoronto.ca/TFTools.php

from RNACompete import secondarystructure, matrix
import sys
import time
import glob
import os
from optparse import OptionParser
from collections import defaultdict
import multiprocessing
import pandas as pd
from Bio import motifs, SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

__version__ = 'v0.2.0'


def getoptions():
    usage = "usage: python %prog [options] sequences.fa"
    desc = "Scan sequence for motif binding sites."
    parser = OptionParser(usage=usage, description=desc, version=__version__)
    parser.add_option('-d', type='string', dest="pwm_dir",
                      default=os.path.dirname(os.path.abspath(__file__)) +
                      "/db/pwms",
                      help="Directory of PWMs [%default]")
    parser.add_option('-p', '--pseudocount', type="float", dest="pseudocount",
                      default=0,
                      help="Pseudocount for normalizing PWM. [%default]")
    #parser.add_option('-r', '--rbpinfo', type='string', dest='rbpinfo',
        #default=os.path.dirname(os.path.abspath(__file__)) +
        #"/db/RBP_Information_all_motifs.txt",
        #help="RBP info for adding meta data to results. [%default]")
    parser.add_option('-t', '--type', type='choice', dest='seqtype',
                      choices=['DNA', 'RNA', 'SS'],
                      default="RNA",
                      help=("Alphabet of PWM (DNA|RNA|SS for "
                            "RNAContextualSecondaryStructure). [%default]"))
    parser.add_option('-m', '--minscore', type='float', dest='minscore',
                      default=6,
                      help="Minimum score for motif hits. [%default]")
    parser.add_option('-s', '--seq', type='string', dest='testseq',
                      default=None,
                      help=("Supply a test sequence to scan. FASTA files will"
                            " be ignored."))
    parser.add_option('-c', type="int", default=8,
                      dest="cores", metavar="CORES",
                      help="Number of processing cores [%default]")
    #parser.add_option('-x', '--excel', action="store_true", dest="excel",
            #default=False,
            #help="Format the RBP_ID column with =HYPERLINK(url) for " +
            #"import into Excel [%default]")
    parser.add_option('-b', action="store_true", default=False,
                      dest="use_background",
                      help=("Perform background calculation from input "
                            "sequences [%default]"))
    parser.add_option('-x', '--debug', action="store_true", default=False,
                      dest="debug",
                      help=("Turn on debug mode (aka disable parallelization) "
                            "[%default]"))
    (opts, args) = parser.parse_args()

    if opts.testseq is None and len(args) < 1:
        parser.error("Missing input FASTA file\n")

    if opts.seqtype == 'SS':
        opts.alphabet = secondarystructure.RNAContextualSecondaryStructure()
    elif opts.seqtype == 'RNA':
        opts.alphabet = IUPAC.IUPACUnambiguousRNA()
    else:
        opts.alphabet = IUPAC.IUPACUnambiguousDNA()

    return (opts, args)


def load_motifs(db, *args):
    """
    Load all motifs from given directory. Will look for *.txt files
    """
    motifs = {}
    print >> sys.stderr, "Loading motifs ",
    tic = time.time()
    for file in glob.glob(db + "/*.txt"):
        try:
            id = os.path.splitext(os.path.basename(file))[0]
            motifs[id] = pwm2pssm(file, *args)
        except:
            print >> sys.stderr, "\nFailed to load motif %s" % file
            print >> sys.stderr, "Check that you are using the correct --type"
            raise
        print >> sys.stderr, "\b.",
        sys.stderr.flush()
    toc = time.time()
    print >> sys.stderr, "done in %0.2f seconds!" % (float(toc - tic))
    print >> sys.stderr, "Found %d motifs" % len(motifs)
    return motifs


def pwm2pssm(file, pseudocount, alphabet, bg=None):
    """
    Convert load PWM and convert it to PSSM (take the log_odds)
    """
    pwm = pd.read_table(file)
    pwm = pwm.drop(pwm.columns[0], 1).to_dict(orient='list')
    pwm = motifs.Motif(alphabet=alphabet, counts=pwm)
    pwm = pwm.counts.normalize(pseudocount)

    # Can optionally add background, but for now assuming uniform probability
    pssm = pwm.log_odds(background=bg)

    #if isinstance(alphabet, secondarystructure.RNAContextualSecondaryStructure):
    pssm = matrix.ExtendedPositionSpecificScoringMatrix(pssm.alphabet, pssm)

    return(pssm)


def collect(x):
    """
    Finalize results into a DataFrame for output
    """

    # Create DataFrame from motif hits
    hits = pd.DataFrame(x, columns=['Motif_ID', 'Start', 'End', 'Sequence',
                                    'LogOdds'])

    # Merge metadata with hits
    return hits.sort_values(['Start', 'Motif_ID'])


def scan(pssm, seq, minscore, motif_id):
    results = []
    for position, score in pssm.search(seq, threshold=minscore, both=False):
        end_position = position + len(pssm.consensus)

        fragment = seq[position:end_position]
        #if isinstance(seq.alphabet, IUPAC.IUPACUnambiguousRNA):
            #fragment = fragment.transcribe()

        values = [motif_id,
                  position + 1, end_position,
                  str(fragment),
                  round(score, 3)]
        results.append(values)
    return results


def scan_all(pssms, seq, opts, bg=None):
    """
    Scan seq for all motifs in pssms
    """
    hits = []
    tasks = []

    if opts.debug:
        for motif_id, pssm in pssms.iteritems():
            results = scan(pssm, seq, opts.minscore, motif_id)
            hits.extend(results)
    else:
        p = multiprocessing.Pool(opts.cores)
        for motif_id, pssm in pssms.iteritems():
            tasks.append((pssm, seq, opts.minscore, motif_id,))
        results = [p.apply_async(scan, t) for t in tasks]

        for r in results:
            hits.extend(r.get())
        p.close()

    # Collect results
    final = collect(hits)
    return final


def _set_seq(seq, alphabet):
    if not isinstance(seq, Seq):
        raise TypeError("Seq object must be supplied")

    if isinstance(alphabet, IUPAC.IUPACUnambiguousRNA):
        # If RNA alphabet is specified and input sequences are in DNA, we need
        # to transcribe them to RNA
        try:
            seq = seq.transcribe()
            seq.alphabet = alphabet
        except:
            print >> sys.stderr, "Unexpected error:", sys.exc_info()[0]
            raise
    return seq


def compute_background(fasta, alphabet):
    """Compute background probabiilities from all input sequences
    """
    print >> sys.stderr, "Calculating background"
    content = defaultdict(int)
    total = 0
    for seqrecord in SeqIO.parse(open(fasta), "fasta"):
        seqobj = _set_seq(seqrecord.seq, alphabet)
        for letter in alphabet.letters:
            content[letter] += seqobj.count(letter)
        total += len(seqobj)
    for letter, count in content.iteritems():
        content[letter] = float(count) / total
    return content


def main():
    (opts, args) = getoptions()
    tic = time.time()
    final = pd.DataFrame()
    count = 0

    # Calculate sequence background from input
    bg = None
    if opts.use_background:
        bg = compute_background(args[0], opts.alphabet)

    # Load PWMs
    pssms = load_motifs(opts.pwm_dir, opts.pseudocount, opts.alphabet, bg)

    if opts.testseq is not None:
        seq = _set_seq(Seq(opts.testseq), opts.alphabet)
        final = scan_all(pssms, seq, opts)
        final['Sequence_ID'] = 'testseq'
        count += 1
    else:
        print >> sys.stderr, "Scanning sequences ",

        results = []
        for seqrecord in SeqIO.parse(open(args[0]), "fasta"):
            seq = _set_seq(seqrecord.seq, opts.alphabet)

            m = scan_all(pssms, seq, opts)
            m['Sequence_ID'] = seqrecord.id

            results.append(m)
            count += 1
        final = pd.concat(results)

    cols = final.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    print final[cols].to_csv(sep="\t", index=False),
    toc = time.time()

    print >> sys.stderr, "done in %0.2f seconds!" % (float(toc - tic))
    print >> sys.stderr, "Processed %d sequences" % count

if __name__ == '__main__':
    main()
