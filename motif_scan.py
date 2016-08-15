#!/usr/bin/env python
#
# Calculates motif scores for all PWMs in a given set of sequences in FASTA
# format
#
# Requires python 2.7+
#
# This tool was motivated by the RNA RBP motif scanning tool from CISBP-RNA:
# http://cisbp-rna.ccbr.utoronto.ca/TFTools.php

from RNACompete import secondarystructure, matrix
import sys
import time
import glob
import os
import re
import argparse
from collections import defaultdict
import multiprocessing
from itertools import izip_longest, izip, repeat
import pandas as pd
from Bio import motifs, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
#sys.path.append(os.path.dirname(os.path.abspath(__file__)))

__version__ = 'v0.4.0'


def getoptions():
    desc = "Scan sequence for motif binding sites."
    parser = argparse.ArgumentParser(description=desc,
                                     version=__version__)
    parser.add_argument('fastafile', metavar='FASTA', nargs='?',
                        help="Input FASTA file")
    parser.add_argument('-d', dest="pwm_dir",
                        default=os.path.dirname(os.path.abspath(__file__)) +
                        "/db/pwms",
                        help="Directory of PWMs [%(default)s]")
    parser.add_argument('-p', '--pseudocount', type=float,
                        dest="pseudocount", default=0,
                        help="Pseudocount for normalizing PWM. [%(default)s]")
    #parser.add_argument('-r', '--rbpinfo', type='string', dest='rbpinfo',
        #default=os.path.dirname(os.path.abspath(__file__)) +
        #"/db/RBP_Information_all_motifs.txt",
        #help="RBP info for adding meta data to results. [%(default)s]")
    parser.add_argument('-t', '--type', dest='seqtype',
                        choices=['DNA', 'RNA', 'SS'], default="RNA",
                        help=("Alphabet of PWM (DNA|RNA|SS for "
                              "RNAContextualSecondaryStructure). "
                              "[%(default)s]"))
    parser.add_argument('-m', '--minscore', type=float, dest='minscore',
                        default=6,
                        help="Minimum score for motif hits. [%(default)s]")
    parser.add_argument('-s', '--seq', dest='testseq', default=None,
                        help=("Supply a test sequence to scan. FASTA files "
                              " will be ignored."))
    parser.add_argument('-c', type=int, default=8,
                        dest="cores", metavar="CORES",
                        help="Number of processing cores [%(default)s]")
    #parser.add_argument('-x', '--excel', action="store_true", dest="excel",
            #default=False,
            #help="Format the RBP_ID column with =HYPERLINK(url) for " +
            #"import into Excel [%(default)s]")
    parser.add_argument('-u', action="store_true", default=False,
                        dest="uniform_background",
                        help=("Use uniform background for calculating "
                              "log-odds [%(default)s]. Default is to compute "
                              "background from input sequences"))
    parser.add_argument('-x', '--debug', action="store_true", default=False,
                        dest="debug",
                        help=("Turn on debug mode  "
                              "(aka disable parallelization) [%(default)s]"))
    args = parser.parse_args()

    if args.testseq is None and args.fastafile is None:
        parser.error("Missing input FASTA file or test sequence\n")

    if args.seqtype == 'SS':
        args.alphabet = secondarystructure.RNAContextualSecondaryStructure()
    elif args.seqtype == 'RNA':
        args.alphabet = IUPAC.IUPACUnambiguousRNA()
    else:
        args.alphabet = IUPAC.IUPACUnambiguousDNA()

    return (args)

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.

    Source: http://biopython.org/wiki/Split_large_file
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


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
        except ValueError:
            print >> sys.stderr, "\nFailed to load motif %s" % file
        except KeyError:
            print >> sys.stderr, "\nFailed to load motif %s" % file
            print >> sys.stderr, "Check that you are using the correct --type"
            raise
        except:
            print "Unexpected error:", sys.exc_info()[0]
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


def scan_all(seqrecord, *args):
    """
    Scan seq for all motifs in pssms
    """
    pssms = args[0]
    opts = args[1]
    hits = []
    tasks = []

    seq = _set_seq(seqrecord, opts.alphabet)
    for motif_id, pssm in pssms.iteritems():
        results = scan(pssm, seq, opts.minscore, motif_id)
        hits.extend(results)

    # Collect results
    final = collect(hits)
    return final


def _set_seq(seqrec, alphabet):
    """Pre-process the SeqRecord by setting the alphabet and performing
    transcription if necessary
    """
    if not isinstance(seqrec, SeqRecord):
        raise TypeError("SeqRecord object must be supplied")

    if isinstance(alphabet, IUPAC.IUPACAmbiguousRNA):
        # If RNA alphabet is specified and input sequences are in DNA, we need
        # to transcribe them to RNA
        try:
            seq = seqrec.seq.transcribe()
            seq.alphabet = alphabet
            seq = seq.upper()
        except:
            raise
    
    ## If strand is specified, reverse-complement the sequence
    #strand_match = re.search(r'strand=([+-])', seqrec.description)
    #if strand_match and strand_match.group(1) == "-":
        #seq = seq.reverse_complement()

    return seq

def _scan_all_star(a_b):
    return scan_all(*a_b)

def compute_background(fasta, alphabet, cores=8):
    """Compute background probabiilities from all input sequences
    """
    print >> sys.stderr, "Calculating background probabilities...",
    content = defaultdict(int)
    total = 0
    for seqrecord in SeqIO.parse(open(fasta), "fasta"):
        seqobj = _set_seq(seqrecord, alphabet)
        for letter in alphabet.letters:
            content[letter] += seqobj.count(letter)
            total += seqobj.count(letter)
    pct_sum = 0

    for letter, count in content.iteritems():
        content[letter] = float(count) / total
        pct_sum += content[letter]
        print >> sys.stderr, "%s: %f" % (letter, content[letter]),

    print >> sys.stderr, ""
    assert (1.0 - pct_sum) < 0.0001, "Background sums to %f" % pct_sum
    return content


def main():
    args = getoptions()
    tic = time.time()
    final = pd.DataFrame()
    count = 0

    # Calculate sequence background from input
    if args.uniform_background or args.testseq is not None:
        print >> sys.stderr, "Using uniform background probabilities"
        bg = None
    else:
        bg = compute_background(args.fastafile, args.alphabet, args.cores)

    # Load PWMs
    pssms = load_motifs(args.pwm_dir, args.pseudocount, args.alphabet, bg)

    if args.testseq is not None:
        #seq = _set_seq(SeqRecord(Seq(args.testseq)), args.alphabet)
        final = scan_all(SeqRecord(Seq(args.testseq)), pssms, args)
        final['Sequence_ID'] = 'testseq'
        final['Description'] = ''
        count += 1
    else:
        print >> sys.stderr, "Scanning sequences ",

        results = []
        seq_iter = SeqIO.parse(open(args.fastafile), "fasta")

        if args.debug:
            for seqrecord in seq_iter:
                hits = scan_all(seqrecord, pssms, args)
                hits['Sequence_ID'] = seqrecord.id
                hits['Description'] = seqrecord.description
                results.append(hits)
                count += 1
        else:
            p = multiprocessing.Pool(args.cores)
            for i, batch in enumerate(batch_iterator(seq_iter, 500)):
                batch_results = p.map(_scan_all_star, 
                                izip(batch, 
                                     repeat(pssms),
                                     repeat(args)
                                    )
                                )

                # Process each result
                for j, hits in enumerate(batch_results):
                    if hits is None: continue
                    hits['Sequence_ID'] = batch[j].id
                    hits['Description'] = batch[j].description
                    count += 1
                    results.append(hits)
            p.close()

        final = pd.concat(results)

    cols = final.columns.tolist()
    cols = cols[-2:] + cols[:-2]
    final[cols].to_csv(sys.stdout, sep="\t", index=False)
    toc = time.time()

    print >> sys.stderr, "done in %0.2f seconds!" % (float(toc - tic))
    print >> sys.stderr, "Processed %d sequences" % count

if __name__ == '__main__':
    main()
