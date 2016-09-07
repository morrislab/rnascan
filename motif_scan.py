#!/usr/bin/env python
#
# Calculates motif scores for all PWMs in a given set of sequences in FASTA
# format
#
# Requires python 2.7+
#
# This tool was motivated by the RNA RBP motif scanning tool from CISBP-RNA:
# http://cisbp-rna.ccbr.utoronto.ca/TFTools.php

import sys
import time
import glob
import fileinput
import os
import warnings
import argparse
import ast
from collections import defaultdict
import multiprocessing
from itertools import izip, repeat
from RNACompete import secondarystructure, matrix
import pandas as pd
from Bio import motifs, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import RNAAlphabet, IUPAC
#sys.path.append(os.path.dirname(os.path.abspath(__file__)))

__version__ = 'v0.5.0'


def getoptions():
    desc = "Scan sequence for motif binding sites. Results sent to STDOUT."
    parser = argparse.ArgumentParser(description=desc, version=__version__)
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
    parser.add_argument('-c', '--cores', type=int, default=8, dest="cores",
                        help="Number of processing cores [%(default)s]")
    #parser.add_argument('-x', '--excel', action="store_true", dest="excel",
            #default=False,
            #help="Format the RBP_ID column with =HYPERLINK(url) for " +
            #"import into Excel [%(default)s]")
    parser.add_argument('-u', '--uniformbg', action="store_true",
                        default=False, dest="uniform_background",
                        help=("Use uniform background for calculating "
                              "log-odds [%(default)s]. Default is to compute "
                              "background from input sequences. This option "
                              "is mutually exclusive with -B."))
    parser.add_argument('-B', '--bgonly', action="store_true", default=False,
                        dest="bgonly",
                        help=("Compute background probabilities from input "
                              "sequences (STDOUT) and exit. Useful for "
                              "getting background probabilities from a "
                              "superset of sequences. Then, these values can "
                              "be subsequently supplied using -b. "
                              "[%(default)s]"))
    parser.add_argument('-b', '--bg', default=None, dest="custom_background",
                        help=("Load file of pre-computed background "
                              "probabilities"))
    parser.add_argument('-x', '--debug', action="store_true", default=False,
                        dest="debug",
                        help=("Turn on debug mode  "
                              "(aka disable parallelization) [%(default)s]"))
    args = parser.parse_args()

    if args.uniform_background is True and args.custom_background is not None:
        parser.error("You cannot set uniform and custom background options "
                     "at the same time\n")

    if args.testseq is None and args.fastafile is None:
        parser.error("Missing input FASTA file or test sequence\n")

    if args.seqtype == 'SS':
        args.alphabet = secondarystructure.RNAContextualSecondaryStructure()
    elif args.seqtype == 'RNA':
        args.alphabet = IUPAC.IUPACUnambiguousRNA()
    else:
        args.alphabet = IUPAC.IUPACUnambiguousDNA()

    return args


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


def load_motifs(dbdir, *args):
    """
    Load all motifs from given directory. Will look for *.txt files
    """
    motifs_set = {}
    print >> sys.stderr, "Loading motifs ",
    tic = time.time()
    for mfile in glob.glob(dbdir + "/*.txt"):
        try:
            motif_id = os.path.splitext(os.path.basename(mfile))[0]
            motifs_set[motif_id] = pwm2pssm(mfile, *args)
        except ValueError:
            print >> sys.stderr, "\nFailed to load motif %s" % mfile
        except KeyError:
            print >> sys.stderr, "\nFailed to load motif %s" % mfile
            print >> sys.stderr, "Check that you are using the correct --type"
            raise
        except:
            print "Unexpected error:", sys.exc_info()[0]
            raise
        print >> sys.stderr, "\b.",
        sys.stderr.flush()
    toc = time.time()
    print >> sys.stderr, "done in %0.2f seconds!" % (float(toc - tic))
    print >> sys.stderr, "Found %d motifs" % len(motifs_set)
    return motifs_set


def pwm2pssm(pwm_file, pseudocount, alphabet, background=None):
    """
    Convert load PWM and convert it to PSSM (take the log_odds)
    """
    pwm = pd.read_table(pwm_file)
    pwm = pwm.drop(pwm.columns[0], 1).to_dict(orient='list')
    pwm = motifs.Motif(alphabet=alphabet, counts=pwm)
    pwm = pwm.counts.normalize(pseudocount)

    # Can optionally add background, but for now assuming uniform probability
    pssm = pwm.log_odds(background=background)

    pssm = matrix.ExtendedPositionSpecificScoringMatrix(pssm.alphabet, pssm)

    return pssm


def preprocess_seq(seqrec, alphabet):
    """Pre-process the SeqRecord by setting the alphabet and performing
    transcription if necessary.

    Return Seq object
    """
    if not isinstance(seqrec, SeqRecord):
        raise TypeError("SeqRecord object must be supplied")

    if isinstance(alphabet, IUPAC.IUPACAmbiguousRNA) and \
        not isinstance(seqrec.seq.alphabet, RNAAlphabet):
        # If RNA alphabet is specified and input sequences are in DNA, we need
        # to transcribe them to RNA
        try:
            seq = seqrec.seq.transcribe()
            seq.alphabet = alphabet
            seq = seq.upper()
        except:
            raise
    else:
        seq = seqrec.seq

    ## If strand is specified, reverse-complement the sequence
    #strand_match = re.search(r'strand=([+-])', seqrec.description)
    #if strand_match and strand_match.group(1) == "-":
        #seq = seq.reverse_complement()

    return seq


def collect(motif_hits):
    """ Finalize results into a DataFrame for output
    """

    # Create DataFrame from motif hits
    hits = pd.DataFrame(motif_hits, columns=['Motif_ID', 'Start', 'End', 'Sequence',
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

    seq = preprocess_seq(seqrecord, opts.alphabet)
    for motif_id, pssm in pssms.iteritems():
        results = scan(pssm, seq, opts.minscore, motif_id)
        hits.extend(results)

    # Collect results
    final = collect(hits)
    return final


def _scan_all_star(a_b):
    return scan_all(*a_b)


def compute_background(fasta, alphabet, cores=8):
    """Compute background probabiilities from all input sequences
    """
    print >> sys.stderr, "Calculating background probabilities..."
    content = defaultdict(int)
    total = 0
    fin = fileinput.input(fasta,
                          openhook=fileinput.hook_compressed)
    for seqrecord in SeqIO.parse(fin, "fasta"):
        seqobj = preprocess_seq(seqrecord, alphabet)
        for letter in alphabet.letters:
            content[letter] += seqobj.count(letter)
            total += seqobj.count(letter)
    pct_sum = 0

    for letter, count in content.iteritems():
        content[letter] = float(count) / total
        if content[letter] < 0.0001:
            warnings.warn("Letter %s has very low content: %0.2f" % (letter, content[letter]), Warning)
        pct_sum += content[letter]

    for letter, value in content.iteritems():
        print >> sys.stderr, "%s: %f" % (letter, value),

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
        if args.custom_background is not None:
            print >> sys.stderr, ("Reading custom background probabilities "                  "from %s" % args.custom_background)
            # load custom background
            # http://stackoverflow.com/a/11027069
            with open(args.custom_background, 'r') as fin:
                bg = fin.read()
                bg = ast.literal_eval(bg)
                print >> sys.stderr, dict(bg)
        else:
            bg = compute_background(args.fastafile, args.alphabet, args.cores)

            if args.bgonly:
                # Print background probabilities and quit
                print >> sys.stderr, "Saving background probabilties in %s"
                print dict(bg)
                sys.exit()


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
        fin = fileinput.input(args.fastafile,
                              openhook=fileinput.hook_compressed)
        seq_iter = SeqIO.parse(fin, "fasta")

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
                batch_results = p.map(_scan_all_star, izip(batch,
                                                           repeat(pssms),
                                                           repeat(args)
                                                          )
                                     )

                # Process each result
                for j, hits in enumerate(batch_results):
                    if hits is None:
                        continue
                    hits['Sequence_ID'] = batch[j].id
                    hits['Description'] = batch[j].description
                    count += 1
                    results.append(hits)
            p.close()

        fin.close()
        final = pd.concat(results)

    cols = final.columns.tolist()
    cols = cols[-2:] + cols[:-2]
    final[cols].to_csv(sys.stdout, sep="\t", index=False)
    toc = time.time()

    print >> sys.stderr, "done in %0.2f seconds!" % (float(toc - tic))
    print >> sys.stderr, "Processed %d sequences" % count

if __name__ == '__main__':
    main()
