#!/usr/bin/env python
# Copyright 2016-2017 by Kevin Ha
"""
Calculates motif scores for all PFMs in a given set of sequences in FASTA
format

Requires python 2.7+

This tool was motivated by the RNA RBP motif scanning tool from CISBP-RNA:
http://cisbp-rna.ccbr.utoronto.ca/TFTools.php
"""

import sys
import time
import fileinput
import os
import os.path
import warnings
import argparse
import ast
from collections import defaultdict
import multiprocessing
import pandas as pd
from itertools import izip, repeat
from BioAddons.Alphabet import ContextualSecondaryStructure
from BioAddons.motifs import matrix
from Bio import motifs, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import RNAAlphabet, IUPAC

__version__ = 'v0.8.0'


def getoptions():
    desc = "Scan sequence for motif binding sites. Results sent to STDOUT."
    parser = argparse.ArgumentParser(description=desc, version=__version__)
    parser.add_argument('fastafiles', metavar='FASTA', nargs='*',
                        help="Input sequence and structure FASTA files")
    pfm_grp = parser.add_argument_group("PFM options")
    pfm_grp.add_argument('-p', '--pfm_seq', dest="pfm_seq", type=str,
                         help="Sequence PFM")
    pfm_grp.add_argument('-q', '--pfm_struct', dest="pfm_struct", type=str,
                         help="Structure PFM")

    parser.add_argument('-C', '--pseudocount', type=float,
                        dest="pseudocount", default=0,
                        help="Pseudocount for normalizing PFM. [%(default)s]")
    parser.add_argument('-m', '--minscore', type=float, dest='minscore',
                        default=6,
                        help="Minimum score for motif hits. [%(default)s]")
    parser.add_argument('-t', '--testseq', dest='testseq', default=None,
                        help=("Supply a test sequence to scan. FASTA files "
                              "will be ignored. Can supply sequence and "
                              "structure as single string separated by "
                              " comma."))
    parser.add_argument('-c', '--cores', type=int, default=8, dest="cores",
                        help="Number of processing cores [%(default)s]")

    bg_grp = parser.add_argument_group('Background frequency options')
    bg_grp.add_argument('-u', '--uniformbg', action="store_true",
                        default=False, dest="uniform_background",
                        help=("Use uniform background for calculating "
                              "log-odds [%(default)s]. Default is to compute "
                              "background from input sequences. This option "
                              "is mutually exclusive with -B."))
    bg_grp.add_argument('-g', '--bgonly', action="store_true", default=False,
                        dest="bgonly",
                        help=("Compute background probabilities from input "
                              "sequences (STDOUT) and exit. Useful for "
                              "getting background probabilities from a "
                              "superset of sequences. Then, these values can "
                              "be subsequently supplied using -b. "
                              "[%(default)s]"))
    bg_grp.add_argument('-b', '--bg_seq', default=None, dest="bg_seq",
                        help=("Load file of pre-computed background "
                              "probabilities for nucleotide sequences"))
    bg_grp.add_argument('-B', '--bg_struct', default=None, dest="bg_struct",
                        help=("Load file of pre-computed background "
                          "probabilities for nucleotide sequences"))
    parser.add_argument('-x', '--debug', action="store_true", default=False,
                        dest="debug",
                        help=("Turn on debug mode  "
                              "(aka disable parallelization) [%(default)s]"))
    args = parser.parse_args()

    if not (args.pfm_seq or args.pfm_struct or args.bgonly):
        parser.error("Must specify PFMs with -p and/or -q")

    if args.uniform_background and (args.bg_seq or args.bg_struct):
        parser.error("You cannot set uniform and custom background options "
                     "at the same time\n")

    nfiles = len(args.fastafiles)

    if nfiles == 2:
        if not (args.pfm_seq or args.pfm_struct):
            parser.error("Missing PFMs")
        args.seq_type == "RNASS"
    else:   # nfiles == 1
        if args.pfm_seq and args.pfm_struct and not args.testseq:
            parser.error("Can't specify two PFMs with one input file")
        elif args.pfm_seq and args.pfm_struct and args.testseq:
            args.seq_type = "RNASS"
        elif args.pfm_seq:
            args.seq_type = "RNA"
        elif args.pfm_struct:
            args.seq_type = "SS"
        else:
            parser.error("Must specify PFMs with -p and/or -q")

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


def load_motif(pfm_file, *args):
    """ Load PFM
    """
    motifs_set = {}
    print >> sys.stderr, "Loading PFM %s" % pfm_file,
    tic = time.time()
    try:
        motif_id = os.path.splitext(os.path.basename(pfm_file))[0]
        motifs_set[motif_id] = pfm2pssm(pfm_file, *args)
    except ValueError:
        print >> sys.stderr, "\nFailed to load motif %s" % pfm_file
    except KeyError:
        print >> sys.stderr, "\nFailed to load motif %s" % pfm_file
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


def pfm2pssm(pfm_file, pseudocount, alphabet, background=None):
    """
    Convert load PFM and convert it to PSSM (take the log_odds)
    """
    pfm = pd.read_table(pfm_file)
    pfm = pfm.drop(pfm.columns[0], 1).to_dict(orient='list')
    pfm = motifs.Motif(alphabet=alphabet, counts=pfm)
    pfm = pfm.counts.normalize(pseudocount)


    pssm = pfm.log_odds(background=background)

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
    columns=['Motif_ID', 'Start', 'End', 'Sequence', 'LogOdds']

    # Create DataFrame from motif hits
    hits = pd.DataFrame(motif_hits, columns=columns)

    # Merge metadata with hits
    return hits.sort_values(['Start', 'Motif_ID'])


def scan(pssm, seq, alphabet, minscore):
    results = []
    (motif_id, pm) = pssm.items()[0]
    for position, score in pm.search(seq, threshold=minscore, both=False):
        end_position = position + len(pm.consensus)

        fragment = seq[position:end_position]
        #if isinstance(seq.alphabet, IUPAC.IUPACUnambiguousRNA):
            #fragment = fragment.transcribe()

        values = [motif_id,
                  position + 1, end_position,
                  str(fragment),
                  round(score, 3)]
        results.append(values)
    return results


def scan_all(seqrecord, pssm, alphabet, *args):
    """ Scan seq for all motifs in pssms
    """

    seq = preprocess_seq(seqrecord, alphabet)
    results = scan(pssm, seq, alphabet, *args)

    columns=['Motif_ID', 'Start', 'End', 'Sequence', 'LogOdds']
    final = pd.DataFrame(results, columns=columns)
    return final.sort_values(['Start', 'Motif_ID'])


def _scan_all_star(a_b):
    return scan_all(*a_b)


def parse_sequences(fasta_file):
    """Load FASTA sequence and return SeqRecord iterator
    """
    fin = fileinput.input(fasta_file, openhook=fileinput.hook_compressed)
    return SeqIO.parse(fin, 'fasta')


def compute_background(fastas, alphabet, verbose=True):
    """ Compute background probabiilities from all input sequences
    """
    print >> sys.stderr, "Calculating background probabilities..."
    content = defaultdict(int)
    total = len(alphabet.letters)       # add psuedocount for each letter
    seq_iter = parse_sequences(fastas)

    for seqrecord in seq_iter:
        seqobj = preprocess_seq(seqrecord, alphabet)
        for letter in alphabet.letters:
            amount = seqobj.count(letter)
            content[letter] += amount
            total += amount
    pct_sum = 0

    for letter, count in content.iteritems():
        content[letter] = (float(count) + 1) / total    # add pseudocount
        if content[letter] <= 0.05:
            warnings.warn("Letter %s has low content: %0.2f"
                          % (letter, content[letter]), Warning)
        pct_sum += content[letter]

    if verbose: print >> sys.stderr, dict(content)
    assert abs(1.0 - pct_sum) < 0.0001, "Background sums to %f" % pct_sum
    return content


def load_background(bg_file, uniform, *args):
    """ Load background probabilities if available, otherwise compute from
    input files or use uniform
    """
    if bg_file:
        print >> sys.stderr, ("Reading custom background probabilities "
            "from %s" % file)
        # load custom background
        # http://stackoverflow.com/a/11027069
        with open(file, 'r') as fin:
            bg = fin.read()
            bg = ast.literal_eval(bg)
            print >> sys.stderr, dict(bg)
    elif not uniform:
        bg = compute_background(args)
    else:
        bg = None
    return bg


def scan_main(fasta_file, pssm, alphabet, bg, args):

    final = pd.DataFrame()
    count = 0

    if isinstance(fasta_file, SeqRecord):
        final = scan_all(fasta_file, pssm, alphabet, args.minscore)
        final['Sequence_ID'] = 'testseq'
        final['Description'] = ''
        count += 1
    else:
        print >> sys.stderr, "Scanning sequences "

        results = []
        seq_iter = parse_sequences(fasta_file)

        if args.debug:
            for seqrecord in seq_iter:
                hits = scan_all(seqrecord, pssm, args)
                hits['Sequence_ID'] = seqrecord.id
                hits['Description'] = seqrecord.description
                results.append(hits)
                count += 1
        else:
            p = multiprocessing.Pool(args.cores)
            for i, batch in enumerate(batch_iterator(seq_iter, 500)):
                batch_results = p.map(_scan_all_star, izip(batch,
                                                           repeat(pssm),
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

        if len(results) != 0:
            final = pd.concat(results)

    print >> sys.stderr, "Processed %d sequences" % count
    cols = final.columns.tolist()
    cols = cols[-2:] + cols[:-2]
    return final[cols]


def combine(seq_results, struct_results):
    """ If scoring sequence and structure together, add up the log odds at
    every position
    """

    # Keys to match by : Sequence_ID, Start, End
    # Motif_ID may not match
    result = pd.merge(seq_results, struct_results,
                      on=['Sequence_ID', 'Start', 'End'])
    result.rename(columns={'Description_x': 'Description.Seq',
                           'Description_y': 'Description.Struct',
                           'Sequence_x': 'Sequence.Seq',
                           'Sequence_y': 'Sequence.Struct',
                           'Motif_ID_x': 'Motif_ID.Seq',
                           'Motif_ID_y': 'Motif_ID.Struct',
                           'LogOdds_x': 'LogOdds.Seq',
                           'LogOdds_y': 'LogOdds.Struct'}, inplace=True)
    result['LogOdds.SeqStruct'] = result['LogOdds.Seq'] + \
        result['LogOdds.Struct']
    return result


def main():
    tic = time.time()
    args = getoptions()
    bg = None

    if args.testseq:
        args.testseq = args.testseq.split(',')

    ## Sequence
    if args.seq_type in ['RNA', 'RNASS']:
        if args.testseq:
            seq_file = SeqRecord(Seq(args.testseq[0]))
        else:
            seq_file = args.fastafiles[0]

        if args.bg_seq and not args.testseq:
            bg = load_background(args.bg_seq,
                                 IUPAC.IUPACUnambiguousRNA(),
                                 args.fastafiles[0],
                                 args.uniform_background)

        if not args.bgonly:
            pssm = load_motif(args.pfm_seq,
                              args.pseudocount,
                              IUPAC.IUPACUnambiguousRNA(),
                              bg)
            seq_results = scan_main(seq_file,
                                    pssm,
                                    IUPAC.IUPACUnambiguousRNA(),
                                    bg, args)

    ## Structure
    if args.seq_type in ['SS', 'RNASS']:
        if args.testseq:
            struct_file = SeqRecord(Seq(args.testseq[1]))
        elif args.seq_type == 'SS':
            struct_file = args.fastafiles[0]
        else:
            struct_file = args.fastafiles[1]

        if args.bg_struct and not args.testseq:
            bg = load_background(args.bg_struct,
                                 ContextualSecondaryStructure(),
                                 struct_file,
                                 args.uniform_background)

        if not args.bgonly:
            pssm = load_motif(args.pfm_struct,
                              args.pseudocount,
                              ContextualSecondaryStructure(),
                              bg)
            struct_results = scan_main(struct_file,
                                       pssm,
                                       ContextualSecondaryStructure(),
                                       bg, args)


    if args.seq_type == 'RNASS':
        combined_results = combine(seq_results, struct_results)
        combined_results.to_csv(sys.stdout, sep="\t", index=False)
    elif args.seq_type == 'RNA':
        seq_results.to_csv(sys.stdout, sep="\t", index=False)
    else:
        struct_results.to_csv(sys.stdout, sep="\t", index=False)

    toc = time.time()

    runtime = float(toc - tic)
    if runtime > 60:
        print >> sys.stderr, "Done in %0.4f minutes!" % (runtime / 60)
    else:
        print >> sys.stderr, "Done in %0.4f seconds!" % (runtime)


if __name__ == '__main__':
    main()
