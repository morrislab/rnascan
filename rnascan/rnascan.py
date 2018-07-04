# Copyright (C) 2014-2015 Kate Cook, 2016-2017 Kevin Ha
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



import sys
import time
import glob
import fileinput
import os
import os.path
import warnings
import argparse
import ast
import re
from collections import defaultdict
import multiprocessing
import pandas as pd
import numpy as np
from itertools import repeat
from .BioAddons.Alphabet import ContextualSecondaryStructure
from .BioAddons.motifs import matrix
from Bio import motifs, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import RNAAlphabet, IUPAC

__version__ = 'v0.10.0'


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

    if not (args.pfm_seq or args.pfm_struct):
        parser.error("Must specify PFMs with -p and/or -q")

    if args.uniform_background and (args.bg_seq or args.bg_struct):
        parser.error("You cannot set uniform and custom background options "
                     "at the same time\n")

    return args


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

###############################################################################
# Sequence functions
###############################################################################
def _guess_seq_type(args):
    """Given arguments, determine the sequence analysis type: RNA, SS, or RNASS
    """
    nfiles = len(args.fastafiles)

    if nfiles == 2:
        if not (args.pfm_seq or args.pfm_struct):
            eprint("Missing PFMs")
            sys.exit(1)
        seq_type = "RNASS"
    else:   # nfiles == 1
        if args.pfm_seq and args.pfm_struct and not args.testseq:
            eprint("Can't specify two PFMs with one input file")
            sys.exit(1)
        elif args.pfm_seq and args.pfm_struct and args.testseq:
            seq_type = "RNASS"
        elif args.pfm_seq:
            seq_type = "RNA"
        elif args.pfm_struct:
            seq_type = "SS"
        else:
            eprint("Must specify PFMs with -p and/or -q")
            sys.exit(1)
    return seq_type


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
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


def parse_sequences(fasta_file):
    """Load FASTA sequence and return SeqRecord iterator
    """
    fin = fileinput.input(fasta_file, openhook=fileinput.hook_compressed)
    return SeqIO.parse(fin, 'fasta')


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


###############################################################################
# PFM functions
###############################################################################
def load_motif(pfm_file, *args):
    """ Load PFM
    """
    motifs_set = {}
    eprint("Loading PFM %s" % pfm_file, end="")
    tic = time.time()
    try:
        motif_id = os.path.splitext(os.path.basename(pfm_file))[0]
        motifs_set[motif_id] = pfm2pssm(pfm_file, *args)
    except ValueError:
        eprint("\nFailed to load motif %s" % pfm_file)
    except KeyError:
        eprint("\nFailed to load motif %s" % pfm_file)
        eprint("Check that you are using the correct --type")
        raise
    except:
        eprint("Unexpected error: %s" % sys.exc_info()[0])
        raise
    eprint("\b.", end="")
    sys.stderr.flush()
    toc = time.time()
    eprint("done in %0.2f seconds!" % (float(toc - tic)))
    eprint("Found %d motifs" % len(motifs_set))
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


###############################################################################
# Motif scan functions
###############################################################################
def scan(pssm, seq, alphabet, minscore):
    """ Core scanning function
    """
    results = []
    (motif_id, pm) = list(pssm.items())[0]
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


def _scan_all(a_b):
    return scan_all(*a_b)


def scan_averaged_structure(struct_file, pssm, minscore):
    """Scan PSSM on an averaged secondary structure model
    """
    struct = pd.read_table(struct_file)
    del struct['PO']
    (motif_id, pm) = list(pssm.items())[0]
    motif_scores = []
    pm = pd.DataFrame(pm)     # Convert dict back to data frame
    N = len(pm.index)
    for i in range(0, len(struct.index) - N + 1):
        score = 0
        for j in range(0, N):
            # Multiply by SSM
            score += np.nan_to_num(np.dot(struct.iloc[i + j, :],
                                          pm.iloc[j, :]))

        # Sum log odds
        if score > minscore:
            motif_scores.append(pd.Series([motif_id, i + 1, i + N, '.', score],
                                          index=['Motif_ID', 'Start', 'End',
                                                 'Sequence',
                                                 'LogOdds']))
    return pd.DataFrame(motif_scores)


def _scan_averaged_structure(a_b):
    return scan_averaged_structure(*a_b)


def _add_sequence_id(df, seq_id, description):
    """ Helper function to add Sequence_ID and Description (df is a reference)
    """
    df['Sequence_ID'] = seq_id
    df['Description'] = description


def scan_main(fasta_file, pssm, alphabet, bg, args):
    """ Main function for handling scanning of PSSM and a sequence/structure
    """
    final = pd.DataFrame()
    count = 0

    if isinstance(fasta_file, SeqRecord):
        final = scan_all(fasta_file, pssm, alphabet, args.minscore)
        _add_sequence_id(final, 'testseq', '')
        count += 1
    else:
        results = []

        if os.path.isdir(fasta_file):
            eprint("Scanning averaged secondary structures ")

            structures = glob.glob(fasta_file + "/structure.*.txt")
            if len(structures) == 0:
                raise IOError("No averaged structure files found")

            if args.debug:
                for struct_file in structures:
                    hits = scan_averaged_structure(struct_file, pssm,
                                                 args.minscore)
                    _add_sequence_id(hits, struct_file, '')
                    results.append(hits)
                    count += 1
            else:
                p = multiprocessing.Pool(args.cores)
                batch_results = p.map(_scan_averaged_structure,
                                      zip(structures, repeat(pssm),
                                           repeat(args.minscore)))
                for j, hits in enumerate(batch_results):
                    if hits is None:
                        continue
                    match = re.search(r'^structure\.(.*)\.txt$',
                                      os.path.basename(structures[j]))
                    _add_sequence_id(hits, match.group(1), '')
                    count += 1
                    results.append(hits)
                p.close()
        else:
            eprint("Scanning sequences ")

            seq_iter = parse_sequences(fasta_file)

            if args.debug:
                for seqrecord in seq_iter:
                    hits = scan_all(seqrecord, pssm, alphabet, args.minscore)
                    _add_sequence_id(hits, seqrecord.id, seqrecord.description)
                    results.append(hits)
                    count += 1
            else:
                p = multiprocessing.Pool(args.cores)
                for i, batch in enumerate(batch_iterator(seq_iter, 2000)):
                    batch_results = p.map(_scan_all, zip(batch,
                                                          repeat(pssm),
                                                          repeat(alphabet),
                                                          repeat(args.minscore)
                                                          )
                                          )

                    # Process each result
                    for j, hits in enumerate(batch_results):
                        if hits is None:
                            continue
                        _add_sequence_id(hits, batch[j].id,
                                         batch[j].description)
                        count += 1
                        results.append(hits)
                p.close()

        if len(results) != 0:
            final = pd.concat(results)

    eprint("Processed %d sequences" % count)
    cols = final.columns.tolist()
    cols = cols[-2:] + cols[:-2]
    return final[cols]


def combine(seq_results, struct_results):
    """ If scoring sequence and structure together, add up the log odds at
    every position
    """
    # Keys to match by : Sequence_ID, Start, End
    # NB: Motif_ID may not match
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


###############################################################################
# Background functions
###############################################################################
def compute_background(fastas, alphabet, verbose=True):
    """ Compute background probabiilities from all input sequences
    """
    eprint("Calculating background probabilities...")
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

    for letter, count in content.items():
        content[letter] = (float(count) + 1) / total    # add pseudocount
        if content[letter] <= 0.05:
            warnings.warn("Letter %s has low content: %0.2f"
                          % (letter, content[letter]), Warning)
        pct_sum += content[letter]

    if verbose: eprint(dict(content))
    assert abs(1.0 - pct_sum) < 0.0001, "Background sums to %f" % pct_sum
    return content


def load_background(bg_file, uniform, *args):
    """ Load background probabilities if available, otherwise compute from
    input files or use uniform
    """
    if bg_file:
        eprint("Reading custom background probabilities from %s" % file)
        # load custom background
        # http://stackoverflow.com/a/11027069
        with open(bg_file, 'r') as fin:
            bg = fin.read()
            bg = ast.literal_eval(bg)
            eprint(dict(bg))
    elif not uniform:
        bg = compute_background(*args)
    else:
        bg = None
    return bg


###############################################################################
# Main
###############################################################################
def main():
    tic = time.time()
    args = getoptions()
    seq_type = _guess_seq_type(args)
    bg = None

    if args.testseq:
        testseq_stack = args.testseq.split(',')[::-1]    # make a stack

    ## Sequence
    if seq_type in ['RNA', 'RNASS']:
        if args.testseq:
            seq_file = SeqRecord(Seq(testseq_stack.pop()))
        else:
            seq_file = args.fastafiles[0]

        if not args.testseq:
            bg = load_background(args.bg_seq,
                                 args.uniform_background,
                                 seq_file,
                                 IUPAC.IUPACUnambiguousRNA(),
                                 not args.bgonly)

        if not args.bgonly:
            pssm = load_motif(args.pfm_seq,
                              args.pseudocount,
                              IUPAC.IUPACUnambiguousRNA(),
                              bg)
            seq_results = scan_main(seq_file,
                                    pssm,
                                    IUPAC.IUPACUnambiguousRNA(),
                                    bg, args)
        else:
            print(dict(bg))
            sys.exit()

    ## Structure
    if seq_type in ['SS', 'RNASS']:
        if args.testseq:
            struct_file = SeqRecord(Seq(testseq_stack.pop()))
        elif seq_type == 'SS':
            struct_file = args.fastafiles[0]
        else:
            struct_file = args.fastafiles[1]

        if not args.testseq:
            bg = load_background(args.bg_struct,
                                 args.uniform_background,
                                 struct_file,
                                 ContextualSecondaryStructure(),
                                 not args.bgonly)

        if not args.bgonly:
            pssm = load_motif(args.pfm_struct,
                              args.pseudocount,
                              ContextualSecondaryStructure(),
                              bg)
            struct_results = scan_main(struct_file,
                                       pssm,
                                       ContextualSecondaryStructure(),
                                       bg, args)
        else:
            print(dict(bg))
            sys.exit()

    if seq_type == 'RNASS':
        combined_results = combine(seq_results, struct_results)
        combined_results.to_csv(sys.stdout, sep="\t", index=False)
    elif seq_type == 'RNA':
        seq_results.to_csv(sys.stdout, sep="\t", index=False)
    else:
        struct_results.to_csv(sys.stdout, sep="\t", index=False)

    toc = time.time()

    runtime = float(toc - tic)
    if runtime > 60:
        eprint("Done in %0.4f minutes!" % (runtime / 60))
    else:
        eprint("Done in %0.4f seconds!" % (runtime))


if __name__ == '__main__':
    main()
