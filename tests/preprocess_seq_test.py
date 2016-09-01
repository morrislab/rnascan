#!/usr/bin/env python

import unittest
import sys
sys.path.append("..")
import motif_scan as ms
#from Bio import motifs, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

class PreprocessSeqTestCase(unittest.TestCase):

    def test_preprocessSeq_1(self):
        '''Test preprocess_seq() on DNA alphabet'''
        seqrec = SeqRecord(Seq('GATTACA', IUPAC.IUPACUnambiguousDNA()))
        target = ms.preprocess_seq(seqrec, IUPAC.IUPACUnambiguousRNA())
        expected = 'GAUUACA'
        self.assertEqual(str(target), expected)


    def test_preprocessSeq_2(self):
        '''Test preprocess_seq() on RNA alphabet'''
        seqrec = SeqRecord(Seq('GAUUACA', IUPAC.IUPACUnambiguousRNA()))
        target = ms.preprocess_seq(seqrec, IUPAC.IUPACUnambiguousRNA())
        expected = 'GAUUACA'
        self.assertEqual(str(target), expected)


    def test_preprocessSeq_3(self):
        '''Test preprocess_seq() on non-DNA/RNA alphabet'''
        seqrec = SeqRecord(Seq('GAUUACA', IUPAC.IUPACUnambiguousRNA()))
        target = ms.preprocess_seq(seqrec, IUPAC.IUPACUnambiguousDNA())
        expected = 'GAUUACA'
        self.assertEqual(str(target), expected)


    def test_preprocessSeq_4(self):
        '''Test preprocess_seq() on non-DNA/RNA alphabet'''
        seqrec = SeqRecord(Seq('GATTACA', IUPAC.IUPACUnambiguousDNA()))
        target = ms.preprocess_seq(seqrec, IUPAC.IUPACUnambiguousDNA())
        expected = 'GATTACA'
        self.assertEqual(str(target), expected)

if __name__ == '__main__':
    print sys.argv[0]
    unittest.main()
