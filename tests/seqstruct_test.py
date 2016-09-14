import unittest
import sys
sys.path.append("..")
from motif_scan.BioAddons.SeqStruct import SeqStruct
from motif_scan.BioAddons.Alphabet import *


class SeqStructTestCase(unittest.TestCase):

    def test_convert_1(self):
        '''Test converting RNA sequence and secondary structure sequence
        to ContextualSequenceSecondaryStructure sequence
        '''
        target = SeqStruct.convert('AG', 'BE')
        expected = 'KA'
        self.assertEqual(target, expected)

    def test_convert_2(self):
        '''Test converting RNA sequence and secondary structure sequence
        to ContextualSequenceSecondaryStructure sequence
        '''
        target = SeqStruct.convert('GCAUG', 'EEEEE')
        expected = 'AVHOA'
        self.assertEqual(target, expected)

    def test_convert_3(self):
        '''Test converting RNA sequence and secondary structure sequence
        to ContextualSequenceSecondaryStructure sequence
        '''
        target = SeqStruct.convert('CCCC', 'MMRR')
        expected = ']][['
        self.assertEqual(target, expected)

    def test_create_seqstruct_1(self):
        '''Test instantiation of a SeqStruct object
        '''
        target = SeqStruct('AG', 'BE')
        self.assertIsInstance(target, SeqStruct)
        self.assertIsInstance(target.alphabet,
                              ContextualSequenceSecondaryStructure)

    def test_create_seqstruct_2(self):
        '''Test instantiation of a SeqStruct object but the input sequences
        have different lengths
        '''
        with self.assertRaises(ValueError):
            SeqStruct('GATTACA', 'KEVIN')

    def test_reverse_convert_1(self):
        '''Test reverse convert of a SeqStruct sequence
        '''
        seqstructobj = SeqStruct('GCAUGAAAAAACC', 'BEBEBHLMTRLMR')
        target = seqstructobj.reverse_convert()
        expected = ('GCAUGAAAAAACC', 'BEBEBHLMTRLMR')
        self.assertEqual(target, expected)

    def test_reverse_convert_2(self):
        '''Test reverse convert of a SeqStruct sequence and return substring
        from beginning of sequence
        '''
        seqstructobj = SeqStruct('GCAUGAAAAAACC', 'BEBEBHLMTRLMR')
        target = seqstructobj.reverse_convert(start=0, end=4)
        expected = ('GCAU', 'BEBE')
        self.assertEqual(target, expected)

    def test_reverse_convert_3(self):
        '''Test reverse convert of a SeqStruct sequence and return substring
        from the end of the sequence
        '''
        seqstructobj = SeqStruct('GCAUGAAAAAACC', 'BEBEBHLMTRLMR')
        target = seqstructobj.reverse_convert(start=10, end=13)
        expected = ('ACC', 'LMR')
        self.assertEqual(target, expected)

    def test_reverse_convert_4(self):
        '''Test reverse convert of a SeqStruct sequence and return substring
        but with an invalid start position
        '''
        seqstructobj = SeqStruct('GCAUGAAAAAACC', 'BEBEBHLMTRLMR')
        with self.assertRaises(ValueError):
            target = seqstructobj.reverse_convert(start=13, end=13)

    def test_reverse_convert_5(self):
        '''Test reverse convert of a SeqStruct sequence and return substring
        but with an end position greater than start position
        '''
        seqstructobj = SeqStruct('GCAUGAAAAAACC', 'BEBEBHLMTRLMR')
        with self.assertRaises(ValueError):
            target = seqstructobj.reverse_convert(start=10, end=5)

if __name__ == '__main__':
    print sys.argv[0]
    unittest.main()
