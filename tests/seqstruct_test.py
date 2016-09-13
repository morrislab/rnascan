import unittest
import sys
sys.path.append("..")
from BioAddons.SeqStruct import SeqStruct
from BioAddons.Alphabet import *


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
        seqstructobj = SeqStruct('GCAUGAAAAAAAAAAAA', 'BEBEBHLMRTLMTRLMR')
        target = seqstructobj.reverse_convert()
        expected = ('GCAUGAAAAAAAAAAAA', 'BEBEBHLMRTLMTRLMR')
        self.assertEqual(target, expected)

if __name__ == '__main__':
    print sys.argv[0]
    unittest.main()
