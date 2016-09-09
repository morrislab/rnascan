import unittest
import sys
sys.path.append("..")
import combine_pfms as cp
import pandas as pd


class CombinePfmsTestCase(unittest.TestCase):

    def test_load_pfm_1(self):
        '''Test loading of PFM and removing first column
        '''
        target = cp.load_pfm('test_seq_pfm.txt')
        self.assertEqual(target.shape[1], 4)

    def test_combine_pfms_1(self):
        '''Test combining two test PFMs
        '''
        seqpfm = cp.load_pfm('test_seq_pfm.txt')
        structpfm = cp.load_pfm('test_struct_pfm.txt')

        target = cp.combine_pfms(seqpfm, structpfm)

        # Test dimensions
        self.assertEqual(target.shape[0], seqpfm.shape[0])
        self.assertEqual(target.shape[1], seqpfm.shape[1]*structpfm.shape[1])

        # Test products
        self.assertAlmostEqual(target.loc[0, 'A'], 0.043247417362, places=4)
        self.assertAlmostEqual(target.loc[3, 'Z'], 0.015526, places=4)
        self.assertAlmostEqual(target.loc[1, '['], 0.026736, places=4)
        self.assertAlmostEqual(target.loc[2, ']'], 0.0750003, places=4)

if __name__ == '__main__':
    print sys.argv[0]
    unittest.main()
