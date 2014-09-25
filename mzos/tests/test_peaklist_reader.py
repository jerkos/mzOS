__author__ = 'Marc'

import unittest
import os.path as op
from mzos.peaklist_reader import PeakListReader
from mzos.exp_design import ExperimentalSettings


class TestPeakListReader(unittest.TestCase):
    FILE = op.normcase("mzos/ressources/peaks_matrix_NEG.tsv")

    def test_read_tsv(self):
        pkl_reader = PeakListReader(self.FILE, ExperimentalSettings(mz_tol_ppm=10.0))
        peakels = pkl_reader.get_peakels()
        self.assertEqual(len(peakels), 3238)
        self.assertTrue(all([isinstance(p.moz, float) for p in peakels]))
        self.assertTrue(all([p.moz != 0.0] for p in peakels))