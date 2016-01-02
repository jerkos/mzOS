from __future__ import absolute_import
import unittest

from mzos.database_finder import DatabaseSearch
from mzos.database_finder import get_moz_bounds
from mzos.feature import Peakel
from mzos.formula import Formula
from mzos.tests import WithHMDBMixin


class TestDatabaseSearch(WithHMDBMixin, unittest.TestCase):
    def setUp(self):
        TestDatabaseSearch.unzip_hmdb()

    def tearDown(self):
        TestDatabaseSearch.remove_hmdb()

    def test_database_search(self):
        mass_fruc_6p = 260.029718526
        name = "Fructose 6-phosphate"
        peakel = Peakel(mass_fruc_6p - 1.007276, 0.0, 0.0, 0.0)
        peakel.charge = 1
        peakel.polarity = -1
        peakels = [peakel]
        db_search = DatabaseSearch('hmdb + lmsd', None)
        db_search.assign_formula(peakels, ['H1'])
        m_names = set()
        for f in peakels:
            for annot in f.annotations:
                m_names.add(annot.metabolite.name)
        self.assertIn(name, m_names)

    def test_moz_bounds(self):
        f2 = Peakel(260.029718526 + 1.007276, 0.0, 0.0, 260.029718526)
        m, m_min, m_max = get_moz_bounds(f2, Formula.from_str('H1'), 10)
        self.assertAlmostEqual(m, 260.029718526, places=2)


