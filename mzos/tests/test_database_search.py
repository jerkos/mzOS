__author__ = 'marc.dubois@omics-services.com'

import unittest

from mzos.scripts.hmdb_sqlite_creator import ELEMENT_PATTERN, add_element, remove_element
# from mzos.database_finder import DatabaseSearch
# from mzos.feature import Peakel
# from mzos.exp_design import ExperimentalSettings


class TestDatabaseSearch(unittest.TestCase):

    def test_add_formula(self):
        f = {x[0]: x[1] for x in ELEMENT_PATTERN.findall('C6H6O12')}
        s = add_element(f, 'C', 12)
        self.assertEqual('C18H6O12', s)

    def test_remove_formula(self):
        f = {x[0]: x[1] for x in ELEMENT_PATTERN.findall('C6H6O12')}
        s = remove_element(f, 'C', 12)
        self.assertEqual('H6O12', s)

        f = {x[0]: x[1] for x in ELEMENT_PATTERN.findall('C6H6O12')}
        s = remove_element(f, 'C', 5)
        self.assertEqual('CH6O12', s)

        f = {x[0]: x[1] for x in ELEMENT_PATTERN.findall('C6H6O12')}
        s = remove_element(f, 'C', 2)
        self.assertEqual('C4H6O12', s)

    # WONT PASS ON TRAVIS
    # def test_database_search(self):
    #     mass_fruc_6p = 260.029718526
    #     name = "Fructose 6-phosphate"
    #     peakel = Peakel(mass_fruc_6p - 1.007276, 0.0, 0.0, 0.0)
    #     peakel.charge = -1
    #     peakels = [peakel]
    #     dbSearch = DatabaseSearch('hmdb', None)
    #     metabolites_by_feature = dbSearch.assign_formula(peakels)
    #     m_names = set()
    #     for f, metabs in metabolites_by_feature.iteritems():
    #         for m in metabs:
    #             m_names.add(m.name)
    #     self.assertIn(name, m_names)

