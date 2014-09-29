__author__ = 'marc.dubois@omics-services.com'

import unittest
import zipfile
import os.path as op
import os
import shutil

from mzos.database_finder import DatabaseSearch
from mzos.feature import Peakel


class TestDatabaseSearch(unittest.TestCase):

    def setUp(self):
        z = zipfile.ZipFile(op.normcase('mzos/ressources/hmdb.zip'))
        self.hmdb_path = z.extract('hmdb.sqlite')
        print("Moving extracted archive...")
        shutil.move(self.hmdb_path, 'mzos/ressources/hmdb.sqlite')
        print("Done")

    def tearDown(self):
        print("removing 'hmdb.sqlite'...")
        try:
            os.remove(op.normcase('mzos/ressources/hmdb.sqlite'))
            print("Done")
        except OSError:
            pass

    def test_database_search(self):
        mass_fruc_6p = 260.029718526
        name = "Fructose 6-phosphate"
        peakel = Peakel(mass_fruc_6p - 1.007276, 0.0, 0.0, 0.0)
        peakel.charge = 1
        peakel.polarity = -1
        peakels = [peakel]
        db_search = DatabaseSearch('hmdb', None)
        db_search.assign_formula(peakels)
        m_names = set()
        for f in peakels:
            for annot in f.annotations:
                m_names.add(annot.metabolite.name)
        self.assertIn(name, m_names)

