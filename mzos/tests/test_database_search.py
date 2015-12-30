import logging
import shutil
import unittest
import zipfile
import os
import os.path as op

from mzos.database_finder import DatabaseSearch
from mzos.feature import Peakel
from mzos.tests import WithHMDBMixin


class TestDatabaseSearch(WithHMDBMixin, unittest.TestCase):

    def setUp(self):
        z = zipfile.ZipFile(op.abspath('mzos/ressources/hmdb.zip'))
        self.hmdb_path = z.extract('hmdb.sqlite')
        logging.info("Moving extracted archive...")
        shutil.move(self.hmdb_path, 'mzos/ressources/hmdb.sqlite')
        logging.info("Done")

    def tearDown(self):
        logging.info("removing 'hmdb.sqlite'...")
        try:
            os.remove(op.normcase('mzos/ressources/hmdb.sqlite'))
            logging.info("Done")
        except OSError:
            logging.error("Unable to remove sqlite file or file does not exist")

    def test_database_search(self):
        mass_fruc_6p = 260.029718526
        name = "Fructose 6-phosphate"
        peakel = Peakel(mass_fruc_6p - 1.007276, 0.0, 0.0, 0.0)
        peakel.charge = 1
        peakel.polarity = -1
        peakels = [peakel]
        db_search = DatabaseSearch('hmdb', None)
        db_search.assign_formula(peakels, ['H1'])
        m_names = set()
        for f in peakels:
            for annot in f.annotations:
                m_names.add(annot.metabolite.name)
        self.assertIn(name, m_names)
