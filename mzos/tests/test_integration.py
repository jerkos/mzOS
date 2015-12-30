import logging
import shutil
import unittest
import zipfile
import os
import os.path as op

from mzos.exp_design import ExperimentalSettings
from mzos.peaklist_reader import PeakListReader
from mzos.annotator import PeakelsAnnotator
from mzos.stats import StatsModel
from mzos.bayesian_inference import BayesianInferer
from mzos.results_exporter import ResultsExporter
from mzos.database_finder import DatabaseSearch
from mzos.tests import WithHMDBMixin


class TestResultExporter(WithHMDBMixin, unittest.TestCase):

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

    def test_result_exporter(self):
        polarity = -1
        exp_settings = ExperimentalSettings(10, polarity, is_dims_experiment=False)

        peakels = PeakListReader(op.normcase("mzos/ressources/peaks_matrix_NEG.tsv"), exp_settings).get_peakels()

        # annotation
        peakels_annotator = PeakelsAnnotator(peakels, exp_settings)
        peakels_annotator.annotate()

        # database finding
        db_search = DatabaseSearch('hmdb', exp_settings)
        adducts_l = ['H1']
        db_search.assign_formula(peakels, adducts_l, exp_settings.mz_tol_ppm)

        # scoring
        # first simplistic
        model = StatsModel(peakels, exp_settings.mz_tol_ppm * 1.5)
        # populate annotations objects
        model.calculate_score()

        # scoring 2
        bi = BayesianInferer(peakels, exp_settings)
        # populate annotations object
        bi.infer_assignment_probabilities()

        exporter = ResultsExporter(op.normcase("annotations.tsv"), peakels)
        exporter.write()
