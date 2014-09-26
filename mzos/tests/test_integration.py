__author__ = 'Marc'

import unittest
import os.path as op
import os
import zipfile
import shutil

from mzos.exp_design import ExperimentalSettings
from mzos.peaklist_reader import PeakListReader
from mzos.annotator import PeakelsAnnotator
from mzos.stats import StatsModel
from mzos.bayesian_inference import BayesianInferer
from mzos.results_exporter import ResultsExporter
from mzos.database_finder import DatabaseSearch


class TestResultExporter(unittest.TestCase):

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

    def test_result_exporter(self):
        polarity = -1
        exp_settings = ExperimentalSettings(10, polarity, is_dims_experiment=False)

        peakels = PeakListReader(op.normcase("mzos/ressources/peaks_matrix_NEG.tsv"), exp_settings).get_peakels()

        ##annotation##
        peakels_annotator = PeakelsAnnotator(peakels, exp_settings)
        best_monos = peakels_annotator.annotate_()

         ##database finding##
        db_search = DatabaseSearch('hmdb', exp_settings)
        nb_metabs, not_found = db_search.assign_formula(peakels, exp_settings.mz_tol_ppm)

        ##scoring
        #first simplistic
        model = StatsModel(peakels, exp_settings.mz_tol_ppm * 1.5)
        #populate annotations objects
        model.calculate_score()

        ##scoring 2##
        bi = BayesianInferer(peakels, exp_settings)
        #populate annotations object
        bi.infer_assignment_probabilities()

        exporter = ResultsExporter(op.normcase("annotations.tsv"), peakels)
        exporter.write()