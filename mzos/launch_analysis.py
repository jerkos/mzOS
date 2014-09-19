# Copyright (C) 2014  omics-services.com
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

__email__ = 'marc.dubois@omics-services.com'

import sys
import argparse
import os
import logging
import time

from peaklist_reader import PeakListReader
from annotator import PeakelsAnnotator
from database_finder import DatabaseSearch
from stats import StatsModel
from exp_design import ExperimentalSettings
from exp_design import IONISATION_MODE
from bayesian_inference import BayesianInferer
from utils import merge
from results_exporter import ResultsExporter


def get_arg_parser():
    """

    @return:
    """
    parser = argparse.ArgumentParser(description="Annotation analysis of xcms peaklist")
    parser.add_argument("--xcms_pkl", type=str)
    parser.add_argument("--polarity", type=int, default=-1)
    parser.add_argument("--mz_tol_ppm", type=float, default=10.0)
    parser.add_argument("--dims", type=bool, default=False)
    parser.add_argument("--output", type=str, default="annotations.tsv")
    return parser


if __name__ == '__main__':
    arg_parser = get_arg_parser()
    arguments = arg_parser.parse_args(sys.argv[1:])
    if arguments.xcms_pkl is None or not arguments.xcms_pkl:
        raise ValueError("Supply a XCMS peaklist.")
    if not os.path.isfile(arguments.xcms_pkl):
        raise ValueError("XCMS peaklist path does not seem to be valid.")

    if arguments.polarity is None:
        raise ValueError("polarity must be '-1' or '+1.'")
    polarity = IONISATION_MODE.NEG if arguments.polarity <= 0 else IONISATION_MODE.POS
    exp_settings = ExperimentalSettings(arguments.mz_tol_ppm, polarity, arguments.dims)

    #clear logging
    logging.getLogger('').handlers = []

    logging.basicConfig(level=logging.INFO)
    t1 = time.clock()

    peakels = PeakListReader(arguments.xcms_pkl, exp_settings).get_peakels()
    logging.info("Peaklist loaded.")

    ##annotation##
    peakels_annotator = PeakelsAnnotator(peakels, exp_settings)
    logging.info("Annotating...")

    best_monos = peakels_annotator.annotate_()
    logging.info("Monoisotopic found: #{}".format(len(best_monos)))

    ##database finding##
    dbSearch = DatabaseSearch('hmdb', exp_settings)
    logging.info("Searching in database...")
    nb_metabs, not_found = dbSearch.assign_formula(peakels, exp_settings.mz_tol_ppm)
    logging.info("Found #{} metabolites, #{} "
                 "elution peak with no metabolite assignments".format(nb_metabs, not_found))

    ##scoring
    #first simplistic
    model = StatsModel(peakels, exp_settings.mz_tol_ppm * 1.5)
    logging.info("Compute score 1....")
    #populate annotations objects
    model.calculate_score()
    logging.info("Done")


    ##scoring 2##
    bi = BayesianInferer(peakels, exp_settings)
    logging.info("Compute score 2...")
    #populate annotations object
    bi.infer_assignment_probabilities()
    logging.info("Done")

    exporter = ResultsExporter(arguments.output, peakels)
    exporter.write()