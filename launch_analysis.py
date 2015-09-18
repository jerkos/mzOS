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
from numba.cuda.cudadrv.driver import require_device_memory

__email__ = 'marc.dubois@omics-services.com'

import sys
import argparse
import os
import logging
import time

from mzos.peaklist_reader import PeakListReader
from mzos.annotator import PeakelsAnnotator
from mzos.database_finder import DatabaseSearch
from mzos.stats import StatsModel
from mzos.exp_design import ExperimentalSettings
from mzos.exp_design import IONISATION_MODE
from mzos.bayesian_inference import BayesianInferer
from mzos.results_exporter import ResultsExporter
from gooey import Gooey


def get_arg_parser():

    parser = argparse.ArgumentParser(prog="mzOS", description="Annotation analysis of xcms peaklist")
    parser.add_argument("-x", "--xcms_pkl", type=str, help="path to the xcms peaklist", required=True)
    parser.add_argument("-p", "--polarity", default='negative', choices=['negative', 'positive'],
                        help='experiment polarity', required=True)
    parser.add_argument("--mz_tol_ppm", type=float, default=10.0, help='mass over charge tolerance', required=False)
    parser.add_argument("--dims", default=False, action='store_true', help='direct infusion MS experiment', required=False)
    parser.add_argument('--db', default='hmdb', choices=['hmdb', 'lmsd', 'hmdb + lmsd'], required=False)
    parser.add_argument("--output", type=str, default="annotations.tsv", required=False)
    parser.add_argument("--bayes", default=True, required=False)
    return parser


# @Gooey(advanced=True, config=True, program_name='mzOS')
def main():
    arg_parser = get_arg_parser()
    arguments = arg_parser.parse_args(sys.argv[1:])
    if arguments.xcms_pkl is None or not arguments.xcms_pkl:
        raise ValueError("Supply a XCMS peaklist.")
    if not os.path.isfile(arguments.xcms_pkl):
        raise ValueError("XCMS peaklist path does not seem to be valid.")

    if arguments.polarity is None:
        raise ValueError("polarity must be '-1' or '+1.'")
    polarity = IONISATION_MODE.NEG if arguments.polarity == 'negative' else IONISATION_MODE.POS
    exp_settings = ExperimentalSettings(arguments.mz_tol_ppm, polarity, arguments.dims)

    # clear logging
    logging.getLogger('').handlers = []

    logging.basicConfig(level=logging.INFO)
    t1 = time.clock()

    peakels = PeakListReader(arguments.xcms_pkl, exp_settings).get_peakels()
    logging.info("Peaklist loaded.")

    # annotation
    peakels_annotator = PeakelsAnnotator(peakels, exp_settings)
    logging.info("Annotating...")

    best_monos = peakels_annotator.annotate()
    logging.info("Monoisotopic found: #{}".format(len(best_monos)))

    # database finding

    db = 'hmdb'
    if arguments.db == 'hmdb':
        pass
    elif arguments.db == 'lmsd':
        db = 'lmsd'
    else:
        logging.warn('Error specifying db, default to hmdb...')
    db_search = DatabaseSearch(db, exp_settings)
    logging.info("Searching in database...")
    adducts_l = ['H1']
    nb_metabs, not_found = db_search.assign_formula(peakels, adducts_l, exp_settings.mz_tol_ppm)
    logging.info("Found #{} metabolites, #{} "
                 "elution peak with no metabolite assignments".format(nb_metabs, not_found))

    # scoring first simplistic
    model = StatsModel(peakels, exp_settings.mz_tol_ppm * 1.5)
    logging.info("Compute score 1....")
    # populate annotations objects
    model.calculate_score()
    logging.info("Done.")

    # scoring bayesian inferer
    bi = BayesianInferer(peakels, exp_settings)
    logging.info("Compute score 2...")
    # populate annotations object
    bi.infer_assignment_probabilities()
    # logging.info("Finished")

    logging.info('Exporting results...')
    exporter = ResultsExporter(arguments.output, sorted(peakels, key=lambda _: _.id))
    exporter.write()
    logging.info("Done.")

if __name__ == '__main__':
    main()