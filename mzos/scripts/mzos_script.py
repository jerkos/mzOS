from __future__ import absolute_import

import sys
import argparse
import os
import logging
import time
import shutil

import yaml

from mzos import ressources
from mzos.peaklist_reader import PeakListReader
from mzos.annotator import PeakelsAnnotator
from mzos.database_finder import DatabaseSearch
from mzos.stats import StatsModel
from mzos.exp_design import ExperimentalSettings
from mzos.exp_design import IONISATION_MODE
from mzos.bayesian_inference import BayesianInferer
from mzos.results_exporter import ResultsExporter


NEG_ADDUCTS = "NEG_ADDUCTS_IMS.csv"
POS_ADDUCTS = "POS_ADDUCTS_IMS.csv"
FRAGMENTS = "FRAGMENTS_IMS.csv"
MZOS_YML = "mzos.yml"


def run_analysis(**kwargs):
    """
    function to run analysis (all the pipeline)
    :return:
    """

    xcms_pkl = kwargs['xcms_pkl']
    kwargs.pop('xcms_pkl')

    polarity = kwargs['polarity']
    kwargs.pop('polarity')

    mz_tol_ppm = kwargs['mz_tol_ppm']
    kwargs.pop('mz_tol_ppm')

    is_dims = kwargs['is_dims']
    kwargs.pop('is_dims')

    db_search = kwargs['db_search']
    kwargs.pop('db_search')

    output = kwargs['output']
    kwargs.pop('output')

    bayes = kwargs['bayes']
    kwargs.pop('bayes')

    if xcms_pkl is None or not xcms_pkl:
        raise ValueError("Supply a XCMS peaklist.")
    if not os.path.isfile(xcms_pkl):
        raise ValueError("XCMS peaklist path does not exist.")

    if polarity is None:
        raise ValueError("polarity must be 'negative' or 'positive'.")
    ionisation_mode = IONISATION_MODE.NEG if polarity == 'negative' else IONISATION_MODE.POS

    frag_conf = kwargs.get('frag_conf')
    neg_adducts_conf = kwargs.get('neg_adducts_conf')
    pos_adducts_conf = kwargs.get('pos_adducts_conf')

    exp_settings = ExperimentalSettings(mz_tol_ppm, ionisation_mode, is_dims,
                                        frag_conf=frag_conf,
                                        neg_adducts_conf=neg_adducts_conf,
                                        pos_adducts_conf=pos_adducts_conf)

    # clear logging
    logging.getLogger('').handlers = []

    logging.basicConfig(level=logging.INFO)
    t1 = time.clock()

    peakels = PeakListReader(xcms_pkl, exp_settings).get_peakels()
    logging.info("Peaklist loaded.")

    # annotation
    peakels_annotator = PeakelsAnnotator(peakels, exp_settings)
    logging.info("Annotating...")

    best_monos = peakels_annotator.annotate()
    logging.info("Monoisotopic found: #{0}".format(len(best_monos)))

    # database finding
    db = []
    for d in db_search:
        if d not in ('hmdb', 'lmsd'):
            logging.warn('Error specifying db (got {}), only hmdb and lmsd are supported...'.format(d))
        else:
            db.append(d)
    db = '+'.join(db)

    search = DatabaseSearch(db, exp_settings)
    logging.info("Searching in database...")
    adducts_l = ['H1']
    nb_metabs, not_found = search.assign_formula(peakels, adducts_l, exp_settings.mz_tol_ppm)
    logging.info("Found #{} metabolites, #{} "
                 "elution peak with no metabolite assignments".format(nb_metabs, not_found))

    # scoring first simplistic
    model = StatsModel(peakels, exp_settings.mz_tol_ppm * 1.5)
    logging.info("Compute score 1....")
    # populate annotations objects
    model.calculate_score()
    logging.info("Done.")

    # scoring bayesian inferer
    if bayes:
        bi = BayesianInferer(peakels, exp_settings)
        logging.info("Compute score 2...")
        # populate annotations object
        bi.infer_assignment_probabilities()
        # logging.info("Finished")

    logging.info('Exporting results...')
    exporter = ResultsExporter(output, sorted(peakels, key=lambda _: _.id))
    exporter.write()
    logging.info("Done.")


def main():
    """
    :return:
    """

    default_opts = {
        'xcms_pkl': 'peaklist.csv',
        'polarity': 'negative',
        'mz_tol_ppm': 5,
        'is_dims': False,
        'db_search': ['hmdb', 'lmsd'],
        'bayes': True,
        'output': 'results.csv'
    }

    current_files = set(os.listdir(os.curdir))
    missing = []

    absolute_path = os.path.normpath(os.path.split(ressources.__file__)[0])

    for f in (NEG_ADDUCTS, POS_ADDUCTS, FRAGMENTS):
        if f not in current_files:

            missing.append(f)
            # copy the current file
            shutil.copy(os.path.join(absolute_path, f), os.curdir)

    # dealing with
    if MZOS_YML not in current_files:
        with open(MZOS_YML, 'w') as mzos_yml:
            mzos_yml.write(yaml.dump(default_opts))
            missing.append(MZOS_YML)

    if missing:
        print("The following configuration files are missing: ")
        for f in missing:
            print(f)
        print("Please modify in order to match your experimental conditions. Exiting...")
        sys.exit()

    # load yml data
    logging.info("loading configuration")
    data = {}
    try:
        data = yaml.load(open(MZOS_YML, 'r').read())
    except yaml.ParserError as e:
        logging.error("Found error in mzos.yml: {}. Please modify the mzos.yml file.".format(e.message()))
    data['frag_conf'] = os.path.abspath(FRAGMENTS)
    data['neg_adducts_conf'] = os.path.abspath(NEG_ADDUCTS)
    data['pos_adducts_conf'] = os.path.abspath(POS_ADDUCTS)

    run_analysis(**data)


if __name__ == '__main__':
    main()
