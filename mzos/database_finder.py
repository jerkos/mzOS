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

import logging
import sqlite3
from collections import defaultdict as ddict
from itertools import izip
import multiprocessing
from collections import namedtuple
from feature import Annotation

Metabolite = namedtuple("Metabolite", "acession, name, formula, inchi, mono_mass, average_mass, "
                                      "description, status, origin, kegg_id, isotopic_pattern_pos, "
                                      "isotopic_pattern_neg")


class IDatabaseSearcher(object):
    """
    Interface of databse search
    """
    def search_moz(self, moz, moz_tol_ppm):
        """
        :param moz:
        :param moz_tol_ppm:
        """
        raise NotImplementedError

    def search_formula(self, molecular_formula):
        """
        :param molecular_formula:
        """
        raise NotImplementedError


def search_metabolites_for(args):
    """
    pickling problem if use inside class
    """
    database, feature, with_tol_ppm = args[0], args[1], args[2]
    mass = feature.get_real_mass()
    tol_da = mass * with_tol_ppm / 1e6
    min_mass, max_mass = mass - tol_da, mass + tol_da

    c = sqlite3.connect(database).cursor()
    metabolites = []
    for row in c.execute('select * from metabolite where mono_mass >=  ? and mono_mass <= ?', (min_mass, max_mass)):
        m = Metabolite._make(row)
        if m.kegg_id is not None:
            metabolites.append(m)  # got warning du to the underscore
    metabolites.sort(key=lambda _: abs(_.mono_mass - mass))
    return metabolites


class DatabaseSearch(IDatabaseSearcher):

    HMDB_FILE = "ressources/hmdb.sqlite"

    def __init__(self, bank, exp_design):
        self.exp_design = exp_design
        self.metabolites_by_feature = {}
        self.bank = 'hmdb' if bank not in ['hmdb, kegg'] else bank
        logging.info("Performing database search in {} {}".format(self.bank, 'v3.5'))

    def assign_formula(self, features, with_tol_ppm=10.0):
        """
        assign molecular formula to features using multiprocessing module
        :param features: list or set ? of features
        :param with_tol_ppm: mz tolerance in order to perform the look up
        :return: dictionary with key: feature, value: list of metabolites
        """
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
        args = [(self.HMDB_FILE, f_, with_tol_ppm) for f_ in features]
        metabs = pool.map(search_metabolites_for, args, chunksize=20)
        pool.close()
        # create Annotation objects
        m_count, not_found = 0, 0
        for f, metabs in izip(features, metabs):
            if not metabs:
                not_found += 1
            else:
                m_count += len(metabs)
            f.annotations = [Annotation(m) for m in metabs]
        return m_count, not_found


if __name__ == '__main__':
    from peaklist_reader import PeakListReader
    from annotator import PeakelsAnnotator
    from stats import StatsModel
    from exp_design import ExperimentalSettings
    import os.path
    import time

    logging.basicConfig(level=logging.INFO)

    t1 = time.clock()

    polarity = -1
    path = os.path.normcase("""../../tests/data/peaks_matrix_NEG.tsv""")

    exp = ExperimentalSettings(10.0, -1)

    peakels = PeakListReader(path, exp).get_peakels()
    logging.info("Data loaded.")

    dbSearch = DatabaseSearch('hmdb', None)
    metabolites_by_feature = dbSearch.assign_formula(peakels)
    logging.info("search done for #{} peakels".format(len(metabolites_by_feature)))
