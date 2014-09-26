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
import os.path as op
from collections import defaultdict as ddict
from itertools import izip
import multiprocessing
from collections import namedtuple
from feature import Annotation


class MolecularEntity(object):
    """
    Molecular entity
    """
    def __init__(self):
        self.name = None

        self.kegg_id = None
        self.hmdb_id = None

        self.description = None
        self.formula = None
        self.inchi = None

        self.mono_mass = None

        self.isotopic_pattern_neg = None
        self.isotopic_pattern_pos = None


class Metabolite(MolecularEntity):
    """
    Metabolite entity
    """
    COLUMNS = "acession,name,formula,inchi,mono_mass,average_mass," \
              "description,status,origin,kegg_id,isotopic_pattern_pos,isotopic_pattern_neg".split(",")

    MAPPING = {'acession': 'hmdb_id'}

    def __init__(self, *args):
        super(MolecularEntity, self).__init__()
        for name, value in dict(izip(Metabolite.COLUMNS, args)).iteritems():
            if name in Metabolite.MAPPING:
                setattr(self, Metabolite.MAPPING[name], value)
            else:
                setattr(self, name, value)

class Lipid(MolecularEntity):
    pass

# Metabolite = namedtuple("Metabolite", "acession, name, formula, inchi, mono_mass, average_mass, "
#                                       "description, status, origin, kegg_id, isotopic_pattern_pos, "
#                                       "isotopic_pattern_neg")
# Lipid = namedtuple("Lipid", 'systematic_name, kegg_id, inchi_key')


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

    @staticmethod
    def get_moz_bounds(feature, mz_tol_ppm):
        """
        :param feature:
        :param mz_tol_ppm:
        """
        mass = feature.get_real_mass()
        tol_da = mass * mz_tol_ppm / 1e6
        return mass, mass - tol_da, mass + tol_da


def search_metabolites_for(args):
    """
    pickling problem if use inside class
    """
    database, feature, with_tol_ppm = args[0], args[1], args[2]
    #mass = feature.get_real_mass()
    # tol_da = mass * with_tol_ppm / 1e6
    # min_mass, max_mass = mass - tol_da, mass + tol_da
    mass, min_mass, max_mass = IDatabaseSearcher.get_moz_bounds(feature, with_tol_ppm)

    conn = sqlite3.connect(database)
    c = conn.cursor()
    metabolites = []
    for row in c.execute('select * from metabolite where mono_mass >=  ? and mono_mass <= ?', (min_mass, max_mass)):
        m = Metabolite(*row)#Metabolite._make(row)
        if m.kegg_id is not None:
            metabolites.append(m)  # got warning du to the underscore
    conn.close()
    metabolites.sort(key=lambda _: abs(_.mono_mass - mass))
    return metabolites


# def search_lipids_for(args):
#     database, feature, with_tol_ppm = args[0], args[1], args[2]
#     mass, min_mass, max_mass = IDatabaseSearcher.get_moz_bounds(feature, with_tol_ppm)
#
#     conn = sqlite3.connect(database)
#     c = conn.cursor()
#     lipids = []
#     for row in c.execute('select SYSTEMATIC_NAME, KEGG_ID, INCHI_KEY from lipids where EXACT_MASS >= ? and EXACT_MASS <= ?', (min_mass, max_mass)):
#         l = Lipid._r


class DatabaseSearch(IDatabaseSearcher):
    """
    :param bank:
    :param exp_design:
    """
    HMDB_FILE = op.normcase("mzos/ressources/hmdb.sqlite")
    LMSD_FILE = op.normcase("ressources/lmsd.sqlite")

    def __init__(self, bank, exp_design):
        self.exp_design = exp_design
        self.metabolites_by_feature = {}
        self.bank = 'hmdb' if bank not in ['hmdb, kegg'] else bank  #self.exp_design.databases
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
