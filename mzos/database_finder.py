import logging
import sqlite3
import os.path as op
from itertools import izip
import multiprocessing

from mzos.feature import Annotation
from mzos.formula import Formula


class MolecularEntity(object):
    """
    Molecular entity
    """
    def __init__(self):
        self.name = None

        self.kegg_id = None
        self.hmdb_id = None
        self.lm_id = None

        self.description = None
        self.formula = None
        self.inchi_key = None

        self.mono_mass = None

        self.isotopic_pattern_neg = None
        self.isotopic_pattern_pos = None


class Metabolite(MolecularEntity):
    """
    Metabolite entity
    """
    COLUMNS = "acession,name,formula,inchi,mono_mass,average_mass," \
              "description,status,origin,kegg_id,isotopic_pattern_pos,isotopic_pattern_neg".split(",")

    MAPPING = {'acession': 'hmdb_id', 'inchi': 'inchi_key'}

    def __init__(self, *args):
        MolecularEntity.__init__(self)
        for name, value in dict(izip(Metabolite.COLUMNS, args)).iteritems():
            if name in Metabolite.MAPPING:
                setattr(self, Metabolite.MAPPING[name], value)
            else:
                setattr(self, name, value)


class Lipid(MolecularEntity):
    def __init__(self, lm_id):
        MolecularEntity.__init__(self)

        self.lm_id = lm_id


def get_moz_bounds(feature, for_adduct, mz_tol_ppm):
    """
    :param for_adduct: Formula
    :param feature:
    :param mz_tol_ppm:
    """
    mass = feature.moz * feature.charge
    ad_mass = for_adduct.mono_mass()
    if feature.polarity < 0:
        mass += ad_mass
    else:
        mass -= ad_mass

    tol_da = mass * mz_tol_ppm / 1e6
    return mass, mass - tol_da, mass + tol_da


def search_metabolites_for(args):
    """
    pickling problem if use inside class
    :param args: database, feature, with_tol_ppm
    """
    database, feature, formula, with_tol_ppm = args[0], args[1], args[2], args[3]
    mass, min_mass, max_mass = get_moz_bounds(feature, formula, with_tol_ppm)

    conn = sqlite3.connect(database)
    c = conn.cursor()
    metabolites = []
    for row in c.execute('select * from metabolite where mono_mass >=  ? and mono_mass <= ?', (min_mass, max_mass)):
        m = Metabolite(*row)  # Metabolite._make(row)  got warning du to the underscore
        if m.kegg_id is not None:
            metabolites.append(m)
    conn.close()
    metabolites.sort(key=lambda _: abs(_.mono_mass - mass))
    return metabolites


def search_lipids_for(args):
    database, feature, formula, with_tol_ppm = args  # args[0], args[1], args[2], args[3]
    mass, min_mass, max_mass = get_moz_bounds(feature, formula, with_tol_ppm)

    conn = sqlite3.connect(database)
    c = conn.cursor()
    lipids = []
    for row in c.execute('select SYSTEMATIC_NAME, FORMULA, EXACT_MASS, LM_ID, KEGG_ID, HMDBID, INCHI_KEY '
                         'from lipids where EXACT_MASS >= ? and EXACT_MASS <= ?', (min_mass, max_mass)):
        name, chem_formula, exact_mass, lm_id, kegg_id, hmdb_id, inchi = row
        l = Lipid(lm_id)
        # warning lmsd names can be empty
        l.name = name or ''
        l.mono_mass = exact_mass or ''
        l.formula = chem_formula or ''
        l.kegg_id = kegg_id or ''
        l.hmdb_id = hmdb_id or ''
        l.inchi_key = inchi or ''
        lipids.append(l)
    conn.close()
    lipids.sort(key=lambda _: abs(_.mono_mass - mass))
    return lipids


class DatabaseSearch(object):
    """
    :param bank:
    :param exp_design:
    """
    HMDB_FILE = op.abspath("mzos/ressources/hmdb.sqlite")
    LMSD_FILE = op.abspath("mzos/ressources/lmsd.sqlite")

    def __init__(self, bank, exp_design):
        self.exp_design = exp_design
        self.metabolites_by_feature = {}
        self.bank = 'hmdb' if bank not in {'hmdb', 'kegg', 'lmsd', 'hmdb + lmsd', 'lmsd + hmdb'} else bank  # self.exp_design.databases
        logging.info("Performing database search in {} {}".format(self.bank, 'v3.5'))

    def assign_formula(self, features, for_adducts, with_tol_ppm=10.0):
        """
        assign molecular formula to features using multiprocessing module

        :param for_adducts: string adducts list
        :param features: list or set ? of features
        :param with_tol_ppm: mz tolerance in order to perform the look up
        :return: dictionary with key: feature, value: list of metabolites
        """
        m_count, not_found = 0, 0
        for for_adduct in for_adducts:
            formula = Formula.from_str(for_adduct)
            logging.info("searching for adducts: {}".format(str(formula)))

            pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

            metabolites = []
            if self.bank == 'hmdb':
                args = [(self.HMDB_FILE, f, formula, with_tol_ppm) for f in features]
                metabolites = pool.map(search_metabolites_for, args, chunksize=20)
            elif self.bank == 'lmsd':
                args = [(self.LMSD_FILE, f, formula, with_tol_ppm) for f in features]
                metabolites = pool.map(search_lipids_for, args, chunksize=20)
            elif self.bank in {'lmsd + hmdb', 'hmdb + lmsd'}:
                logging.info('Searching in LMSD...')
                args_lmsd = [(self.LMSD_FILE, f, formula, with_tol_ppm) for f in features]
                metabs_lmsd = pool.map(search_lipids_for, args_lmsd, chunksize=20)

                logging.info('Searching in HMDB...')
                args_hmdb = [(self.HMDB_FILE, f, formula, with_tol_ppm) for f in features]
                metabs_hmdb = pool.map(search_metabolites_for, args_hmdb, chunksize=20)

                # merge the 2 results set
                for lmsd_met, hmdb_met in izip(metabs_lmsd, metabs_hmdb):
                    metabolites.append(lmsd_met + hmdb_met)

            # create Annotation objects
            for f, metabs in izip(features, metabolites):
                if not metabs:
                    not_found += 1
                else:
                    m_count += len(metabs)
                for_adducts_str = '[M{}{}]='.format('-' if f.polarity > 0 else '+', formula)
                f.annotations += [Annotation(m, for_adducts_str) for m in metabs]

            pool.close()
            try:
                pool.terminate()
            except OSError:
                pass
        return m_count, not_found
