from __future__ import absolute_import
from __future__ import print_function
import sys
import os.path as op
from xml.etree import cElementTree
import sqlite3
import glob
import logging
import multiprocessing
import time

try:
    from mzos.formula import Formula
    from mzos.utils import get_theo_ip
except ImportError:
    logging.warn('print called hmdb creator from main')
    abspath = op.abspath(op.normcase('../../'))
    sys.path.append(abspath)
    from mzos.formula import Formula
    from mzos.utils import get_theo_ip


def build_sqlite3_file(sqlite_filepath):
    """
    :param sqlite_filepath:
    :return:
    """
    conn = sqlite3.connect(sqlite_filepath)
    c = conn.cursor()
    try:
        c.execute('drop table metabolite;')
    except sqlite3.OperationalError:
        pass
    c.execute('create table metabolite (acession text primary key, name text, formula text, inchi text, inchikey text,'
              'mono_mass double, average_mass double, description text, status text, origin text, kegg_id text, '
              'isotopic_pattern_pos text, isotopic_pattern_neg text)')
    return conn, c


def parse_metabolite_card(args):
    """
    @return:
    """
    filepath = args[0]
    emass_path = args[1]
    context = cElementTree.iterparse(filepath, events=('end',))
    metab = list()
    try:
        for action, elem in context:
            if elem.tag == 'metabolite' and action == 'end':
                f = elem.find('chemical_formula').text
                metab.append(elem.find('accession').text)
                metab.append(elem.find('name').text)
                metab.append(f)
                metab.append(elem.find('inchi').text)
                metab.append(elem.find('inchikey').text[9:])
                try:
                    metab.append(float(elem.find('monisotopic_moleculate_weight').text))
                except (ValueError, TypeError):
                    logging.info("Failed to parse monisotopic_moleculate_weight in {0}".format(filepath))
                    metab.append(0.0)
                try:
                    metab.append(float(elem.find('average_molecular_weight').text))
                except (ValueError, TypeError):
                    logging.info("Failed to parse average_molecular_weight in {0}".format(filepath))
                    metab.append(0.0)
                metab.append(elem.find('description').text)
                # float(elem.find('description').text)
                metab.append("")
                # float(elem.find('description').text)
                metab.append("")
                # float(elem.find('description').text)
                # isotopic pattern elem.find('kegg_id').text
                metab.append(elem.find('kegg_id').text)

                formula = Formula.from_str(f)
                # add one H to have positive
                f1 = formula.add('H', new_obj=True)
                metab.append(get_theo_ip(emass_path, str(f1), min_rel_int=1.0))
                #
                f2 = formula.remove('H', new_obj=True)
                metab.append(get_theo_ip(emass_path, str(f2), min_rel_int=1.0))
    except (WindowsError, IOError, ValueError, TypeError) as e:
        logging.warn("Error parsing metabolite card : {0} with following exception : \n {1}".format(filepath, e))
        return None
    return tuple(metab)


def build_library(sqlite_filepath, cards_directory, emass_path, nb_procs=multiprocessing.cpu_count()):
    """
    :param sqlite_filepath:
    :param cards_directory:
    :param nb_procs:
    """
    files = glob.glob(cards_directory + "\\*.xml")
    if not files:
        print("No metabolites card found")
    files = [(f, emass_path) for f in files]
    t1 = time.clock()
    conn, c = build_sqlite3_file(sqlite_filepath)
    p = multiprocessing.Pool(processes=nb_procs)
    metabs = p.map(parse_metabolite_card, files, chunksize=20)
    p.close()
    c.executemany('insert into metabolite values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                  [m for m in metabs if m is not None])
    conn.commit()
    c.execute('create index mass_index on metabolite(mono_mass)')
    c.execute('create index inchi_key_idx on metabolite(inchikey)')
    conn.close()
    logging.info("Finished, elpased time {0} seconds".format(time.clock() - t1))


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    # from mzos.formula import Formula
    build_library('hmdb_test.sqlite',
                  op.normcase('../ressources/hmdb_metabolites'),
                  op.normcase("../../third_party/emass/emass.exe -i ../../third_party/emass/ISOTOPE.DAT"))
