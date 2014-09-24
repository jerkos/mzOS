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

from xml.etree import cElementTree
import sqlite3
import glob
import subprocess
import logging
import multiprocessing
import time
import re

ELEMENT_PATTERN = re.compile(r'''
            ([A-Z][a-z]{0,2})
            ([\-]?[\d]*)
''', re.X)


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
    c.execute('create table metabolite (acession text primary key, name text, formula text, inchi text, '
              'mono_mass double, average_mass double, description text, status text, origin text, kegg_id text, '
              'isotopic_pattern_pos text, isotopic_pattern_neg text)')
    return conn, c


def parse_metabolite_card(filepath):
    """
    @param filepath:  str
    @return:
    """
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
                try:
                    metab.append(float(elem.find('monisotopic_moleculate_weight').text))
                except (ValueError, TypeError):
                    logging.info("Failed to parse monisotopic_moleculate_weight in {}".format(filepath))
                    metab.append(0.0)
                try:
                    metab.append(float(elem.find('average_molecular_weight').text))
                except (ValueError, TypeError):
                    logging.info("Failed to parse average_molecular_weight in {}".format(filepath))
                    metab.append(0.0)
                metab.append(elem.find('description').text)
                #float(elem.find('description').text)
                metab.append("")
                #float(elem.find('description').text)
                metab.append("")
                #float(elem.find('description').text)
                metab.append(elem.find('kegg_id').text)
                metab.append(get_theo_ip(f, min_rel_int=1.0, polarity=1))  # isotopic pattern elem.find('kegg_id').text
                metab.append(get_theo_ip(f, min_rel_int=1.0, polarity=-1))
    except Exception as e:
        logging.warn("Error parsing metabolite card : {} with following exception : \n {}".format(filepath, e.message))
        return None
    return tuple(metab)


def build_library(sqlite_filepath, cards_directory, nb_procs=multiprocessing.cpu_count()):
    """
    :param sqlite_filepath:
    :param cards_directory:
    :param nb_procs:
    """
    files = glob.glob(cards_directory + "\\*.xml")
    if not files:
        print("No metabolites card found")
    t1 = time.clock()
    conn, c = build_sqlite3_file(sqlite_filepath)
    p = multiprocessing.Pool(processes=nb_procs)
    metabs = p.map(parse_metabolite_card, files, chunksize=20)
    c.executemany('insert into metabolite values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                  [m for m in metabs if m is not None])
    conn.commit()
    c.execute('create index mass_index on metabolite(mono_mass)')
    conn.close()
    logging.info("Finished, elpased time {} seconds".format(time.clock() - t1))


def add_element(f, element, n):
    """
    :param f: dict key element, value number of element
    :param n: add or remove n element
    :return:
    """
    assert(isinstance(element, str))
    if element in f:
        nb_h = int(f[element]) if f[element] != '' else 1
        f[element] = str(nb_h + n)
    else:
        f[element] = str(n)
    return "".join(["".join((k, v)) for k, v in sorted(f.iteritems(), key=lambda _: _[0])])


def remove_element(f, element, n):
    """
    :param f:
    :param element:
    :param n:
    :return:
    """
    assert(isinstance(element, str))
    if element in f:
        nb_h = int(f[element]) if f[element] else 1
        if nb_h - n <= 0:
            del f[element]
        elif nb_h - n == 1:
            f[element] = ''
        else:
            f[element] = str(nb_h - n)
    return "".join(["".join((k, v)) for k, v in sorted(f.iteritems(), key=lambda _: _[0])])


def get_theo_ip(formula, min_rel_int=5.0, polarity=1):
    """
    # todo ask wich adducts to pass in parameter
    formula is a string meaning compound
    :param formula:
    :param min_rel_int:
    :param polarity:
    """
    if not isinstance(formula, str):
        raise Exception("[generate theoritical isotopic pattern]"
                        " formula parameter must be a string:%s" % repr(formula))
    if not formula:
        logging.warning("Formula seems to be empty.")
        return

    f = {x[0]: x[1] for x in ELEMENT_PATTERN.findall(formula)}
    if polarity == 1:
        formula = add_element(f, 'H', 1)
    else:
        formula = remove_element(f, 'H', 1)

    p = subprocess.Popen("emass/emass.exe", stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = p.communicate(input=formula)
    if not out:
        logging.warn("Error computing isotopic pattern with formula: {0}.Skip it".format(formula))
        return

    iso = ""
    try:
        iso = repr(filter(lambda x: x[1] > min_rel_int,
                          [(lambda x: (float(x[0]), float(x[1])))(l.rstrip().split(" "))
                           for l in out.split('\n')[1:-1]]))
    except Exception:
        logging.warn("Error parsing isotopic pattern.Skip it")
        return

    return iso

if __name__ == '__main__':
    import sys
    if not len(sys.argv) > 1:
        logging.error("""Usage: hmdb_sqlite_creator path/to/hmdb/directory""")
    hmdb_cards_dir = sys.argv[1]
    logging.basicConfig(level=logging.INFO)
    build_library('hmdb.sqlite', 'hmdb_metabolites')
