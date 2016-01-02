from __future__ import absolute_import
import subprocess
import logging
import os.path as op

EMASS_PATH = op.abspath('mzos/third_party/emass/emass -i mzos/third_party/emass/ISOTOPE.DAT')


def get_theo_ip(f, min_rel_int=5.0):
    """
    # todo ask wich adducts to pass in parameter
    formula is a string meaning compound
    :param f:
    :param min_rel_int:
    """
    p = subprocess.Popen(EMASS_PATH, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    out, err = p.communicate(input=f)
    if not out:
        logging.warn("Error computing isotopic pattern with formula: {0}.Skip it".format(f))
        return

    try:
        iso = repr([x for x in [(lambda x: (float(x[0]), float(x[1])))(l.rstrip().split(" "))
                           for l in out.split('\n')[1:-1]] if x[1] > min_rel_int])
    except IndexError:
        logging.warn("Error parsing isotopic pattern.Skip it")
        return

    return iso


def calculate_mass_diff_da(feature, moz_metabolite, include_isotopes=False):
    if not include_isotopes:
        return abs(feature.get_real_mass() - moz_metabolite)
    # TODO kind of rmsd but on masses ?
    return 0
