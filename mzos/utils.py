import subprocess
import logging


def get_theo_ip(emass_path, f, min_rel_int=5.0):
    """
    # todo ask wich adducts to pass in parameter
    formula is a string meaning compound
    :param f:
    :param emass_path:
    :param min_rel_int:
    """
    p = subprocess.Popen(emass_path, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    out, err = p.communicate(input=f)
    if not out:
        logging.warn("Error computing isotopic pattern with formula: {0}.Skip it".format(f))
        return

    try:
        iso = repr(filter(lambda x: x[1] > min_rel_int,
                          [(lambda x: (float(x[0]), float(x[1])))(l.rstrip().split(" "))
                           for l in out.split('\n')[1:-1]]))
    except IndexError:
        logging.warn("Error parsing isotopic pattern.Skip it")
        return

    return iso


def calculate_mass_diff_da(feature, moz_metabolite, include_isotopes=False):
    if not include_isotopes:
        return abs(feature.get_real_mass() - moz_metabolite)
    # TODO kind of rmsd but on masses ?
    return 0
