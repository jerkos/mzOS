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

import subprocess
import logging


def get_theo_ip(emass_path, f, min_rel_int=5.0):
    """
    # todo ask wich adducts to pass in parameter
    formula is a string meaning compound
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