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

from math import sqrt
from collections import defaultdict as ddict
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


# def calculate_rmsd(feature, theo_ip, moz_tol_ppm):
#         """
#         caculate the root mean square deviation between observed and
#         theortical isotopic pattern
#         @param feature:
#         @param theo_ip:
#         @return: root mean square deviation
#         """
#         def get_nearest_peakel(mz, mz_tol_ppm, peaks):
#             p = peaks.sort(key=lambda x: abs(x.moz - moz))[0]
#             if abs(p.moz - mz) < mz * mz_tol_ppm / 1e6:
#                 return p
#             return None
#
#         if not feature.isotopes:
#             return float('nan')
#
#         max_rel_int = sorted(theo_ip, key=lambda x: x[1])[-1][1]
#         max_real_int = sorted(feature.get_isotopic_pattern_as_peakel(), key=lambda x: x.area)[-1].area
#         rmsd = 0
#         for idx, (moz, rel_int) in enumerate(theo_ip):
#             p = get_nearest_peakel(moz, moz_tol_ppm, feature.get_isotopic_pattern_as_peakel())
#             if p is not None:
#                 rmsd += (p.area * max_rel_int / max_real_int) ** 2
#             else:
#                 #could do something like the first quartile of the distribution of all intensities
#                 pass
#         return sqrt(rmsd)


def calculate_mass_diff_da(feature, moz_metabolite, include_isotopes=False):
    if not include_isotopes:
        return abs(feature.get_real_mass() - moz_metabolite)
    #TODO kind of rmsd but on masses ?
    return 0

#
#
# def merge(d1, d2):
#     d3 = {}
#     for f in d1.keys():
#         m1 = d1[f]
#         m2 = d2.get(f, None)
#         if m2 is None:
#             raise ValueError("feature not found in second score")
#         # if m2 is None:
#         #     d3[f] = []
#         #     for m, s1 in m1:
#         #         d3[f].append((m, s1, float('nan')))
#         #     continue
#         d3[f] = []
#         for (m, s1) in m1:
#             for m_, s2 in m2:
#                 if m == m_:
#                     d3[f].append((m, s1, s2))
#                     break
#     return d3
#
#
# def inverse_dict(d):
#     new_d = ddict(set)
#     for mo, sons in d.items():
#         for s in sons:
#             new_d[s].append(mo)
#     return new_d
#
#
# def group_by(l, callable_):
#     if not isinstance(callable_, callable):
#         raise TypeError("The second argument must be a callable")
#     d = dict(list)
#     for c in l:
#         d[callable_(c)].append(c)
#     return d