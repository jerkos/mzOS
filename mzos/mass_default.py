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

import os.path as op
import sqlite3
import re
import math

ELEMENT_PATTERN = re.compile(r'''
            ([A-Z][a-z]{0,2})
            ([\-]?[\d]*)
''', re.X)

db = op.normcase("ressources/hmdb.sqlite")
request = "select formula, mono_mass from metabolite"


def is_good_formula(formula, allowed_atoms={"C", "H", "N", "O", "P", "S"}):
    """
    :param formula:
    :param allowed_atoms:
    :return:
    """
    atoms = [x[0] for x in ELEMENT_PATTERN.findall(formula)]
    if not all([a in allowed_atoms for a in atoms]):
        return False
    return True

if __name__ == "__main__":
    c = sqlite3.connect(db)
    cursor = c.cursor()
    metabolites = []
    for row in cursor.execute(request):
        if is_good_formula(row[0]):
            mass = float(row[1])
            metabolites.append((mass, mass - math.floor(mass)))
    with open("ressources/default_mass_C_H_N_O_P_S.tsv", 'w') as fh:
        for m in metabolites:
            fh.write("{}\t{}\n".format(str(m[0]), str(m[1])))


