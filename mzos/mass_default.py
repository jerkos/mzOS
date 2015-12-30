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


def is_good_formula(formula, allowed_atoms=frozenset({"C", "H", "N", "O", "P", "S"})):
    """
    :param formula:
    :param allowed_atoms:
    :return:
    """
    atoms = [x[0] for x in ELEMENT_PATTERN.findall(formula)]
    if not all(a in allowed_atoms for a in atoms):
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
            fh.write("{0}\t{1}\n".format(str(m[0]), str(m[1])))


