__author__ = 'Marc'

import re
import subprocess
import logging
import os.path as op
from collections import Counter


class Formula(dict):
    """Modelcular formula"""

    ELEMENT_PATTERN = re.compile(r'''
            ([A-Z][a-z]{0,2})
            ([\-]?[\d]*)
    ''', re.X)

    EMASS_PATH = op.normcase('third_party/emass/emass.exe -i third_party/emass/ISOTOPE.DAT')

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)

    @staticmethod
    def from_str(fstr):
        """ utility function, create formula from a string
        :param fstr:
        """
        assert(isinstance(fstr, str))
        v = Formula.ELEMENT_PATTERN.findall(fstr)
        if not v:
            logging.warn("provided string does not match a molecular formula")
            return None

        return Formula({elem: (int(nb) if nb else 1) for elem, nb in v})

    def __str__(self):
        return "".join(["".join((k, (str(v) if v > 1 else '')))
                        for k, v in sorted(self.iteritems(), key=lambda _: _[0])])

    @staticmethod
    def _check_input(f):
        """
        :param f:
        :return:
        """
        if isinstance(f, str):
            fd = Formula.from_str(f)
        elif isinstance(f, dict) or isinstance(f, Counter):
            fd = Formula(f)
        else:
            if not isinstance(f, Formula):
                raise TypeError('provided argument must be str, dict, or formula')
            fd = f
        return fd

    def add(self, f, new_obj=False):
        """
        :param new_obj:
        :param f: dict key element, value number of element
        :param n: add or remove n element
        :return:
        """
        fd = self._check_input(f)
        #set the working formula obj
        #shallow copy if needed
        wd = Formula(self) if new_obj else self

        for elem, nb in fd.iteritems():
            wd[elem] = wd.get(elem, 0) + nb
        return wd

    def remove(self, f, new_obj=False):
        """

        :param new_obj:
        :param f:
        :param element:
        :param n:
        :return:
        """
        fd = self._check_input(f)

        wd = Formula(self) if new_obj else self

        for elem, nb in fd.iteritems():
            if elem in wd:
                nb_e = wd[elem] - nb
                if nb_e <= 0:
                    del wd[elem]
                else:
                    wd[elem] = nb_e
        return wd

    # def get_theo_ip(self, min_rel_int=5.0):
    #     """
    #     # todo ask wich adducts to pass in parameter
    #     formula is a string meaning compound
    #     :param min_rel_int:
    #     """
    #     p = subprocess.Popen(Formula.EMASS_PATH, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    #
    #     out, err = p.communicate(input=str(self))
    #     if not out:
    #         logging.warn("Error computing isotopic pattern with formula: {0}.Skip it".format(str(self)))
    #         return
    #
    #     try:
    #         iso = repr(filter(lambda x: x[1] > min_rel_int,
    #                           [(lambda x: (float(x[0]), float(x[1])))(l.rstrip().split(" "))
    #                            for l in out.split('\n')[1:-1]]))
    #     except IndexError:
    #         logging.warn("Error parsing isotopic pattern.Skip it")
    #         return
    #
    #     return iso
