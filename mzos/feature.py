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

from itertools import count
from itertools import groupby
from collections import defaultdict as ddict
import logging
import numpy as np


class Peak(object):
    """mz peak object"""
    def __init__(self, mz, intensity):
        self.moz = mz
        self.intensity = intensity


class Attribution(object):
    """
    Object corresponding to one possible attribution
    An attribution corresponds to `kind of peak`, i.e
    one of the following tag: isotope, adduct, or monoisotope.
    tag corresponds to the peak kind identification
    """
    def __init__(self, attribution, parent_id, charge):

        assert(isinstance(attribution, str))
        self.attribution = attribution
        self.parent_id = parent_id
        self.charge = charge

    def __str__(self):
        return "{} of {} for charge={}".format(self.attribution, self.parent_id, self.charge)


class Annotation(object):
    """
    Annotation is the result of a matching metabolite
    to one detected feature.
    """
    ADDUCTS = {'+Na': '[M+Na]',
               '+H': '[M+H]',
               '-H': '[M-H]',
               '-Cl': '[M-Cl]'
               }

    def __init__(self,
                 metabolite,
                 for_adduct,
                 score_isos=0.0,
                 score_network=0.0):

        self.metabolite = metabolite
        self.score_isos = score_isos
        self.score_network = score_network
        self.for_adduct = for_adduct

    def has_br_atom(self):
        """could be introduced in the model """
        return 'Br' in self.metabolite.formula

    def has_s_atom(self):
        """could be introduced in the model """
        return 'S' in self.metabolite.formula


class Peakel(object):
    """
    peakel or elution peak
    """

    ADDUCTS_MASS = {'H': 1.007276}

    _ids = count(1)

    def __init__(self, moz, mozmin, mozmax, rt, rtmin=0, rtmax=0):
        """
        moz :calculated moz
        rt : rt of the centroid
        area : raw AUC
        rtmin : minRt observed
        rtmax : maxrt               
        """
        # used in __hash__
        self.id = Peakel._ids.next()

        # all seen areas
        self.area_by_sample_name = {}
        self.area = self.get_median_area()

        # not set yet
        self.polarity = 0

        self.moz = moz
        self.mozmin = mozmin
        self.mozmax = mozmax
        self.rt = rt
        self.rtmin = rtmin
        self.rtmax = rtmax

        # by default
        self.charge = 1

        # could be in a feature object
        self.isotopes = set()
        self.adducts = set()

        # chemical formula
        self.main_annotation = ""
        self.annotations = []

        # None here means monoisotopic
        self.main_attribution = None
        self.attributions = set()

        self.ip_score_isotopes = set()

        # in concordance with the main tag
        self.is_isotope = False
        self.is_adduct_or_fragment = False

        ##
        self.peaks = []

    def get_metabolites(self):
        """return metabolites from annotations"""
        return [a.metabolite for a in self.annotations]

    def _build_branch_str_v2(self, parent, charge, son, isos_adds):
        """
        @param parent:
        @param son:
        @param isos_adds:
        @return:
        """
        s = ""
        if son is not self:
            attrib = ""
            try:
                attrib += son.get_attributions_by_parent_id()[parent.id][0].attribution
            except KeyError:  # see index error
                attrib += son.main_attribution.attribution

            s += str(son.id) + "=" + attrib

        isos = set([si for si in son.isotopes if si.get_attributions_by_parent_id()[son.id][0].charge == charge])
        n_ = isos.union(son.adducts)
        nb_isos, nb_adducts = len(isos), len(son.adducts)

        if n_:
            if self is not son:
                s += ":"
            for new_son in n_:
                s += "("
                # Peakel._build_branch_str(son, charge, new_son, isos_adds.union(n_))
                m, q, r = son._build_branch_str_v2(son, charge, new_son, isos_adds.union(n_))
                s += m
                nb_isos += q
                nb_adducts += r
                s += ");"
        return s, nb_isos, nb_adducts

    def get_top_down_attribution_tree(self):
        """
        :return:
        """
        output, n_isos, n_adducts = self._build_branch_str_v2(None,
                                                              self.charge,
                                                              self,
                                                              self.isotopes.union(self.adducts))
        output = output if not output.endswith(";") else output[:-1]
        output = output.replace(";)", ")")
        output = output.replace(";;", ";")
        return output, n_isos, n_adducts

    def get_bottom_up_attribution_tree(self, feature_by_id):
        """
        @param feature_by_id: dictionnary key:peakel, value: id
        @return:
        """
        # "{} of {} ".format(self.main_attribution.attribution, self.main_attribution.parent_id)
        s = str(self.main_attribution)
        p = feature_by_id[self.main_attribution.parent_id]
        while p.main_attribution is not None:
            # s += "{} of {} ".format(p.main_attribution.attribution, p.main_attribution.parent_id)
            s += " of " + str(p.main_attribution)
            p = feature_by_id[p.main_attribution.parent_id]
        return s

    @staticmethod
    def get_others_bottom_up_attribution_tree(attribution, feature_by_id):
        """
        @param attribution:
        @param feature_by_id: dictionnary key:peakel, value: id
        @return:
        """
        # s = "{} of {} ".format(attribution.attribution, attribution.parent_id)
        s = str(attribution)
        p = feature_by_id[attribution.parent_id]
        while p.main_attribution is not None:
            # s += "{} of {} ".format(p.main_attribution.attribution, p.main_attribution.parent_id)
            s += " of " + str(p.main_attribution)
            p = feature_by_id[p.main_attribution.parent_id]
        return s

    # def fill_peaks(self):
    #     """
    #     mzdb provided ?
    #     """
    #     pass

    # ----------------attributions stuffs
    def get_attributions_by(self, callable_):
        """
        @param callable_:
        @return:
        """
        return {k: list(v) for k, v in groupby(self.attributions, callable_)}

    def get_attributions_by_parent_id(self):
        """
        :return:
        """
        return {k: list(v) for k, v in groupby(self.attributions, lambda x: x.parent_id)}

    def get_attributions_by_charge(self):
        """
        :return:
        """
        return {k: list(v) for k, v in groupby(self.attributions, lambda x: x.charge)}

    def remove_attribution_with_parent(self, parent_id):
        """
        @param parent_id:
        @return:
        """
        attr = None
        for a in self.attributions:
            if a.parent_id == parent_id:
                attr = a
                break
        if attr is not None:
            self.attributions.remove(attr)
            if self.main_attribution == attr:
                self.main_attribution = None
                if self.attributions:
                    self.main_attribution = list(self.attributions)[0]  # the next one ?

    def add_attribution(self, attrib):
        """
        @param attrib:
        @return:
        """
        self.attributions.add(attrib)
        if self.main_attribution is None:
            self.main_attribution = attrib

    def set_main_attribution(self, attrib):
        """
        @param attrib:
        @return:
        """
        # TODO add or not main attribution to attributions set
        # in that case remove unecessary code
        m = None
        # save previous attribution
        if self.main_attribution is not None:
            m = self.main_attribution

        self.main_attribution = attrib
        # FIX
        # add the new main attribution to all attributions set
        self.attributions.add(attrib)
        if m is not None and m not in self.attributions:
            self.attributions.add(m)

    def get_areas(self):
        """
        :return:
        """
        return self.area_by_sample_name.values()

    def get_median_area(self):
        """
        :return:
        """
        return np.median(self.get_areas())

    def corr_intensity_against(self, peakel):
        """
        @param peakel:
        @return:
        """
        values = [peakel.area_by_sample_name[k] for k in self.area_by_sample_name.keys()]
        return np.corrcoef(self.area_by_sample_name.values(), values)[1, 0]

    def corr_shape_against(self, peakel):
        """
        :param peakel:
        :return:
        """
        pass

    def get_isotopic_pattern(self):
        """
        :return:
        """
        isos = [(iso.moz, iso.area) for iso in self.isotopes]
        isos.append((self.moz, self.area))
        isos.sort(key=lambda x: x[1])
        return isos

    def get_isotopic_pattern_as_peakel(self):
        """
        :return:
        """
        isos = list(self.isotopes)[:]
        isos.insert(0, self)
        return isos

    def get_real_mass(self, adducts=None):
        m = self.moz * self.charge
        charge_mass = self.charge * (adducts or Peakel.ADDUCTS_MASS['H'])
        return m + charge_mass if self.polarity < 0 else m - charge_mass


class Feature(object):
    """
    :param mono_mz:
    :param rt:
    :param peakels:
    """
    def __init__(self, mono_mz, rt, peakels=None):
        self.moz = mono_mz
        self.rt = rt
        self.isotopes = peakels or []
        self.adducts = set()


class PeakelIndex(object):
    """
    :param peakels:
    :param scan_ids:
    :param bin_size:
    """
    mozgetter = staticmethod(lambda x: x.moz)

    def __init__(self, peakels, scan_ids=None, bin_size=1):
        """
        peakels: list of peakels taken from the output of xcms
        scanid
        bin_size : size of one bin in mz dimension
        """
        self.scan_ids = scan_ids
        self._index = ddict(list)
        self.inv_bin_size = 1.0 / bin_size

        self.sorted_peaks = sorted(peakels, key=lambda x: x.moz)
        self.min_moz, self.max_moz = self.sorted_peaks[0].moz, self.sorted_peaks[-1].moz

        for p in peakels:
            bin_ = int(p.moz * self.inv_bin_size)
            self._index[bin_].append(p)

    def empty(self):
        """
        check emptyness of the index        
        """
        return len(self._index)

    def get_nearest_peakel(self, moz, mz_tol_ppm):
        """
        :param moz:
        :param mz_tol_ppm:
        could return a None value
        """
        if moz < self.min_moz or moz > self.max_moz:
            return None

        bin_ = int(moz * self.inv_bin_size)
        tol_da = (mz_tol_ppm * moz) / 1e6
        peaks = []
        for i in xrange(bin_ - 1, bin_ + 2):
            peaks.extend(self._index.get(i, []))
        if not peaks:
            return None
        peaks.sort(key=lambda x: abs(x.moz - moz))

        return peaks[0] if (abs(peaks[0].moz - moz) < tol_da) else None
