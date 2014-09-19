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


from cmath import isnan
from math import sqrt, log10
from feature import PeakelIndex
from utils import calculate_mass_diff_da
from collections import defaultdict as ddict
import numpy as np

mass_H = 1.00794
mass_electron = 0  # negligeable 9E-31


class StatsModel(object):
    """
    assign score to metabolites that match in mz to specified feature
    """

    @staticmethod
    def model(x):
        """
        @param x:
        @return:
        """
        return x ** 2

    @staticmethod
    def transform_score(x):
        """
        @param x:
        @return:
        """
        return -10 * log10(x)

    #name metrics, weight allowed
    def __init__(self, features, moz_tol_ppm):
        self.features = features
        self.moz_tol_ppm = moz_tol_ppm
        self.metrics = [("isotopic_pattern_rmsd", 2.0), ("mass_difference", 1.0)]

        #sorted_features_by_area = sorted(features, key=lambda _: _.area)

        feature_by_index = ddict(list)
        #assuming they are sorted by mass
        for f in self.features:
            feature_by_index[1].append(f)
            for (i, isotope) in enumerate(f.isotopes):
                feature_by_index[i + 2].append(isotope)
        #retrive min max intensity value by intensity index
        self.min_max_values_by_isotopes_index = dict()

        for (index, values) in feature_by_index.iteritems():
            sorted_values = sorted(values, key=lambda __: __.area)
            self.min_max_values_by_isotopes_index[index] = (sorted_values[0].area, sorted_values[-1].area)

    def _calculate_worst_cases(self, feature, theoritical_isotopes):
        """
        @param: feature, the feature to evaluate
        @param: theortical isotopic pattern
        """
        peakel_index = PeakelIndex(feature.get_isotopic_pattern_as_peakel())
        rmsd, mass_diff = 0, 0
        for idx, iso in enumerate(theoritical_isotopes):
            # find isotopes
            p = peakel_index.get_nearest_peakel(iso[0], self.moz_tol_ppm)
            if p is not None:
                try:
                    #TODO why get a key error here ?
                    mn, mx = self.min_max_values_by_isotopes_index[idx + 1]  #we start counting at one
                    rmsd += (mx - mn) ** 2
                except KeyError:
                    pass

        return sqrt(rmsd), feature.moz * self.moz_tol_ppm / 1e6, peakel_index

    def calculate_metabolites_score(self, feature):
        """
        feature: peakel instance with several isotopes
        metabolites: list of metabolites
        """
        for annot in feature.annotations:

            m = annot.metabolite
            # worst case
            ip = m.isotopic_pattern_pos if feature.polarity == 1 else m.isotopic_pattern_neg
            isotopic_pattern = [(float(a), float(b)) for a, b in eval(ip)]
            worst_rmsd, worst_mass_diff, peakel_index = self._calculate_worst_cases(feature, isotopic_pattern)

            # interpol_worst_rmsd, interpol_worst_mass_diff = 1.0, 1.0
            #as we interpolate al line y = x we use directly the result and pass
            #it to the model
            mass_diff = calculate_mass_diff_da(feature, m.mono_mass)
            interpol_mass_diff = mass_diff / worst_mass_diff
            ponderated_mass_diff = self.transform_score(self.model(interpol_mass_diff)) * self.metrics[1][1]

            rmsd = self._calculate_rmsd_2(feature, peakel_index, isotopic_pattern)
            if isnan(rmsd) or rmsd == 0.0 or worst_rmsd == 0.0:
                #metab_with_score.append((m, ponderated_mass_diff))
                annot.score_isos = ponderated_mass_diff
                continue

            interpol_rmsd = rmsd / worst_rmsd
            ponderated_rmsd = self.transform_score(self.model(interpol_rmsd)) * self.metrics[0][1]

            final_score = (ponderated_mass_diff + ponderated_rmsd) / (self.metrics[0][1] + self.metrics[1][1])

            annot.score_isos = final_score

    def calculate_score(self):
        """
        todo use multiprocessing
        :return: None
        """
        for f in self.features:
            self.calculate_metabolites_score(f)

    def _calculate_rmsd_2(self, feature, peakel_idx, theo_ip, method="mean"):
        """

        @param feature:
        @param peakel_idx:
        @param theo_ip:
        @param method: str
        @return:
        """
        if not feature.isotopes:
            return float('nan')
        max_rel_int = max(theo_ip, key=lambda x: x[1])[1]
        isotopic_pattern = feature.get_isotopic_pattern_as_peakel()

        sample_rmsd = []
        for sample in feature.area_by_sample_name.keys():
            rmsd = 0.0
            max_real_int = max(isotopic_pattern,
                               key=lambda x: x.area_by_sample_name[sample]).area_by_sample_name[sample]

            if not max_real_int:
                continue

            for idx, (moz, rel_int) in enumerate(theo_ip):
                p = peakel_idx.get_nearest_peakel(moz, self.moz_tol_ppm)
                if p is not None:
                    feature.ip_score_isotopes.add(p)
                    area = p.area_by_sample_name[sample]
                    ## fixme: could be penalized ?
                    if not area:
                        continue
                    #rmsd += ((area * max_rel_int / max_real_int) - rel_int) ** 2
                    rmsd += ((area / max_real_int * 100) * max_rel_int - rel_int) ** 2
                else:
                    #could do something like the first quartile of the distribution of all intensities
                    pass

            sample_rmsd.append(sqrt(rmsd))

        return np.mean(sample_rmsd) if method == "mean" else np.median(sample_rmsd)

    def _calculate_rmsd(self, feature, peakel_idx, theo_ip):
        """
        @param feature: Peakel
        caculate the root mean square deviation between observed and
        theortical isotopic pattern
        """
        if not feature.isotopes:
            return float('nan')

        max_rel_int = sorted(theo_ip, key=lambda x: x[1])[-1][1]

        max_real_int = sorted(feature.get_isotopic_pattern_as_peakel(), key=lambda x: x.area)[-1].area
        rmsd = 0
        for idx, (moz, rel_int) in enumerate(theo_ip):
            p = peakel_idx.get_nearest_peakel(moz, self.moz_tol_ppm)
            if p is not None:
                rmsd += (p.area * max_rel_int / max_real_int) ** 2
                feature.ip_score_isotopes.add(p)
            else:
                #could do something like the first quartile of the distribution of all intensities
                pass
        return sqrt(rmsd)





