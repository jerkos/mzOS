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

import logging
from collections import defaultdict as ddict
from itertools import chain

from mzos.feature import PeakelIndex, Attribution
from mzos.peakel_clusterer import PeakelClusterer



class PeakelsAnnotator(object):
    """

    Main function to annotates elution peak
    =======================================

    Try to distinguish isotopic pattern from extracted mass trace.
    Basically, we use for now XCMS as a backend. It uses a SVM model
    taken from FeatureFidnerMetabo developped in OpenMS

    We also provide a method for finding adducts and fragments

    """

    # isotopes to look for
    ISOTOPES = {1.003355: "Isotope C13",
                0.997035: "Isotope N15",
                1.995796: "Isotope S34"}  # 1.997953: "Br79-81"}

    def __init__(self, peakels, exp_settings):
        self.peakels = peakels
        self.exp_settings = exp_settings
        self.index = PeakelIndex(peakels)

        self.peakel_clusterer = None
        self.adducts_or_fragments = self.exp_settings.get_mass_to_check()

    def set_peakels(self, peakels):
        """
        :param peakels:
        :return:
        """
        self.peakels = peakels

    def get_nearest_peakel(self, moz, mz_tol_ppm):
        """
        proxy
        :param mz_tol_ppm:
        :param moz:
        """
        return self.index.get_nearest_peakel(moz, mz_tol_ppm)

    @staticmethod
    def _get_theoritical_isotope_mass(idx, mz, charge):
        """
        mean_theo = 1.000857 * j + 0.001091, sigma_theo = 0.0016633 * j - 0.0004571

        taken form the model FeatureFinderMetabo openms

        :param idx:
        :param mz:
        :param charge:
        :return:
        """
        theo_mass = mz + ((1.000857 * idx + 0.001091) / charge)
        d = {x[1]: abs(theo_mass - (mz + x[0])) for x in PeakelsAnnotator.ISOTOPES.items()}
        return theo_mass, min(d)

    def _look_for_isotopes(self,
                           peakel,
                           mos_by_iso,
                           max_charge,
                           max_gap,
                           max_isotopes_nb,
                           error_rt):
        """
        :param peakel: Peakel object peakel to consider
        :param mos_by_iso: dict key:iso, value:set
        :param max_charge: int max charge to check
        :param max_isotopes_nb: int
        :param error_rt: float
        :return: dict
        """
        moz_tol_ppm = self.exp_settings.mz_tol_ppm

        half_error_rt = error_rt * 0.5

        # will hold set of isotopes
        result_by_charge = {}

        # iterate over possible charges
        for charge in xrange(1, max_charge + 1):

            isotopes = set()
            gap = 0

            # allow uniquely intensity desc
            desc_iso_intensity = False

            # assign last iso to considered peakel
            last_iso = peakel

            # assume that supposed monoisotopic elution peak has the highest intensity,
            # but may not be true in reality. The measured mass over charge is known to be
            # more accurate with higher intensities
            ref_moz = peakel.moz

            # iterate over  number of isotopes
            for j in xrange(1, max_isotopes_nb + 1):
                # generate all possible masses for this isotopes index
                mass_to_check, isotope_tag = PeakelsAnnotator._get_theoritical_isotope_mass(j, ref_moz, charge)

                # iterate over those generated masses
                peak = self.get_nearest_peakel(mass_to_check, moz_tol_ppm)

                # error rt check even we are inside a rt cluster ? pass it into a flag
                if peak is not None and abs(peak.rt - peakel.rt) < half_error_rt:

                    if desc_iso_intensity and peak.area > last_iso.area:
                        break

                    isotopes.add(peak)

                    # add it as a possible parent
                    mos_by_iso[peak].add(peakel)

                    # add a new tag
                    peak.attributions.add(Attribution(isotope_tag, peakel.id, charge))

                    if peak.area > last_iso.area:
                        desc_iso_intensity = True
                        ref_moz = peak.moz

                    last_iso = peak
                # if no peaks has been found this isotope index increase the gap
                else:
                    gap += 1
                    if gap >= max_gap:
                        break

            # if isotopes not empty save it in the dictionnary
            if isotopes:
                result_by_charge[charge] = isotopes
        return result_by_charge

    @staticmethod
    def _remove_peakel_from_wrong_parents(peakel, wrong_parents, mos_by_iso, isotopes_clustered):
        """
        allow to remove an ambiguity, i.e. removing a peakel from parent' s isotopes.

        :param peakel:
        :param wrong_parents:
        :param mos_by_iso:
        :param isotopes_clustered:
        :return:
        """
        for wrong_parent in wrong_parents:
            ordered_isos = sorted(list(wrong_parent.isotopes), key=lambda _: _.moz)
            if peakel in ordered_isos:
                # remove its tag
                peakel.remove_attribution_with_parent(wrong_parent.id)

                peakel_index = ordered_isos.index(peakel)
                for p in ordered_isos[peakel_index + 1:]:
                    p.remove_attribution_with_parent(wrong_parent.id)

                    # if this peakel has only one parent, can safely remove it from isotopes set
                    if len(mos_by_iso[p]) == 1:
                        try:
                            isotopes_clustered.remove(p)
                        except KeyError:
                            pass
                # remove isos
                ordered_isos = ordered_isos[:peakel_index]  # exclusive
                wrong_parent.isotopes = set(ordered_isos)

    @staticmethod
    def _promote_to_mo(peakel, isotopes_clustered):
        """
        function used as a callback inside `_find_isotopes`

        :param peakel:
        :param isotopes_clustered:
        :return:
        """

        # selected isos are all different
        # remove isotopes from isotopes, keep all the possibilities
        if peakel.main_attribution is not None:
            peakel.attributions.add(peakel.main_attribution)
            peakel.main_attribution = None
        try:
            isotopes_clustered.remove(peakel)
        except KeyError:
            pass
        return 1

    @staticmethod
    def _set_isotopes(peakel, isotopes):
        """
        function used as a callback inside ``find_isotopes``

        :param peakel:
        :param isotopes:
        :return:
        """
        peakel.isotopes = isotopes

    def _find_isotopes(self, rt_cluster, error_rt=6.0, moz_tol_ppm=10, max_charge=2, max_isotopes_nb=5, max_gap=0):
        """
        a mo can appear locally in an isotope set of another mo, forming an ambiguity
        If we have some clue that an elution peak could be an mo we force it to be in the mo set

        :param rt_cluster
        :param error_rt
        :param moz_tol_ppm
        :param max_charge
        :param max_isotopes_nb
        :param max_gap
        :return two sets: mos, and isotopes
        """
        several_parents_conflicts = 0

        isotopes_clustered = set()

        # link to a direct parent
        mos_by_iso = ddict(set)

        for peakel in rt_cluster:

            # if this considered peakel is a previously detected isotope
            detected_as_iso = True if peakel in mos_by_iso.keys() else False

            result_by_charge = self._look_for_isotopes(peakel,
                                                       mos_by_iso,
                                                       max_charge,
                                                       max_gap,
                                                       max_isotopes_nb,
                                                       error_rt)

            # nothing found
            if not result_by_charge:
                continue

            # best result is the one with the longest
            # what to do when there is draw match ?
            best_charge_result = max(result_by_charge.keys(), key=lambda y: len(result_by_charge[y]))

            # select best isotopes
            selected_isos = result_by_charge[best_charge_result]

            # add annoations to selected isos
            for iso in selected_isos:
                iso.main_attribution = iso.get_attributions_by_charge()[best_charge_result][0]

            if detected_as_iso:
                # if this considered peakel is a previously detected isotope
                # if detected isotopes peak are subset of parent detected isotopes
                selected_isos_including_himself = set(selected_isos)
                selected_isos_including_himself.add(peakel)

                parents = mos_by_iso[peakel]

                n_parents = len(parents)

                # keeping trace of this king of problem
                if n_parents > 1:
                    several_parents_conflicts += 1

                # key: parent_id, value: isotopes, lambda, args
                isos_by_parent = {}

                for mo_parent in parents:
                    # ensure that parents is not in isotopes set

                    if mo_parent.charge == best_charge_result:

                        if selected_isos_including_himself.issubset(mo_parent.isotopes):
                            # no problem, do nothing
                            isos_by_parent[mo_parent] = (selected_isos, None, tuple())

                        else:
                            # an intersection exists, bring back the difference to the parent
                            diff = mo_parent.isotopes.difference(selected_isos_including_himself)
                            isos = mo_parent.isotopes.union(diff)
                            isos_by_parent[mo_parent] = (isos, PeakelsAnnotator._set_isotopes, (mo_parent, isos))
                    else:
                        # if the charge is different, we promote
                        isos_by_parent[mo_parent] = (selected_isos,
                                                     PeakelsAnnotator._promote_to_mo,
                                                     (peakel, isotopes_clustered))
                # end for

                # get the one with the max length
                max_key_len = max(isos_by_parent.keys(), key=lambda l: len(isos_by_parent[l][0]))
                max_len_value = len(isos_by_parent[max_key_len][0])
                best_parents = filter(lambda z: len(isos_by_parent[z][0]) == max_len_value, isos_by_parent.keys())

                if not best_parents:
                    pass
                else:
                    # set main tag and remove from wrong parents
                    best_parent = best_parents[0]
                    best_isos, callback, args = isos_by_parent[best_parent]

                    i = None
                    if callback is not None:
                        i = callback(*args)

                    # we can remove ata least all others wrong parents
                    # if len(best_parents) == 1:
                    PeakelsAnnotator._remove_peakel_from_wrong_parents(peakel,
                                                                       # {best_parent}),
                                                                       parents.difference(best_parents),
                                                                       mos_by_iso,
                                                                       isotopes_clustered)
                    if i is None:
                        peakel.set_main_attribution(peakel.get_attributions_by_parent_id()[best_parent.id][0])

            # end if detected as iso

            peakel.isotopes = selected_isos
            peakel.charge = best_charge_result
            isotopes_clustered = isotopes_clustered.union(selected_isos)

        return list(set(rt_cluster).difference(isotopes_clustered)), isotopes_clustered

    def _find_adducts_and_fragments_in_cluster(self, cluster, mz_tol_ppm=10, max_charge=2):
        """
        :param cluster : set of peakels grouped by retention time
        :param mz_tol_ppm: float, mass tolerance
        :param max_charge: int maximum charge to consider

        could return in case mo not found ?
        """

        # avoid to modify the model using a reference to the parent
        parents_by_son = ddict(list)

        # indexing of the cluster
        index = PeakelIndex(cluster)

        for peakel in cluster:
            for charge in xrange(1, max_charge + 1):
                #
                for adds, attr in self.adducts_or_fragments:
                    mz = peakel.moz / 1 + adds[0]
                    master_peak = index.get_nearest_peakel(mz, mz_tol_ppm)
                    # master peak not found
                    if master_peak is None:
                        pass
                    else:
                        # add possible match son lead to parents
                        attribution = Attribution(attr, master_peak.id, 1)
                        parents_by_son[peakel].append((master_peak, attribution))

        # reverse the dictionary
        adducts_by_mo = ddict(list)
        for add, possible_mos in parents_by_son.iteritems():
            for possible_mo, attrib in possible_mos:
                adducts_by_mo[possible_mo].append((add, attrib))

        # mos = adducts_by_mo.keys()
        if adducts_by_mo:
            # in case there is only one mo per cluster
            best_mos_as_tuple = sorted(list(adducts_by_mo.items()), key=lambda x: len(x[1]))
            best_mos_as_tuple.reverse()
            # best_mo, best_mo_frags = best_mos_as_tuple[0][0], best_mos_as_tuple[0][1]

            best_mos, fragsset = set(), set()

            for mo, frags in best_mos_as_tuple:  # adducts_by_mo.iteritems():
                for frag, attrib in frags:
                    if frag not in fragsset:
                        # frag.main_attribution = attrib
                        frag.set_main_attribution(attrib)
                        mo.adducts.add(frag)
                        best_mos.add(mo)

                        fragsset.add(frag)
                        fragsset.add(mo)

            # todo need more work here
            # best_mo = max(mos, key=lambda y: len(adducts_by_mo[y]))

            # for frag, attrib in best_mo_frags:  #best_mo.adducts:
            #    frag.main_attribution = attrib  #Attribution(annot, best_mo.id, best_mo.charge)

            # best_mo.adducts.union(set(best_mo_frags))

            # best_mos = {best_mo}
            # fragsset = best_mo.adducts.union(best_mo)
            return list(best_mos)  # [best_mo]
        else:
            # for the moment will return the entire list
            return list(cluster)

    def find_adducts_and_fragments(self, clusters):
        """
        todo could be done using multiprocessing or executor
        wrapper for each clusters
        :param clusters:
        """

        return list(chain.from_iterable([self._find_adducts_and_fragments_in_cluster(x) for x in clusters]))

    def annotate(self, error_rt=6.0,
                 max_charge=2,
                 max_isotopes_nb=3,
                 max_gap=0,
                 distance_corr_shape=PeakelClusterer.DEFAULT_SHAPE_CORR,
                 distance_corr_intensity=PeakelClusterer.DEFAULT_INT_CORR):
        """
        Wrapper function around clustering and fragments/adducts
        @param error_rt:
        @param max_charge:
        @param max_isotopes_nb:
        @param max_gap:
        @param distance_corr_shape:
        @param distance_corr_intensity:
        @return:
        """
        self.peakel_clusterer = PeakelClusterer(self.peakels)
        rt_clusters = self.peakel_clusterer.clusterize_by_rt(error_rt=error_rt)
        logging.info('{} rt clusters found'.format(len(rt_clusters)))

        less_isotopes = []  # will be list of list
        for rt_cluster in rt_clusters:
            less_isotopes.append(self._find_isotopes(rt_cluster, error_rt,
                                                     self.exp_settings.mz_tol_ppm,
                                                     max_charge, max_isotopes_nb, max_gap)[0])

        curated_clusters = self.peakel_clusterer._check_update_corrs(less_isotopes,
                                                                     distance_corr_shape,
                                                                     distance_corr_intensity)
        best_mos = self.find_adducts_and_fragments(curated_clusters)
        return best_mos
