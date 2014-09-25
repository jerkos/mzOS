from __future__ import division

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
import sys
from collections import defaultdict as ddict
from collections import Counter
import random as rdm
import os.path as op
import math
import cPickle

#from mzos.scripts.get_kegg_reactions import load_reactions
from scipy.stats import norm
import numpy as np

import reac
sys.modules['reac'] = reac


def _sample_metabolite(args):  #feature, mz_tol_ppm, probs_by_metab_id, counter):
        """
        Sample metabolites from previous assigned metabolites (assigned_features)
        assigned
        """
        assigned_features = set()

        feature, mz_tol_ppm, probs_by_metab, counter = args[0], args[1], args[2], args[3]

        feature_mass = feature.get_real_mass()

        # retrieve metabolites from that feature
        metabs = feature.get_metabolites()

        #no metabolites
        if not metabs:
            return

        #one metabolite
        if len(metabs) == 1:
            counter[metabs[0].kegg_id] += 1
            assigned_features.add(metabs[0].kegg_id)
            return

        #calc prior probabilities: N(m, mz_tol) * n_assigned / n
        probs = [BayesianInferer._calc_prob(m.mono_mass, feature_mass,
                 mz_tol_ppm) * probs_by_metab[m.kegg_id] for m in metabs]

        #normalize probabilities, which sum is almost one
        norm_probs = BayesianInferer._norm_prob(probs)

        #sample metabolites using numpy multinomial function
        sample = np.random.multinomial(1, norm_probs, size=1)[0]

        # find index
        sampled_metab_id = metabs[np.where(sample == 1)[0][0]].kegg_id

        # add to the assigned features
        assigned_features.add(sampled_metab_id)

        #update counters
        counter[sampled_metab_id] += 1

        return assigned_features, counter


class BayesianInferer(object):
    """
    class implementing  a MCMC sampler in order
    to infer metabolite assignment to a peakel knowing
    informations of the network. The prior probabilities
    are normal and are ponderated by a coefficient resulting
    og the analysis of the network
    """
    def __init__(self, features, experiment):
        """
        features: list of peakels instances
        experiment: experimentdesign object
        """
        self.features = features
        self.experiment = experiment

        logging.info("Loading reaction...")
        self.reactions = self.load_reactions()
        logging.info("#{} reactions loaded".format(len(self.reactions)))
        self.assigned_compounds = set()

    @staticmethod
    def load_reactions():
        """
        can raise IOError
        :return:
        """
        with open(op.normcase("mzos/ressources/reaction.reac"), 'rb') as f:
            reactions = cPickle.load(f)
        return reactions

    def _assign_without_ambiguity(self):
        """
        assign a formula to mass if there is only one possible assignment
        """
        logging.info("First assignment pass...")
        rmsds = []
        c = 0
        for f in self.features:
            metabs = f.get_metabolites()
            if len(metabs) == 1:
                m = metabs[0]
                self.assigned_compounds.add(m.kegg_id)
                c += 1
                #rmsds.append(calculate_rmsd(f, m.isotopic_pattern_pos))

        # fit distribution
        #self.mu_i, self.sigma_i = norm.fit(rmsds)  # could be a poisson distribution
        logging.info("Assign #{} features".format(c))

    @staticmethod
    def _norm_prob(probs):
        """
        Normalize the probabilities by the sum of all
        probabilities. Return a list of normalized probability
        """
        s = sum(probs)
        return [p / s for p in probs]
        #s = np.sum(probs)
        #probs /= s

    @staticmethod
    def _calc_prob(m_obs, m_theo, mz_tol_ppm):
        """
        Return the probability using mass information of
        a metabolite and its theoritical mass and a mass
        tolerance
        """
        n = norm.cdf(m_obs, m_theo, m_theo * mz_tol_ppm / 1e6)
        return (1 - n) * 2 if n > 0.5 else n * 2

    def _calc_network_prob(self, assigned_features, prob_by_metab_by_feature):
        """
        compute network probabilities and update a dictionnary
        for the next step of the MCMC
        """
        for f in self.features:
            metabs = f.get_metabolites()

            probs = []

            sum_prob_assignment = 0.0

            for metab in metabs:
                as_r, as_p = self.reactions[metab.kegg_id].get()
                u = as_r.union(as_p)
                intersec_with_already_identified = assigned_features.intersection(u)
                n = len(intersec_with_already_identified)
                probs.append(n)
                sum_prob_assignment += n
            if sum_prob_assignment:
                prob_by_metab_by_feature[f].update({m.kegg_id: prob / sum_prob_assignment
                                                    for m, prob in zip(metabs, probs)})
            else:
                prob_by_metab_by_feature[f].update({m.kegg_id: 1.0 for m in metabs})

        return prob_by_metab_by_feature

    def _init_probs(self):
        """
        Init assignment and then probabilities
        Do not count into dedicated counter sample
        In order to obtain a faster convergence, we use
        the prior probabilities to sample at the beginning
        """
        assigned_features = set()

        for f in self.features:

            f_mass = f.get_real_mass()
            metabs = f.get_metabolites()
            n = len(metabs)
            if not n:
                continue

            probs = [BayesianInferer._calc_prob(m.mono_mass, f_mass, self.experiment.mz_tol_ppm) for m in metabs]

            norm_probs = BayesianInferer._norm_prob(probs)

            sample = np.random.multinomial(1, norm_probs, size=1)[0]

            sampled_metab_id = metabs[np.where(sample == 1)[0][0]].kegg_id

            # add to the assigned features
            assigned_features.add(sampled_metab_id)

        # dict containing network probability as value
        # and metab.kegg_id as key
        #we do a first pass to know all prob du to previous assignment
        return self._calc_network_prob(assigned_features, ddict(dict))

    def _sample_metabolite_from_feature(self, feature, assigned_features,
                                        probs_by_metab_id_by_feature, counter_by_feature):
        """
        Sample metabolites from previous assigned metabolites (assigned_features)
        assigned
        """
        feature_mass = feature.get_real_mass()

        # retrieve metabolites from that feature
        metabs = feature.get_metabolites()

        #no metabolites
        if not metabs:
            return

        #one metabolite
        if len(metabs) == 1:
            counter_by_feature[feature][metabs[0].kegg_id] += 1
            assigned_features.add(metabs[0].kegg_id)
            return

        #calc prior probabilities: N(m, mz_tol) * n_assigned / n
        probs = [BayesianInferer._calc_prob(m.mono_mass, feature_mass,
                 self.experiment.mz_tol_ppm) * probs_by_metab_id_by_feature[feature][m.kegg_id] for m in metabs]

        #normalize probabilities, which sum is almost one
        norm_probs = BayesianInferer._norm_prob(probs)

        #sample metabolites using numpy multinomial function
        sample = np.random.multinomial(1, norm_probs, size=1)[0]

        # find index
        sampled_metab_id = metabs[np.where(sample == 1)[0][0]].kegg_id

        # add to the assigned features
        assigned_features.add(sampled_metab_id)

        #update counters
        counter_by_feature[feature][sampled_metab_id] += 1

    def _sample(self, assigned_features, probs_by_metab_id_by_f, counter_by_feature):
        """
        """
        #using multiprocessing it is much slower
        # pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
        # args = [(f, self.experiment.mz_tol_ppm, probs_by_metab_id_by_f[f], counter_by_feature[f]) for f in self.features]
        # assigned_mets_and_counters = pool.map(_sample_metabolite, args, chunksize=500)
        # pool.close()
        #
        # for f, t in zip(self.features, assigned_mets_and_counters):
        #     if t is not None:
        #         assigned_features = assigned_features.union(t[0])
        #         counter_by_feature[f].update(t[1])

        for feature in self.features:
            self._sample_metabolite_from_feature(feature, assigned_features,
                                                 probs_by_metab_id_by_f,
                                                 counter_by_feature)

    def _update_probs(self, assigned_compound, prob_by_metab_by_feature):
        self._calc_network_prob(assigned_compound, prob_by_metab_by_feature)

    @staticmethod
    def _calc_posterior_probs(counter_by_feature, prob_by_metab_by_feature, n_samples, n_burning_samples):
        """
        Reduce step, compute frequencies
        """
        #n_s = float(n_samples - n_burning_samples)
        for f, counter in counter_by_feature.iteritems():
            metabs = f.get_metabolites()
            for m in metabs:
                counts = counter[m.kegg_id]
                prob_by_metab_by_feature[f][m.kegg_id] = counts / n_samples  #(counts - n_burning_samples ) / float(n_s)

    def infer_assignment_probabilities(self, n_samples=200, n_burning_sample=10):
        """
        Main function
        :param n_samples:
        :param n_burning_sample:
        """
        assert n_samples > n_burning_sample, "n_burning_samples seems to be > to n_samples"
        #assign obvious
        self._assign_without_ambiguity()

        #g et counter for each feature in order to store sampling

        # mapping between annotation and its metabolite
        # may be useful to retrieve probabilities and assign
        # it the posterior prbability
        annotation_by_metab_id_by_feature = ddict(dict)

        # dict for counting the sampling result in each iteration
        # map a feature and collections.Counter
        counter_by_feature = ddict(Counter)

        # fill the 2 previous declared dicts
        for f in self.features:
            annotations = f.annotations
            #counter = Counter()
            for a in annotations:
                m = a.metabolite
                annotation_by_metab_id_by_feature[f][m.kegg_id] = a
                # we are counting the string
                #counter[m.kegg_id]  #= 0
                counter_by_feature[f][m.kegg_id] = 0  #= counter

        # initalization of prior probabilities
        logging.info("Init probabilities...")
        prob_by_metab_by_feature = self._init_probs()

        for i in xrange(n_samples):
            assigned_compounds = self.assigned_compounds
            rdm.shuffle(self.features)
            #for f in self.features:  #rdm.shuffle(self.features):
            self._sample(assigned_compounds, prob_by_metab_by_feature, counter_by_feature)
            self._update_probs(assigned_compounds, prob_by_metab_by_feature)
            sys.stdout.write("Progression:[*** " + str(int(round(float(i) / n_samples * 100))) + "% ***]" + chr(13))

        logging.info("Computing posterior probabilities...")
        # posterior probabilities calculation
        self._calc_posterior_probs(counter_by_feature,
                                   prob_by_metab_by_feature,
                                   n_samples, n_burning_sample)


        for f, prob_by_metab in prob_by_metab_by_feature.iteritems():
            for metab_id, prob in prob_by_metab.iteritems():
                if math.isnan(prob):
                    raise ValueError("nan probability")
                annotation_by_metab_id_by_feature[f][metab_id].score_network = prob

        logging.info("Done.")