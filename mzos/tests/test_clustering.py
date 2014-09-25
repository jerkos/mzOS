from mzos.database_finder import Metabolite

__author__ = 'marc.dubois@omics-services.com'

import unittest

import scipy as sp
import numpy as np

from mzos.feature import Peakel, Annotation, Attribution
from mzos.clustering import clusterize_hierarchical, clusterize_basic, clusterize_dbscan
from mzos.peakel_clusterer import PeakelClusterer
from mzos.exp_design import ExperimentalSettings

class TestClustering(unittest.TestCase):

    def setUp(self):
        self.f1 = Peakel(1256.52, 0.0, 0.0, 1256.52)
        self.f1.area_by_sample_name = {'a': 102564, 'b': 130156, 'c': 150000, 'd': 10000}
        self.f1.charge = 1
        self.f1.polarity = -1

        self.f2 = Peakel(1258.52, 0.0, 0.0, 1258.52)
        self.f2.area_by_sample_name = {'a': 102564 / 10.0, 'b': 130156 / 10.0, 'c': 150000 / 10.0, 'd': 10000 / 10.0}

        self.f3 = Peakel(1261.52, 0.0, 0.0, 1261.52)
        self.f3.area_by_sample_name = {'a': 102564 / 2.0, 'b': 130156 / 2.3, 'c': 150000 / 2.25, 'd': 10000 / 1.89}

        self.f4 = Peakel(1262.52, 0.0, 0.0, 1262.52)
        self.f4.area_by_sample_name = {'a': 102564 * 3.0, 'b': 130156 / 1.8, 'c': 150000 / 4.1, 'd': 10000 * 0.5}

        self.f5 = Peakel(1274.52, 0.0, 0.0, 1274.52)
        self.f5.area_by_sample_name = {'a': 102564 * 12.23, 'b': 130156 * 30.0, 'c': 150000 * 33.0, 'd': 10000 * 44.0}

        self.f6 = Peakel(1275.52, 0.0, 0.0, 1275.52)
        self.f6.area_by_sample_name = {'a': 102564 * 44.0, 'b': 130156 * 60.0, 'c': 150000 * 66.0, 'd': 10000 * 88.0}

        self.f7 = Peakel(1281.52, 0.0, 0.0, 1281.52)
        self.f7.area_by_sample_name = {'a': 102564 * 6.0, 'b': 130156 * 4.56, 'c': 150000 / 78.0, 'd': 10000 / 1236.0}

        self.features = [self.f1, self.f2, self.f3, self.f4, self.f5, self.f6, self.f7]

        self.f1.isotopes.add(self.f2)
        self.f2.set_main_attribution(Attribution('isotope c13', self.f1.id, 1))
        self.f1.adducts.add(self.f3)
        self.f3.set_main_attribution(Attribution('[M+Na+]', self.f1, 1))
        t = ("acession", "name", "formula", "inchi", "mono_mass", "average_mass", "description", "status", "origin",
             "kegg_id", "isotopic_pattern_pos", "isotopic_pattern_neg")
        self.f1.annotations.append(Annotation(metabolite=Metabolite._make(t)))

    def test_clusterize_basic(self):
        clusters = clusterize_basic(self.features, PeakelClusterer.BASIC_RT_CALLABLE, 6.0)
        print("len clusters basic: {}".format(len(clusters)))
        self.assertGreaterEqual(4, len(clusters))

    def test_clusterize_hierarchical(self):
        rts = [[f.rt] for f in self.features]
        matrix_dist = sp.spatial.distance.pdist(np.array(rts))  # euclidean distance
        clusters = clusterize_hierarchical(self.features, matrix_dist, 3.0)
        print("len clusters hierarchical: {}".format(len(clusters)))
        self.assertGreaterEqual(4, len(clusters))

    def test_clusterize_dbscan_rt(self):
        clusters = clusterize_dbscan([[x.rt] for x in self.features], self.features, eps=3.0, min_samples=1)
        print("len clusters dbscan: {}".format(len(clusters)))
        self.assertGreaterEqual(4, len(clusters))

    def test_clusterize_hierarchical_int(self):
        ints = [f.area_by_sample_name.values() for f in self.features]
        matrix_dist = sp.spatial.distance.pdist(np.array(ints), metric="correlation")  # euclidean distance
        clusters = clusterize_hierarchical(self.features, matrix_dist, 0.1, clip=True)
        print("len clusters hierarchical int: {}".format(len(clusters)))
        self.assertGreaterEqual(4, len(clusters))

    def test_peakel_clusterer_main(self):
        peakel_clusterer = PeakelClusterer(self.features, rt_clust_method=3, corr_int_method=2)
        clusters = peakel_clusterer.clusterize(error_rt=6.0)
        self.assertGreaterEqual(4, len(clusters))

    def test_peakel_cluster_basic(self):
        peakel_clusterer = PeakelClusterer(self.features, rt_clust_method=1, corr_int_method=2)
        clusters = peakel_clusterer.clusterize(error_rt=6.0)
        self.assertGreaterEqual(4, len(clusters))

    def test_peakel_cluster_hier(self):
        peakel_clusterer = PeakelClusterer(self.features, rt_clust_method=3, corr_int_method=2)
        clusters = peakel_clusterer.clusterize(error_rt=6.0)
        self.assertGreaterEqual(4, len(clusters))

    def test_error_1(self):
        #peakel_clusterer = PeakelClusterer(self.features, )
        self.assertRaises(ValueError, PeakelClusterer, self.features, rt_clust_method=3,
                                           corr_int_method=2, corr_shape_method=2)

    def test_error_2(self):
        peakel_clusterer = PeakelClusterer(self.features, rt_clust_method=3)
        self.assertEqual(peakel_clusterer.corr_int_method, 2)
        self.assertIsNone(peakel_clusterer.corr_shape_method)

    def test_feature(self):
        ms_name = [m.name for m in self.f1.get_metabolites()]
        self.assertIn('name', ms_name)

        d = self.f2.get_attributions_by_parent_id()
        self.assertIn(self.f1.id, d.keys())

        s, nb_isos, nb_adducts = self.f1.get_top_down_attribution_tree()
        self.assertEqual(nb_isos, 1)
        self.assertEqual(nb_adducts, 1)
        print s
        self.assertIn('({}=isotope c13)'.format(self.f2.id), s)
        self.assertIn('({}=[M+Na+])'.format(self.f3.id), s)

        ip = self.f1.get_isotopic_pattern()
        self.assertEqual(len(ip), 2)
        self.assertTrue(all([len(x) == 2 for x in ip]))

        ip2 = self.f1.get_isotopic_pattern_as_peakel()
        self.assertIn(self.f1, ip2)

        real_mass = self.f1.get_real_mass()
        self.assertAlmostEqual(real_mass, self.f1.moz + 1.007276)

    def test_exp_design(self):
        exp = ExperimentalSettings(mz_tol_ppm=10.0)
        samples = ['sample1', 'sample2', 'sample3', 'sample4']
        group1 = exp.create_group("Blanks", samples[:2])
        group2 = exp.create_group("Treated", samples[2:])
        self.assertEqual(group1, exp.get_group("Blanks"))

        self.assertIs(group2, exp.get_group_of('sample4'))

        self.assertEqual('Treated', exp.get_group_id_of('sample4'))
        self.assertIsNone(exp.get_group_id_of('sample5'))


if __name__ == '__main__':
    unittest.main()