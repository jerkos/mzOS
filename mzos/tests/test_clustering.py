__author__ = 'marc.dubois@omics-services.com'

import unittest

import scipy as sp
import numpy as np

from mzos.feature import Peakel
from mzos.clustering import clusterize_hierarchical, clusterize_basic, clusterize_dbscan
from mzos.peakel_clusterer import PeakelClusterer


class TestClustering(unittest.TestCase):

    def setUp(self):
        f1 = Peakel(1256.52, 0.0, 0.0, 1256.52)
        f1.area_by_sample_name = {'a': 102564, 'b': 130156, 'c': 150000, 'd': 10000}

        f2 = Peakel(1258.52, 0.0, 0.0, 1258.52)
        f2.area_by_sample_name = {'a': 102564 / 10.0, 'b': 130156 / 10.0, 'c': 150000 / 10.0, 'd': 10000 / 10.0}

        f3 = Peakel(1261.52, 0.0, 0.0, 1261.52)
        f3.area_by_sample_name = {'a': 102564 / 2.0, 'b': 130156 / 2.3, 'c': 150000 / 2.25, 'd': 10000 / 1.89}

        f4 = Peakel(1262.52, 0.0, 0.0, 1262.52)
        f4.area_by_sample_name = {'a': 102564 * 3.0, 'b': 130156 / 1.8, 'c': 150000 / 4.1, 'd': 10000 * 0.5}

        f5 = Peakel(1274.52, 0.0, 0.0, 1274.52)
        f5.area_by_sample_name = {'a': 102564 * 12.23, 'b': 130156 * 30.0, 'c': 150000 * 33.0, 'd': 10000 * 44.0}

        f6 = Peakel(1275.52, 0.0, 0.0, 1275.52)
        f6.area_by_sample_name = {'a': 102564 * 44.0, 'b': 130156 * 60.0, 'c': 150000 * 66.0, 'd': 10000 * 88.0}

        f7 = Peakel(1281.52, 0.0, 0.0, 1281.52)
        f7.area_by_sample_name = {'a': 102564 * 6.0, 'b': 130156 * 4.56, 'c': 150000 / 78.0, 'd': 10000 / 1236.0}

        self.features = [f1, f2, f3, f4, f5, f6, f7]


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
        for c in clusters:
            for p in c:
                print p.moz
            print "\n"
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

if __name__ == '__main__':
    unittest.main()