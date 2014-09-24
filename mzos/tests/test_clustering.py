__author__ = 'marc.dubois@omics-services.com'

import unittest
import logging

import scipy as sp
import numpy as np

from mzos.feature import Feature
from mzos.clustering import clusterize_hierarchical, clusterize_basic, clusterize_dbscan
from mzos.peakel_clusterer import PeakelClusterer


class TestClustering(unittest.TestCase):

    def setUp(self):
        self.features = []
        self.features.append(Feature(1256.52, 1256.52))
        self.features.append(Feature(1258.52, 1258.52))
        self.features.append(Feature(1261.52, 1261.52))
        self.features.append(Feature(1262.52, 1262.52))
        self.features.append(Feature(1274.52, 1274.52))
        self.features.append(Feature(1275.52, 1275.52))
        self.features.append(Feature(1281.52, 1281.52))

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

    def test_clusterize_dbscan(self):
        clusters = clusterize_dbscan(self.features, eps=3.0, min_samples=1)
        print("len clusters dbscan: {}".format(len(clusters)))
        self.assertGreaterEqual(4, len(clusters))


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    unittest.main()