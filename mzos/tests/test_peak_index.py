__author__ = 'Marc'

from mzos.feature import Feature, PeakelIndex
import unittest


class RtClustererTest(unittest.TestCase):

    def setUp(self):
        self.features = []
        self.features.append(Feature(1256.52, 1256.52))
        self.features.append(Feature(1258.52, 1258.52))
        self.features.append(Feature(1261.52, 1261.52))
        self.features.append(Feature(1262.52, 1262.52))
        self.features.append(Feature(1274.52, 1274.52))
        self.features.append(Feature(1275.52, 1275.52))
        self.features.append(Feature(1281.52, 1281.52))

    def test_nearest_peak(self):
        findex = PeakelIndex(self.features)
        ppm = 1261.52 * 10 / 1e6
        p = findex.get_nearest_peakel(1261.52, 10)
        print p.moz

    def test_fail(self):
        findex = PeakelIndex(self.features)
        ppm = 1261.52 * 10 / 1e6
        val = 1261.52 - ppm
        p = findex.get_nearest_peakel(val, 10)
        self.assertIsNone(p)

if __name__ == '__main__':
    unittest.main()