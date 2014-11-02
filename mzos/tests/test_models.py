from mzos.database_finder import Metabolite

__author__ = 'marc.dubois@omics-services.com'

import unittest
import os.path as op
import logging

import scipy as sp
import numpy as np

from mzos.feature import Peakel, Annotation, Attribution, PeakelIndex
from mzos.clustering import clusterize_hierarchical, clusterize_basic, clusterize_dbscan
from mzos.peakel_clusterer import PeakelClusterer
from mzos.exp_design import ExperimentalSettings
from mzos.formula import Formula


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
        self.f3.set_main_attribution(Attribution('[M+Na+]', self.f1.id, 1))
        t = ("acession", "name", "formula", "inchi", "mono_mass", "average_mass", "description", "status", "origin",
             "kegg_id", "isotopic_pattern_pos", "isotopic_pattern_neg")
        self.f1.annotations.append(Annotation(metabolite=Metabolite(*t), for_adduct='H2'))

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
        for c in clusters:
            for f in c:
                print f.rt
            print '\n'

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
        metabolites = self.f1.get_metabolites()
        print metabolites[0].__dict__
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

        f_by_id = {f.id: f for f in self.features}
        s = self.f2.get_bottom_up_attribution_tree(f_by_id)
        self.assertEqual(s, 'isotope c13 of {} for charge={}'.format(self.f1.id, self.f2.main_attribution.charge))

        a = Attribution('isotope s34', self.f3.id, charge=2)
        self.f2.add_attribution(a)
        s = Peakel.get_others_bottom_up_attribution_tree(a, f_by_id)
        self.assertEqual(s, 'isotope s34 of {} for charge={} of [M+Na+] of {} for charge=1'.format(self.f3.id, 2, self.f1.id))

        self.f2.get_attributions_by(lambda _: _.attribution)

        self.f2.get_attributions_by_charge()

        self.f2.get_attributions_by_parent_id()

        self.f2.remove_attribution_with_parent(self.f3.id)

        real_mass = self.f1.get_real_mass()
        self.assertAlmostEqual(real_mass, self.f1.moz + 1.007276)

    def test_nearest_peak(self):
        findex = PeakelIndex(self.features)
        ppm = 1261.52 * 10 / 1e6
        p = findex.get_nearest_peakel(1261.52, 10)
        self.assertAlmostEqual(p.moz, 1261.52)

    def test_fail(self):
        findex = PeakelIndex(self.features)
        ppm = 1261.52 * 10 / 1e6
        val = 1261.52 - ppm
        p = findex.get_nearest_peakel(val, 10)
        self.assertIsNone(p)

    def test_exp_design(self):
        exp = ExperimentalSettings(mz_tol_ppm=10.0)
        samples = ['sample1', 'sample2', 'sample3', 'sample4']
        group1 = exp.create_group("Blanks", samples[:2])
        group2 = exp.create_group("Treated", samples[2:])
        self.assertEqual(group1, exp.get_group("Blanks"))

        self.assertIs(group2, exp.get_group_of('sample4'))

        self.assertEqual('Treated', exp.get_group_id_of('sample4'))
        self.assertIsNone(exp.get_group_id_of('sample5'))

    def test_fromula_obj(self):
        fstr = 'C6H12O6'
        f = Formula.from_str(fstr)
        self.assertIsNotNone(f)

        self.assertEqual(f['C'], 6)
        self.assertEqual(f['H'], 12)
        self.assertEqual(f['O'], 6)

        self.assertEqual(str(f), fstr)

        f.remove('C5H3O12')
        self.assertEqual(str(f), 'CH9')

        d = {'C': 148, 'H': 122}
        f.add(d)
        self.assertEqual(str(f), 'C149H131')

        f2 = Formula(d)
        f.remove(f2)
        self.assertEqual(str(f), 'CH9')

        self.assertRaises(TypeError, f.add, ['C', 12, 'H', 32])

        f3 = Formula.from_str('C4H12')
        f4 = Formula(f3)
        self.assertIsNot(f3, f4)

        self.assertEqual(str(f4), 'C4H12')

        f5 = f3.remove('CH', new_obj=True)

        self.assertEqual(str(f3), 'C4H12')
        self.assertIsNot(f5, f3)
        self.assertEqual(str(f5), 'C3H11')

        #test + -
        i = Formula.from_str('C4H4O4')
        i += 'C4H4O4'
        self.assertEqual(str(i), 'C8H8O8')

        j = i + 'C2H2O2'
        self.assertIsNot(j, i)
        self.assertEqual(str(j), 'C10H10O10')
        self.assertEqual(str(i), 'C8H8O8')

        #not exits anymore print f5.get_theo_ip()

    def test_script_hmdb(self):
        from mzos.scripts.hmdb_sqlite_creator import build_library
        build_library(op.normcase('hmdb_test.sqlite'),
                      op.normcase('mzos/ressources/hmdb_metabolites'),
                      op.normcase(Formula.EMASS_PATH))
        self.assertTrue(op.exists('hmdb_test.sqlite'))

    def test_dbscan_clustering_for_alignment(self):
        f1 = Peakel(1256.52, 0.0, 0.0, 100.0)
        f2 = Peakel(1258.52, 0.0, 0.0, 500.52)
        f3 = Peakel(1257.52, 0.0, 0.0, 101.52)
        f4 = Peakel(1600.52, 0.0, 0.0, 99.52)
        f7 = Peakel(1600.86, 0.0, 0.0, 107.12)
        f5 = Peakel(1600.52, 0.0, 0.0, 3205.52)
        f6 = Peakel(1456.52, 0.0, 0.0, 600.52)

        peakels = [f1, f2, f3, f4, f5, f6, f7]
        peakels_by_sample = {'a': {f1, f2, f4}, 'b': {f3, f5, f6}}

        sample_by_peakel = {f1: 'a', f2: 'a', f4: 'a', f3: 'b', f5: 'b', f6: 'b', f7: 'b'}
        values = [[x.moz, x.rt] for x in peakels]

        clusters = clusterize_dbscan(values, peakels, eps=5, min_samples=1)

        for c in clusters:
            for p in c:
                print p.moz, p.rt, sample_by_peakel[p]
            print '\n'


    @classmethod
    def tearDownClass(cls):
        import os
        try:
            os.remove('hmdb_test.sqlite')
        except (WindowsError, IOError):
            pass