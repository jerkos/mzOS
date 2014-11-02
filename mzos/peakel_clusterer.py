import logging

import scipy as sp
import numpy as np
from sklearn.metrics.pairwise import pairwise_distances

from clustering import clusterize_basic, clusterize_hierarchical, clusterize_dbscan


class PeakelClusterer(object):
    """
    PeakelClusterer class
    provides several to clusterize a list of elution peaks using Rt and cross-abundances

    Methods:
    -------

    DBSCAN: scikit learn implementation
    SIMPLE: clustering.py
    H_CLUST: hierarchical clustering

    """
    CLUST_METHOD = {"basic": 1, "hierarchical": 2, "dbscan": 3}
    REV_CLUST_METHOD = {v: k for k, v in CLUST_METHOD.iteritems()}

    DEFAULT_SHAPE_CORR = 0.4
    DEFAULT_INT_CORR = 0.4

    #aggregation condition
    BASIC_RT_CALLABLE = staticmethod(lambda x, y, z: True if abs(x.rt - z.rt) <= z else False)
       
    #callable for correlation distance calculation
    BASIC_CORR_SHAPE_CALLABLE = staticmethod(lambda x, y, z: True if x.corr_shape_against(y) >= z else False)
    BASIC_CORR_INT_CALLABLE = staticmethod(lambda x, y, z: True if x.corr_intensity_against(y) >= z else False)
    
    def __init__(self, peakels, **kw):
        """
        peakels   
        kw : string to use for clustering method
        """
        self.peakels = peakels
        self.rt_method = kw.get('rt_clust_method', 3)  # dbscan clustering by default
        
        self.corr_shape_method = kw.get('corr_shape_method')
        self.corr_int_method = kw.get('corr_int_method')

        # no correlation method provided, default intensity method hierarchical clustering
        if not self.corr_shape_method and not self.corr_int_method:
            self.corr_int_method = self.CLUST_METHOD['hierarchical']
            #raise ValueError("no correlation method found")

        # the 2 correlations method are provided, just need one
        if self.corr_int_method and self.corr_shape_method:
            raise ValueError("both correlation method are defined")

        corr_method_used = ""
        if self.corr_shape_method:
            corr_method_used += " ".join(["correlation shape", self.REV_CLUST_METHOD[self.corr_shape_method]])
        else:
            corr_method_used += " ".join(["correlation intensity", self.REV_CLUST_METHOD[self.corr_int_method]])
        #logger
        logging.info("rt clustering_method used: {}".format(self.REV_CLUST_METHOD[self.rt_method]))
        logging.info("correlation clustering method used: {}".format(corr_method_used))

    def set_peakels(self, peakels):
        """
        :param peakels:
        """
        self.peakels = peakels
    
    def clusterize_by_rt(self, error_rt):
        """
        PUBLIC function        
        Provide a basic clustering home made
        :param error_rt:
        return: list of clusters(as set)
        
        """
        if self.rt_method == 1:
            logging.info("Basic clustering with rt_error:{}".format(error_rt))
            if not isinstance(error_rt, float):
                raise TypeError("[clusterize]: args[0] is not a float")
            return clusterize_basic(self.peakels, self.BASIC_RT_CALLABLE, error_rt)

        elif self.rt_method == 2:
            rts = [[x.rt] for x in self.peakels]
            matrix_dist = sp.spatial.distance.pdist(np.array(rts))  #metric = eclidean by default
            return clusterize_hierarchical(self.peakels, matrix_dist, "", error_rt).values()

        elif self.rt_method == 3:
            logging.info('DB SCAN clustering with error_rt:{}'.format(error_rt))
            rts = [[x.rt] for x in self.peakels]
            clusters = clusterize_dbscan(rts, self.peakels, eps=0.35)
            # with open('clusters.txt', 'w') as f:
            #     for c in clusters:
            #         for fe in c:
            #             f.write(str(fe.rt) + '\n')
            #         f.write('\n')
            return clusters  # eps=error_rt / 2.0, min_samples=1)

        else:
            raise ValueError("wrong clustering technique !")

    @staticmethod
    def _split_rt_cluster(clust_list):
        clust_list.sort(key=lambda x: len(x))  # sort in place
        main_peakels_set = clust_list[-1]  # longest sets of peakels the most representative
        to_split_from_cluster = clust_list[:-2]
        return main_peakels_set, to_split_from_cluster
        
    def _check_update_corr_shape_in_rt_cluster(self, rt_cluster, distance_corr=DEFAULT_SHAPE_CORR):
        """
        PRIVATE function
        calculate corral        
        """
        clust_list = None
        if self.corr_shape_method == 1:
            clust_list = clusterize_basic(rt_cluster, self.BASIC_CORR_SHAPE_CALLABLE, distance_corr)
        
        elif self.corr_shape_method == 2:
            ints = map(lambda x: map(lambda y: y.intensity, x.peaks) if len(x.peaks) else [0], rt_cluster)
            matrix_dist = sp.spatial.distance.pdist(np.array(ints), metric='correlation')
            clust_list = clusterize_hierarchical(rt_cluster, matrix_dist, distance_corr, clip=True)
        
        return self._split_rt_cluster(clust_list)

    def _check_update_corr_intensity_in_rt_cluster(self, rt_cluster, distance_corr=DEFAULT_INT_CORR):
        """
        Private function
        """
        if len(rt_cluster) == 1:
            return []  #rt_cluster, []
        
        #clust_list = None
        if self.corr_int_method == 1:
            clust_list = clusterize_basic(rt_cluster, self.BASIC_CORR_INT_CALLABLE, distance_corr)
        
        elif self.corr_int_method == 2:
            ints = [x.area_by_sample_name.values() for x in rt_cluster]  #
            #matrix_dist = sp.spatial.distance.pdist(np.array(ints), metric='correlation')
            #ude by default all cores on the machine
            matrix_dist = pairwise_distances(np.array(ints), metric='correlation', n_jobs=-1)
            clust_list = clusterize_hierarchical(rt_cluster, matrix_dist, distance_corr, clip=True)
        else:
            raise ValueError("dbscan not supported for intensities correlation clustering")
        
        return clust_list  #self._split_rt_cluster(clust_list)

    def check_update_corrs(self, rt_clusters, corr_shape_dist, corr_int_dist):
        """
        Private function
        :param rt_clusters:
        :param corr_shape_dist:
        :param corr_int_dist:
        """
        new_curated_rt_clusters = []
        if self.corr_shape_method:
            for rt_cluster in rt_clusters:  # curated_rt_clusters):
                new_curated_rt_clusters += self._check_update_corr_shape_in_rt_cluster(rt_cluster, corr_shape_dist)
        else:
            for rt_cluster in rt_clusters:
                new_curated_rt_clusters += self._check_update_corr_intensity_in_rt_cluster(rt_cluster, corr_int_dist)

        return new_curated_rt_clusters

    def clusterize(self, error_rt=10.0,
                   corr_shape_dist=DEFAULT_SHAPE_CORR,
                   corr_int_dist=DEFAULT_INT_CORR):
        """
        :param corr_shape_dist:
        :param corr_int_dist:
        :param error_rt:
        """
        
        rt_clusters = self.clusterize_by_rt(error_rt)
        logging.info("Rt clustering done, nb clusters: {}".format(len(rt_clusters)))
        if rt_clusters is None:
            logging.error("could not make rt clusters...")
            return []
            
        curated = self.check_update_corrs(rt_clusters, corr_shape_dist, corr_int_dist)
        logging.info("Intensity clustering done, nb clusters:{}".format(len(curated)))
        return curated

    # def clusterize_(self, error_rt=10.0,
    #                 distance_corr_shape=DEFAULT_SHAPE_CORR,
    #                 distance_corr_intensity=DEFAULT_INT_CORR):
    #     pass
