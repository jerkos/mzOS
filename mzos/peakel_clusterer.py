import logging

from clustering import clusterize_basic, clusterize_hierarchical

sp_np = False
try:
    import scipy as sp
    import numpy as np
    sp_np = True
except ImportError:
    print('no numpy, or scipy. Fallback to basic clustering')
    

CLUST_METHOD = {"b_clust": 1, "h_clust": 2}
#==============================================================================
# 
#==============================================================================


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
    default_shape_corr_distance = 0.4
    default_intensity_corr_distance = 0.4  #0.3
    
    #for logging only    
    which_clust_method = staticmethod(lambda x: "basic_clustering" if x == 1 else "hierarchical_clustering mlpy based")
    
    #aggregation condition
    rt_clust_callable = staticmethod(lambda x: True if abs(x[0].rt - x[1].rt) <= x[2] else False)
       
    #callable for correlation distance calculation
    corr_clust_shape_callable = staticmethod(lambda x: True if x[0].corr_shape_against(x[1]) >= x[2] else False)
    corr_clust_intensity_callable = staticmethod(lambda x: True if x[0].corr_intensity_against(x[1]) >= x[2] else False)
    
    def __init__(self, peakels, **kw):
        """
        peakels   
        kw : string to use for clustering method
        """
        self.peakels = peakels
        self.rt_clust_method = kw.get('rt_clust_method', 2)  # basic clustering by default
        
        self.corr_clust_method_shape = kw.get('corr_clust_method_shape', 2)  # mlpy based hierarchical clustering
        self.corr_clust_method_intensity = kw.get('corr_clust_method_intensity', 2)
        
        #logger
        logging.info("rt_clustering_method used:%s" % self.which_clust_method(self.rt_clust_method))
        logging.info("correlation_clustering_method used:%s" %
                     self.which_clust_method(self.corr_clust_method_intensity))

    def set_peakels(self, peakels):
        """
        """        
        self.peakels = peakels
    
    def clusterize_by_rt(self, error_rt):
        """
        PUBLIC function        
        Provide a basic clustering home made
        return: list of clusters(as set) 
        
        """
        if not isinstance(error_rt, float):
            raise TypeError("[clusterize]: args[0] is not a float")
        
        if self.rt_clust_method == 1:
            return clusterize_basic(self.peakels, self.rt_clust_callable, error_rt)  

        elif self.rt_clust_method == 2:
            from sklearn.cluster import DBSCAN
            from collections import defaultdict as ddict
            rt = [[x.rt] for x in self.peakels]
            db = DBSCAN(eps=0.2, min_samples=1).fit(np.array(rt))
            labels = db.labels_
            #print("len labels: {}, len peakels: {}".format(len(labels), len(self.peakels)))
            clust_by_id = ddict(set)
            for i, label in enumerate(labels):
                clust_by_id[label].add(self.peakels[i])
            return clust_by_id.values()
            # rt = [[x.rt] for x in self.peakels]
            # matrix_dist = sp.spatial.distance.pdist(np.array(rt))
            # clust_by_id = clusterize_hierarchical(self.peakels, matrix_dist, "median", error_rt)
            # return clust_by_id.values()
        else:
            raise ValueError("wrong clustering technique !")

    @staticmethod
    def _split_rt_cluster(clust_list):
        clust_list.sort(key=lambda x: len(x))  # sort in place
        main_peakels_set = clust_list[-1]  # longest sets of peakels the most representative
        to_split_from_cluster = clust_list[:-2]
        return main_peakels_set, to_split_from_cluster
        
    def _check_update_corr_shape_in_rt_cluster(self, rt_cluster, distance_corr=default_shape_corr_distance):
        """
        PRIVATE function
        calculate corral        
        """
        clust_list = None
        if self.corr_clust_method_shape == 1:
            clust_list = clusterize_basic(rt_cluster, self.corr_clust_shape_callable, distance_corr)
        
        elif self.corr_clust_method_shape == 2:
            ints = map(lambda x: map(lambda y: y.intensity, x.peaks) if len(x.peaks) else [0], rt_cluster)
            matrix_dist = sp.spatial.distance.pdist(np.array(ints), metric='correlation')
            clust_list = clusterize_hierarchical(rt_cluster, matrix_dist, 'complete', distance_corr).values()
        
        return self._split_rt_cluster(clust_list)

    def _check_update_corr_intensity_in_rt_cluster(self, rt_cluster, distance_corr=default_intensity_corr_distance):
        """
        Private function
        """
        if len(rt_cluster) == 1:
            return rt_cluster, []
        
        clust_list = None
        if self.corr_clust_method_intensity == 1:
            clust_list = clusterize_basic(rt_cluster, self.corr_clust_intensity_callable, distance_corr)
        
        elif self.corr_clust_method_intensity == 2:
            ints = [x.area_by_sample_name.values() for x in rt_cluster]  #
            matrix_dist = sp.spatial.distance.pdist(np.array(ints), metric='correlation')
            clust_list = clusterize_hierarchical(rt_cluster, matrix_dist, 'median', distance_corr).values()
        
        return self._split_rt_cluster(clust_list)

    def check_update_corrs(self, rt_clusters, distance_corr_shape, distance_corr_intensity):
        """
        Private function        
        
        """
#        def split_now_if_needed(rt_cluster, callable_, distance_corr, curated_clusters):
#            rt_cluster, to_split = callable_(rt_cluster, distance_corr)            
#            #first add to split            
#            if not to_split: #to split is empty 
#                for rt_cluster_ in to_split:
#                    curated_clusters.append(rt_cluster_)
#            if not rt_cluster:
#                raise Exception("Check correlation error")
#            curated_clusters.append(rt_cluster)
#            return curated_clusters
        ##########################################################################
#        curated_rt_clusters = []
#        if abs(distance_corr_shape) > 1: 
#            logging.info("shape_correlation skipped because wrong distance_correlation: abs(dist) > 1")
#            curated_rt_clusters = rt_clusters
#        elif all( [ not len(x.peaks) for x in self.peakels ] ):
#            logging.info("shape_correlation skipped because of lacks of peak")
#            curated_rt_clusters = rt_clusters
#        else:
#             for index, rt_cluster in enumerate(rt_clusters):
#                split_now_if_needed(rt_cluster, 
#                                    self._check_update_corr_shape_in_rt_cluster, 
#                                    distance_corr_shape, 
#                                    curated_rt_clusters)
#                                    
#        logging.info("Shape clustering done !")
#        
#        if not distance_corr_intensity:
#            logging.info("shape_correlation skipped because wrong distance_correlation: abs(dist) > 1")
#            return curated_rt_clusters
#            
#        elif all( [not len(x.area_by_sample_name.values()) for x in self.peakels] ):
#            logging.info("intensity_correlation skipped because of lacks of intensities")
#            return curated_rt_clusters
#        
#        else:
        new_curated_rt_clusters = []
        for index, rt_cluster in enumerate(rt_clusters):  # curated_rt_clusters):
            # split_now_if_needed(rt_cluster, self._check_update_corr_intensity_in_rt_cluster,
            # distance_corr_intensity, new_curated_rt_clusters)
            rt_cluster, to_split = self._check_update_corr_intensity_in_rt_cluster(rt_cluster, distance_corr_intensity)
            #first add to split            
            #if to_split: #to split is empty 
            for rt_cluster_ in to_split:
                new_curated_rt_clusters.append(rt_cluster_)
            #if not rt_cluster:
            #    raise Exception("Check correlation error")
            new_curated_rt_clusters.append(rt_cluster)
        return new_curated_rt_clusters

    def clusterize(self, error_rt=10.0,
                   distance_corr_shape=default_shape_corr_distance,
                   distance_corr_intensity=default_intensity_corr_distance):
        """
        """
        
        rt_clusters = self.clusterize_by_rt(error_rt)
        logging.info("Rt clustering done, nb clusters:%d" % len(rt_clusters))
        if rt_clusters is None:
            logging.error("could not make rt clusters...")
            return []
            
        curated = self.check_update_corrs(rt_clusters, distance_corr_shape, distance_corr_intensity)
        logging.info("Intensity clustering done, nb clusters:%d" % len(curated))
        return curated

    def clusterize_(self, error_rt=10.0,
                    distance_corr_shape=default_shape_corr_distance,
                    distance_corr_intensity=default_intensity_corr_distance):
        pass
