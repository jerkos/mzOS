from collections import defaultdict as ddict

from sklearn.cluster import DBSCAN
from scipy.cluster.hierarchy import linkage, fcluster
import numpy as np


def clusterize_basic(peakels, dist_func, *args):
    """
    Provide a basic clustering based on provided functions, kind of agglomerative
    clustering with mean retention time values for each cluster
    :param peakels:
    :param dist_func: is a callable
    :param args:
    return: list of clusters (as set)
    """
    half_width = args[0] * 0.5
    rt_clusters = []
    peakels_clustered = set()
    for i, peakel in enumerate(peakels):
        if peakel in peakels_clustered:
            continue
        l = set()
        l.add(peakel)
        peakels_clustered.add(peakel)

        for j, peakel_ in enumerate(peakels):
            if peakel_ not in peakels_clustered and i != j:
                if abs(np.mean([p.rt for p in l]) - peakel_.rt) < half_width:
                    l.add(peakel_)
                    peakels_clustered.add(peakel_)

        rt_clusters.append(list(l))
    return rt_clusters
    

def clusterize_hierarchical(peakels, matrix_dist, cut, clip=False):
    """

    :param clip:
    :param peakels:
    :param matrix_dist:
    :param method:
    :param cut:
    """
    # having negative value in the matrix distance
    # leading to a valueerror
    # clip i order to prevent negative value in the matrix distance
    if clip:
        np.clip(matrix_dist, 0, 1, matrix_dist)
    k = linkage(matrix_dist, method='complete')

    # dist = maxdists(k)
    # fit = norm.fit(dist)
    # cut = np.percentile(dist, 10.0)  #norm.ppf(5.0, loc=fit[0], scale=fit[1])

    k2 = fcluster(k, cut, criterion='distance')  # , criterion='distance')
    clust_by_id = ddict(list)
    for i, v in enumerate(k2):
        clust_by_id[v].append(peakels[i])
    return clust_by_id.values()


def clusterize_dbscan(values, peakels, eps=0.2, min_samples=1):
    """
    :param values:
    :param peakels:
    :param eps:
    :param min_samples:
    :return:
    """
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(np.array(values))
    labels = db.labels_
    clust_by_id = ddict(list)
    for i, label in enumerate(labels):
        clust_by_id[label].append(peakels[i])
    return clust_by_id.values()