from __future__ import absolute_import
import csv
import logging

import numpy as np

from mzos.feature import Peakel
from six.moves import range


class PeakListReader(object):
    """
    peaklist reader
    """
    KEYS = ['mz', 'mzmin', 'mzmax', 'rt', 'rtmin', 'rtmax', 'npeaks']
    
    def __init__(self, filepath, exp_design):
        """
        @param filepath:
        @param exp_design @type ExperimentalSettings
        """
        self.peaklist_filepath = filepath
        self.exp_design = exp_design
        self.directories = self._find_directories()  #

    def get_peakels(self):
        """return peakels objects """
        reader = csv.DictReader(open(self.peaklist_filepath, 'rb'), delimiter="\t")
        return [self._to_peakel_obj(row) for row in reader]
    
    def _find_directories(self):
        import os
        curr_dir = os.path.dirname(self.peaklist_filepath)
        # find dir
        gen = os.walk(curr_dir)
        all_dirs = None
        try:
            all_dirs = next(gen)
        except StopIteration:
            pass
        if all_dirs is None or not all_dirs[1]:
            return []
        group_directories = [1]
        for i in range(len(group_directories)):
            c_dir, dirs, files = next(gen)
            self.exp_design.create_group(c_dir, files)

        return group_directories

    def _to_peakel_obj(self, d):
        """ convert csv data to peakel objects  """
        p = Peakel(float(d[PeakListReader.KEYS[0]]),
                   float(d[PeakListReader.KEYS[1]]),
                   float(d[PeakListReader.KEYS[2]]),
                   float(d[PeakListReader.KEYS[3]]),
                   float(d[PeakListReader.KEYS[4]]),
                   float(d[PeakListReader.KEYS[5]]))

        #  set the right polarity
        polarity = None
        if "Mode" in list(d.keys()):
            polarity = 1 if d["Mode"] == 'Positif' else -1
        else:
            p.polarity = self.exp_design.polarity

        #  remove keys
        for k in (PeakListReader.KEYS + self.directories + ["", "BIO", "mzmed", "rt.minutes", "Var", "Blc.Ext", "BLC",
                                                            "Mode", "Correlation_Dilution_Log", "NOT_M.QC", "NOT_M.Blc",
                                                            "NOT_QC.Blc", "NOT_CV..", "NOT_CV", "NOT_Correl", "Correl",
                                                            "NOT_BIO.Blc", "NOT_nom", "rt.min", "Negatifs"]):
            try:
                del d[k]
            except KeyError:
                pass

        #  assign area of each sample
        p.area_by_sample_name.update({a: float(b) for a, b in list(d.items())})
        p.area = np.median(list(p.area_by_sample_name.values()))

        if not p.area:
            p.area = np.mean(list(p.area_by_sample_name.values()))
            logging.debug("an elution peak has the median of "
                          "its area equals to 0 ! Using mean instead: {}.".format(p.area))
        return p