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

import csv
import logging
from feature import Peakel
import numpy as np


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
            all_dirs = gen.next()
        except StopIteration:
            pass
        if all_dirs is None or not all_dirs[1]:
            return []
        group_directories = [1]
        for i in range(len(group_directories)):
            c_dir, dirs, files = gen.next()
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
        if "Mode" in d.keys():
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
        p.area_by_sample_name.update({a: float(b) for a, b in d.items()})
        p.area = np.median(p.area_by_sample_name.values())

        if not p.area:
            p.area = np.mean(p.area_by_sample_name.values())
            logging.debug("an elution peak has the median of "
                          "its area equals to 0 ! Using mean instead: {}.".format(p.area))
        return p