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

import glob
from xml.etree.cElementTree import iterparse
#import time
import base64
import struct
from itertools import chain
import numpy as np
import re


class Scan(dict):
    """

    """

    def __init__(self, mz_array, int_array, metadata={}):
        dict.__init__(self, metadata)
        self.mz_array = mz_array
        self.int_array = int_array


def decode_spectrum(string64):
    """

    """
    data = base64.b64decode(string64)
    count = len(data) / struct.calcsize('!d')
    data = struct.unpack('!' + 'd' * count, data[0:len(data)])
    # data = np.array(data)
    # print data
    # data.shape = (-1, 2)
    # return data
    mz_array, int_array = [], []
    for i in xrange(0, len(data) - 1, 2):
        mz_array.append(data[i])
        int_array.append(data[i + 1])
    return np.array(mz_array), np.array(int_array)


def encode_spectrum(mz_array, int_array):
    """

    """
    t = list(chain.from_iterable(zip(mz_array, int_array)))
    scanbyte = struct.pack('!' + 'd' * len(t), *t)
    return base64.b64encode(scanbyte)


def load_spectra(filepath):
    #------------------------------------------------------------------------------
    def _get_infos():
        """
        case we have a mzxml file, get run length and number of scans
        sentinel MAX
        """
        p = re.compile('\s+<msRun\sscanCount="(\d+)"\sstartTime="PT(\d+\.\d+|\d+)S"\sendTime="PT(\d+\.\d+)S">')
        q = re.compile('<mzXML\sxmlns="(.+)"')  #to get the namespace
        MAX = 10  #ten lines max the skip
        with open(filepath) as fd:
            line = fd.readline()
            i = 0
            while q.match(line) is None and i < MAX:
                line = fd.readline()
                i += 1
            prefix = "{" + q.match(line).group(1) + "}"
            line = fd.readline()
            i = 0
            while p.match(line) is None and i < MAX:
                line = fd.readline()
                i += 1
            scan_count = int(p.match(line).group(1))
            max_time = float(p.match(line).group(3))
            min_time = float(p.match(line).group(2))
            return max_time, min_time, scan_count, prefix

    def _get_header():
        """
        get the header of file until "scan" balise starts
        useful for creating new xml files

        """
        MAX = 35
        p = re.compile('\s+<scan\snum="1"\n')
        header = ""
        with open(filepath) as fd:
            line = fd.readline()
            i = 0
            while p.match(line) is None and i < MAX:
                header += line
                line = fd.readline()
                i += 1
            if i == 34:
                print("header parsing failed")
        return header

    #----------------------------------------------------------------------
    max_time, min_time, scan_count, prefix = _get_infos()
    header = _get_header()
    #t = time.clock()
    scans = []

    context = iterparse(filepath, events=('end',))

    for action, elem in context:
        if elem.tag == "".join([prefix, "scan"]) and action == 'end':
            scan = {'scanNumber': int(elem.attrib.get('num', 0)),
                    'msLevel': int(elem.attrib.get('msLevel', 1)),
                    'peaksCount': int(elem.attrib.get('peaksCount', 0)),
                    'centroided': elem.attrib.get('centroided', '0'),
                    'scanType': elem.attrib.get('scanType', ""),
                    'retentionTime': float(elem.attrib['retentionTime'].strip('PTS')),
                    'basePeakMz': float(elem.attrib.get('basePeakMz', 0)),
                    'basePeakIntensity': float(elem.attrib.get('basePeakIntensity', 0)),
                    'totIonCurrent': float(elem.attrib.get('totIonCurrent', 0.0)),
                    'polarity': elem.attrib.get('polarity', "")}
            if scan['msLevel'] == 1:
                #get data points
                for e in elem.getchildren():
                    if e.tag == "".join([prefix, "peaks"]) and action == 'end':
                        scan['byteOrder'] = e.attrib.get('byteOrder', 'network')
                        scan['precision'] = int(e.attrib.get('precision', '64'))
                        scan['compressionType'] = e.attrib.get('compressionType', 'none')
                        mz_array, int_array = decode_spectrum(e.text)
                        scans.append(Scan(mz_array, int_array, scan))
                        #treat points
                        break
    del context
    return header, scans


def merge_spectra(scans, size_bin=0.002, method="mean"):
    """

    """
    if len(scans) == 1:
        return scans[0]

    if not all([isinstance(x, Scan) for x in scans]):
        raise TypeError(["merged spectra: wrong data types"])

    minmz = np.array([o.mz_array.min() for o in scans]).min()
    maxmz = np.array([o.mz_array.max() for o in scans]).max()

    diff_mz = maxmz - minmz
    length = int(round(diff_mz / size_bin))

    xbinning = np.array([minmz + i * size_bin for i in xrange(length + 1)])
    intensity = np.zeros(length + 1)

    if method == "mean":
        for e in scans:
            bin_ = ((e.mz_array - minmz) / size_bin).astype(int)
            intensity[bin_] += e.int_array  #/ e['totIonCurrent'])
        intensity /= float(len(scans))
    elif method == "tic":
        for e in scans:
            bin_ = ((e.mz_array - minmz) / size_bin).astype(int)
            intensity[bin_] += e.int_array / e['totIonCurrent']

    scan = {'scanNumber': 1,
            'msLevel': 1,
            'peaksCount': len(xbinning),
            'centroided': '0',
            'scanType': 'Full',
            'retentionTime': np.array([e['retentionTime'] for e in scans]).mean(),
            #'basePeakMz': xbinning[np.where(intensity == intensity.max())[1]],
            'basePeakIntensity': intensity.max(),
            'totIonCurrent': intensity.sum(),
            'polarity': '-',
            'byteOrder': 'network',
            'precision': 64,
            'compressionType': 'none'}

    return Scan(xbinning, intensity, scan)


def write_spectrum(header, scan):
    """
    """
    if not isinstance(header, str):
        raise TypeError("writeSpectrum")
    string = """    <scan num="1"
      scanType="Full"
      centroided="0"
      msLevel="1"
      peaksCount="{}"
      polarity="-"
      retentionTime="PT{}S"
      basePeakIntensity="{}"
      totIonCurrent="{}"
      msInstrumentID="1">
    <peaks compressionType="none"
            compressionLen="0"
            precision="64"
            byteOrder="network"
            contentType="m/z-int">{}</peaks>
    </scan>
  </msRun>
</mzXML>""".format(scan['peaksCount'],
                   scan['retentionTime'],
                   scan['basePeakIntensity'],
                   scan['totIonCurrent'],
                   encode_spectrum(scan.mz_array, scan.int_array))
    header += string
    return header


def create_file((filepath, mintime, maxtime)):
    """

    """
    if not isinstance(mintime, float) or not isinstance(maxtime, float):
        raise TypeError("min time and maxtime must be floating numbers")
    header, scans = load_spectra(filepath)
    good_time_scans = [scan for scan in scans]  #[scan for scan in scans if mintime < scan['retentionTime'] < maxtime]
    new_scan = None
    if not good_time_scans:
        print("No scan defined in range: {}, {} in file {}".format(mintime, maxtime, filepath))
        return new_scan
    elif len(good_time_scans) == 1:
        new_scan = good_time_scans[0]
    else:
        new_scan = merge_spectra(good_time_scans, size_bin=0.0005)

    header = write_spectrum(header, new_scan)
    with open(filepath.split(".")[0] + "modified.mzXML", 'w') as f:
        f.write(header)
    return new_scan

if __name__ == '__main__':
    #import multiprocessing
    files = glob.glob("tests/200-1000/*.mzXML")
    filesplusargs = [(f, 0.4, 1.6) for f in files]
    #print filesplusargs
    for t in filesplusargs:
        print "Working on {}".format(t[0])
        create_file(t)
    # p = multiprocessing.Pool(processes=4)
    # r = p.map(create_file, filesplusargs, chunksize=2)
