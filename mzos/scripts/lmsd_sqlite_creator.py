__author__ = 'marc.dubois@omics-services.com'

import pandas
import os.path as op
import sqlite3

lmsd_tsv = op.normcase("../ressources/LMSD_20130306_All.tsv")


def write_sqlite():
    """
    :return:
    """
    df = pandas.read_csv(lmsd_tsv, sep="\t", error_bad_lines=False)
    conn = sqlite3.connect(op.normcase("../ressources/lmsd.sqlite"))
    df.to_sql("lipids", conn)

if __name__ == "__main__":
    write_sqlite()