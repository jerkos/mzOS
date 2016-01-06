from __future__ import absolute_import

import pandas
import os.path as op
import sqlite3


def main():
    lmsd_tsv = op.abspath("mzos/ressources/LMSD_20130306_All.tsv")
    df = pandas.read_csv(lmsd_tsv, sep="\t", error_bad_lines=False)
    conn = sqlite3.connect(op.abspath("mzos/ressources/lmsd.sqlite"))
    df.to_sql("lipids", conn)


if __name__ == "__main__":
    main()
