from __future__ import absolute_import
import zipfile
import os.path as op
import os
import shutil
import logging


class WithHMDBMixin(object):

    @staticmethod
    def unzip_hmdb():
        """
        Utility to unzip hmdb for test purposes
        :param self:
        :return:
        """
        z = zipfile.ZipFile(op.abspath('mzos/ressources/hmdb.zip'))
        hmdb_path = z.extract('hmdb.sqlite')
        logging.info("Moving extracted archive...")
        shutil.move(hmdb_path, 'mzos/ressources/hmdb.sqlite')
        logging.info("Done")

    @staticmethod
    def remove_hmdb():
        logging.info("removing 'hmdb.sqlite'...")
        try:
            os.remove(op.abspath('mzos/ressources/hmdb.sqlite'))
            logging.info("Done")
        except OSError:
            logging.error("Unable to remove sqlite file or file does not exist")