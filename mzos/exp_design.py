from collections import defaultdict as ddict
import os.path as op


def enum(**enums):
    """#enumeration
    #backward compatible
    :param enums:
    """
    return type('Enum', (), enums)


IONISATION_MODE = enum(NEG=-1, POS=1)


#==============================================================================
# Experiment
#==============================================================================
class ExperimentalSettings(object):

    """
    :param mz_tol_ppm:
    :param ionisation_mode:
    :param is_dims_experiment:
    """
    ADDUCTS_POS = op.normcase("ressources/POS_ADDUCTS_IMS.csv")
    ADDUCTS_NEG = op.normcase("ressources/NEG_ADDUCTS_IMS.csv")
    FRAGMENTS = op.normcase("ressources/FRAGMENTS_IMS.csv")

    ISOS = "ressources/"

    def __init__(self, mz_tol_ppm, ionisation_mode=None, is_dims_experiment=False, databases={'hmdb'}):
        self.samples = set()

        self.polarity = ionisation_mode  # warning is an ENUM
        self.mz_tol_ppm = mz_tol_ppm
        self.is_dims_exp = is_dims_experiment
        self.databases = databases

        self.group_by_id = ddict(set)
        self.group_by_sample = {}

        #setting isos file, same for both polarity
        self.isos_file = ExperimentalSettings.ISOS

        #setting good frags_file
        self.frags_file, self.adducts_file = ExperimentalSettings.FRAGMENTS, None
        self.adducts_file = ExperimentalSettings.ADDUCTS_NEG \
            if ionisation_mode == IONISATION_MODE.NEG else ExperimentalSettings.ADDUCTS_POS

    def get_frags(self):
        """
        :return:
        """
        lines = list()
        with open(self.frags_file) as f:
            lines += [l.split(",") for l in f.readlines()[1:]]
        return [((float(l[3]), 1), l[0]) for l in lines]

    def get_adducts(self):
        """
        :return:
        """
        lines = list()
        with open(self.adducts_file) as f:
            lines += [l.split(",") for l in f.readlines()[1:]]
        return [((float(l[3]), 1), l[0]) for l in lines]

    def get_mass_to_check(self):
        """
        :return:
        """
        if self.is_dims_exp:
            return self.get_frags()
        return self.get_adducts() + self.get_frags()

    def create_group(self, id_, samples):
        """
        :param id_:
        :param samples:
        :return:
        """
        group = Group(id_, samples)
        for s in list(samples):
            self.group_by_sample[s] = group
        self.group_by_id[id_] = group
        self.samples.union(set(samples))
        return group

    def get_group(self, id_):
        """
        :param id_:
        :return:
        """
        return self.group_by_id.get(id_, None)

    def get_group_of(self, sample):
        """
        :param sample:
        :return: return group or None
        """
        return self.group_by_sample.get(sample, None)

    def get_group_id_of(self, sample):
        """
        :param sample:
        :return:
        """
        group = self.get_group_of(sample)
        if group is None:
            return None
        return group.name_id


class Group (list):
    """
    :param name_id:
    :param samples:
    :param description:
    """

    def __init__(self, name_id, samples=[], description=""):
        list.__init__(samples)

        self.description = description
        self.name_id = name_id