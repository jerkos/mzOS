from __future__ import absolute_import
from __future__ import print_function

import os.path as op
import os
import logging
import six.moves.cPickle
from collections import defaultdict as ddict

from bioservices import KEGGParser

from mzos.reac import RPrecord


kegg_parser = KEGGParser(verbose=False)


def get_compounds(reaction_id):
    """
    :param reaction_id:
    :return:
    """
    print("treating #reaction_id: {0}".format(reaction_id))
    r = kegg_parser.get(reaction_id)
    reaction = kegg_parser.parse(r)
    reactants, products = reaction["equation"].split("=")

    def get_cpd_ids(string):
        """parse cp ids from kegg
        :param string:
        """
        return [x for x in string.split(" ") if x.startswith("C")]

    reactants_ids = get_cpd_ids(reactants)
    products_ids = get_cpd_ids(products)

    return reactants_ids, products_ids


def get_kegg_reactions():
    """
    :return:
    """
    import multiprocessing
    rp_record_by_id = ddict(RPrecord)

    reac_ids = kegg_parser.reactionIds
    logging.info("# reacids: {0}".format(len(reac_ids)))

    p = multiprocessing.Pool(processes=multiprocessing.cpu_count())

    t = p.map(get_compounds, reac_ids, chunksize=10)

    for reactants_ids, product_ids in t:
        for id__ in reactants_ids:
            for id_ in product_ids:
                rp_record_by_id[id__].as_r.add(id_)
                rp_record_by_id[id_].as_p.add(id__)
    return rp_record_by_id


if __name__ == "__main__":
    d = get_kegg_reactions()
    logging.info("writing reation data to file")

    if not op.exists(op.normcase('ressources')):
        logging.info("mkdir ressources")
        os.mkdir("ressources")

    with open('ressources/reaction.r', 'wb') as output:
        # output = open("reaction.r", 'wb')
        six.moves.cPickle.dump(d, output)