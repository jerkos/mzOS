from __future__ import absolute_import
from __future__ import print_function

import logging
from collections import defaultdict as ddict

from bioservices import KEGGParser

import json

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
        return [x.strip() for x in string.split(" ") if x.startswith("C")]

    reactants_ids = get_cpd_ids(reactants)
    products_ids = get_cpd_ids(products)

    return reactants_ids, products_ids


def get_kegg_reactions():
    """
    :return:
    """
    import multiprocessing
    rp_record_by_id = ddict(lambda: ddict(set))

    reac_ids = kegg_parser.reactionIds
    print("# reacids: {0}".format(len(reac_ids)))

    p = multiprocessing.Pool(processes=multiprocessing.cpu_count())

    t = p.map(get_compounds, reac_ids, chunksize=20)

    for reactants_ids, product_ids in t:
        for id__ in reactants_ids:
            for id_ in product_ids:
                rp_record_by_id[id__]['as_r'].add(id_)
                rp_record_by_id[id_]['as_p'].add(id__)

    # transform value to list
    for key in rp_record_by_id:
        v_r = rp_record_by_id[key]['as_r']
        rp_record_by_id[key]['as_r'] = list(v_r)
        v_p = rp_record_by_id[key]['as_p']
        rp_record_by_id[key]['as_p'] = list(v_p)
    print("len rp record: {}".format(len(rp_record_by_id)))
    return rp_record_by_id


def main():
    d = get_kegg_reactions()
    logging.info("writing reation data to file")

    with open('reactions.json', 'wb') as output:
        json.dump(d, output)


if __name__ == "__main__":
    main()
