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

import os.path as op
import os
import logging
import cPickle
from collections import defaultdict as ddict

from bioservices import KEGGParser  #KeggParser

from mzos.reac import RPrecord


kegg_parser = KEGGParser(verbose=False)


def get_compounds(reaction_id):
    """
    :param reaction_id:
    :return:
    """
    print "treating #reaction_id: {}".format(reaction_id)
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
    logging.info("# reacids: {}".format(len(reac_ids)))

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
        #output = open("reaction.r", 'wb')
        cPickle.dump(d, output)