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

from bioservices import KeggParser
from collections import defaultdict as ddict
from reac import RPrecord
import cPickle


kegg_parser = KeggParser(verbose=False)


def get_compounds(reaction_id):
    print "treating #reaction_id: {}".format(reaction_id)
    r = kegg_parser.get(reaction_id)
    reaction = kegg_parser.parse(r)
    reactants, products = reaction["equation"].split("=")

    def get_cpd_ids(string):
        return [x for x in string.split(" ") if x.startswith("C")]

    reactants_ids = get_cpd_ids(reactants)
    products_ids = get_cpd_ids(products)

    return reactants_ids, products_ids


def get_kegg_reactions():
    import multiprocessing
    rp_record_by_id = ddict(RPrecord)

    reac_ids = kegg_parser.reactionIds
    print "# reacids: {}".format(len(reac_ids))

    p = multiprocessing.Pool(processes=6)

    t = p.map(get_compounds, reac_ids, chunksize=10)

    # for reac_id in reac_ids:
    #     print "reaction id: {}".format(reac_id)
    #     reactants_ids, product_ids = get_compounds(kegg, reac_id)
    #     for id in reactants_ids:
    #         for id_ in product_ids:
    #             rp_record_by_id[id].as_r.add(id_)
    #             rp_record_by_id[id_].as_p.add(id)
    for reactants_ids, product_ids in t:
        for id in reactants_ids:
            for id_ in product_ids:
                rp_record_by_id[id].as_r.add(id_)
                rp_record_by_id[id_].as_p.add(id)
    return rp_record_by_id


def load_reactions():
    return cPickle.load(open("ressources/reaction.reac"))
#
# if __name__ == "__main__":
#     d = get_kegg_reactions()
#     output = open("reaction.r", 'wb')
#     cPickle.dump(d, output)
#     #d = cPickle.load(open("reaction.reac", 'rb'))
#     #print "load #cpds:{}".format(len(d.keys()))
#     #print len(d['C00002'].as_r)
#     # my ask in parallel