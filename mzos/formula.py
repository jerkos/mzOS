__author__ = 'Marc'

import re
import subprocess
import logging
import os.path as op
from collections import Counter


class Element:
    """Element object definition.
        name: (str) name
        symbol: (str) symbol
        atomic_nb: (int) atomic number
        isotopes: (dict) dict of isotopes {mass number:(mass, abundance),...}
        valence: (int)
    """

    def __init__(self, name, symbol, atomic_nb, isotopes, valence=None):

        self.name = name
        self.symbol = symbol
        self.atomic_nb = int(atomic_nb)
        self.isotopes = isotopes
        self.valence = valence

        # init masses
        mass_mo = 0
        mass_av = 0
        max_abundance = 0
        for isotop in self.isotopes.values():
            mass_av += isotop[0] * isotop[1]
            if max_abundance < isotop[1]:
                mass_mo = isotop[0]
                max_abundance = isotop[1]
        if mass_mo == 0 or mass_av == 0:
            mass_mo = isotopes[0][0]
            mass_av = isotopes[0][0]

        self.mass = (mass_mo, mass_av)


ELEMENTS = {
    'Ac': Element(name='Actinium', symbol='Ac', atomic_nb=89, isotopes={227: (227.02774700000001, 1.0)}, valence=3),
    'Ag': Element(name='Silver', symbol='Ag', atomic_nb=47, isotopes={107: (106.90509299999999, 0.51839000000000002),
                                                                      109: (108.90475600000001, 0.48160999999999998)},
                  valence=1),
    'Al': Element(name='Aluminium', symbol='Al', atomic_nb=13, isotopes={27: (26.981538440000001, 1.0)}, valence=3),
    'Am': Element(name='Americium', symbol='Am', atomic_nb=95,
                  isotopes={241: (241.05682289999999, 0.0), 243: (243.06137269999999, 1.0)}, valence=3),
    'Ar': Element(name='Argon', symbol='Ar', atomic_nb=18, isotopes={40: (39.962383123000002, 0.99600299999999997),
                                                                     36: (35.967546280000001, 0.0033649999999999999),
                                                                     38: (37.962732199999998, 0.00063199999999999997)},
                  valence=0),
    'As': Element(name='Arsenic', symbol='As', atomic_nb=33, isotopes={75: (74.921596399999999, 1.0)}, valence=3),
    'At': Element(name='Astatine', symbol='At', atomic_nb=85,
                  isotopes={210: (209.98713100000001, 0.0), 211: (210.987481, 1.0)}, valence=1),
    'Au': Element(name='Gold', symbol='Au', atomic_nb=79, isotopes={197: (196.96655200000001, 1.0)}, valence=1),
    'B': Element(name='Boron', symbol='B', atomic_nb=5,
                 isotopes={10: (10.012937000000001, 0.19900000000000001), 11: (11.0093055, 0.80100000000000005)},
                 valence=3),
    'Ba': Element(name='Barium', symbol='Ba', atomic_nb=56,
                  isotopes={130: (129.90630999999999, 0.00106), 132: (131.905056, 0.00101),
                            134: (133.90450300000001, 0.024170000000000001),
                            135: (134.90568300000001, 0.065920000000000006),
                            136: (135.90457000000001, 0.078539999999999999), 137: (136.905821, 0.11232),
                            138: (137.90524099999999, 0.71697999999999995)}, valence=2),
    'Be': Element(name='Beryllium', symbol='Be', atomic_nb=4, isotopes={9: (9.0121821000000004, 1.0)}, valence=2),
    'Bh': Element(name='Bohrium', symbol='Bh', atomic_nb=107, isotopes={264: (264.12473, 1.0)}, valence=0),
    'Bi': Element(name='Bismuth', symbol='Bi', atomic_nb=83, isotopes={209: (208.98038299999999, 1.0)}, valence=5),
    'Bk': Element(name='Berkelium', symbol='Bk', atomic_nb=97,
                  isotopes={249: (249.07498000000001, 0.0), 247: (247.07029900000001, 1.0)}, valence=3),
    'Br': Element(name='Bromine', symbol='Br', atomic_nb=35, isotopes={81: (80.916291000000001, 0.49309999999999998),
                                                                       79: (78.918337600000001, 0.50690000000000002)},
                  valence=1),
    'C': Element(name='Carbon', symbol='C', atomic_nb=6,
                 isotopes={12: (12.0, 0.98929999999999996), 13: (13.0033548378, 0.010699999999999999),
                           14: (14.003241987999999, 0.0)}, valence=4),
    'Ca': Element(name='Calcium', symbol='Ca', atomic_nb=20, isotopes={40: (39.962591199999999, 0.96940999999999999),
                                                                       42: (41.958618299999998, 0.0064700000000000001),
                                                                       43: (42.958766799999999, 0.0013500000000000001),
                                                                       44: (43.9554811, 0.02086),
                                                                       46: (45.953692799999999, 4.0000000000000003e-05),
                                                                       48: (47.952534, 0.0018699999999999999)},
                  valence=2),
    'Cd': Element(name='Cadmium', symbol='Cd', atomic_nb=48,
                  isotopes={106: (105.906458, 0.012500000000000001), 108: (107.904183, 0.0088999999999999999),
                            110: (109.903006, 0.1249), 111: (110.90418200000001, 0.128),
                            112: (111.9027572, 0.24129999999999999), 113: (112.9044009, 0.1222),
                            114: (113.90335810000001, 0.2873), 116: (115.90475499999999, 0.074899999999999994)},
                  valence=2),
    'Ce': Element(name='Cerium', symbol='Ce', atomic_nb=58,
                  isotopes={136: (135.90714, 0.0018500000000000001), 138: (137.90598600000001, 0.0025100000000000001),
                            140: (139.90543400000001, 0.88449999999999995), 142: (141.90924000000001, 0.11114)},
                  valence=4),
    'Cf': Element(name='Californium', symbol='Cf', atomic_nb=98,
                  isotopes={249: (249.07484700000001, 0.0), 250: (250.07640000000001, 0.0),
                            251: (251.07957999999999, 1.0), 252: (252.08161999999999, 0.0)}, valence=3),
    'Cl': Element(name='Chlorine', symbol='Cl', atomic_nb=17,
                  isotopes={35: (34.96885271, 0.75780000000000003), 37: (36.9659026, 0.2422)}, valence=1),
    'Cm': Element(name='Curium', symbol='Cm', atomic_nb=96,
                  isotopes={243: (243.0613822, 0.0), 244: (244.06274629999999, 0.0), 245: (245.06548559999999, 0.0),
                            246: (246.06721759999999, 0.0), 247: (247.070347, 1.0), 248: (248.07234199999999, 0.0)},
                  valence=3),
    'Co': Element(name='Cobalt', symbol='Co', atomic_nb=27, isotopes={59: (58.933200200000002, 1.0)}, valence=2),
    'Cr': Element(name='Chromium', symbol='Cr', atomic_nb=24, isotopes={50: (49.946049600000002, 0.043450000000000003),
                                                                        52: (51.940511899999997, 0.83789000000000002),
                                                                        53: (52.9406538, 0.095009999999999997),
                                                                        54: (53.938884899999998, 0.023650000000000001)},
                  valence=3),
    'Cs': Element(name='Caesium', symbol='Cs', atomic_nb=55, isotopes={133: (132.90544700000001, 1.0)}, valence=1),
    'Cu': Element(name='Copper', symbol='Cu', atomic_nb=29, isotopes={65: (64.927793699999995, 0.30830000000000002),
                                                                      63: (62.929601099999999, 0.69169999999999998)},
                  valence=1),
    'Db': Element(name='Dubnium', symbol='Db', atomic_nb=105, isotopes={262: (262.11415, 1.0)}, valence=0),
    'Dy': Element(name='Dysprosium', symbol='Dy', atomic_nb=66,
                  isotopes={160: (159.925194, 0.023400000000000001), 161: (160.92693, 0.18909999999999999),
                            162: (161.926795, 0.25509999999999999), 163: (162.92872800000001, 0.249),
                            164: (163.929171, 0.28179999999999999), 156: (155.92427799999999, 0.00059999999999999995),
                            158: (157.92440500000001, 0.001)}, valence=3),
    'Er': Element(name='Erbium', symbol='Er', atomic_nb=68,
                  isotopes={162: (161.928775, 0.0014), 164: (163.92919699999999, 0.0161),
                            166: (165.93029000000001, 0.33610000000000001), 167: (166.93204499999999, 0.2293),
                            168: (167.932368, 0.26779999999999998), 170: (169.93546000000001, 0.14929999999999999)},
                  valence=3),
    'Es': Element(name='Einsteinium', symbol='Es', atomic_nb=99, isotopes={252: (252.08296999999999, 1.0)}, valence=3),
    'Eu': Element(name='Europium', symbol='Eu', atomic_nb=63, isotopes={153: (152.92122599999999, 0.52190000000000003),
                                                                        151: (150.91984600000001, 0.47810000000000002)},
                  valence=2),
    'F': Element(name='Fluorine', symbol='F', atomic_nb=9, isotopes={19: (18.998403199999998, 1.0)}, valence=1),
    'Fe': Element(name='Iron', symbol='Fe', atomic_nb=26,
                  isotopes={56: (55.934942100000001, 0.91754000000000002), 57: (56.9353987, 0.021190000000000001),
                            58: (57.933280500000002, 0.00282), 54: (53.939614800000001, 0.058450000000000002)},
                  valence=2),
    'Fm': Element(name='Fermium', symbol='Fm', atomic_nb=100, isotopes={257: (257.095099, 1.0)}, valence=3),
    'Fr': Element(name='Francium', symbol='Fr', atomic_nb=87, isotopes={223: (223.0197307, 1.0)}, valence=1),
    'Ga': Element(name='Gallium', symbol='Ga', atomic_nb=31,
                  isotopes={69: (68.925580999999994, 0.60107999999999995), 71: (70.924705000000003, 0.39892)},
                  valence=3),
    'Gd': Element(name='Gadolinium', symbol='Gd', atomic_nb=64,
                  isotopes={160: (159.92705100000001, 0.21859999999999999), 152: (151.91978800000001, 0.002),
                            154: (153.920862, 0.0218), 155: (154.922619, 0.14799999999999999),
                            156: (155.92212000000001, 0.20469999999999999), 157: (156.923957, 0.1565),
                            158: (157.92410100000001, 0.24840000000000001)}, valence=3),
    'Ge': Element(name='Germanium', symbol='Ge', atomic_nb=32, isotopes={72: (71.922076200000006, 0.27539999999999998),
                                                                         73: (72.923459399999999, 0.077299999999999994),
                                                                         74: (73.9211782, 0.36280000000000001),
                                                                         76: (75.921402700000002, 0.076100000000000001),
                                                                         70: (69.924250400000005, 0.2084)}, valence=4),
    'H': Element(name='Hydrogen', symbol='H', atomic_nb=1,
                 isotopes={1: (1.0078250321, 0.99988500000000002), 2: (2.0141017780000001, 0.000115),
                           3: (3.0160492675000001, 0.0)}, valence=1),
    'He': Element(name='Helium', symbol='He', atomic_nb=2,
                  isotopes={3: (3.0160293096999999, 1.37e-06), 4: (4.0026032496999999, 0.99999863)}, valence=0),
    'Hf': Element(name='Hafnium', symbol='Hf', atomic_nb=72, isotopes={174: (173.94004000000001, 0.0016000000000000001),
                                                                       176: (175.94140179999999, 0.052600000000000001),
                                                                       177: (176.94322, 0.186),
                                                                       178: (177.9436977, 0.27279999999999999),
                                                                       179: (178.9458151, 0.13619999999999999),
                                                                       180: (179.94654879999999, 0.3508)}, valence=4),
    'Hg': Element(name='Mercury', symbol='Hg', atomic_nb=80,
                  isotopes={196: (195.96581499999999, 0.0015), 198: (197.96675200000001, 0.099699999999999997),
                            199: (198.96826200000001, 0.16869999999999999), 200: (199.968309, 0.23100000000000001),
                            201: (200.97028499999999, 0.1318), 202: (201.97062600000001, 0.29859999999999998),
                            204: (203.97347600000001, 0.068699999999999997)}, valence=2),
    'Ho': Element(name='Holmium', symbol='Ho', atomic_nb=67, isotopes={165: (164.930319, 1.0)}, valence=3),
    'I': Element(name='Iodine', symbol='I', atomic_nb=53, isotopes={127: (126.90446799999999, 1.0)}, valence=1),
    'In': Element(name='Indium', symbol='In', atomic_nb=49,
                  isotopes={113: (112.904061, 0.042900000000000001), 115: (114.90387800000001, 0.95709999999999995)},
                  valence=3),
    'Ir': Element(name='Iridium', symbol='Ir', atomic_nb=77,
                  isotopes={193: (192.96292399999999, 0.627), 191: (190.96059099999999, 0.373)}, valence=3),
    'K': Element(name='Potassium', symbol='K', atomic_nb=19,
                 isotopes={40: (39.963998670000002, 0.000117), 41: (40.96182597, 0.067302000000000001),
                           39: (38.963706899999998, 0.93258099999999999)}, valence=1),
    'Kr': Element(name='Krypton', symbol='Kr', atomic_nb=36, isotopes={78: (77.920385999999993, 0.0035000000000000001),
                                                                       80: (79.916377999999995, 0.022800000000000001),
                                                                       82: (81.913484600000004, 0.1158),
                                                                       83: (82.914135999999999, 0.1149),
                                                                       84: (83.911507, 0.56999999999999995),
                                                                       86: (85.910610300000002, 0.17299999999999999)},
                  valence=0),
    'La': Element(name='Lanthanum', symbol='La', atomic_nb=57,
                  isotopes={138: (137.907107, 0.00089999999999999998), 139: (138.90634800000001, 0.99909999999999999)},
                  valence=3),
    'Li': Element(name='Lithium', symbol='Li', atomic_nb=3, isotopes={6: (6.0151222999999998, 0.075899999999999995),
                                                                      7: (7.0160039999999997, 0.92410000000000003)},
                  valence=1),
    'Lr': Element(name='Lawrencium', symbol='Lr', atomic_nb=103, isotopes={262: (262.10969, 1.0)}, valence=3),
    'Lu': Element(name='Lutetium', symbol='Lu', atomic_nb=71,
                  isotopes={176: (175.9426824, 0.025899999999999999), 175: (174.9407679, 0.97409999999999997)},
                  valence=3),
    'Md': Element(name='Mendelevium', symbol='Md', atomic_nb=101,
                  isotopes={256: (256.09404999999998, 0.0), 258: (258.09842500000002, 1.0)}, valence=3),
    'Mg': Element(name='Magnesium', symbol='Mg', atomic_nb=12, isotopes={24: (23.985041899999999, 0.78990000000000005),
                                                                         25: (24.985837020000002, 0.10000000000000001),
                                                                         26: (25.982593040000001, 0.1101)}, valence=2),
    'Mn': Element(name='Manganese', symbol='Mn', atomic_nb=25, isotopes={55: (54.938049599999999, 1.0)}, valence=2),
    'Mo': Element(name='Molybdenum', symbol='Mo', atomic_nb=42,
                  isotopes={96: (95.904678899999993, 0.1668), 97: (96.906020999999996, 0.095500000000000002),
                            98: (97.905407800000006, 0.24129999999999999), 100: (99.907477, 0.096299999999999997),
                            92: (91.906809999999993, 0.1484), 94: (93.905087600000002, 0.092499999999999999),
                            95: (94.905841499999994, 0.15920000000000001)}, valence=6),
    'Mt': Element(name='Meitnerium', symbol='Mt', atomic_nb=109, isotopes={268: (268.13882000000001, 1.0)}, valence=0),
    'N': Element(name='Nitrogen', symbol='N', atomic_nb=7,
                 isotopes={14: (14.0030740052, 0.99631999999999998), 15: (15.000108898400001, 0.0036800000000000001)},
                 valence=3),
    'Na': Element(name='Sodium', symbol='Na', atomic_nb=11, isotopes={23: (22.989769670000001, 1.0)}, valence=1),
    'Nb': Element(name='Niobium', symbol='Nb', atomic_nb=41, isotopes={93: (92.906377500000005, 1.0)}, valence=5),
    'Nd': Element(name='Neodymium', symbol='Nd', atomic_nb=60,
                  isotopes={142: (141.90771899999999, 0.27200000000000002), 143: (142.90980999999999, 0.122),
                            144: (143.91008299999999, 0.23799999999999999),
                            145: (144.91256899999999, 0.083000000000000004),
                            146: (145.91311200000001, 0.17199999999999999), 148: (147.916889, 0.057000000000000002),
                            150: (149.92088699999999, 0.056000000000000001)}, valence=3),
    'Ne': Element(name='Neon', symbol='Ne', atomic_nb=10, isotopes={20: (19.992440175900001, 0.90480000000000005),
                                                                    21: (20.993846739999999, 0.0027000000000000001),
                                                                    22: (21.991385510000001, 0.092499999999999999)},
                  valence=0),
    'Ni': Element(name='Nickel', symbol='Ni', atomic_nb=28, isotopes={64: (63.927969599999997, 0.0092560000000000003),
                                                                      58: (57.935347899999996, 0.68076899999999996),
                                                                      60: (59.930790600000002, 0.26223099999999999),
                                                                      61: (60.9310604, 0.011398999999999999),
                                                                      62: (61.928348800000002, 0.036345000000000002)},
                  valence=2),
    'No': Element(name='Nobelium', symbol='No', atomic_nb=102, isotopes={259: (259.10102000000001, 1.0)}, valence=2),
    'Np': Element(name='Neptunium', symbol='Np', atomic_nb=93,
                  isotopes={237: (237.04816729999999, 1.0), 239: (239.05293140000001, 0.0)}, valence=3),
    'O': Element(name='Oxygen', symbol='O', atomic_nb=8,
                 isotopes={16: (15.9949146221, 0.99756999999999996), 17: (16.999131500000001, 0.00038000000000000002),
                           18: (17.999160400000001, 0.0020500000000000002)}, valence=2),
    'Os': Element(name='Osmium', symbol='Os', atomic_nb=76,
                  isotopes={192: (191.961479, 0.4078), 184: (183.95249100000001, 0.00020000000000000001),
                            186: (185.95383799999999, 0.015900000000000001),
                            187: (186.95574790000001, 0.019599999999999999),
                            188: (187.95583600000001, 0.13239999999999999), 189: (188.95814490000001, 0.1615),
                            190: (189.95844500000001, 0.2626)}, valence=3),
    'P': Element(name='Phosphorus', symbol='P', atomic_nb=15, isotopes={31: (30.973761509999999, 1.0)}, valence=3),
    'Pa': Element(name='Protactinium', symbol='Pa', atomic_nb=91, isotopes={231: (231.0358789, 1.0)}, valence=4),
    'Pb': Element(name='Lead', symbol='Pb', atomic_nb=82,
                  isotopes={208: (207.97663600000001, 0.52400000000000002), 204: (203.973029, 0.014),
                            206: (205.97444899999999, 0.24099999999999999), 207: (206.97588099999999, 0.221)},
                  valence=4),
    'Pd': Element(name='Palladium', symbol='Pd', atomic_nb=46,
                  isotopes={102: (101.905608, 0.010200000000000001), 104: (103.90403499999999, 0.1114),
                            105: (104.905084, 0.2233), 106: (105.90348299999999, 0.27329999999999999),
                            108: (107.90389399999999, 0.2646), 110: (109.905152, 0.1172)}, valence=2),
    'Pm': Element(name='Promethium', symbol='Pm', atomic_nb=61,
                  isotopes={145: (144.912744, 1.0), 147: (146.91513399999999, 0.0)}, valence=3),
    'Po': Element(name='Polonium', symbol='Po', atomic_nb=84, isotopes={209: (208.982416, 1.0), 210: (209.982857, 0.0)},
                  valence=2),
    'Pr': Element(name='Praseodymium', symbol='Pr', atomic_nb=59, isotopes={141: (140.90764799999999, 1.0)}, valence=3),
    'Pt': Element(name='Platinum', symbol='Pt', atomic_nb=78,
                  isotopes={192: (191.96103500000001, 0.0078200000000000006),
                            194: (193.96266399999999, 0.32967000000000002),
                            195: (194.96477400000001, 0.33832000000000001), 196: (195.964935, 0.25241999999999998),
                            198: (197.96787599999999, 0.071629999999999999),
                            190: (189.95993000000001, 0.00013999999999999999)}, valence=2),
    'Pu': Element(name='Plutonium', symbol='Pu', atomic_nb=94,
                  isotopes={238: (238.04955340000001, 0.0), 239: (239.0521565, 0.0), 240: (240.0538075, 0.0),
                            241: (241.05684529999999, 0.0), 242: (242.05873679999999, 0.0), 244: (244.064198, 1.0)},
                  valence=3),
    'Ra': Element(name='Radium', symbol='Ra', atomic_nb=88,
                  isotopes={224: (224.02020200000001, 0.0), 226: (226.02540260000001, 1.0),
                            228: (228.03106410000001, 0.0), 223: (223.018497, 0.0)}, valence=2),
    'Rb': Element(name='Rubidium', symbol='Rb', atomic_nb=37, isotopes={85: (84.911789299999995, 0.72170000000000001),
                                                                        87: (86.909183499999997, 0.27829999999999999)},
                  valence=1),
    'Re': Element(name='Rhenium', symbol='Re', atomic_nb=75,
                  isotopes={185: (184.95295569999999, 0.374), 187: (186.9557508, 0.626)}, valence=4),
    'Rf': Element(name='Rutherfordium', symbol='Rf', atomic_nb=104, isotopes={261: (261.10874999999999, 1.0)},
                  valence=0),
    'Rh': Element(name='Rhodium', symbol='Rh', atomic_nb=45, isotopes={103: (102.90550399999999, 1.0)}, valence=3),
    'Rn': Element(name='Radon', symbol='Rn', atomic_nb=86,
                  isotopes={211: (210.99058500000001, 0.0), 220: (220.01138409999999, 0.0),
                            222: (222.01757050000001, 1.0)}, valence=0),
    'Ru': Element(name='Ruthenium', symbol='Ru', atomic_nb=44, isotopes={96: (95.907597999999993, 0.055399999999999998),
                                                                         98: (97.905287000000001, 0.018700000000000001),
                                                                         99: (98.9059393, 0.12759999999999999),
                                                                         100: (99.904219699999999, 0.126),
                                                                         101: (100.9055822, 0.1706),
                                                                         102: (101.9043495, 0.3155),
                                                                         104: (103.90543, 0.1862)}, valence=3),
    'S': Element(name='Sulfur', symbol='S', atomic_nb=16,
                 isotopes={32: (31.972070689999999, 0.94930000000000003), 33: (32.971458499999997, 0.0076),
                           34: (33.967866829999998, 0.042900000000000001),
                           36: (35.967080879999997, 0.00020000000000000001)}, valence=2),
    'Sb': Element(name='Antimony', symbol='Sb', atomic_nb=51,
                  isotopes={121: (120.903818, 0.57210000000000005), 123: (122.90421569999999, 0.4279)}, valence=5),
    'Sc': Element(name='Scandium', symbol='Sc', atomic_nb=21, isotopes={45: (44.955910199999998, 1.0)}, valence=3),
    'Se': Element(name='Selenium', symbol='Se', atomic_nb=34, isotopes={74: (73.922476599999996, 0.0088999999999999999),
                                                                        76: (75.919214100000005, 0.093700000000000006),
                                                                        77: (76.919914599999998, 0.076300000000000007),
                                                                        78: (77.917309500000002, 0.23769999999999999),
                                                                        80: (79.916521799999998, 0.49609999999999999),
                                                                        82: (81.916700000000006, 0.087300000000000003)},
                  valence=2),
    'Sg': Element(name='Seaborgium', symbol='Sg', atomic_nb=106, isotopes={266: (266.12193000000002, 1.0)}, valence=0),
    'Si': Element(name='Silicon', symbol='Si', atomic_nb=14, isotopes={28: (27.976926532699999, 0.92229700000000003),
                                                                       29: (28.976494720000002, 0.046831999999999999),
                                                                       30: (29.973770219999999, 0.030872)}, valence=4),
    'Sm': Element(name='Samarium', symbol='Sm', atomic_nb=62, isotopes={144: (143.91199499999999, 0.030700000000000002),
                                                                        147: (146.91489300000001, 0.14990000000000001),
                                                                        148: (147.914818, 0.1124),
                                                                        149: (148.91718, 0.13819999999999999),
                                                                        150: (149.917271, 0.073800000000000004),
                                                                        152: (151.91972799999999, 0.26750000000000002),
                                                                        154: (153.92220499999999, 0.22750000000000001)},
                  valence=2),
    'Sn': Element(name='Tin', symbol='Sn', atomic_nb=50,
                  isotopes={112: (111.904821, 0.0097000000000000003), 114: (113.902782, 0.0066),
                            115: (114.903346, 0.0033999999999999998), 116: (115.90174399999999, 0.1454),
                            117: (116.90295399999999, 0.076799999999999993), 118: (117.901606, 0.2422),
                            119: (118.90330899999999, 0.085900000000000004), 120: (119.9021966, 0.32579999999999998),
                            122: (121.9034401, 0.046300000000000001), 124: (123.9052746, 0.0579)}, valence=4),
    'Sr': Element(name='Strontium', symbol='Sr', atomic_nb=38, isotopes={88: (87.905614299999996, 0.82579999999999998),
                                                                         84: (
                                                                         83.913425000000004, 0.0055999999999999999),
                                                                         86: (85.909262400000003, 0.098599999999999993),
                                                                         87: (
                                                                         86.908879299999995, 0.070000000000000007)},
                  valence=2),
    'Ta': Element(name='Tantalum', symbol='Ta', atomic_nb=73,
                  isotopes={180: (179.94746599999999, 0.00012), 181: (180.94799599999999, 0.99987999999999999)},
                  valence=5),
    'Tb': Element(name='Terbium', symbol='Tb', atomic_nb=65, isotopes={159: (158.925343, 1.0)}, valence=3),
    'Tc': Element(name='Technetium', symbol='Tc', atomic_nb=43,
                  isotopes={97: (96.906364999999994, 0.0), 98: (97.907216000000005, 1.0),
                            99: (98.906254599999997, 0.0)}, valence=4),
    'Te': Element(name='Tellurium', symbol='Te', atomic_nb=52,
                  isotopes={128: (127.9044614, 0.31740000000000002), 130: (129.90622279999999, 0.34079999999999999),
                            120: (119.90402, 0.00089999999999999998), 122: (121.90304709999999, 0.025499999999999998),
                            123: (122.904273, 0.0088999999999999999), 124: (123.90281950000001, 0.047399999999999998),
                            125: (124.90442470000001, 0.070699999999999999), 126: (125.9033055, 0.18840000000000001)},
                  valence=2),
    'Th': Element(name='Thorium', symbol='Th', atomic_nb=90,
                  isotopes={232: (232.0380504, 1.0), 230: (230.0331266, 2.3203809999999998)}, valence=4),
    'Ti': Element(name='Titanium', symbol='Ti', atomic_nb=22,
                  isotopes={48: (47.9479471, 0.73719999999999997), 49: (48.947870799999997, 0.054100000000000002),
                            50: (49.944792100000001, 0.051799999999999999), 46: (45.9526295, 0.082500000000000004),
                            47: (46.951763800000002, 0.074399999999999994)}, valence=4),
    'Tl': Element(name='Thallium', symbol='Tl', atomic_nb=81,
                  isotopes={203: (202.972329, 0.29524), 205: (204.974412, 0.70476000000000005)}, valence=3),
    'Tm': Element(name='Thulium', symbol='Tm', atomic_nb=69, isotopes={169: (168.934211, 1.0)}, valence=3),
    'U': Element(name='Uranium', symbol='U', atomic_nb=92, isotopes={233: (233.03962799999999, 2.3802891000000002),
                                                                     234: (234.04094559999999, 5.5000000000000002e-05),
                                                                     235: (235.0439231, 0.0071999999999999998),
                                                                     236: (236.0455619, 0.0),
                                                                     238: (238.05078259999999, 0.99274499999999999)},
                 valence=3),
    'V': Element(name='Vanadium', symbol='V', atomic_nb=23, isotopes={50: (49.947162800000001, 0.0025000000000000001),
                                                                      51: (50.943963699999998, 0.99750000000000005)},
                 valence=5),
    'W': Element(name='Tungsten', symbol='W', atomic_nb=74,
                 isotopes={184: (183.95093259999999, 0.30640000000000001), 186: (185.954362, 0.2843),
                           180: (179.94670600000001, 0.0011999999999999999), 182: (181.948206, 0.26500000000000001),
                           183: (182.95022449999999, 0.1431)}, valence=6),
    'Xe': Element(name='Xenon', symbol='Xe', atomic_nb=54, isotopes={128: (127.90353039999999, 0.019199999999999998),
                                                                     129: (128.90477949999999, 0.26440000000000002),
                                                                     130: (129.90350789999999, 0.040800000000000003),
                                                                     131: (130.9050819, 0.21179999999999999),
                                                                     132: (131.9041545, 0.26889999999999997),
                                                                     134: (133.9053945, 0.10440000000000001),
                                                                     136: (135.90722, 0.088700000000000001),
                                                                     124: (123.9058958, 0.00089999999999999998),
                                                                     126: (125.904269, 0.00089999999999999998)},
                  valence=0),
    'Y': Element(name='Yttrium', symbol='Y', atomic_nb=39, isotopes={89: (88.905847899999998, 1.0)}, valence=3),
    'Yb': Element(name='Ytterbium', symbol='Yb', atomic_nb=70,
                  isotopes={168: (167.93389400000001, 0.0012999999999999999), 170: (169.93475900000001, 0.0304),
                            171: (170.93632199999999, 0.14280000000000001),
                            172: (171.93637770000001, 0.21829999999999999), 173: (172.93820679999999, 0.1613),
                            174: (173.9388581, 0.31830000000000003), 176: (175.94256799999999, 0.12759999999999999)},
                  valence=2),
    'Zn': Element(name='Zinc', symbol='Zn', atomic_nb=30, isotopes={64: (63.929146600000003, 0.48630000000000001),
                                                                    66: (65.926036800000006, 0.27900000000000003),
                                                                    67: (66.927130899999995, 0.041000000000000002),
                                                                    68: (67.924847600000007, 0.1875),
                                                                    70: (69.925325000000001, 0.0061999999999999998)},
                  valence=2),
    'Zr': Element(name='Zirconium', symbol='Zr', atomic_nb=40, isotopes={96: (95.908276000000001, 0.028000000000000001),
                                                                         90: (89.904703699999999, 0.51449999999999996),
                                                                         91: (90.905645000000007, 0.11219999999999999),
                                                                         92: (91.905040099999994, 0.17150000000000001),
                                                                         94: (93.906315800000002, 0.17380000000000001)},
                  valence=4),
}


class Formula(dict):
    """Modelcular formula"""

    ELEMENT_PATTERN = re.compile(r'''
            ([A-Z][a-z]{0,2})
            ([\-]?[\d]*)
    ''', re.X)

    EMASS_PATH = op.normcase('third_party/emass/emass.exe -i third_party/emass/ISOTOPE.DAT')

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)

    @staticmethod
    def from_str(fstr):
        """ utility function, create formula from a string
        :param fstr:
        """
        assert (isinstance(fstr, str))
        v = Formula.ELEMENT_PATTERN.findall(fstr)
        if not v:
            logging.warn("provided string does not match a molecular formula")
            return None

        return Formula({elem: (int(nb) if nb else 1) for elem, nb in v})

    def __str__(self):
        return "".join(["".join((k, (str(v) if v > 1 else '')))
                        for k, v in sorted(self.iteritems(), key=lambda _: _[0])])

    def mono_mass(self):
        """return the mono mass of the formula"""
        return sum(ELEMENTS[elem].mass[0] * nb for elem, nb in self.iteritems())

    @staticmethod
    def _check_input(f):
        """
        :param f:
        :return:
        """
        if isinstance(f, str):
            fd = Formula.from_str(f)
        elif isinstance(f, dict) or isinstance(f, Counter):
            fd = Formula(f)
        else:
            if not isinstance(f, Formula):
                raise TypeError('provided argument must be str, dict, or formula')
            fd = f
        return fd

    def __add__(self, other):
        return self.add(other, new_obj=True)

    def __iadd__(self, other):
        return self.add(other, new_obj=False)

    def __sub__(self, other):
        return self.remove(other, new_obj=True)

    def __isub__(self, other):
        return self.remove(other, new_obj=False)

    def add(self, f, new_obj=False):
        """
        :param new_obj:
        :param f: dict key Element, value number of Element
        :param n: add or remove n Element
        :return:
        """
        fd = self._check_input(f)
        # set the working formula obj
        #shallow copy if needed
        wd = Formula(self) if new_obj else self

        for elem, nb in fd.iteritems():
            wd[elem] = wd.get(elem, 0) + nb
        return wd

    def remove(self, f, new_obj=False):
        """

        :param new_obj:
        :param f:
        :param Element:
        :param n:
        :return:
        """
        fd = self._check_input(f)

        wd = Formula(self) if new_obj else self

        for elem, nb in fd.iteritems():
            if elem in wd:
                nb_e = wd[elem] - nb
                if nb_e <= 0:
                    del wd[elem]
                else:
                    wd[elem] = nb_e
        return wd

        # def get_theo_ip(self, min_rel_int=5.0):
        # """
        #     # todo ask wich adducts to pass in parameter
        #     formula is a string meaning compound
        #     :param min_rel_int:
        #     """
        #     p = subprocess.Popen(Formula.EMASS_PATH, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        #
        #     out, err = p.communicate(input=str(self))
        #     if not out:
        #         logging.warn("Error computing isotopic pattern with formula: {0}.Skip it".format(str(self)))
        #         return
        #
        #     try:
        #         iso = repr(filter(lambda x: x[1] > min_rel_int,
        #                           [(lambda x: (float(x[0]), float(x[1])))(l.rstrip().split(" "))
        #                            for l in out.split('\n')[1:-1]]))
        #     except IndexError:
        #         logging.warn("Error parsing isotopic pattern.Skip it")
        #         return
        #
        #     return iso
