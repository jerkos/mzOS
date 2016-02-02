# -*- coding: utf-8 -*-

import json
import os

REAC_FILE = os.path.normpath('mzos/ressources/reactions.json')
REACTIONS = json.load(open(REAC_FILE, 'rb'))
