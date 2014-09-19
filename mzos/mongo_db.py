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

from pymongo.connection import Connection
import humongolus as orm
import humongolus.field as field


class Annotation(orm.EmbeddedDocument):
    """
    metabolite name with its 2 scores
    """
    annotation = field.Char()
    score1 = field.Float()
    score2 = field.Float()


class Abundance(orm.EmbeddedDocument):
    """
    Custom abundances class
    """
    sample = field.Char()
    abundance = field.Float()


class Feature(orm.Document):
    """
    Feature class
    """
    _db = "os_mongo"
    _collection = "features"
    feature_id = field.AutoIncrement(collection="experiment")
    experiment_id = field.Integer(required=True)
    mass = field.Float(required=True)
    rt = field.Float(required=True)
    abundances = orm.List(type=Abundance)
    main_attribution = field.Char()
    annotations = orm.List(type=Annotation)


class Experiment(orm.Document):
    """
    Base class for various experiments
    """
    _db = "os_mongo"
    _collection = "experiments"
    experiment_id = field.AutoIncrement(collection="experiment")
    organization = field.Char(required=True)
    title = field.Char(required=True)
    description = field.Char(required=False)
    date = field.Date(required=True)


class MetabolomicsExperiment(Experiment):
    """
    Metabolomics experiments
    """
    software = field.Char(required=True)
    version = field.Char(required=True)
    parameters = field.Char()
    filenames = orm.List(type=str)
