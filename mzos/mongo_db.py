from __future__ import absolute_import
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
