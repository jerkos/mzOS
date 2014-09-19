__author__ = 'Marc'


class RPrecord(object):
    def __init__(self):
        """
        @type self:
        """
        self.as_r = set()
        self.as_p = set()

    def get(self):
        return self.as_r, self.as_p