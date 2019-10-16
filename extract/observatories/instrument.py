class Instrument(object):
    """ Instrument class. """

    def __init__(self, _instrument):
        self.instrument = _instrument
        self.latitude = None
        self.longitude = None
        self.elevation = None
        self.blocked_regions = None

    def __repr__(self):
        return 'Instrument object.'
