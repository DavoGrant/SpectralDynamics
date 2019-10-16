from extract.observatories.instrument import Instrument


class GeminiTelescope(object):
    """ Gemini telescope class. """

    def __init__(self):
        self.telescope = 'Gemini'

    def __repr__(self):
        return 'Gemini library of instruments.'

    @property
    def gmos(self):
        """ GMOS instrument. """
        # Instantiate instrument object.
        instrument = Instrument('gemini_gmos')
        instrument.latitude = None
        instrument.longitude = None
        instrument.elevation = None
        instrument.blocked_regions = {
            'Blue_one': (3662.0, 4003.0),
            'Blue_two_I': (4015.0, 4188.0),
            'Blue_two_II': (4220.0, 4610.0),
            'Blue_three_I': (4685.0, 4827.0),
            'Blue_three_II': (4830.0, 4984.0),
            'Blue_four': (4987.7, 5603.0),
            'GMOS_artifact': (5658.0, 5661.0),
            'He_Na': (5864.0, 5932.0),
            'S_Fe': (6130.0, 6396.0),
            'H_Alpha': (6495.0, 6634.0),
            'He_6678': (6658.0, 6700.0),
            'He_7065': (7020.0, 7100.0),
            'Unknown1': (7115.0, 7160.0),
            'Unknown2': (7364.0, 7397.0),
            'telluric_6400': (6857.0, 6933.0),
            'telluric_7200': (7162.0, 7345.0),
            'telluric_7600': (7581.0, 7700.0)
        }

        return instrument


GeminiTelescopeInstruments = GeminiTelescope()
