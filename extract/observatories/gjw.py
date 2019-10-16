from extract.observatories.instrument import Instrument


class GJWNetwork(object):
    """ GJW telescope network class. """

    def __init__(self):
        self.telescope = 'GJWNetwork'

    def __repr__(self):
        return 'GJW library of instruments.'

    @property
    def gjw_cl_spec_red(self):
        """ GJW-CL Red Optical Spectrograph instrument. """
        # Instantiate instrument object.
        instrument = Instrument('gjw_cl_spec_red')
        instrument.latitude = None
        instrument.longitude = None
        instrument.elevation = None
        instrument.blocked_regions = {
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

    @property
    def gjw_sa_spec_red(self):
        """ GJW-SA Red Optical Spectrograph instrument. """
        # Instantiate instrument object.
        instrument = Instrument('gjw_sa_spec_red')
        instrument.latitude = None
        instrument.longitude = None
        instrument.elevation = None
        instrument.blocked_regions = {
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

    @property
    def gjw_in_spec_red(self):
        """ GJW-IN Red Optical Spectrograph instrument. """
        # Instantiate instrument object.
        instrument = Instrument('gjw_in_spec_red')
        instrument.latitude = None
        instrument.longitude = None
        instrument.elevation = None
        instrument.blocked_regions = {
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

    @property
    def gjw_wa_spec_red(self):
        """ GJW-WA Red Optical Spectrograph instrument. """
        # Instantiate instrument object.
        instrument = Instrument('gjw_wa_spec_red')
        instrument.latitude = None
        instrument.longitude = None
        instrument.elevation = None
        instrument.blocked_regions = {
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

    @property
    def gjw_oz_spec_red(self):
        """ GJW-OZ Red Optical Spectrograph instrument. """
        # Instantiate instrument object.
        instrument = Instrument('gjw_oz_spec_red')
        instrument.latitude = None
        instrument.longitude = None
        instrument.elevation = None
        instrument.blocked_regions = {
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


GJWNetworkInstruments = GJWNetwork()
