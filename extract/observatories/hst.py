from extract.observatories.instrument import Instrument


class HubbleSpaceTelescope(object):
    """ Hubble space telescope class. """

    def __init__(self):
        self.telescope = 'HST'

    def __repr__(self):
        return 'Hubble space telescope library of instruments.'

    @property
    def stis(self):
        """ STIS instrument. """
        # Instantiate instrument object.
        instrument = Instrument('hst_stis')
        instrument.latitude = None
        instrument.longitude = None
        instrument.elevation = None
        instrument.blocked_regions = {
            'H_Zeta': (3878.0, 3909.0),
            'H_Epsilon': (3957.0, 3987.0),
            'H_Delta': (4089.0, 4135.0),
            'Blue_one': (4165.0, 4185.0),
            'Blue_two': (4225.0, 4253.0),
            'H_Beta': (4833.0, 4888.0),
            'Blue_three': (4913.0, 4934.0),
            'Blue_four': (4995.0, 5032.0),
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


HubbleSpaceTelescopeInstruments = HubbleSpaceTelescope()
