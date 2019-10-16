class Template(object):
    """ Spectral template class. """

    def __init__(self, _target, _sky_coord, _feature, _version, _solver):
        # Args.
        self.target = _target
        self.sky_coord = _sky_coord
        self.feature = _feature
        self.template_version = _version
        self.solver = _solver

        # Templates.
        self.window = None
        self.n_components = None
        self.component_number = None
        self.guess = None
        self.priors = None

    def __repr__(self):
        return 'Spectral template for {} {} feature between {}A ' \
               'and {}A: template_version={}, solver={}.'.format(
                self.target, self.feature, self.window[0],
                self.window[1], self.template_version, self.solver)
