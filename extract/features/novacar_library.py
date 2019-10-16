from extract.features.template import Template


class NovaCarSpectralTemplateLibrary(object):
    """ Nova Carinae's spectral template library class.

    Instantiate an object specifying the spectral feature,
    template version and solver to fully define which fitting
    routine to implement. Properties return a template object
    with all required attributes for parsing to line fitting.
    Most wavelengths in defined in air.

    Properties
    ----------

        h_i_6563 : H-alpha at 6562.8A line.

    """

    def __init__(self):
        self.verbose = True
        self.target = 'NovaCar2018'
        self.sky_coord = '10 36 13.71 -59 35 55.1'

    def __repr__(self):
        return 'Nova Carinae spectral template library.'

    def h_i_6563(self, version, solver):
        """ H I at ~6562.8A. """
        # Instantiate template object.
        template = Template(self.target, self.sky_coord,
                            'h_i_6563', version, solver)

        if version == 'BHM' and solver == 'BHM':
            template.window = (6520.0, 6607.0)
            template.guess = 0.5  # Fraction of max.

        elif version == 'Simple' and solver == 'CF':
            template.window = (6520.0, 6607.0)
            template.n_components = 1
            template.component_number = [1]
            template.guess = [6557, 10, 0]
            template.priors = ((6550, 0, -20),
                               (6570, 30, 20))

        else:
            # Version and solver combination not defined.
            raise ValueError('Unknown version and solver.')

        return template


novacar_templates = NovaCarSpectralTemplateLibrary()
