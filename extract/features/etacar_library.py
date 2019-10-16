from extract.features.template import Template


class EtaCarSpectralTemplateLibrary(object):
    """ Eta Carinae's spectral template library class.

    Instantiate an object specifying the spectral feature,
    template version and solver to fully define which fitting
    routine to implement. Properties return a template object
    with all required attributes for parsing to line fitting.
    Most wavelengths in defined in air.

    Properties
    ----------

        h_i_3750 : H-kappa at 3750.2A emission line.

        h_i_3771 : H-iota at 3770.6A emission line.

        h_i_3798 : H-theta at 3797.9A emission line.

        h_i_3835 : H-eta at 3835.4A emission line.

        h_i_3889 : H-zeta at 3889.0A emission line.

        h_i_3970 : H-epsilon at 3970.1A emission line.

        h_i_4101 : H-delta at 4101.7A emission line.

        h_i_4340 : H-gamma at 4340.5A emission line.

        h_i_4861 : H-beta at 4861.3A emission line.

        h_i_6563 : H-alpha at 6562.8A emission line.

    """

    def __init__(self):
        self.verbose = True
        self.target = 'EtaCar'
        self.sky_coord = '10 45 03.5362075818 -59 41 04.053436648'

    def __repr__(self):
        return 'Eta Carinae spectral template library.'

    def h_i_3750(self, version, solver):
        """ H I kappa at ~3750.2A. """
        # Instantiate template object.
        template = Template(self.target, self.sky_coord,
                            'h_i_3750', version, solver)

        if version == 'BHM' and solver == 'BHM':
            template.window = (3740.0, 3758.0)
            template.guess = 0.5  # Fraction of max.

        elif version == 'Harvard_GMOS_C1' and solver == 'CF':
            template.window = (3740.0, 3758.0)
            template.n_components = 6
            template.component_number = [1, 2, 3, 4, -1, -2]
            template.guess = [3749, 0.7, 0.3,
                              3751, 0.7, 0.2,
                              3753, 0.7, 0.2,
                              3748, 0.7, 0.2,
                              3743.5, 0.5, -0.05,
                              3745.5, 0.5, -0.08]
            template.priors = ((3747, 0, 0,
                                3749, 0, 0,
                                3751, 0, 0,
                                3746.5, 0, 0,
                                3741.5, 0, -0.4,
                                3743.5, 0, -0.4),
                               (3751, 1.2, 0.7,
                                3753, 1.2, 0.7,
                                3755, 1.2, 0.7,
                                3750, 1, 0.4,
                                3745, 1.5, 0,
                                3746.3, 1.2, 0))

        elif version == 'Harvard_GMOS_C2' and solver == 'CF':
            template.window = (3740.0, 3758.0)
            template.n_components = 6
            template.component_number = [1, 2, 3, 4, -1, -2]
            template.guess = [3750, 0.7, 0.3,
                              3751, 0.7, 0.2,
                              3753, 0.7, 0.2,
                              3749, 0.7, 0.2,
                              3743.5, 0.5, -0.05,
                              3745.5, 0.5, -0.08]
            template.priors = ((3749.5, 0, 0,
                                3750, 0, 0,
                                3751, 0, 0,
                                3748.5, 0, 0,
                                3741.5, 0, -0.4,
                                3743.5, 0, -0.4),
                               (3752, 1.2, 0.7,
                                3753, 1.2, 0.7,
                                3755, 1.2, 0.7,
                                3751, 1, 0.4,
                                3745, 1.5, 0,
                                3746.3, 1.2, 0))

        else:
            # Version and solver combination not defined.
            raise ValueError('Unknown version and solver.')

        return template

    def h_i_3771(self, version, solver):
        """ H I iota at ~3770.6A. """
        # Instantiate template object.
        template = Template(self.target, self.sky_coord,
                            'h_i_3771', version, solver)

        if version == 'BHM' and solver == 'BHM':
            template.window = (3759.0, 3779.0)
            template.guess = 0.5  # Fraction of max.

        elif version == 'Harvard_GMOS_C1' and solver == 'CF':
            template.window = (3759.0, 3781.0)
            template.n_components = 8
            template.component_number = [1, 2, 3, 4, 20, 22, -1, -2]
            template.guess = [3769, 0.7, 0.3,
                              3771, 0.7, 0.2,
                              3773, 0.7, 0.2,
                              3768, 0.7, 0.2,
                              3760.7, 0.4, 0.03,
                              3763.5, 0.4, 0.03,
                              3762, 0.5, -0.05,
                              3765.5, 0.5, -0.08]
            template.priors = ((3767, 0, 0,
                                3769, 0, 0,
                                3771, 0, 0,
                                3766.5, 0, 0,
                                3760.3, 0, 0,
                                3763.1, 0, 0,
                                3761.5, 0, -0.4,
                                3763.5, 0, -0.4),
                               (3771, 1.2, 0.7,
                                3773, 1.2, 0.7,
                                3775, 1.2, 0.7,
                                3770.2, 1, 0.4,
                                3761.1, 0.5, 0.08,
                                3764, 0.5, 0.08,
                                3765.5, 1.5, 0,
                                3768, 1.5, 0))

        elif version == 'Harvard_GMOS_C2' and solver == 'CF':
            template.window = (3759.0, 3781.0)
            template.n_components = 6
            template.component_number = [1, 2, 3, 4, -1, -2]
            template.guess = [3769, 0.7, 0.3,
                              3771, 0.7, 0.2,
                              3773, 0.7, 0.2,
                              3768, 0.7, 0.2,
                              3763.5, 0.5, -0.05,
                              3765.5, 0.5, -0.08]
            template.priors = ((3767, 0, 0,
                                3769, 0, 0,
                                3771, 0, 0,
                                3766.5, 0, 0,
                                3761.5, 0, -0.4,
                                3763.5, 0, -0.4),
                               (3771, 1.2, 0.7,
                                3773, 1.2, 0.7,
                                3775, 1.2, 0.7,
                                3770, 1, 0.4,
                                3765.5, 1.5, 0,
                                3768, 1.5, 0))

        elif version == 'Harvard_GMOS_C3' and solver == 'CF':
            template.window = (3759.0, 3781.0)
            template.n_components = 8
            template.component_number = [1, 2, 3, 4, 20, 22, -1, -2]
            template.guess = [3769, 0.7, 0.3,
                              3771, 0.7, 0.2,
                              3773, 0.7, 0.2,
                              3768, 0.7, 0.2,
                              3760.7, 0.4, 0.03,
                              3763.5, 0.4, 0.03,
                              3764, 0.5, -0.05,
                              3765.5, 0.5, -0.08]
            template.priors = ((3767, 0, 0,
                                3769, 0, 0,
                                3771, 0, 0,
                                3766.5, 0, 0,
                                3760.3, 0, 0,
                                3763.1, 0, 0,
                                3761.5, 0, -0.5,
                                3763.5, 0, -0.5),
                               (3771, 1.2, 0.7,
                                3773, 1.2, 0.7,
                                3775, 1.2, 0.7,
                                3770.2, 1, 0.4,
                                3761.1, 0.5, 0.08,
                                3764, 0.5, 0.08,
                                3767, 1.5, 0,
                                3768, 1.5, 0))

        else:
            # Version and solver combination not defined.
            raise ValueError('Unknown version and solver.')

        return template

    def h_i_3798(self, version, solver):
        """ H I theta at ~3797.9A. """
        # Instantiate template object.
        template = Template(self.target, self.sky_coord,
                            'h_i_3798', version, solver)

        if version == 'BHM' and solver == 'BHM':
            template.window = (3784.0, 3810.0)
            template.guess = 0.5  # Fraction of max.

        elif version == 'Harvard_GMOS_C1' and solver == 'CF':
            template.window = (3784.0, 3810.0)
            template.n_components = 10
            template.component_number = [1, 2, 3, 4, 7, 8, 21, -1, -3, -2]
            template.guess = [3797, 0.7, 0.3,
                              3799, 0.7, 0.2,
                              3801, 0.7, 0.2,
                              3795, 0.7, 0.2,
                              3803, 0.7, 0.2,
                              3805, 0.7, 0.2,
                              3806, 0.5, 0.05,
                              3788, 0.5, -0.02,
                              3791.5, 0.5, -0.05,
                              3793.5, 0.5, -0.08]
            template.priors = ((3795, 0, 0,
                                3797, 0, 0,
                                3799, 0, 0,
                                3794, 0, 0,
                                3802.5, 0, 0,
                                3803.5, 0, 0,
                                3805, 0, 0,
                                3786, 0, -0.2,
                                3789, 0, -0.4,
                                3791, 0, -0.4),
                               (3799, 1.4, 0.7,
                                3801, 1.4, 0.7,
                                3802.5, 1.4, 0.7,
                                3797, 1.4, 0.4,
                                3805, 1.4, 0.7,
                                3806, 1.4, 0.7,
                                3807, 0.8, 0.2,
                                3791, 1.5, 0,
                                3793, 1.5, 0,
                                3794, 1.5, 0))

        else:
            # Version and solver combination not defined.
            raise ValueError('Unknown version and solver.')

        return template

    def h_i_3835(self, version, solver):
        """ H I eta at ~3835.4A. """
        # Instantiate template object.
        template = Template(self.target, self.sky_coord,
                            'h_i_3835', version, solver)

        if version == 'BHM' and solver == 'BHM':
            template.window = (3824.0, 3844.0)
            template.guess = 0.5  # Fraction of max.

        elif version == 'Harvard_GMOS_C1' and solver == 'CF':
            template.window = (3824.0, 3844.0)
            template.n_components = 8
            template.component_number = [1, 2, 3, 4, 21, -1, -3, -2]
            template.guess = [3835, 1, 0.7,
                              3836.5, 1, 0.5,
                              3838, 1, 0.5,
                              3833, 1, 0.3,
                              3824, 0.7, 0.05,
                              3827.5, 0.5, -0.05,
                              3830, 0.5, -0.1,
                              3831, 0.5, -0.1]
            template.priors = ((3833, 0, 0,
                                3835, 0, 0,
                                3836, 0, 0,
                                3831, 0, 0,
                                3823, 0, 0,
                                3825, 0, -0.3,
                                3827, 0, -0.3,
                                3828, 0, -0.5),
                               (3837, 1.5, 1,
                                3838, 1.5, 1,
                                3841, 1.5, 1,
                                3835, 1.5, 0.5,
                                3825, 1, 0.2,
                                3830, 2, 0,
                                3831, 2, 0,
                                3832, 2, 0))

        elif version == 'Harvard_GMOS_S1' and solver == 'CF':
            template.window = (3824.0, 3844.0)
            template.n_components = 8
            template.component_number = [1, 2, 3, 4, 20, -1, -3, -2]
            template.guess = [3835, 1, 0.7,
                              3836.5, 1, 0.5,
                              3838, 1, 0.5,
                              3833.5, 0.5, 0.1,
                              3824, 0.7, 0.05,
                              3827.5, 0.5, -0.05,
                              3830, 0.5, -0.1,
                              3831, 0.5, -0.1]
            template.priors = ((3833, 0, 0,
                                3835, 0, 0,
                                3836, 0, 0,
                                3831, 0, 0,
                                3823, 0, 0,
                                3825, 0, -0.3,
                                3827, 0, -0.3,
                                3828, 0, -0.5),
                               (3837, 1.5, 1,
                                3838, 1.5, 1,
                                3841, 1.5, 1,
                                3835, 0.9, 0.5,
                                3825, 1, 0.2,
                                3830, 2, 0,
                                3831, 2, 0,
                                3832, 2, 0))

        else:
            # Version and solver combination not defined.
            raise ValueError('Unknown version and solver.')

        return template

    def h_i_3889(self, version, solver):
        """ H I zeta at ~3889.0A. """
        # Instantiate template object.
        template = Template(self.target, self.sky_coord,
                            'h_i_3889', version, solver)

        if version == 'BHM' and solver == 'BHM':
            template.window = (3870, 3906)
            template.guess = 0.5  # Fraction of max.

        elif version == 'Harvard_GMOS_C1' and solver == 'CF':
            template.window = (3870, 3906)
            template.n_components = 12
            template.component_number = [1, 2, 3, 4, 5, 7, 8, -5, -1, -4, -2, -3]
            template.guess = [3887, 1, 1,
                              3889, 1, 0.7,
                              3891, 1, 0.7,
                              3885, 1, 0.4,
                              3893, 3, 0.1,
                              3899.4, 0.7, 0.03,
                              3902.2, 0.7, 0.03,
                              3876, 0.5, -0.1,
                              3879, 0.5, -0.2,
                              3881, 0.5, -0.2,
                              3883, 0.5, -0.2,
                              3886.2, 0.2, -0.3]
            template.priors = ((3885, 0, 0,
                                3887, 0, 0,
                                3889, 0, 0,
                                3884, 0, 0,
                                3887, 0, 0,
                                3898.5, 0, 0,
                                3901, 0, 0,
                                3874, 0, -0.3,
                                3877, 0, -0.8,
                                3879, 0, -0.8,
                                3881, 0, -0.8,
                                3885.7, 0, -0.5),
                               (3889, 3, 3,
                                3891, 2, 3,
                                3893, 3, 3,
                                3887, 3, 1,
                                3898, 10, 0.5,
                                3900.5, 1, 0.2,
                                3903.5, 1, 0.2,
                                3878, 2, 0,
                                3881, 2, 0,
                                3883, 1.25, 0,
                                3885, 1.25, 0,
                                3886.7, 0.5, -0.2))

        elif version == 'Harvard_GMOS_C2' and solver == 'CF':
            template.window = (3870, 3906)
            template.n_components = 11
            template.component_number = [1, 2, 3, 5, 7, 8, -5, -1, -4, -2, -3]
            template.guess = [3887, 1, 1,
                              3889, 1, 0.7,
                              3891, 1, 0.7,
                              3893, 3, 0.1,
                              3899.4, 0.7, 0.03,
                              3902.2, 0.7, 0.03,
                              3876, 0.5, -0.1,
                              3879, 0.5, -0.2,
                              3881, 0.5, -0.2,
                              3883, 0.5, -0.2,
                              3886.2, 0.2, -0.3]
            template.priors = ((3885, 0, 0,
                                3887, 0, 0,
                                3889, 0, 0,
                                3887, 0, 0,
                                3898.5, 0, 0,
                                3901, 0, 0,
                                3874, 0, -0.3,
                                3877, 0, -0.8,
                                3879, 0, -0.8,
                                3881, 0, -0.8,
                                3885.7, 0, -0.5),
                               (3889, 3, 3,
                                3891, 3, 3,
                                3893, 3, 3,
                                3898, 10, 0.5,
                                3900.5, 1, 0.2,
                                3903.5, 1, 0.2,
                                3878, 2, 0,
                                3881, 2, 0,
                                3883, 1.25, 0,
                                3885, 1.25, 0,
                                3886.7, 0.5, 0))

        else:
            # Version and solver combination not defined.
            raise ValueError('Unknown version and solver.')

        return template

    def h_i_3970(self, version, solver):
        """ H I epsilon at ~3970.1A. """
        # Instantiate template object.
        template = Template(self.target, self.sky_coord,
                            'h_i_3970', version, solver)

        if version == 'BHM' and solver == 'BHM':
            template.window = (3953.0, 3988.0)
            template.guess = 0.5  # Fraction of max.

        elif version == 'Harvard_GMOS_C1' and solver == 'CF':
            template.window = (3953.0, 3988.0)
            template.n_components = 10
            template.component_number = [1, 2, 3, 4, 5, 6, -1, -4, -2, -3]
            template.guess = [3969, 1, 1.5,
                              3971, 1, 1,
                              3973, 1, 1,
                              3966, 1, 1,
                              3971, 6, 0.3,
                              3981, 0.7, 0.02,
                              3959.5, 0.5, -0.2,
                              3960.5, 0.5, -0.2,
                              3962, 0.5, -0.2,
                              3967.9, 0.2, -0.3]
            template.priors = ((3967, 0, 0,
                                3969.5, 0, 0,
                                3971, 0, 0,
                                3965, 0, 0,
                                3967, 0, 0,
                                3980, 0, 0,
                                3957, 0, -0.5,
                                3958.5, 0, -0.5,
                                3960, 0, -0.5,
                                3967.3, 0, -0.5),
                               (3972, 3, 3,
                                3973, 3, 3,
                                3975, 3, 3,
                                3968, 3, 3,
                                3976, 10, 0.6,
                                3982, 1, 0.07,
                                3960, 2, 0,
                                3962, 2, 0,
                                3964, 2, 0,
                                3968.3, 0.5, 0))

        elif version == 'Harvard_GMOS_C2' and solver == 'CF':
            template.window = (3953.0, 3988.0)
            template.n_components = 10
            template.component_number = [1, 2, 3, 4, 5, 6, -1, -4, -2, -3]
            template.guess = [3969, 1, 1.5,
                              3971, 1, 1,
                              3973, 1, 1,
                              3966, 1, 1,
                              3971, 6, 0.3,
                              3981, 0.7, 0.02,
                              3959.5, 0.5, -0.2,
                              3963, 0.5, -0.2,
                              3964, 0.5, -0.2,
                              3967.9, 0.3, -0.3]
            template.priors = ((3967, 0, 0,
                                3969.5, 0, 0,
                                3971, 0, 0,
                                3965, 0, 0,
                                3967, 0, 0,
                                3980, 0, 0,
                                3957, 0, -0.5,
                                3958.5, 0, -0.5,
                                3960, 0, -0.5,
                                3967.3, 0.2, -0.5),
                               (3972, 3, 3,
                                3973, 3, 3,
                                3975, 3, 3,
                                3968, 3, 3,
                                3976, 10, 0.6,
                                3982, 1, 0.07,
                                3960, 2, 0,
                                3964, 2, 0,
                                3966, 2, 0,
                                3968.7, 0.5, 0))

        else:
            # Version and solver combination not defined.
            raise ValueError('Unknown version and solver.')

        return template

    def h_i_4101(self, version, solver):
        """ H I delta at ~4101.7A. """
        # Instantiate template object.
        template = Template(self.target, self.sky_coord,
                            'h_i_4101', version, solver)

        if version == 'BHM' and solver == 'BHM':
            template.window = (4087.0, 4111.0)
            template.guess = 0.5  # Fraction of max.

        elif version == 'Harvard_GMOS_A1' and solver == 'CF':
            template.window = (4087.0, 4111.0)
            template.n_components = 6
            template.component_number = [1, 2, 3, 4, 5, 6]
            template.guess = [4100, 1.5, 1.5,
                              4103, 1.5, 1.5,
                              4105, 1.5, 0.5,
                              4097, 1.5, 0.5,
                              4108, 8, 0.5,
                              4101.5, 0.4, 0.2]
            template.priors = ((4098, 0, 0,
                                4101, 0, 0,
                                4103, 0, 0,
                                4096, 0, 0,
                                4106, 7, 0.2,
                                4101, 0, 0),
                               (4102, 3, 5,
                                4105, 3, 5,
                                4107, 3, 5,
                                4099, 3, 5,
                                4112, 15, 1,
                                4102, 0.6, 0.5))

        elif version == 'Harvard_GMOS_P1' and solver == 'CF':
            template.window = (4087.0, 4111.0)
            template.n_components = 8
            template.component_number = [1, 2, 3, 4, 5, 6, -1, -2]
            template.guess = [4100, 1.5, 1.5,
                              4103, 1.5, 1.5,
                              4105, 1.5, 0.5,
                              4099, 1.5, 0.5,
                              4108, 8, 0.5,
                              4101.5, 0.4, 0.2,
                              4092.5, 1, -0.1,
                              4095, 0.5, -0.2]
            template.priors = ((4098, 0, 0,
                                4101, 0, 0,
                                4103, 0, 0,
                                4096, 0, 0,
                                4106, 7, 0.2,
                                4101, 0, 0,
                                4091, 0, -0.6,
                                4093, 0, -0.6),
                               (4102, 3, 5,
                                4105, 3, 5,
                                4107, 3, 5,
                                4099.5, 3, 5,
                                4112, 15, 1,
                                4102, 0.6, 0.5,
                                4094, 2, 0,
                                4097, 1.5, 0))

        elif version == 'Harvard_GMOS_P2' and solver == 'CF':
            template.window = (4087.0, 4111.0)
            template.n_components = 8
            template.component_number = [1, 2, 3, 4, 5, 6, -1, -2]
            template.guess = [4100, 1.5, 1.5,
                              4103, 1.5, 1.5,
                              4105, 1.5, 0.5,
                              4099, 1.5, 0.5,
                              4108, 8, 0.5,
                              4101.5, 0.4, 0.2,
                              4095.5, 0.3, -0.2,
                              4097.5, 0.3, -0.2]
            template.priors = ((4098, 0, 0,
                                4101, 0, 0,
                                4103, 0, 0,
                                4096, 0, 0,
                                4106, 7, 0.2,
                                4101, 0, 0,
                                4091, 0, -0.3,
                                4093, 0, -0.3),
                               (4102, 3, 5,
                                4105, 3, 5,
                                4107, 3, 5,
                                4099.5, 3, 5,
                                4112, 15, 1,
                                4102, 0.6, 0.5,
                                4096, 1.5, 0,
                                4098, 1.5, 0))

        elif version == 'Harvard_GMOS_S1' and solver == 'CF':
            template.window = (4087.0, 4111.0)
            template.n_components = 7
            template.component_number = [1, 2, 3, 4, 5, 6, -2]
            template.guess = [4100, 1.5, 1.5,
                              4103, 1.5, 1.5,
                              4105, 1.5, 0.5,
                              4097, 1.5, 0.5,
                              4108, 8, 0.5,
                              4101.5, 0.4, 0.2,
                              4095, 0.5, -0.1]
            template.priors = ((4098, 0, 0,
                                4101, 0, 0,
                                4103, 0, 0,
                                4096, 0, 0,
                                4106, 7, 0.2,
                                4101, 0, 0,
                                4093, 0, -0.2),
                               (4102, 3, 5,
                                4105, 3, 5,
                                4107, 3, 5,
                                4099, 3, 5,
                                4112, 15, 1,
                                4102, 0.6, 0.5,
                                4097, 1, 0))

        elif version == 'Harvard_GMOS_S2' and solver == 'CF':
            template.window = (4087.0, 4111.0)
            template.n_components = 8
            template.component_number = [1, 2, 3, 4, 5, 6, -1, -2]
            template.guess = [4100, 1.5, 1.5,
                              4103, 1.5, 1.5,
                              4105, 1.5, 0.5,
                              4099, 1.5, 0.5,
                              4108, 8, 0.5,
                              4101.5, 0.4, 0.2,
                              4092.5, 1, -0.1,
                              4094, 0.5, -0.1]
            template.priors = ((4098, 0, 0,
                                4101, 0, 0,
                                4103, 0, 0,
                                4096, 0, 0,
                                4106, 7, 0.2,
                                4101, 0, 0,
                                4091, 0, -0.6,
                                4093, 0, -0.1),
                               (4102, 3, 5,
                                4105, 3, 5,
                                4107, 3, 5,
                                4099.5, 3, 5,
                                4112, 15, 1,
                                4102, 0.6, 0.5,
                                4094, 2, 0,
                                4097, 1.5, 0))

        else:
            # Version and solver combination not defined.
            raise ValueError('Unknown version and solver.')

        return template

    def h_i_4340(self, version, solver):
        """ H I gamma at ~4340.5A. """
        # Instantiate template object.
        template = Template(self.target, self.sky_coord,
                            'h_i_4340', version, solver)

        if version == 'BHM' and solver == 'BHM':
            template.window = (4322.0, 4349.0)
            template.guess = 0.5  # Fraction of max.

        elif version == 'Harvard_GMOS_A1' and solver == 'CF':
            template.window = (4322.0, 4349.0)
            template.n_components = 9
            template.component_number = [1, 2, 3, 4, 55, 5, 6, 20, 8]
            template.guess = [4338, 1, 3,
                              4340, 1, 3,
                              4343, 1, 3,
                              4334, 0.7, 3,
                              4320, 5, 0.2,
                              4385, 25, 3,
                              4340, 0.4, 0.2,
                              4325.5, 0.5, 0.5,
                              4345, 1, 0.5]
            template.priors = ((4336, 0, 0,
                                4338, 0, 0,
                                4341, 0, 0,
                                4333, 0, 0,
                                4315, 3, 0.1,
                                4370, 20, 1.8,
                                4339.5, 0, 0,
                                4324.5, 0, 0,
                                4344, 0, 0),
                               (4341, 2, 10,
                                4343, 2, 10,
                                4345, 2, 10,
                                4340, 1, 10,
                                4325, 7, 0.5,
                                4395, 30, 4,
                                4340.5, 0.7, 0.5,
                                4326.5, 1, 1,
                                4347, 2, 1))

        elif version == 'Harvard_GMOS_A2' and solver == 'CF':
            template.window = (4322.0, 4349.0)
            template.n_components = 8
            template.component_number = [1, 2, 3, 4, 55, 5, 20, 8]
            template.guess = [4338, 1, 3,
                              4340, 1, 3,
                              4343, 1, 3,
                              4334, 0.7, 3,
                              4320, 5, 0.2,
                              4385, 25, 3,
                              4325.5, 0.5, 0.5,
                              4345, 1, 0.5]
            template.priors = ((4336, 0, 0,
                                4338, 0, 0,
                                4341, 0, 0,
                                4333, 0, 0,
                                4315, 3, 0.1,
                                4370, 20, 1.8,
                                4324.5, 0, 0,
                                4344, 0, 0),
                               (4341, 2, 10,
                                4343, 2, 10,
                                4344, 1.5, 10,
                                4340, 1, 10,
                                4325, 7, 0.5,
                                4395, 30, 4,
                                4326.5, 1, 1,
                                4347, 2, 1))

        elif version == 'Harvard_GMOS_P1' and solver == 'CF':
            template.window = (4322.0, 4349.0)
            template.n_components = 10
            template.component_number = [1, 2, 3, 4, 55, 5, 20, 8, -1, -2]
            template.guess = [4338, 1, 3,
                              4340, 1, 3,
                              4343, 1, 3,
                              4335, 0.7, 3,
                              4320, 5, 0.2,
                              4385, 25, 3,
                              4325.5, 0.5, 0.5,
                              4345, 1, 0.5,
                              4330, 1, -0.1,
                              4333, 1, -0.2]
            template.priors = ((4336, 0, 0,
                                4338, 0, 0,
                                4341, 0, 0,
                                4333, 0, 0,
                                4315, 3, 0.1,
                                4370, 20, 1.8,
                                4324.5, 0, 0,
                                4344, 0, 0,
                                4328, 0, -1,
                                4331, 0, -1),
                               (4341, 2, 10,
                                4343, 2, 10,
                                4344, 1.5, 10,
                                4340, 1, 10,
                                4325, 7, 0.5,
                                4395, 30, 4,
                                4326.5, 1, 1,
                                4347, 2, 1,
                                4332, 2, 0,
                                4333.5, 1.5, 0))

        elif version == 'Harvard_GMOS_P2' and solver == 'CF':
            template.window = (4322.0, 4349.0)
            template.n_components = 10
            template.component_number = [1, 2, 3, 4, 55, 5, 20, 8, -1, -2]
            template.guess = [4338, 1, 3,
                              4340, 1, 3,
                              4343, 1, 3,
                              4335, 0.7, 3,
                              4320, 5, 0.2,
                              4385, 25, 3,
                              4325.5, 0.5, 0.5,
                              4345, 1, 0.5,
                              4330, 1, -0.1,
                              4333, 1, -0.2]
            template.priors = ((4336, 0, 0,
                                4338, 0, 0,
                                4341, 0, 0,
                                4333, 0, 0,
                                4315, 3, 0.1,
                                4370, 20, 1.8,
                                4324.5, 0, 0,
                                4344, 0, 0,
                                4328, 0, -1,
                                4331, 0, -1),
                               (4341, 2, 10,
                                4343, 2, 10,
                                4344, 1.5, 10,
                                4340, 1, 10,
                                4325, 7, 0.5,
                                4395, 30, 4,
                                4326.5, 1, 1,
                                4347, 2, 1,
                                4332, 2, 0,
                                4333.5, 1.5, 0))

        elif version == 'Harvard_GMOS_P3' and solver == 'CF':
            template.window = (4322.0, 4349.0)
            template.n_components = 9
            template.component_number = [1, 2, 3, 4, 55, 5, 20, 8, -2]
            template.guess = [4338, 1, 3,
                              4340, 1, 3,
                              4343, 1, 3,
                              4336, 0.7, 3,
                              4320, 5, 0.2,
                              4385, 25, 3,
                              4325.5, 0.5, 0.5,
                              4345, 1, 0.5,
                              4334, 1, -0.2]
            template.priors = ((4336, 0, 0,
                                4338, 0, 0,
                                4341, 0, 0,
                                4333, 0, 0,
                                4315, 3, 0.1,
                                4370, 20, 1.8,
                                4324.5, 0, 0,
                                4344, 0, 0,
                                4331, 0, -1),
                               (4341, 2, 10,
                                4343, 2, 10,
                                4344, 1.5, 10,
                                4340, 1, 10,
                                4325, 7, 0.5,
                                4395, 30, 4,
                                4326.5, 1, 1,
                                4347, 2, 1,
                                4335, 2, 0))

        else:
            # Version and solver combination not defined.
            raise ValueError('Unknown version and solver.')

        return template

    def h_i_4861(self, version, solver):
        """ H I beta at ~4861.3A. """
        # Instantiate template object.
        template = Template(self.target, self.sky_coord,
                            'h_i_4861', version, solver)

        if version == 'BHM' and solver == 'BHM':
            template.window = (4838.0, 4881.5)
            template.guess = 0.5  # Fraction of max.

        elif version == 'Harvard_GMOS_A1' and solver == 'CF':
            template.window = (4838.0, 4881.5)
            template.n_components = 10
            template.component_number = [1, 2, 3, 4, 5, 6, 20, 21, 7, 8]
            template.guess = [4859, 1.5, 7,
                              4862, 1.5, 7,
                              4865, 1.5, 2,
                              4856, 1.5, 2,
                              4865, 10, 2,
                              4861, 0.5, 1,
                              4847.5, 0.5, 0.3,
                              4873.5, 0.5, 0.3,
                              4867, 1, 0.3,
                              4852, 1, 0.3]
            template.priors = ((4857, 0, 0,
                                4860, 0, 0,
                                4863, 0, 0,
                                4854, 0, 0,
                                4860, 0, 0,
                                4860.5, 0, 0,
                                4846, 0, 0,
                                4872, 0, 0,
                                4865, 0, 0,
                                4850, 0, 0),
                               (4861, 3, 20,
                                4864, 3, 20,
                                4867, 2.5, 10,
                                4858, 3, 10,
                                4870, 20, 5,
                                4861.5, 0.7, 1.5,
                                4849, 1, 1,
                                4875, 1, 1,
                                4871, 1.5, 1,
                                4854, 1.5, 1))

        elif version == 'Harvard_GMOS_P1' and solver == 'CF':
            template.window = (4838.0, 4881.5)
            template.n_components = 11
            template.component_number = [1, 2, 3, 4, 5, 6, 20, 21, 7, -1, -2]
            template.guess = [4859, 1.5, 7,
                              4862, 1.5, 7,
                              4865, 1.5, 2,
                              4856, 1.5, 2,
                              4865, 10, 2,
                              4861, 0.5, .3,
                              4847.5, 0.5, 0.3,
                              4873.5, 0.5, 0.3,
                              4867, 1, 0.3,
                              4852, 1, -0.2,
                              4853, 1, -0.5]
            template.priors = ((4857, 0, 0,
                                4860, 0, 0,
                                4863, 0, 0,
                                4854, 0, 0,
                                4860, 0, 0,
                                4860.4, 0, 0,
                                4846, 0, 0,
                                4872, 0, 0,
                                4865, 0, 0,
                                4849, 0, -0.5,
                                4851, 0, -0.8),
                               (4861, 3, 20,
                                4864, 3, 20,
                                4867, 2.5, 10,
                                4858, 3, 10,
                                4870, 20, 5,
                                4861.5, 0.7, 1.5,
                                4849, 1, 1,
                                4875, 1, 1,
                                4871, 1.5, 1,
                                4853, 2, 0,
                                4855, 1.5, 0))

        else:
            # Version and solver combination not defined.
            raise ValueError('Unknown version and solver.')

        return template

    def h_i_6563(self, version, solver):
        """ H I at ~6562.8A. """
        # Instantiate template object.
        template = Template(self.target, self.sky_coord,
                            'h_i_6563', version, solver)

        if version == 'BHM' and solver == 'BHM':
            template.window = (6520.0, 6607.0)
            template.guess = 0.5  # Fraction of max.

        elif version == 'Harvard_GMOS_A1' and solver == 'CF':
            template.window = (6520.0, 6607.0)
            template.n_components = 14
            template.component_number = [1, 2, 3, 4, 5, 6, 20, 21, 7, 8, 9, 10, 11, 12]
            template.guess = [6562, 3, 15,
                              6565, 3, 15,
                              6570, 3, 15,
                              6555, 1, 2,
                              6565, 15, 5,
                              6562, 0.5, 10,
                              6546.5, 1, 0.5,
                              6582, 1, 2,
                              6540, 2, 0.5,
                              6587, 2, 0.5,
                              6599, 2, 0.5,
                              6550, 1, 0.5,
                              6576, 1, 0.5,
                              6580, 1, 0.5]
            template.priors = ((6557.6, 0, 0,
                                6563, 0, 0,
                                6567, 0, 0,
                                6553, 0, 0,
                                6560, 0, 0,
                                6560, 0, 0,
                                6546, 0, 0,
                                6581, 0, 0,
                                6538, 0, 0,
                                6582, 0, 0,
                                6592, 0, 0,
                                6548, 0, 0,
                                6574, 0, 0,
                                6578, 0, 0),
                               (6570, 5, 60,
                                6571, 5, 60,
                                6575, 5, 60,
                                6558, 3, 20,
                                6570, 30, 10,
                                6564, 1, 30,
                                6548, 1.4, 1,
                                6583, 1.4, 2,
                                6543, 3, 1,
                                6590, 3, 1,
                                6603, 3, 1,
                                6552, 2, 3,
                                6578, 2, 2,
                                6582, 2, 1))

        elif version == 'Harvard_GMOS_T1' and solver == 'CF':
            template.window = (6520.0, 6607.0)
            template.n_components = 15
            template.component_number = [1, 2, 3, 4, 5, 6, 20, 21, 7, 8, 9, 10, 11, 12, -3]
            template.guess = [6562, 3, 15,
                              6565, 3, 15,
                              6570, 3, 15,
                              6555, 1, 2,
                              6565, 15, 5,
                              6562, 0.5, 10,
                              6546.5, 1, 0.5,
                              6582, 1, 2,
                              6540, 2, 0.5,
                              6587, 2, 0.5,
                              6599, 2, 0.5,
                              6550, 1, 0.5,
                              6576, 1, 0.5,
                              6580, 1, 0.5,
                              6559, 0.5, -2]
            template.priors = ((6557.6, 0, 0,
                                6563, 0, 0,
                                6567, 0, 0,
                                6553, 0, 0,
                                6560, 0, 0,
                                6560, 0, 0,
                                6546, 0, 0,
                                6581, 0, 0,
                                6538, 0, 0,
                                6582, 0, 0,
                                6592, 0, 0,
                                6548, 0, 0,
                                6574, 0, 0,
                                6578, 0, 0,
                                6558, 0.2, -4),
                               (6570, 5, 60,
                                6571, 5, 60,
                                6575, 5, 60,
                                6558, 3, 20,
                                6570, 30, 10,
                                6564, 1, 30,
                                6548, 1.4, 1,
                                6583, 1.4, 2,
                                6543, 3, 1,
                                6590, 3, 1,
                                6603, 3, 1,
                                6552, 2, 3,
                                6578, 2, 1,
                                6582, 2, 1,
                                6560, 1, -1))

        elif version == 'Harvard_GMOS_P1' and solver == 'CF':
            template.window = (6520.0, 6607.0)
            template.n_components = 16
            template.component_number = [1, 2, 3, 4, 5, 6, 20, 21, 7, 8, 9, 11, 12, -1, -2, -3]
            template.guess = [6562, 3, 15,
                              6565, 3, 15,
                              6570, 3, 15,
                              6555, 1, 2,
                              6565, 15, 5,
                              6562, 0.5, 10,
                              6546.5, 1, 0.5,
                              6582, 1, 2,
                              6540, 2, 0.5,
                              6587, 2, 0.5,
                              6599, 2, 0.5,
                              6576, 1, 0.5,
                              6580, 1, 0.5,
                              6548, 1, -0.5,
                              6552, 1, -0.5,
                              6559, 0.5, -2]
            template.priors = ((6557.6, 0, 0,
                                6563, 0, 0,
                                6567, 0, 0,
                                6553, 0, 0,
                                6560, 0, 0,
                                6560, 0, 0,
                                6546, 0, 0,
                                6581, 0, 0,
                                6538, 0, 0,
                                6582, 0, 0,
                                6592, 0, 0,
                                6574, 0, 0,
                                6578, 0, 0,
                                6545, 0, -3,
                                6550, 0, -3,
                                6558, 0.2, -4),
                               (6570, 5, 60,
                                6571, 5, 60,
                                6575, 5, 60,
                                6558, 3, 20,
                                6570, 30, 10,
                                6564, 1, 30,
                                6548, 1.4, 1,
                                6583, 1.4, 2,
                                6543, 3, 1,
                                6590, 3, 1,
                                6603, 3, 1,
                                6578, 2, 1,
                                6582, 2, 1,
                                6554, 3, 0,
                                6556, 2, 0,
                                6560, 1, -1))

        elif version == 'Harvard_GMOS_P2' and solver == 'CF':
            template.window = (6520.0, 6607.0)
            template.n_components = 15
            template.component_number = [1, 2, 3, 4, 5, 6, 21, 7, 8, 9, 11, 12, -1, -2, -3]
            template.guess = [6562, 3, 15,
                              6565, 3, 15,
                              6570, 3, 15,
                              6555, 1, 2,
                              6565, 15, 5,
                              6562, 0.5, 10,
                              6582, 1, 2,
                              6540, 2, 0.5,
                              6587, 2, 0.5,
                              6599, 2, 0.5,
                              6576, 1, 0.5,
                              6580, 1, 0.5,
                              6548, 1, -0.5,
                              6552, 1, -0.5,
                              6559, 0.5, -2]
            template.priors = ((6557.6, 0, 0,
                                6563, 0, 0,
                                6567, 0, 0,
                                6553, 0, 0,
                                6560, 0, 0,
                                6560, 0, 0,
                                6581, 0, 0,
                                6538, 0, 0,
                                6582, 0, 0,
                                6592, 0, 0,
                                6574, 0, 0,
                                6578, 0, 0,
                                6545, 0, -3,
                                6550, 0, -3,
                                6558, 0.2, -4),
                               (6570, 5, 60,
                                6571, 5, 60,
                                6575, 5, 60,
                                6558, 3, 20,
                                6570, 30, 10,
                                6564, 1, 30,
                                6583, 1.4, 2,
                                6543, 3, 1,
                                6590, 3, 1,
                                6603, 3, 1,
                                6578, 2, 1,
                                6582, 2, 1,
                                6554, 3, 0,
                                6556, 2, 0,
                                6560, 1, -1))

        else:
            # Version and solver combination not defined.
            raise ValueError('Unknown version and solver.')

        return template


etacar_templates = EtaCarSpectralTemplateLibrary()
