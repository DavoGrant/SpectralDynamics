import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import gaussian_filter1d

from extract.helpers import find_central_dispersion_axis_y_pixel, \
    std_convoluted, calculate_velocity


class StatisticalFits(object):
    """ Statistical fitting class.

    Methods
    -------

        reduce_gmos_spectrum_and_errors : calculate wavelength, flux
            and errors from 2D to 1D GMOS spectrum.

        reduce_stis_spectrum_and_errors : calculate wavelength, flux
            and errors from 2D to 1D STIS spectrum.

        reduce_uves_spectrum_and_errors : calculate wavelength, flux
            and errors from 2D to 1D UVES spectrum.

        determine_radial_velocity_errors : post fitting calculation
            of errors.

        information_criterion : evaluate model likelihood vs model
            complexity - AIC, AICc, BIC.

        residuals_resemble_gauss_distribution : test if fit residuals
            resemble a Gaussian distribution.


    """

    def __init__(self):
        self.verbose = False

    def __repr__(self):
        return 'Statistical fitting object {}.'

    def reduce_gmos_spectrum_and_errors(self, spectrum_2d, error_mode='poisson',
                                        dynamic_spatial_slit_width=False,
                                        arcseconds_wanted=None):
        """ Calculate wavelength, flux and errors from 2D to 1D GMOS spectrum. """
        # Get wavelength map details from header.
        ref_pixel = spectrum_2d.header['CRPIX1']
        lambda_ref_pixel = spectrum_2d.header['CRVAL1']
        lambda_per_pixel = spectrum_2d.header['CD1_1']
        arcsecs_per_pixel = spectrum_2d.header['CD2_2'] * 3600

        if dynamic_spatial_slit_width:
            # Find central dispersion axis pixel.
            central_dispersion_axis_pixel = find_central_dispersion_axis_y_pixel(
                spectrum_2d.data[:, int(spectrum_2d.data.shape[1] / 2)])

            # Bin flux vertically on ccd within specified width of dispersion axis.
            vertical_pixel_amplitude = int(arcseconds_wanted / arcsecs_per_pixel / 2)
            counts_2d = spectrum_2d.data
            counts_1d = counts_2d[central_dispersion_axis_pixel
                                  - vertical_pixel_amplitude:
                                  central_dispersion_axis_pixel
                                  + vertical_pixel_amplitude, :].sum(axis=0)
            err_2d_kernel = std_convoluted(counts_2d, N=1)
            err_1d_kernel = np.sqrt(np.power(
                err_2d_kernel[central_dispersion_axis_pixel - vertical_pixel_amplitude:
                              central_dispersion_axis_pixel + vertical_pixel_amplitude, :],
                2).sum(axis=0))

        else:
            # Bin flux vertically on total ccd.
            counts_2d = spectrum_2d.data
            counts_1d = counts_2d[0:-1, :].sum(axis=0)
            err_2d_kernel = std_convoluted(counts_2d, N=1)
            err_1d_kernel = np.sqrt(np.power(err_2d_kernel[0:-1, :], 2).sum(axis=0))

        # Calibrate wavelengths.
        pixels = np.arange(0, len(counts_1d))
        wavelengths = lambda_ref_pixel + (lambda_per_pixel * (pixels - ref_pixel))

        # Bin flux.
        flux = pd.Series(counts_1d)

        # Determine flux error.
        if error_mode == 'poisson':
            flux_error = np.sqrt(pd.Series(counts_1d))
        elif error_mode == 'sd_kernel':
            flux_error = pd.Series(err_1d_kernel)
        elif error_mode == 'gauss_smooth':
            flux_error = abs(flux - gaussian_filter1d(flux, 2))
        else:
            raise ValueError('Error mode in spectrum reduction not recognised.')

        return pd.Series(wavelengths), flux, flux_error

    def reduce_stis_spectrum_and_errors(self, spectrum_2d, errors_2d,
                                        dynamic_spatial_slit_width=False,
                                        arcseconds_wanted=None):
        """ Calculate wavelength, flux and errors from 2D to 1D STIS spectrum. """
        # Get wavelength map details from header.
        ref_pixel = spectrum_2d.header['CRPIX1']
        lambda_ref_pixel = spectrum_2d.header['CRVAL1']
        lambda_per_pixel = spectrum_2d.header['CD1_1']
        arcsecs_per_pixel = spectrum_2d.header['CD2_2'] * 3600
        central_dispersion_axis_pixel = int(spectrum_2d.header['CRPIX2'])

        if dynamic_spatial_slit_width:
            # Bin flux vertically on ccd within specified width of dispersion axis.
            vertical_pixel_amplitude = int(arcseconds_wanted / arcsecs_per_pixel / 2)
            flux_2d = spectrum_2d.data
            flux_1d = flux_2d[central_dispersion_axis_pixel - vertical_pixel_amplitude:
                              central_dispersion_axis_pixel + vertical_pixel_amplitude, :].sum(axis=0)
            err_2d = errors_2d.data
            err_1d = np.sqrt(np.power(
                err_2d[central_dispersion_axis_pixel - vertical_pixel_amplitude:
                       central_dispersion_axis_pixel + vertical_pixel_amplitude, :], 2).sum(axis=0))

        else:
            # Bin flux vertically on total ccd.
            flux_2d = spectrum_2d.data
            flux_1d = flux_2d[0:-1, :].sum(axis=0)
            err_2d = errors_2d.data
            err_1d = np.sqrt(np.power(err_2d[0:-1, :], 2).sum(axis=0))

        # Calibrate wavelengths.
        pixels = np.arange(0, len(flux_1d))
        wavelengths = lambda_ref_pixel + (lambda_per_pixel * (pixels - ref_pixel))

        # Bin flux.
        flux = pd.Series(flux_1d)

        # Flux error.
        flux_error = pd.Series(err_1d)

        return pd.Series(wavelengths), flux, flux_error

    def reduce_uves_spectrum_and_errors(self, spectrum_2d, error_mode='poisson',
                                        dynamic_spatial_slit_width=False,
                                        arcseconds_wanted=None):
        """ Calculate wavelength, flux and errors from 2D to 1D UVES spectrum. """
        # Get wavelength map details from header.
        ref_pixel = spectrum_2d.header['CRPIX1']
        lambda_ref_pixel = spectrum_2d.header['CRVAL1']
        lambda_per_pixel = spectrum_2d.header['CDELT1']
        arcsecs_per_pixel = spectrum_2d.header['CDELT2']

        if dynamic_spatial_slit_width:
            # Find central dispersion axis pixel.
            central_dispersion_axis_pixel = find_central_dispersion_axis_y_pixel(
                spectrum_2d.data[:, int(spectrum_2d.data.shape[1] / 2)])

            # Bin flux vertically on ccd within specified width of dispersion axis.
            vertical_pixel_amplitude = int(arcseconds_wanted / arcsecs_per_pixel / 2)
            counts_2d = spectrum_2d.data
            counts_1d = counts_2d[central_dispersion_axis_pixel
                                  - vertical_pixel_amplitude:
                                  central_dispersion_axis_pixel
                                  + vertical_pixel_amplitude, :].sum(axis=0)
            err_2d_kernel = std_convoluted(counts_2d, N=1)
            err_1d_kernel = np.sqrt(np.power(
                err_2d_kernel[central_dispersion_axis_pixel - vertical_pixel_amplitude:
                              central_dispersion_axis_pixel + vertical_pixel_amplitude, :],
                2).sum(axis=0))

        else:
            # Bin flux vertically on total ccd.
            counts_2d = spectrum_2d.data
            counts_1d = counts_2d[0:-1, :].sum(axis=0)
            err_2d_kernel = std_convoluted(counts_2d, N=1)
            err_1d_kernel = np.sqrt(np.power(err_2d_kernel[0:-1, :], 2).sum(axis=0))

        # Calibrate wavelengths.
        pixels = np.arange(0, len(counts_1d))
        wavelengths = lambda_ref_pixel + (lambda_per_pixel * (pixels - ref_pixel))

        # Bin flux.
        flux = pd.Series(counts_1d)

        # Determine flux error.
        if error_mode == 'poisson':
            flux_error = np.sqrt(pd.Series(counts_1d))
        elif error_mode == 'sd_kernel':
            flux_error = pd.Series(err_1d_kernel)
        elif error_mode == 'gauss_smooth':
            flux_error = abs(flux - gaussian_filter1d(flux, 2))
        else:
            raise ValueError('Error mode in spectrum reduction not recognised.')

        return pd.Series(wavelengths), flux, flux_error

    def determine_radial_velocity_errors(self, _spectral_fits, _velocity_error,
                                         _lab_wavelength=None):
        """ Determine radial velocity errors. """
        if _velocity_error == 'fitting_covariance':
            _spectral_fits['VelocityError'] = \
                calculate_velocity(
                    _spectral_fits['Mean'] + _spectral_fits['MeanError'],
                    _lab_wavelength) - _spectral_fits['Velocity']

        elif _velocity_error == 'intra_night_spread':
            # Bin over 1 day intervals.
            interval = 1
            bins = np.arange(int(_spectral_fits['JD'].min()),
                             int(_spectral_fits['JD'].max()) + 2 * interval,
                             interval)
            _spectral_fits['JDBinned'] = pd.cut(_spectral_fits['JD'], bins)

            # Calc sd within each day bin.
            daily_std_velocity = _spectral_fits.groupby(
                ['JDBinned'])[['Velocity']].std()
            daily_std_velocity = daily_std_velocity['Velocity'].dropna(axis=0)

            # Set daily sd back on all observations within that day.
            _spectral_fits = pd.merge(_spectral_fits, daily_std_velocity,
                                      on=['JDBinned'], how='outer',
                                      suffixes=('', 'Error'))

            # For days with only one observation, use median sd.
            _spectral_fits['VelocityError'].fillna(
                _spectral_fits['VelocityError'].median(), inplace=True)

        elif _velocity_error == 'pooled_intra_night_variance':
            # Bin over 1 day intervals.
            interval = 1
            bins = np.arange(int(_spectral_fits['JD'].min()),
                             int(_spectral_fits['JD'].max()) + 2 * interval,
                             interval)
            _spectral_fits['JDBinned'] = pd.cut(_spectral_fits['JD'], bins)

            # Calc n_observations and sd within each day bin.
            daily_std_velocity = _spectral_fits.groupby(
                ['JDBinned'])[['Velocity']].std()['Velocity']
            daily_count_velocity = _spectral_fits.groupby(
                ['JDBinned'])[['Velocity']].count()['Velocity']

            # Keep only nights with > 1 observations.
            valid_std_nights = daily_std_velocity.ix[
                daily_count_velocity.index[daily_count_velocity > 1].to_list()].copy()
            valid_count_nights = daily_count_velocity.ix[
                daily_count_velocity.index[daily_count_velocity > 1].to_list()].copy()

            # Calc pooled variance.
            pooled_var = ((valid_count_nights - 1) * (valid_std_nights ** 2)).sum() \
                         / (valid_count_nights - 1).sum()

            # Calc pooled standard deviation.
            pooled_std = pooled_var ** 0.5

            _spectral_fits['VelocityError'] = pooled_std

        else:
            _spectral_fits['VelocityError'] = _velocity_error

        return _spectral_fits['VelocityError']

    def information_criterion(self, data, n_params, n_dp, test='AICc'):
        """ Information criterion. """
        # Akaike information criterion (AIC).
        if test == 'AIC':
            data[test] = (2 * n_params) - (2 * (data['LogLikelihood']))

        # Corrected Akaike information criterion (AIC).
        elif test == 'AICc':
            data[test] = (2 * n_params) - (2 * (data['LogLikelihood'])) \
                         + ((2 * (n_params ** 2)) + (2 * n_params)
                            / (n_dp - n_params - 1))

        # Bayesian information criterion (BIC).
        elif test == 'BIC':
            data[test] = (np.log(n_dp) * n_params) \
                         - (2 * (-0.5 * data['ChiSquared']))

        return data

    def residuals_resemble_gauss_distribution(self, residuals, test='NormalTest',
                                              alpha=1e-3):
        """ Check normalised residuals show a Gaussian distribution. """
        # Normalise.
        norm_residuals = residuals / residuals.sum()

        # Stat test.
        # Null hypothesis: residuals come from a normal distribution.
        if test == 'NormalTest':
            k2, p = stats.normaltest(norm_residuals)
        elif test == 'KSTest':
            loc, scale = stats.norm.fit(norm_residuals)
            D, p = stats.kstest(norm_residuals, 'norm', args=(loc, scale))
        else:
            raise ValueError('Stat test not recognised. Must be '
                             'one of NormalTest, KSTest.')

        if self.verbose:
            print('p = {:g}'.format(p))
            if p < alpha:
                print('The null hypothesis can be rejected')
                print('Residuals do not come from a normal distribution.')
            else:
                print('The null hypothesis cannot be rejected')
                print('Residuals resemble a normal distribution.')

        return p
