import os
import datetime
import sqlite3
import pickle
from multiprocessing import Pool
import emcee
import corner
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from astropy.io import fits
from barycorrpy import get_BC_vel
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import constants
from scipy.integrate import simps
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.optimize import curve_fit, differential_evolution
from scipy.signal import savgol_filter, medfilt
from scipy.stats import chisquare
from sklearn.mixture import BayesianGaussianMixture

from extract.helpers import construct_gaussian, \
    construct_weighted_normalised_gaussian, rebin_with_spectres, \
    calculate_velocity, multi_component_weighted_mean
from extract.core.stats import StatisticalFits


class SpectralFits(object):
    """ Spectral fits class.

    Attributes
    ----------

        active_spectrum : pandas series, active spectrum's meta data.

        active_data : pandas dataframe, spectral data for active
            spectrum.

        active_feature_data : pandas dataframe, spectral data for active
            spectrum within the window of the current template.

        spectral_fits : pandas dataframe, successful spectral fits and
            meta data for writing to database table.

    Methods
    -------

        add : Keras influenced add a fitting routine layer. Template
            object is appended to list of fits to be made per spectrum.

        compile : Keras influenced ready added routines and run options.

        fit : Run the added templates over every spectrum passed to the
            method.

        view : View fits saved in the database overlaid on original spectra.

        save : Save fits to pickle for use in the next stage of analysis.

    """

    def __init__(self, _instrument):
        # Setup options.
        self._helio_correction = False
        self._continuum_normalisation = False
        self._re_bin = None
        self._refine_continuum = False
        self._bad_jds = []
        self._diagnostics = False
        self._draw = False
        self._db_path = None

        # Runtime attributes.
        self.instrument = _instrument
        print(self)
        self._templates = []
        self.active_spectrum = pd.Series()
        self.active_data = pd.DataFrame()
        self.active_feature_data = pd.DataFrame()
        self._priors = None
        self._popt = None
        self._pcov = None
        self._chisq = None
        self._log_likelihood = None
        self._chisq_statistic = None
        self.guass_res = []
        self.guass_res_p_vals = []
        self.spectral_fits = pd.DataFrame()
        self.table_names = []
        self.rcs_edge_size = (15, 15)

        # Clock.
        self._start_time = datetime.datetime.utcnow()
        self.image_number = 0

    def __repr__(self):
        return 'Spectral fits object optimised for ' \
               'instrument={}.'.format(self.instrument.instrument)

    def add(self, template):
        """ Keras influenced add fitting routine layer. """
        self._templates.append(template)

    def compile(self, helio_correction=False, continuum_normalisation=False,
                re_bin=None, refine_continuum=False, bad_jds=None):
        """ Keras influenced ready added routines and run options. """
        self._helio_correction = helio_correction
        self._continuum_normalisation = continuum_normalisation
        self._re_bin = re_bin
        self._refine_continuum = refine_continuum
        if len(bad_jds) > 0:
            self._bad_jds = [round(float(j), 5) for j in bad_jds.tolist()]

        print('Readying the following templates:')
        for template in self._templates:
            print(template)

    def fit(self, dataset, diagnostics=False, draw=False, db_path=None):
        """ Run the added and compiled fitting routines. """
        self._diagnostics = diagnostics
        self._draw = draw
        self._db_path = db_path
        if not os.path.exists(self._db_path):
            raise FileExistsError('Database path not found.')

        print('Preparing to fit with options: \nHeliocentric correction={}, '
              'Continuum Normalisation={}, Continuum Refinement={}.'.format(
                self._helio_correction, self._continuum_normalisation,
                self._refine_continuum))

        for spec_idx, spectrum in dataset.iterrows():
            for temp_idx, template in enumerate(self._templates):

                try:
                    # Setup and pre-processing.
                    self._set_active_spectrum(spectrum)
                    if self._helio_correction:
                        self._apply_heliocentric_correction(template.sky_coord)
                    if self._continuum_normalisation:
                        self._apply_continuum_normalisation()

                    # Fitting
                    self._fit_spectral_feature(template)
                    print('Fitted spectrum {}/{} with template {}/'
                          '{}.'.format(spec_idx + 1, len(dataset),
                                       temp_idx + 1, len(self._templates)))

                except Exception as e:
                    # If failure, jump to next spectrum.
                    print('Failure: {}'.format(e))
                    continue

        end_time = datetime.datetime.utcnow()
        print('Fitting completed in {}.'.format(str(end_time - self._start_time)))

        # Write fits to db.
        self._saving_fits()

    def view(self, db_path, table, start_at=0, specific_jd=None,
             bad_jds=None, comp_number=0, save=False, fold=False,
             save_loc='Fit_comp{}_{}.jpg',
             custom_ax=None, custom_comp=[]):
        """ View fits made overlaid on spectra. """
        query = 'SELECT * FROM {} '.format(table)

        # Select from db.
        connection = sqlite3.connect(db_path)
        self.spectral_fits = pd.read_sql_query(query, connection)
        connection.close()

        # Float and sort by JD.
        self.spectral_fits['JD'] = self.spectral_fits['JD'].astype('float')
        self.spectral_fits = self.spectral_fits.sort_values(
            by=['JD'], ascending=True)

        # Apply clean jds.
        self.spectral_fits['RoundedJD'] = self.spectral_fits['JD'].apply(
            lambda x: round(x, 5))
        for bad_jd in bad_jds:
            self.spectral_fits.drop(
                self.spectral_fits[self.spectral_fits['RoundedJD'] == round(
                    float(bad_jd), 5)].index, inplace=True)

        for spec_path, spec_fits in self.spectral_fits.groupby(['FilePath']):

            # Start at.
            if spec_fits['JD'].iloc[0] < start_at:
                continue

            # Specific JD.
            if specific_jd:
                if not spec_fits['JD'].iloc[0] == specific_jd:
                    continue

            try:
                print('JD={}'.format(spec_fits['JD'].iloc[0]))
                # Setup and pre-processing.
                self._set_active_spectrum(spec_fits.iloc[0])
                if spec_fits['HelioCorrection'].iloc[0] == 1:
                    self._apply_heliocentric_correction()
                if spec_fits['ContinuumNormalisation'].iloc[0] == 1:
                    self._apply_continuum_normalisation()

                # Reset and extract window.
                self.active_feature_data = pd.DataFrame()
                self.active_feature_data = self.active_data.loc[
                    (self.active_data['Wavelength']
                     > spec_fits['WindowBlue'].iloc[0])
                    & (self.active_data['Wavelength']
                       < spec_fits['WindowRed'].iloc[0])].copy()
                self.active_feature_data.reset_index(inplace=True)

                # Re-binning.
                if spec_fits['ReBin'].iloc[0] is not None:
                    new_x_vals, re_sampled, re_sampled_errs = rebin_with_spectres(
                        spec_fits['ReBin'], self.active_feature_data['Wavelength'],
                        self.active_feature_data['Flux'],
                        self.active_feature_data['FluxError'])
                    self.active_feature_data = pd.DataFrame()
                    self.active_feature_data['Wavelength'] = new_x_vals
                    self.active_feature_data['Flux'] = re_sampled
                    self.active_feature_data['FluxError'] = re_sampled_errs

                # Continuum refinement.
                if spec_fits['RefineContinuum'].iloc[0] == 1:
                    self.active_feature_data['Flux'] = \
                        self._refine_continuum_subtraction()

                # Viewing
                self._chisq = spec_fits['ChiSquared'].iloc[0]
                self._popt = np.stack((
                    spec_fits['Mean'].values, spec_fits['StandardDeviation'].values,
                    spec_fits['Amplitude'].values), axis=-1).flatten()
                spec_fits.reset_index(drop=True, inplace=True)
                ix_custom_comp = spec_fits.index[spec_fits['ComponentNumber'].isin(
                    custom_comp)].tolist()
                self._interim_fit_plotting(
                    track_fit=comp_number, save_option=save, save_loc=save_loc,
                    custom_ax=custom_ax, custom_comp_ix=ix_custom_comp)

            except RuntimeError as e:
                # If failure, jump to next spectrum.
                print('Failure: {}'.format(e))

    @staticmethod
    def save(db_path, table, comp_number=None, lab_wavelength=None,
             velocity_error='Auto', pickle_path=None):
        """ Save fits for next step in analysis. """
        # SQL command.
        query = 'SELECT * FROM {} '.format(table)

        # Select from db.
        connection = sqlite3.connect(db_path)
        spectral_fits = pd.read_sql_query(query, connection)
        connection.close()

        if isinstance(comp_number, list):
            # Weighted mean.
            tracking_fits = spectral_fits.loc[
                spectral_fits['ComponentNumber'].isin(comp_number)].copy()
            w_mean, w_mean_error = multi_component_weighted_mean(
                tracking_fits, comp_number)

            # Reduce dataset to one row per spectrum.
            spectral_fits = spectral_fits.groupby(['JD']).first().reset_index()
            spectral_fits.drop(columns=[
                'Mean', 'MeanError',
                'StandardDeviation', 'StandardDeviationError',
                'Amplitude', 'AmplitudeError'], inplace=True)
            spectral_fits['Mean'] = w_mean
            spectral_fits['MeanError'] = w_mean_error
            spectral_fits['ComponentNumber'] = 0
            comp_number = 0

        elif comp_number is None:
            # No components i.e bhm method.
            spectral_fits['ComponentNumber'] = 0
            comp_number = 0

        spectral_fits = spectral_fits.loc[
            spectral_fits['ComponentNumber'] == comp_number].copy()

        # Float and sort by JD.
        spectral_fits['JD'] = spectral_fits['JD'].astype('float')
        spectral_fits = spectral_fits.sort_values(
            by=['JD'], ascending=True)

        # Pre-calculate velocity.
        spectral_fits['Velocity'] = calculate_velocity(
            spectral_fits['Mean'], lab_wavelength)

        # Pre-calculate velocity error.
        statistical_fits_object = StatisticalFits()
        spectral_fits['VelocityError'] = \
            statistical_fits_object.determine_radial_velocity_errors(
                spectral_fits, velocity_error, lab_wavelength)

        # Save.
        print('Pickling {} fits as velocities to {}'.format(
            table, pickle_path))
        pickle.dump(spectral_fits, open(pickle_path, 'wb'))

    def _set_active_spectrum(self, _spectrum_series):
        """ Set active spectrum. """
        self.active_spectrum = _spectrum_series

        # Check if in list of bad jds.
        if round(self.active_spectrum['JD'], 5) in self._bad_jds:
            raise ValueError('Bad JD detected, skipping to next spectrum.')

        # Read in .fits from disk to pandas.
        if self.active_spectrum['Datastream'] == 'GlobalJetWatch':

            dat = Table.read(self.active_spectrum['FilePath'], format='fits')
            self.active_data = dat.to_pandas()
            # self.active_data.drop(columns=['uncorrected_wavelength', 'x'], inplace=True)
            self.active_data.rename(index=str, columns={
                'wavelength': 'Wavelength',
                'intensity': 'Flux',
                'error': 'FluxError'}, inplace=True)

        elif self.active_spectrum['Datastream'] == 'GMOS':

            hdul = fits.open(self.active_spectrum['FilePath'])
            hdu_1 = hdul[1]  # Counts.

            statistical_fits_object = StatisticalFits()
            w, f, fe = statistical_fits_object.reduce_gmos_spectrum_and_errors(
                spectrum_2d=hdu_1, error_mode='poisson',
                dynamic_spatial_slit_width=False,
                arcseconds_wanted=3)
            hdul.close()

            self.active_data = pd.DataFrame()
            self.active_data['Wavelength'] = w
            self.active_data['Flux'] = f
            self.active_data['FluxError'] = fe

        elif self.active_spectrum['Datastream'] == 'STIS':

            hdul = fits.open(self.active_spectrum['FilePath'])
            hdu_1 = hdul[1]  # Flux.
            hdu_3 = hdul[3]  # Errors.

            statistical_fits_object = StatisticalFits()
            w, f, fe = statistical_fits_object.reduce_stis_spectrum_and_errors(
                spectrum_2d=hdu_1, errors_2d=hdu_3,
                dynamic_spatial_slit_width=False,
                arcseconds_wanted=3)
            hdul.close()

            self.active_data = pd.DataFrame()
            self.active_data['Wavelength'] = w
            self.active_data['Flux'] = f
            self.active_data['FluxError'] = fe

        elif self.active_spectrum['Datastream'] == 'UVES':

            hdul = fits.open(self.active_spectrum['FilePath'])
            hdu_1 = hdul[1]  # Counts.

            statistical_fits_object = StatisticalFits()
            w, f, fe = statistical_fits_object.reduce_uves_spectrum_and_errors(
                spectrum_2d=hdu_1, error_mode='poisson',
                dynamic_spatial_slit_width=False,
                arcseconds_wanted=3)
            hdul.close()

            self.active_data = pd.DataFrame()
            self.active_data['Wavelength'] = w
            self.active_data['Flux'] = f
            self.active_data['FluxError'] = fe

    def _apply_heliocentric_correction(self, _sky_coord):
        """ Apply heliocentric correction for the active spectrum. """
        # Archive original wavelength values.
        self.active_data.rename(index=str, columns={
            'Wavelength': 'WavelengthOriginal'}, inplace=True)

        # Get sky coordinates.
        c = SkyCoord(_sky_coord, unit=(u.hourangle, u.deg))

        # Calculate heliocentric velocity correction.
        # Up to 3 m/s precision.
        velocity = get_BC_vel(
            Time(float(self.active_spectrum['JD']), format='jd', scale='utc'),
            ra=c.ra.degree, dec=c.dec.degree,
            lat=self.instrument.latitude,
            longi=self.instrument.longitude,
            alt=self.instrument.elevation)

        # Correct wavelengths.
        # lambda_source = lambda_obs * sqrt((1-beta)/(1+beta))
        beta = -velocity[0][0] / constants.value('speed of light in vacuum')
        self.active_data['Wavelength'] = self.active_data['WavelengthOriginal'] \
                                         * np.sqrt((1. - beta) / (1. + beta))

    def _apply_continuum_normalisation(self):
        """ Apply continuum normalisation. """
        # Archive original wavelength values.
        self.active_data.rename(index=str, columns={
            'Flux': 'FluxOriginal',
            'FluxError': 'FluxErrorOriginal'}, inplace=True)

        # Median filter kernel size and savgol window length - assert odd.
        window_length = 2001
        medfilt_kernel_size = 301

        interpolated_regions = self.active_data.copy()
        for region_key, bounds in self.instrument.blocked_regions.items():

            try:
                # Find unblocked region.
                unblocked_region = interpolated_regions.loc[
                    (interpolated_regions['Wavelength'] < bounds[0])
                    | (interpolated_regions['Wavelength'] > bounds[1])]

                # Interpolate flux in blocked region.
                interpolation_function = interp1d(
                    unblocked_region['Wavelength'],
                    unblocked_region['FluxOriginal'],
                    kind="linear", fill_value='extrapolate')

                # Not so readable but way faster than df.apply.
                interpolated_regions.loc[
                    (interpolated_regions['Wavelength'] >= bounds[0])
                    & (interpolated_regions['Wavelength'] <= bounds[1]),
                    'FluxOriginal'] = interpolation_function(
                    interpolated_regions['Wavelength'].loc[
                    (interpolated_regions['Wavelength'] >= bounds[0])
                    & (interpolated_regions['Wavelength'] <= bounds[1])])

            except ValueError as e:
                # Blocked region out of spectral range, onto next region.
                continue

        # Apply median filter.
        # interpolated_regions['Medfilt_Flux'] = medfilt(
        #     interpolated_regions['FluxOriginal'],
        #     kernel_size=medfilt_kernel_size)

        # Apply Savitzky-Golay filter.
        continuum_fit = savgol_filter(
                interpolated_regions['FluxOriginal'],
                window_length, 3, mode="mirror")

        # Subtract continuum - N.B. only normalise OR subtract.
        # self.active_data['Flux'] = \
        #     self.active_data['FluxOriginal'] - continuum_fit

        # Normalise by continuum.
        self.active_data['Flux'] = \
            (self.active_data['FluxOriginal'] / continuum_fit) - 1
        self.active_data['FluxError'] = \
            self.active_data['FluxErrorOriginal'] / continuum_fit

        # Draw.
        if self._draw:
            fig = plt.figure(figsize=(10, 6))
            ax1 = fig.add_subplot(1, 1, 1)
            ax1.plot(self.active_data['Wavelength'],
                     self.active_data['FluxOriginal'],
                     label='Observed', linewidth=1)
            ax1.plot(self.active_data['Wavelength'],
                     self.active_data['FluxErrorOriginal'],
                     label='Error', linewidth=1)
            ax1.plot(self.active_data['Wavelength'],
                     interpolated_regions['FluxOriginal'],
                     label='Interpolated Blocked Regions',
                     linewidth=1)
            # ax1.plot(self.active_data['Wavelength'],
            #          interpolated_regions['Medfilt_Flux'],
            #          label='Interpolated and Median Filtered')
            ax1.plot(self.active_data['Wavelength'],
                     continuum_fit,
                     label='Savgol Filter',
                     linewidth=1)
            ax1.plot(self.active_data['Wavelength'],
                     self.active_data['Flux'],
                     label='Normalised Result',
                     linewidth=1)
            ax1.set_title('Continuum Normalisation: JD={}'.format(
                self.active_spectrum['JD']), fontsize=16)
            ax1.set_ylabel('Fluxes', fontsize=12)
            ax1.set_xlabel('Wavelength / $\AA$', fontsize=12)
            plt.legend()
            plt.tight_layout()
            plt.show()

    def _fit_spectral_feature(self, _template):
        """ Fit template to spectral feature. """
        # Reset and extract window.
        self.active_feature_data = pd.DataFrame()
        self.active_feature_data = self.active_data.loc[
            (self.active_data['Wavelength'] > _template.window[0])
            & (self.active_data['Wavelength'] < _template.window[1])].copy()
        self.active_feature_data.reset_index(inplace=True)

        # Sanity check selected data range.
        try:
            assert self.active_feature_data[
                       'Wavelength'].min() > _template.window[0]
            assert self.active_feature_data[
                       'Wavelength'].max() < _template.window[1]
        except AssertionError as err:
            raise IndexError('Wavelength range required for fitting not available.')

        # Check selected range includes spectral data.
        tolerance = 1
        if (self.active_feature_data['Wavelength'].min()
                - _template.window[0] > tolerance) \
                or (_template.window[1]
                    - self.active_feature_data['Wavelength'].max() > tolerance):
            print('Spectral data does not cover template region.')
            return

        # Re-binning.
        if self._re_bin is not None:
            new_x_vals_orig, re_sampled_orig, re_sampled_errs_orig = rebin_with_spectres(
                self._re_bin, self.active_feature_data['Wavelength'],
                self.active_feature_data['FluxOriginal'],
                self.active_feature_data['FluxErrorOriginal'])
            new_x_vals, re_sampled, re_sampled_errs = rebin_with_spectres(
                self._re_bin, self.active_feature_data['Wavelength'],
                self.active_feature_data['Flux'],
                self.active_feature_data['FluxError'])
            if self._draw:
                self._draw_re_binning(new_x_vals, re_sampled)
            self.active_feature_data = pd.DataFrame()
            self.active_feature_data['Wavelength'] = new_x_vals
            self.active_feature_data['FluxOriginal'] = re_sampled_orig
            self.active_feature_data['FluxErrorOriginal'] = re_sampled_errs_orig
            self.active_feature_data['Flux'] = re_sampled
            self.active_feature_data['FluxError'] = re_sampled_errs

            # Update binning meta data.
            self.active_spectrum['Binning'] = self._re_bin

        # Continuum refinement.
        if self._refine_continuum:
            self.active_feature_data['Flux'] = \
                self._refine_continuum_subtraction()

        if _template.solver == 'CF':
            # Scipy's curve fit algorithm.
            self._curve_fit(_template)
            self._interim_fit_storing(_template)
            if self._draw:
                self._interim_fit_plotting()

        elif _template.solver == 'DE':
            # Scipy's differential evolution algorithm.
            self._differential_evolution_fit(_template)
            self._interim_fit_storing(_template)
            if self._draw:
                self._interim_fit_plotting()

        elif _template.solver == 'MCMC':
            # Emcee's Markov Chain Monte Carlo.
            self._mcmc_fit(_template)
            self._interim_fit_storing(_template)
            if self._draw:
                self._interim_fit_plotting()

        elif _template.solver == 'BGM':
            # SKLearn's Bayesian Gaussian mixture algorithm.
            self._bayesian_gaussian_mixture_fit(_template)
            self._interim_fit_storing(_template)
            if self._draw:
                self._bgm_plotting()

        elif _template.solver[-4:] == 'Peak':
            # Peak value fit.
            self._peak_fit(_template.solver[:-4])
            self._interim_point_storing(_template)
            if self._draw:
                self._interim_point_plotting()

        elif _template.solver[-4:] == 'BHM':
            # Peak value fit.
            self._bisector_half_max_fit(_template)
            self._interim_point_storing(_template)
            if self._draw:
                self._interim_point_plotting()

        else:
            # Unknown fitting method specified.
            raise ValueError('Unknown fitting method specified.')

    def _curve_fit(self, __template):
        """ Fit spectral feature by curve fit. """
        # Fit.
        popt, pcov = curve_fit(
            self._multi_component_model,
            self.active_feature_data['Wavelength'],
            self.active_feature_data['Flux'],
            sigma=self.active_feature_data['FluxError'],
            p0=__template.guess, bounds=__template.priors)

        # Metrics.
        self._popt = popt
        self._pcov = pcov
        self._chisq, self._log_likelihood, self._chisq_statistic = \
            self._chi_squared_of_fit(solver='CF')

    @staticmethod
    def _multi_component_model(wavelengths, *args):
        """ Multi-component Gaussian model. """
        comp_params = np.split(np.array(args), len(args) / 3)
        total_components = []
        for model_params in comp_params:

            # Construct component.
            comp = construct_gaussian(
                wavelengths, model_params[0],
                model_params[1], model_params[2])

            # Add to total components list.
            total_components.append(comp.tolist())

        return np.asarray(total_components).sum(axis=0)

    def _differential_evolution_fit(self, __template):
        """ Fit spectral feature by differential evolution. """
        # Minimisation.
        res = differential_evolution(
            self._multi_component_chi_squared, __template.priors,
            args=(self.active_feature_data['Wavelength'],
                  self.active_feature_data['Flux']),
            strategy='best1bin', maxiter=10,
            popsize=10, tol=0.01, mutation=(0.5, 1), recombination=0.7,
            seed=5, callback=None, disp=False, polish=True,
            init='latinhypercube')

        # Metrics.
        self._popt = res.x
        self._pcov = None
        self._chisq = res.fun

    @staticmethod
    def _multi_component_chi_squared(coefs, *args):
        """ Chi-squared of multi-component Gaussian model. """
        wavelengths = args[0]
        observed = args[1]

        comp_params = np.split(np.array(coefs), len(coefs) / 3)
        total_components = []
        for model_params in comp_params:

            # Construct component.
            comp = construct_gaussian(
                wavelengths, model_params[0],
                model_params[1], model_params[2])

            # Add to total components list.
            total_components.append(comp.tolist())

        model = np.asarray(total_components).sum(axis=0)

        # Chi-squared.
        _chisq, p = chisquare(
            f_obs=observed + 100,
            f_exp=model + 100)

        return _chisq

    def _mcmc_fit(self, __template):
        """ Fit spectral feature by MCMC. """
        # Configure sampling.
        walkers = 256
        burn = 300
        run = 100

        # Set priors.
        self.priors = __template.priors

        # Set starting position as Gaussian ball around initial guess.
        initial_position = [__template.guess + 1e-3 * np.random.randn(len(
            __template.guess)) for i in range(walkers)]

        # Multi-processing on all cores.
        with Pool() as pool:

            # Affine invariant MCMC ensemble sampler.
            print('Burn in.')
            sampler = emcee.EnsembleSampler(
                walkers, len(__template.guess), self._mcmc_lnprob, pool=pool)
            pos, prob, state = sampler.run_mcmc(initial_position, burn)

            print('Walking.')
            sampler.reset()
            sampler.run_mcmc(pos, run)
            print('Finished MCMC.')

        # Get result dimensions.
        walkers, steps, pdims = sampler.chain.shape

        # Get results.
        mcmcm_samples = sampler.flatchain
        mean_acceptance_fraction = np.mean(sampler.acceptance_fraction)

        # Optimized parameter results taken to be 50th percentile.
        # Errors, 16th - 84th percentiles.
        optimized_parameters_median = np.percentile(mcmcm_samples, 50, axis=0)
        optimized_parameters_16pc = np.percentile(mcmcm_samples, 16, axis=0)
        optimized_parameters_84pc = np.percentile(mcmcm_samples, 84, axis=0)
        mcmc_results = pd.DataFrame(
            np.column_stack((
                optimized_parameters_median,
                optimized_parameters_16pc,
                optimized_parameters_84pc)),
            columns=['Median', '16th_percentile', '84th_percentile'])

        if self._draw:
            # Histograms and covariance.
            corner.corner(mcmcm_samples, bins=30,
                          smooth=0.8, quantiles=(0.16, 0.84))
            plt.show()

        # Metrics.
        self._popt = mcmc_results['Median'].values
        self._pcov = mcmc_results[['16th_percentile', '84th_percentile']].values
        self._pcov = None
        self._chisq, self._log_likelihood, self._chisq_statistic = \
            self._chi_squared_of_fit(solver='MCMC')

    def _mcmc_lnprob(self, theta):
        """ Logarithmic posterior or probability given the parameters. """
        # Prior.
        lp = self._mcmc_lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf

        # Likelihood.
        ll = self._mcmc_lnlike(theta)
        if not np.isfinite(ll):
            return -np.inf

        # Posterior.
        prob = lp + ll

        return prob

    def _mcmc_lnprior(self, theta):
        """ Logarithmic prior - up to a constant. """
        for idx, param in enumerate(theta):
            # Check prior.
            if not self.priors[0][idx] < param < self.priors[1][idx]:
                return -np.inf

        return 0.0

    def _mcmc_lnlike(self, theta):
        """ Logarithmic likelihood. """
        # Build model for current params.
        model_flux = self._multi_component_model(
            self.active_feature_data['Wavelength'], *theta)

        # Log likelihood: Chi-squared + ln(2 * pi * sigma^2).
        ln_like = -0.5 * ((((model_flux - self.active_feature_data['Flux']) ** 2)
                          / (self.active_feature_data['FluxError'] ** 2))
                          + np.log(2 * np.pi * (self.active_feature_data['FluxError'] ** 2))).sum()

        return ln_like

    def _bayesian_gaussian_mixture_fit(self, __template):
        """ Fit spectral feature by Bayesian Gaussian Mixture. """
        # Pre-process dataset ready for sci-kit learn's BGM.
        x = self._bgm_pre_processing()

        # Fit - N.B from the docs:
        # If 'warm_start' is True, the solution of the last fitting is used as
        # initialization for the next call of fit(). This can speed up
        # convergence when fit is called several time on similar problems.
        bgm = BayesianGaussianMixture(
            n_components=__template.n_components, tol=1e-3, reg_covar=1e-3,
            n_init=1, init_params='kmeans', max_iter=100, warm_start=False)
        bgm.fit(x)

        # Metrics.
        self._popt = np.column_stack(
            (bgm.means_.flatten(), bgm.covariances_.flatten(),
             bgm.weights_)).flatten()
        self._pcov = None
        self._chisq, self._log_likelihood, self._chisq_statistic = \
            self._chi_squared_of_fit(solver='BGM')

    def _bgm_pre_processing(self):
        """ Pre-processing of spectral data for BGM algorithm. """
        # Refactor data into frequency dataset.
        x = np.array([])
        for wavelength, flux in zip(
                self.active_feature_data['Wavelength'].values,
                self.active_feature_data['Flux'].values):
            if flux < 0:
                # Todo - if negative flux we have problems.
                x = np.append(x, np.zeros(0) + wavelength)
            else:
                x = np.append(x, np.zeros(int(flux)) + wavelength)
        return np.array([x]).reshape(len(x), 1)

    def _peak_fit(self, rule):
        """ Peak value fit. """
        # Max.
        max_index = self.active_feature_data['Flux'].idxmax()
        max_data_point = self.active_feature_data.iloc[max_index]
        if rule == 'Max':
            self._popt = max_data_point['Wavelength']
            return

        # Min.
        min_index = self.active_feature_data['Flux'].idxmin()
        min_data_point = self.active_feature_data.iloc[min_index]
        if rule == 'Min':
            self._popt = min_data_point['Wavelength']
            return

        # Upper peak zero derivative.
        if rule == 'MaxDyDx':
            self._popt = self._zero_derivative_point(max_data_point)
            return

        # Lower peak zero derivative.
        if rule == 'MinDyDx':
            self._popt = self._zero_derivative_point(min_data_point)
            return

    def _zero_derivative_point(self, helper_data_point):
        """ Lower peak zero derivative. """
        # Calc gradient.
        gradient = pd.Series(
            np.gradient(self.active_feature_data['Flux']))

        # Three points where gradient closest to zero.
        zero_grad_points = self.active_feature_data[
            'Wavelength'].iloc[gradient.abs().argsort()[:3]]

        # Select point closest to helper as most likely.
        zero_grad_wavelength = zero_grad_points.iloc[
            (zero_grad_points - helper_data_point[
                'Wavelength']).abs().argsort()]

        return zero_grad_wavelength.iloc[0]

    def _bisector_half_max_fit(self, __template):
        """ Bisector half max value fit. """
        spline = UnivariateSpline(self.active_feature_data['Wavelength'],
                                  self.active_feature_data['Flux']
                                  - (self.active_feature_data['Flux'].max()
                                     * __template.guess), s=0)

        self._popt = (spline.roots().max() + spline.roots().min()) / 2

    def _draw_re_binning(self, _new_x_vals, _re_sampled):
        """ Draw re-binnning. """
        fig = plt.figure(figsize=(10, 6))
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.plot(self.active_feature_data['Wavelength'],
                 self.active_feature_data['Flux'],
                 color='black', label='Original')
        ax1.plot(_new_x_vals, _re_sampled,
                 color='red', label='Bin {}'.format(self._re_bin))
        ax1.set_title('Re-binning', fontsize=16)
        ax1.set_ylabel('Relative Flux', fontsize=12)
        ax1.set_xlabel('Wavelength / $\AA$', fontsize=12)
        plt.legend()
        plt.show()

    def _refine_continuum_subtraction(self):
        """ Subtract any remaining continuum from poor savgol. """
        active_copy = self.active_feature_data.copy()

        # Extract featureless, only noise regions at window edges
        finessing_continuum = active_copy.drop(
            active_copy.index[self.rcs_edge_size[0]:-self.rcs_edge_size[1]])

        # Apply median filter.
        finessing_continuum['Flux'] = medfilt(
            finessing_continuum['Flux'], kernel_size=15)

        # Interpolate flux in gap regions.
        interpolation_function = interp1d(
            finessing_continuum['Wavelength'],
            finessing_continuum['Flux'],
            kind="linear")

        # Subtract remaining continuum.
        active_copy['Flux_Finesse'] = \
            active_copy['Flux'] - interpolation_function(
                active_copy['Wavelength'])

        if self._draw:
            fig = plt.figure(figsize=(10, 6))
            ax1 = fig.add_subplot(1, 1, 1)
            ax1.plot(active_copy['Wavelength'],
                     active_copy['Flux'],
                     color='black')
            ax1.plot(active_copy['Wavelength'],
                     interpolation_function(active_copy['Wavelength']),
                     color='orange')
            ax1.plot(active_copy['Wavelength'],
                     active_copy['Flux_Finesse'],
                     color='blue')
            ax1.axvspan(active_copy['Wavelength'].min(),
                        active_copy['Wavelength'].min()
                        + (self.rcs_edge_size[0] * self.active_spectrum['Binning']),
                        alpha=0.2, color='#000000')
            ax1.axvspan(active_copy['Wavelength'].max()
                        - (self.rcs_edge_size[1] * self.active_spectrum['Binning']),
                        active_copy['Wavelength'].max(),
                        alpha=0.2, color='#000000')
            ax1.set_title('Continuum Refinement', fontsize=16)
            ax1.set_ylabel('Normalised Flux', fontsize=12)
            ax1.set_xlabel('Wavelength / $\AA$', fontsize=12)
            plt.tight_layout()
            plt.show()

        return active_copy['Flux_Finesse']

    def _chi_squared_of_fit(self, solver='CF'):
        """ Chi-squared of fit. """
        # Chi-squared.
        comp_params = np.split(np.array(self._popt), len(self._popt) / 3)
        total_components = []

        if solver == 'BGM':
            for model_params in comp_params:
                comp = construct_weighted_normalised_gaussian(
                    self.active_feature_data['Wavelength'], model_params[0],
                    model_params[1], model_params[2])
                total_components.append(comp.tolist())

            area = simps(
                self.active_feature_data['Flux'],
                self.active_feature_data['Wavelength'], 1)

            chisq = None
            log_likelihood = None
            chisq_statistic, p = chisquare(
                f_obs=(self.active_feature_data['Flux'] / area) + 100,
                f_exp=np.asarray(total_components).sum(axis=0) + 100)

        else:
            for model_params in comp_params:
                comp = construct_gaussian(
                    self.active_feature_data['Wavelength'], model_params[0],
                    model_params[1], model_params[2])
                total_components.append(comp.tolist())

            # Use original bin counts for chi-squared testing.
            continuum = self.active_feature_data['FluxOriginal'] \
                        / (self.active_feature_data['Flux'] + 1)
            chisq = (((self.active_feature_data['FluxOriginal']
                     - ((np.asarray(total_components).sum(axis=0) + 1)
                        * continuum)) ** 2)
                     / (self.active_feature_data['FluxErrorOriginal'] ** 2)).sum()

            # Log likelihood.
            log_likelihood = -0.5 * ((((self.active_feature_data['FluxOriginal']
                                        - ((np.asarray(total_components).sum(axis=0) + 1)
                                           * continuum)) ** 2)
                                      / (self.active_feature_data['FluxErrorOriginal'] ** 2))
                                     + np.log(2 * np.pi * (
                                     self.active_feature_data['FluxErrorOriginal'] ** 2))).sum()

            # Chi-squared statistics.
            chisq_statistic, p = chisquare(
                f_obs=self.active_feature_data['FluxOriginal'] + 100,
                f_exp=((np.asarray(total_components).sum(
                    axis=0) + 1) * continuum) + 100)

        return chisq, log_likelihood, chisq_statistic

    def _interim_fit_plotting(self, track_fit=None, save_option=False,
                              save_loc=None, custom_ax=None,
                              custom_comp_ix=None):
        """ Plotting of fits. """
        # Plt component parameters.
        if track_fit is None:
            # During fitting only show current fit.
            fig = plt.figure(figsize=(10, 6))
            ax1 = fig.add_subplot(1, 1, 1)

        elif isinstance(track_fit, list):
            # Weighted combination of components.
            fig = plt.figure(figsize=(13, 12))
            gs = gridspec.GridSpec(7, 6, figure=fig)
            ax1 = plt.subplot(gs[0:3, :])
            ax2 = plt.subplot(gs[3:5, 0:3])
            ax3 = plt.subplot(gs[3:5, 3:6])
            ax4 = plt.subplot(gs[5:7, 0:2])
            ax5 = plt.subplot(gs[5:7, 2:4])
            ax6 = plt.subplot(gs[5:7, 4:6])
            plt.subplots_adjust(wspace=0.70)

            # Select components.
            tracking_fits = self.spectral_fits.loc[
                self.spectral_fits['ComponentNumber'].isin(track_fit)].copy()
            single_comp_df = tracking_fits.drop_duplicates(
                'JD', keep='first').copy().reset_index()
            single_comp_df['markersize'] = single_comp_df['FilePath'].apply(
                lambda f: '5' if f == self.active_spectrum['FilePath'] else '1')

            # Calc weighted mean by area.
            single_comp_df['w_mean'], single_comp_df['w_mean_error'] = \
                multi_component_weighted_mean(tracking_fits, track_fit)

            # Draw.
            for i, fit in single_comp_df.iterrows():
                markers, caps, bars = ax2.errorbar(
                    fit['JD'], fit['w_mean'],
                    yerr=None, color='#000000', fmt='o',
                    markersize=fit['markersize'],
                    elinewidth=1, capsize=2, capthick=1)
                [bar.set_alpha(0.2) for bar in bars]
                [cap.set_alpha(0.2) for cap in caps]

            ax2.set_xlabel('JD')
            ax2.set_ylabel('Weighted Mean')

            # Iterate components.
            for i, c in enumerate(track_fit):

                # Plot Gaussian params.
                g_param_list = ['Mean', 'StandardDeviation',
                                'Amplitude', 'EquivalentWidth']
                axes = [ax3, ax4, ax5, ax6]
                colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
                self._plot_gaussian_component_params(axes, g_param_list, c,
                                                     colour=colours[i])

        elif track_fit == 'PubPlot':
            ax1 = custom_ax
            comp_params = np.split(np.array(self._popt), len(self._popt) / 3)
            num_emission_components = sum(c[2] > 0 for c in comp_params)
            num_absorption_components = sum(c[2] < 0 for c in comp_params)
            ax1.text(0.03, 0.80, 'Emission x{}'.format(num_emission_components),
                     transform=ax1.transAxes, fontsize=12)
            ax1.text(0.03, 0.70, 'Absorption x{}'.format(num_absorption_components),
                     transform=ax1.transAxes, fontsize=12)

        else:
            # One component selected for plotting.
            fig = plt.figure(figsize=(10, 10))
            gs = gridspec.GridSpec(7, 2, figure=fig)
            ax1 = plt.subplot(gs[0:3, :])
            ax2 = plt.subplot(gs[3:5, 0])
            ax3 = plt.subplot(gs[3:5, 1])
            ax4 = plt.subplot(gs[5:7, 0])
            ax5 = plt.subplot(gs[5:7, 1])
            plt.subplots_adjust(wspace=0.22)

            # Plot Gaussian params.
            g_param_list = ['Mean', 'StandardDeviation',
                            'Amplitude', 'EquivalentWidth']
            axes = [ax2, ax3, ax4, ax5]
            self._plot_gaussian_component_params(axes, g_param_list, track_fit)

        # Plot observed spectrum.
        markers, caps, bars = ax1.errorbar(
            self.active_feature_data['Wavelength'],
            self.active_feature_data['Flux'],
            yerr=self.active_feature_data['FluxError'],
            color='#000000', fmt='o', markersize=1, elinewidth=1, capsize=1,
            capthick=1, alpha=1, label='Observations')
        [bar.set_alpha(0.4) for bar in bars]
        [cap.set_alpha(0.4) for cap in caps]

        # Plot model components.
        comp_params = np.split(np.array(self._popt), len(self._popt) / 3)
        total_components = []
        num = 0
        for model_params in comp_params:

            # Construct component and plot.
            comp = construct_gaussian(
                self.active_feature_data['Wavelength'], model_params[0],
                model_params[1], model_params[2])

            if track_fit == 'PubPlot':
                # Publication specs.
                if num in custom_comp_ix:
                    ax1.plot(self.active_feature_data['Wavelength'],
                             comp, color='#F28B00', linewidth=1, label='_nolegend_')
                else:
                    ax1.plot(self.active_feature_data['Wavelength'],
                             comp, color='#002663', linewidth=1, label='_nolegend_')
            else:
                # Generic specs.
                if num == 0:
                    ax1.plot(self.active_feature_data['Wavelength'],
                             comp, color='#FFD200', linewidth=1, label='Components')
                else:
                    ax1.plot(self.active_feature_data['Wavelength'],
                             comp, color='#FFD200', linewidth=1, label='_nolegend_')

            total_components.append(comp.tolist())
            num += 1

        # Plot summed components.
        total_feature = np.asarray(total_components).sum(axis=0)
        ax1.plot(self.active_feature_data['Wavelength'],
                 total_feature, color='#00B9E4', linewidth=1.5, label='Total Fit')

        # Plot residuals.
        divider = make_axes_locatable(ax1)
        ax1_res = divider.append_axes("bottom", size='30%', pad=0)
        ax1.get_shared_x_axes().join(ax1, ax1_res)
        ax1_res.scatter(self.active_feature_data['Wavelength'],
                        self.active_feature_data['Flux']
                        - total_feature, color='#000000', s=1, alpha=0.6)

        # Stats.
        sf = StatisticalFits()
        _residuals = self.active_feature_data['Flux'] - total_feature
        _p = sf.residuals_resemble_gauss_distribution(
            _residuals, test='KSTest', alpha=0.0027)

        # Draw.
        ax1.legend(fontsize=10)
        ax1.set_xticks([])
        ax1.set_title(
            '{} component fit: chi-squared (reduced) = {} ({}), '
            'Gauss residuals p = {}'.format(
                len(comp_params), round(self._chisq, 5),
                round(self._chisq /
                      (len(self.active_feature_data['Wavelength'])
                       - len(comp_params)), 5), _p), fontsize=12)
        ax1.set_ylabel('Normalised Flux', fontsize=10)
        ax1_res.set_xlabel('Wavelength / $\AA$', fontsize=10)
        ax1_res.set_ylabel('Residuals', fontsize=10)
        ax1_res.set_ylim(-total_feature.max() / 20, total_feature.max() / 20)
        plt.tight_layout()

        # Action.
        if track_fit == 'PubPlot':
            return
        if save_option == 'Image':
            plt.savefig(save_loc.format(track_fit, self.image_number))
            self.image_number += 1
        elif save_option == 'Stats':
            self.guass_res.append(_residuals.values)
            self.guass_res_p_vals.append(_p)
        else:
            plt.show()

        # Close figure.
        plt.clf()
        plt.close('all')

    def _plot_gaussian_component_params(self, ax_list, param_list, comp_s,
                                        colour='#000000'):
        """ Plot Gaussian component params. """
        # Select component.
        tracking_fits = self.spectral_fits.loc[
            self.spectral_fits['ComponentNumber'] == comp_s].copy()
        tracking_fits['markersize'] = tracking_fits['FilePath'].apply(
            lambda f: '5' if f == self.active_spectrum['FilePath'] else '1')

        # Calc EW.
        tracking_fits['EquivalentWidth'] = tracking_fits['Amplitude'] \
                                           * tracking_fits['StandardDeviation'] \
                                           * np.sqrt(2 * np.pi)

        # Draw.
        for axis, param in zip(ax_list, param_list):
            for i, fit in tracking_fits.iterrows():
                markers, caps, bars = axis.errorbar(
                    fit['JD'], fit[param],
                    yerr=None, color=colour, fmt='o',
                    markersize=fit['markersize'],
                    elinewidth=1, capsize=2, capthick=1)
                [bar.set_alpha(0.2) for bar in bars]
                [cap.set_alpha(0.2) for cap in caps]

            axis.set_xlabel('JD')
            axis.set_ylabel(param)

    def _bgm_plotting(self):
        """ Plotting of BGM model. """
        # Integrate area for normalisation.
        area = simps(
            self.active_feature_data['Flux'],
            self.active_feature_data['Wavelength'], 1)

        # Setup figure.
        fig = plt.figure(figsize=(10, 6))
        ax1 = fig.add_subplot(1, 1, 1)

        # Plot observed feature.
        ax1.plot(self.active_feature_data['Wavelength'],
                 self.active_feature_data['Flux'] / area,
                 color='blue', linewidth=0.8)

        # Plot model components.
        comp_params = np.split(np.array(self._popt), len(self._popt) / 3)
        total_components = []
        for model_params in comp_params:

            # Construct component and plot.
            comp = construct_weighted_normalised_gaussian(
                self.active_feature_data['Wavelength'], model_params[0],
                model_params[1], model_params[2])
            ax1.plot(self.active_feature_data['Wavelength'],
                     comp, color='orange', linewidth=0.8)
            total_components.append(comp.tolist())

        # Plot summed components.
        total_feature = np.asarray(total_components).sum(axis=0)
        ax1.plot(self.active_feature_data['Wavelength'],
                 total_feature, color='red', linewidth=0.5)

        # Show.
        plt.show()

        # Close figure.
        plt.clf()
        plt.close('all')

    def _interim_point_plotting(self):
        """ Plotting of point found. """
        # Setup figure.
        fig = plt.figure(figsize=(10, 6))
        ax1 = fig.add_subplot(1, 1, 1)

        # Plot observed feature.
        ax1.plot(self.active_feature_data['Wavelength'],
                 self.active_feature_data['Flux'],
                 color='blue', linewidth=0.8)

        # Plot point.
        ax1.axvline(x=self._popt, linewidth=0.5, color='#440154')

        # Draw and show.
        ax1.set_title('Line at {}'.format(round(self._popt, 3)), fontsize=16)
        ax1.set_ylabel('Normalised Flux', fontsize=12)
        ax1.set_xlabel('Wavelength / $\AA$', fontsize=12)
        plt.show()

        # Close figure.
        plt.clf()
        plt.close('all')

    def _interim_fit_storing(self, __template):
        """ Store model fits. """
        # Get current fits optimal params.
        optimal_params = np.split(np.array(self._popt), len(self._popt) / 3)
        optimal_params_df = pd.DataFrame(
            optimal_params, columns=['Mean', 'StandardDeviation', 'Amplitude'])

        # Get current fit param errors.
        sd_params_df = pd.DataFrame()
        if self._pcov is not None:
            standard_deviations = np.sqrt(np.diag(self._pcov))
            sd_params = np.split(np.array(standard_deviations),
                                 len(standard_deviations) / 3)
            sd_params_df = pd.DataFrame(
                sd_params, columns=['MeanError', 'StandardDeviationError',
                                    'AmplitudeError'])

        # Build df in db table format.
        comp_params_df = pd.concat([optimal_params_df, sd_params_df], axis=1)
        comp_params_df['ChiSquared'] = self._chisq
        comp_params_df['LogLikelihood'] = self._log_likelihood
        comp_params_df['ChiSquaredStatistic'] = self._chisq_statistic

        # Add component label/number to each row.
        comp_params_df.insert(loc=0, column='ComponentNumber',
                              value=pd.Series(__template.component_number))

        # Add meta data to each row.
        comp_params_df['SpectralFeature'] = __template.feature
        comp_params_df['TemplateVersion'] = __template.template_version
        comp_params_df['Solver'] = __template.solver
        comp_params_df['WindowBlue'] = __template.window[0]
        comp_params_df['WindowRed'] = __template.window[1]
        comp_params_df['HelioCorrection'] = self._helio_correction
        comp_params_df['ContinuumNormalisation'] = self._continuum_normalisation
        comp_params_df['ReBin'] = self._re_bin
        comp_params_df['RefineContinuum'] = self._refine_continuum
        for index, meta_data in self.active_spectrum.iteritems():
            comp_params_df[index] = meta_data

        # Append to spectral fits storage df.
        # Index specifies component number of fit.
        self.spectral_fits = pd.concat(
            [self.spectral_fits, comp_params_df], ignore_index=False)

        if self._diagnostics:
            print(self.active_spectrum['JD'])
            print(comp_params_df[['ComponentNumber', 'Mean',
                                  'StandardDeviation', 'Amplitude']])
            try:
                track_fit = [1, 2, 3, 4]
                tracking_fits = pd.DataFrame()
                for c in track_fit:
                    tf = self.spectral_fits.loc[
                        self.spectral_fits['ComponentNumber'] == c].copy()
                    tracking_fits = pd.concat(
                        [tracking_fits, tf], ignore_index=False)
                w_mean, w_mean_error = multi_component_weighted_mean(
                    tracking_fits, track_fit)
                print(w_mean)
            except IndexError as err:
                print('Warning: diagnostics attempting weighted '
                      'mean of unknown components.')

    def _interim_point_storing(self, __template):
        """ Store model points. """
        # Get current points optimal param.
        comp_params_df = pd.DataFrame()
        comp_params_df['Mean'] = np.array([self._popt])

        # Add meta data to each row.
        comp_params_df['SpectralFeature'] = __template.feature
        comp_params_df['TemplateVersion'] = __template.template_version
        comp_params_df['Solver'] = __template.solver
        for index, meta_data in self.active_spectrum.iteritems():
            comp_params_df[index] = meta_data

        # Append to spectral fits storage df.
        # Index specifies component number of fit.
        self.spectral_fits = pd.concat(
            [self.spectral_fits, comp_params_df], ignore_index=False)
        self.spectral_fits.index.name = 'ComponentNumber'

    def _saving_fits(self):
        """ Save fits to database. """
        for temp_idx, template in enumerate(self._templates):

            # Locate all fits for given template.
            template_fits = self.spectral_fits.loc[
                (self.spectral_fits[
                       'SpectralFeature'] == template.feature)
                & (self.spectral_fits[
                       'TemplateVersion'] == template.template_version)
                & (self.spectral_fits[
                       'Solver'] == template.solver)].copy()

            # Table name.
            start_time = str(self._start_time).replace('-', '_').replace(
                ':', '_').replace(' ', '_').split('.')[0]
            table_name = '{}_{}_fit_by_{}_at_{}'.format(
                template.feature, template.template_version,
                template.solver, start_time)
            self.table_names.append(table_name)

            # Write template's fits to table.
            print('Writing fits for template {} to database. '
                  'Table name={}'.format(temp_idx + 1, table_name))
            connection = sqlite3.connect(self._db_path)
            template_fits.to_sql(table_name, connection,
                                 if_exists='append', index=False)
            connection.close()
