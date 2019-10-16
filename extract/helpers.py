import sys
import sqlite3
import datetime
import numpy as np
import pandas as pd
from scipy import constants, signal, interpolate
from astropy.io import fits
from astropy.time import Time
from PyAstronomy import pyasl


def unpack_gjw_fits_file_name_to_series(file_path, spectra_object):
    """ Unpack GJW fits file name into pd series. """
    # Get file name:
    file_name = file_path.split('/')[-1]

    # Extract required info.
    datastream = spectra_object.datastream
    data_release = spectra_object.data_release
    target = file_name.split('-')[0]
    assert target == spectra_object.target
    dimension = spectra_object.dimension
    pre_processing = file_name.split('-')[2]
    binning = float(file_name.split('_')[2][3:])
    exposure = float(file_name.split('_')[1])
    reduction = file_name.split('_')[3].split('.')[0]
    observatory = file_name.split('_')[0][-6:]
    jd = float(file_name.split('-')[1])
    date = file_name.split('-')[3]

    # Make pandas series to return.
    extracted_series = pd.Series(
        [datastream, data_release, target, dimension,
         pre_processing, binning, exposure, reduction,
         observatory, jd, date, file_path])

    return extracted_series


def unpack_gmos_fits_file_to_series(file_path, spectra_object):
    """ Unpack GMOS fits file into pd series. """
    # Open fits file and select scidata from hdu list.
    hdul = fits.open(file_path)
    hdu_0 = hdul[0]
    hdu_1 = hdul[1]

    # Extract required info from header.
    datastream = spectra_object.datastream
    data_release = spectra_object.data_release
    target = hdu_1.header['OBJECT'].replace(' ', '').lower()
    assert target == spectra_object.target.lower()
    dimension = spectra_object.dimension
    pre_processing = 'Single'
    binning = hdu_1.header['CD1_1']
    exposure = hdu_0.header['EXPTIME']
    reduction = 'UMNArchive'
    observatory = hdu_0.header['OBSERVAT']
    # MJD to JD conversion.
    jd = hdu_0.header['TEXPSTRT'] + 2400000.5
    date = hdu_0.header['DATE-OBS']

    # Close fits file.
    hdul.close()

    # Make pandas series to return.
    extracted_series = pd.Series(
        [datastream, data_release, target, dimension,
         pre_processing, binning, exposure, reduction,
         observatory, jd, date, file_path])

    return extracted_series


def unpack_stis_fits_file_to_series(file_path, spectra_object):
    """ Unpack STIS fits file into pd series. """
    # Open fits file and select scidata from hdu list.
    hdul = fits.open(file_path)
    hdu_0 = hdul[0]
    hdu_1 = hdul[1]

    # Extract required info from header.
    datastream = spectra_object.datastream
    data_release = spectra_object.data_release
    target = hdu_0.header['TARGNAME'].replace('-', '')[0:6].lower()
    assert target == spectra_object.target.lower()
    dimension = spectra_object.dimension
    pre_processing = 'Single'
    binning = hdu_1.header['CD1_1']
    exposure = hdu_1.header['EXPTIME']
    reduction = 'UMNArchive'
    observatory = hdu_0.header['TELESCOP']
    # MJD to JD conversion.
    jd = hdu_1.header['EXPSTART'] + 2400000.5
    date = hdu_1.header['DATE-OBS']

    # Close fits file.
    hdul.close()

    # Make pandas series to return.
    extracted_series = pd.Series(
        [datastream, data_release, target, dimension,
         pre_processing, binning, exposure, reduction,
         observatory, jd, date, file_path])

    return extracted_series


def unpack_uves_fits_file_to_series(file_path, spectra_object):
    """ Unpack UVES fits file into pd series. """
    # Open fits file and select scidata from hdu list.
    hdul = fits.open(file_path)
    hdu_0 = hdul[0]
    hdu_1 = hdul[1]

    # Extract required info from header.
    datastream = spectra_object.datastream
    data_release = spectra_object.data_release
    target = hdu_0.header['TARGNAME'].replace('-', '')[0:6].lower()
    assert target == spectra_object.target.lower()
    dimension = spectra_object.dimension
    pre_processing = 'Single'
    binning = hdu_1.header['CDELT1']
    exposure = hdu_1.header['EXPTIME']
    reduction = 'UMNArchive'
    observatory = hdu_1.header['TELESCOP']
    # MJD to JD conversion.
    jd = hdu_1.header['MJD-OBS'] + 2400000.5
    date = hdu_1.header['DATE-OBS']

    # Close fits file.
    hdul.close()

    # Make pandas series to return.
    extracted_series = pd.Series(
        [datastream, data_release, target, dimension,
         pre_processing, binning, exposure, reduction,
         observatory, jd, date, file_path])

    return extracted_series


def construct_gaussian(x_vals, mu, sigma, a, base=0):
    """ Construct Gaussian. """
    y = base + (a * np.exp(-(x_vals - mu) ** 2 / (2. * sigma ** 2)))
    return y


def construct_weighted_normalised_gaussian(x_vals, mu, sigma, weight, base=0):
    """ Construct weighted normalised Gaussian. """
    y = base + weight * ((1 / (sigma * (2 * np.pi)**0.5))
                         * np.exp(-(x_vals - mu) ** 2 / (2. * sigma ** 2)))
    return y


def calculate_velocity(obs, lab):
    """ Calculate velocity between lab and observed wavelengths. """
    # Calculate velocity.
    velocity = constants.value('speed of light in vacuum') \
               * ((obs - lab) / lab)

    # Velocity in km per second.
    velocity /= 1000

    return velocity


def calculate_velocity_error(wavelength_error, lab):
    """ Calculate velocity error. """
    # Calculate velocity error.
    velocity_err = wavelength_error \
                   * constants.value('speed of light in vacuum') \
                   / lab / 1000

    return velocity_err


def calculate_wavelegnth(velocity, lab):
    """ Calculate observed wavelength for lab wavelength. """
    # Velocity in m/s per second.
    velocity *= 1000

    # Calculate wavelength.
    wavelength = lab * ((velocity /
                         constants.value('speed of light in vacuum')) + 1)

    return wavelength


def find_central_dispersion_axis_y_pixel(counts):
    """ Auto detect central dispersion axis position. """
    return np.argmax(counts)


def std_convoluted(image, N):
    """ Standard deviation kernel convolution. """
    im = np.array(image, dtype=float)
    im2 = im**2
    ones = np.ones(im.shape)

    kernel = np.ones((2*N+1, 2*N+1))
    s = signal.convolve2d(im, kernel, mode="same")
    s2 = signal.convolve2d(im2, kernel, mode="same")
    ns = signal.convolve2d(ones, kernel, mode="same")

    return np.sqrt((s2 - s**2 / ns) / ns)


def rebin_with_spectres(binning_wanted, wavelength, flux, flux_errors):
    """ Re-bin spectrum with third party func. """
    new_x_vals = np.arange(min(wavelength), max(wavelength),
                           float(binning_wanted))[2:-2]
    re_sampled, re_sampled_errs = spectres(
        wavelength.values, flux.values, new_x_vals, flux_errors.values)

    return new_x_vals, re_sampled, re_sampled_errs


def make_bins(wavelengths, make_rhs="False"):
    bin_widths = np.zeros(wavelengths.shape[0])

    # This option makes the final entry in the left hand sides array the right hand side of the final bin
    if make_rhs == "True":
        bin_lhs = np.zeros(wavelengths.shape[0] + 1)
        # The first lhs position is assumed to be as far from the first central wavelength as the rhs of the first bin.
        bin_lhs[0] = wavelengths[0] - (wavelengths[1] - wavelengths[0]) / 2
        bin_widths[-1] = (wavelengths[-1] - wavelengths[-2])
        bin_lhs[-1] = wavelengths[-1] + (wavelengths[-1] - wavelengths[-2]) / 2
        bin_lhs[1:-1] = (wavelengths[1:] + wavelengths[:-1]) / 2
        bin_widths[:-1] = bin_lhs[1:-1] - bin_lhs[:-2]

    # Otherwise just return the lhs positions of each bin
    else:
        bin_lhs = np.zeros(wavelengths.shape[0])
        bin_lhs[0] = wavelengths[0] - (wavelengths[1] - wavelengths[0]) / 2
        bin_widths[-1] = (wavelengths[-1] - wavelengths[-2])
        bin_lhs[1:] = (wavelengths[1:] + wavelengths[:-1]) / 2
        bin_widths[:-1] = bin_lhs[1:] - bin_lhs[:-1]

    return bin_lhs, bin_widths


def spectres(spec_wavs, spec_fluxes, resampling, spec_errs=None):
    # Generate arrays of left hand side positions and widths for the old and new bins
    filter_lhs, filter_widths = make_bins(resampling, make_rhs="True")
    spec_lhs, spec_widths = make_bins(spec_wavs)

    # Check that the range of wavelengths to be resampled onto falls within the initial sampling region
    if filter_lhs[0] < spec_lhs[0] or filter_lhs[-1] > spec_lhs[-1]:
        print("Spec_lhs, filter_lhs, filter_rhs, spec_rhs ", spec_lhs[0], filter_lhs[0], filter_lhs[-1], spec_lhs[-1])
        print(
            "spectres was passed a spectrum which did not cover the full wavelength range of the specified filter curve.")

    # Generate output arrays to be populated
    if spec_fluxes.ndim == 1:
        resampled = np.zeros((resampling.shape[0]))

    elif spec_fluxes.ndim == 2:
        resampled = np.zeros((len(resampling), spec_fluxes.shape[1]))

    if spec_errs is not None:
        if spec_errs.shape != spec_fluxes.shape:
            sys.exit("If specified, spec_errs must be the same shape as spec_fluxes.")
        else:
            resampled_errs = np.copy(resampled)

    start = 0
    stop = 0

    # Calculate the new spectral flux and uncertainty values, loop over the new bins
    for j in range(len(filter_lhs) - 1):

        # Find the first old bin which is partially covered by the new bin
        while spec_lhs[start + 1] <= filter_lhs[j]:
            start += 1

        # Find the last old bin which is partially covered by the new bin
        while spec_lhs[stop + 1] < filter_lhs[j + 1]:
            stop += 1

        if spec_fluxes.ndim == 1:

            # If the new bin falls entirely within one old bin the are the same the new flux and new error are the same as for that bin
            if stop == start:

                resampled[j] = spec_fluxes[start]
                if spec_errs is not None:
                    resampled_errs[j] = spec_errs[start]

            # Otherwise multiply the first and last old bin widths by P_ij, all the ones in between have P_ij = 1
            else:

                start_factor = (spec_lhs[start + 1] - filter_lhs[j]) / (spec_lhs[start + 1] - spec_lhs[start])
                end_factor = (filter_lhs[j + 1] - spec_lhs[stop]) / (spec_lhs[stop + 1] - spec_lhs[stop])

                spec_widths[start] *= start_factor
                spec_widths[stop] *= end_factor

                # Populate the resampled spectrum and uncertainty arrays
                resampled[j] = np.sum(spec_widths[start:stop + 1] * spec_fluxes[start:stop + 1]) / np.sum(
                    spec_widths[start:stop + 1])

                if spec_errs is not None:
                    resampled_errs[j] = np.sqrt(
                        np.sum((spec_widths[start:stop + 1] * spec_errs[start:stop + 1]) ** 2)) / np.sum(
                        spec_widths[start:stop + 1])

                # Put back the old bin widths to their initial values for later use
                spec_widths[start] /= start_factor
                spec_widths[stop] /= end_factor


        # The same as above, except operates on each row of the array, resampling all of the input models
        elif spec_fluxes.ndim == 2:

            if stop == start:

                resampled[j, :] = spec_fluxes[start, :]
                if spec_errs is not None:
                    resampled_errs[j, :] = spec_errs[start, :]

            else:

                start_factor = (spec_lhs[start + 1] - filter_lhs[j]) / (spec_lhs[start + 1] - spec_lhs[start])
                end_factor = (filter_lhs[j + 1] - spec_lhs[stop]) / (spec_lhs[stop + 1] - spec_lhs[stop])

                spec_widths[start] *= start_factor
                spec_widths[stop] *= end_factor

                resampled[j, :] = np.sum(
                    np.expand_dims(spec_widths[start:stop + 1], axis=1) * spec_fluxes[start:stop + 1, :],
                    axis=0) / np.sum(spec_widths[start:stop + 1])

                if spec_errs is not None:
                    resampled_errs[j, :] = np.sqrt(
                        np.sum((np.expand_dims(spec_widths[start:stop + 1], axis=1) * spec_errs[start:stop + 1]) ** 2,
                               axis=0)) / np.sum(spec_widths[start:stop + 1])

                spec_widths[start] /= start_factor
                spec_widths[stop] /= end_factor

    # If errors were supplied return the resampled spectrum and error arrays
    if spec_errs is not None:
        return resampled, resampled_errs

    # Otherwise just return the resampled spectrum array
    else:
        return resampled


def modal_bin_from_samples(_chain, param, n_bins):
    """ Get modal bin value form sample chain. """
    hist, bin_edges = np.histogram(
        _chain[param], bins=n_bins, range=(
            _chain[param].min(), _chain[param].max()))
    return (bin_edges[hist.argmax(axis=0)]
            + bin_edges[hist.argmax(axis=0)]) / 2


def fetch_bad_jds(db_path=None, fit=None, comp=None):
    """ Fetch bad spectra julian days. """
    if db_path is None or fit is None or comp is None:
        return []

    query = 'SELECT * FROM {} '.format('Bad_JDs')

    # Select from db.
    connection = sqlite3.connect(db_path)
    bad_jds_tbl = pd.read_sql_query(query, connection)
    connection.close()

    # Get fit and component wanted.
    bad_jds_tbl = bad_jds_tbl.loc[
        (bad_jds_tbl['Fit'] == fit)
        & ((bad_jds_tbl['Component'] == comp)
           | (bad_jds_tbl['Component'] == 'Any'))].copy()

    return bad_jds_tbl['JD']


def fetch_thresholds(db_path=None, fit=None, comp=None):
    """ Fetch thresholds. """
    if db_path is None or fit is None or comp is None:
        return []

    query = 'SELECT * FROM {} '.format('Thresholds')

    # Select from db.
    connection = sqlite3.connect(db_path)
    thresholds_tbl = pd.read_sql_query(query, connection)
    connection.close()

    # Get fit and component wanted.
    thresholds_tbl = thresholds_tbl.loc[
        (thresholds_tbl['Fit'] == fit) &
        (thresholds_tbl['Component'] == comp)].copy()

    return thresholds_tbl.values


def calculate_vacuum_wavelength_rydberg_formula(n_1, n_2, Z=1, N=0):
    """ Calculate vacuum wavelengths in angstroms. """
    assert n_1 < n_2

    # Adjusted Rydberg.
    R_m = constants.value('Rydberg constant') / \
          (1 + (constants.value('electron mass') /
                ((Z * constants.value('proton mass'))
                 + (N * constants.value('neutron mass')))))

    # Rydberg formula.
    lambda_vac = 1 / (R_m * (Z ** 2) * ((1 / (n_1 ** 2)) - (1 / (n_2 ** 2))))

    return lambda_vac * 1e10


def convert_vacuum_to_air_wavelengths(vac):
    """ Convert vacuum to air wavelength (Ciddor 1996). """
    air = pyasl.vactoair2(vac)
    return air


def multi_component_weighted_mean(df, components):
    """ Weighted mean/centroid of multiple components. """
    data = df.copy()
    res = []
    err = []
    for jd, day_fits in data.groupby(['JD']):

        day_comps = day_fits.loc[day_fits['ComponentNumber'].isin(components)]
        means = day_comps['Mean']
        mean_erros = day_comps['MeanError']
        weights = abs(day_comps['Amplitude'] * day_comps['StandardDeviation']
                      * np.sqrt(2 * np.pi))
        quadrature_error = (((mean_erros ** 2) * (weights ** 2)).sum() ** 0.5) \
                           / weights.sum()

        res.append(np.average(means, weights=weights))
        err.append(quadrature_error)

    return pd.Series(res), pd.Series(err)


def multi_component_area_total(df, components):
    """ Total area of multiple components. """
    data = df.copy()
    res = []
    for jd, day_fits in data.groupby(['JD']):

        day_comps = day_fits.loc[day_fits['ComponentNumber'].isin(components)]
        areas = day_comps['Amplitude'] * day_comps['StandardDeviation'] \
                * np.sqrt(2 * np.pi)

        res.append(areas.sum())

    return pd.Series(res)


def convert_jd_to_phase(jds, lit_p, lit_T0):
    """ Convert JD to phase. """
    return (jds - lit_T0) / lit_p


def construct_spectral_profile_grid(spec_fits, comp_subset, wavelengths,
                                    normalise=False):
    """ Construct 2D grid of time-series spectral profiles. """
    grid = []
    for jd, epoch_fits in spec_fits.groupby(['JD']):
        total_components = []
        fit_subset = epoch_fits.loc[
            epoch_fits['ComponentNumber'].isin(comp_subset)]
        if len(fit_subset) == 0:
            total_feature = np.zeros(len(wavelengths))
        else:
            for i, fit in fit_subset.iterrows():

                comp = construct_gaussian(wavelengths, fit['Mean'],
                                          fit['StandardDeviation'],
                                          fit['Amplitude'])
                total_components.append(comp.tolist())

            total_feature = np.asarray(total_components).sum(axis=0)

            if normalise:
                # total_feature /= total_feature.max()
                total_feature /= np.linalg.norm(total_feature, ord=1)

        grid.append(total_feature.tolist())

    return np.asarray(grid)
