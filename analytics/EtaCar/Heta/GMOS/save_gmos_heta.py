import os

from config import REPO_PATH, DB_PATH
from extract.core.lines import SpectralFits
from extract.helpers import convert_vacuum_to_air_wavelengths, \
    calculate_vacuum_wavelength_rydberg_formula


# Pickle combination alg fits.
SpectralFits.save(
    db_path=DB_PATH,
    table='Heta_RegimeConcat',
    comp_number=[1, 2, 3, 4],
    lab_wavelength=convert_vacuum_to_air_wavelengths(
        calculate_vacuum_wavelength_rydberg_formula(2, 9)),
    velocity_error='pooled_intra_night_variance',
    pickle_path=os.path.join(REPO_PATH, 'Pickles', 'EtaCar_v3',
                             'fits_gmos_heta_combi_[1,2,3,4].p'))

# Pickle benchmark bisector half max fits.
SpectralFits.save(
    db_path=DB_PATH,
    table='Heta_BHM',
    comp_number=None,
    lab_wavelength=convert_vacuum_to_air_wavelengths(
        calculate_vacuum_wavelength_rydberg_formula(2, 9)),
    velocity_error=10.0,
    pickle_path=os.path.join(REPO_PATH, 'Pickles', 'EtaCar_v3',
                             'fits_gmos_heta_bhm.p'))
