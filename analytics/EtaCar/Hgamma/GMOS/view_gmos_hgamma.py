from config import DB_PATH
from extract.core.lines import SpectralFits
from extract.observatories.gemini import GeminiTelescopeInstruments
from extract.helpers import fetch_bad_jds


# Define viewing routines.
viewer = SpectralFits(GeminiTelescopeInstruments.gmos)
viewer.view(db_path=DB_PATH,
            table='Hgamma_RegimeConcat',
            start_at=0, comp_number=[1, 2, 3, 4],
            bad_jds=fetch_bad_jds(
                db_path=DB_PATH, fit='h_gamma_master', comp='Any'),
            fold=False, save=False,
            save_loc='/Users/grantd/ScreenGrabs/SpectralDynamics/'
                     'Videos/Hgamma/Hgamma_GMOS_Fit_comp{}_{}.jpg')
