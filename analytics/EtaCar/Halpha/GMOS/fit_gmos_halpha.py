import sqlite3
import pandas as pd

from config import RESULTS_PATH, DB_PATH
from extract.core.spectra import SpectralHandler
from extract.core.lines import SpectralFits
from extract.observatories.gemini import GeminiTelescopeInstruments
from extract.features.etacar_library import etacar_templates
from extract.helpers import fetch_bad_jds


# Initialise dataset from disk.
handler = SpectralHandler()
handler.select_fits_dataset(
    dataset_path=RESULTS_PATH, datastream='GMOS',
    data_release='Archival', target='EtaCar', dimension='2D')

# Fitting config.
concat_table = 'Halpha_RegimeConcat1'
jd_tuples = [{'label': 'in-apastron', 'template': 'Harvard_GMOS_A1',
              'solver': 'CF', 'jd_rule': (2454826.790, 'before')},

             {'label': 'in-transition', 'template': 'Harvard_GMOS_A1',
              'solver': 'CF', 'jd_rule': ((2454826.790, 2454840.667), 'between')},

             {'label': 'early-periastron', 'template': 'Harvard_GMOS_P1',
              'solver': 'CF', 'jd_rule': ((2454840.667, 2454846.678), 'between')},

             {'label': 'late-periastron', 'template': 'Harvard_GMOS_P2',
              'solver': 'CF', 'jd_rule': ((2454846.678, 2454949.535), 'between')},

             {'label': 'out-transition', 'template': 'Harvard_GMOS_T1',
              'solver': 'CF', 'jd_rule': ((2454949.535, 2455204.803), 'between')},

             {'label': 'out-apastron', 'template': 'Harvard_GMOS_A1',
              'solver': 'CF', 'jd_rule': (2455204.803, 'after')},

             {'label': 'benchmark', 'template': 'BHM',
              'solver': 'BHM', 'jd_rule': (None, 'exact')}]

# Iterate different fitting regimes.
res_tables = []
for regime in jd_tuples:
    print('\nStarting new fitting regime={}\n'.format(regime['label']))

    # Select subset of the dataset.
    handler.select_fits_subset(
        pre_processing=None, binning=None, exposure=(None, 'max'),
        reduction=None, observatory=None, jd=regime['jd_rule'])

    # Define fitting routines.
    fitter = SpectralFits(GeminiTelescopeInstruments.gmos)
    fitter.add(etacar_templates.h_i_6563(version=regime['template'],
                                         solver=regime['solver']))

    # Ready fitting routine and pre-processing options.
    fitter.compile(helio_correction=False, continuum_normalisation=True,
                   re_bin=None, refine_continuum=True,
                   bad_jds=fetch_bad_jds(db_path=DB_PATH,
                                         fit='h_alpha_master',
                                         comp='Any'))

    # Execute fitting routines.
    fitter.fit(handler.data_subset, diagnostics=False, draw=False, db_path=DB_PATH)

    # Store consecutive regimes.
    if not regime['label'] == 'benchmark':
        res_tables.extend(fitter.table_names)

# Join regimes.
print('Collating fits into table={}'.format(concat_table))
connection = sqlite3.connect(DB_PATH)
for t in res_tables:
    query = 'SELECT * FROM {} '.format(t)
    data_table = pd.read_sql_query(query, connection)
    data_table.to_sql(concat_table, connection, if_exists='append', index=False)
connection.close()
