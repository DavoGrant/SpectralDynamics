import sqlite3
import pandas as pd

from config import RESULTS_PATH, DB_PATH
from extract.core.spectra import SpectralHandler
from extract.core.lines import SpectralFits
from extract.observatories.gjw import GJWNetworkInstruments
from extract.features.novacar_library import novacar_templates
from extract.helpers import fetch_bad_jds


# Initialise dataset from disk.
handler = SpectralHandler()
handler.select_fits_dataset(
    dataset_path=RESULTS_PATH, datastream='GlobalJetWatch',
    data_release='BraveNewWorld', target='NovaCar2018', dimension='1D')

# Fitting config.
concat_table = 'Halpha_all_observatories'
jd_tuples = [{'label': 'GJW-CL', 'template': 'Simple',
              'solver': 'CF', 'jd_rule': ((2458203, 2458393), 'between'),
              'telescope': GJWNetworkInstruments.gjw_cl_spec_red},
             {'label': 'GJW-SA', 'template': 'Simple',
              'solver': 'CF', 'jd_rule': ((2458203, 2458393), 'between'),
              'telescope': GJWNetworkInstruments.gjw_sa_spec_red},
             {'label': 'GJW-WA', 'template': 'Simple',
              'solver': 'CF', 'jd_rule': ((2458203, 2458393), 'between'),
              'telescope': GJWNetworkInstruments.gjw_wa_spec_red},
             {'label': 'GJW-OZ', 'template': 'Simple',
              'solver': 'CF', 'jd_rule': ((2458203, 2458393), 'between'),
              'telescope': GJWNetworkInstruments.gjw_oz_spec_red}]

# Iterate different fitting regimes.
res_tables = []
for regime in jd_tuples:
    print('\nStarting new fitting regime={}\n'.format(regime['label']))

    # Select subset of the dataset.
    handler.select_fits_subset(
        pre_processing=None, binning=None, exposure=(None, 'max'),
        reduction='endeavour', observatory=regime['label'], jd=regime['jd_rule'])

    # Define fitting routines.
    fitter = SpectralFits(regime['telescope'])
    fitter.add(novacar_templates.h_i_6563(version=regime['template'],
                                          solver=regime['solver']))

    # Ready fitting routine and pre-processing options.
    fitter.compile(helio_correction=False, continuum_normalisation=True,
                   re_bin=None, refine_continuum=True,
                   bad_jds=fetch_bad_jds())

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
