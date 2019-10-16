import os
import pandas as pd

from extract.helpers import unpack_gjw_fits_file_name_to_series, \
    unpack_gmos_fits_file_to_series, unpack_stis_fits_file_to_series, \
    unpack_uves_fits_file_to_series


class SpectralHandler(object):
    """ Spectral data handler class.

    Attributes
    ----------

        all_data : pandas dataframe, all fits files within the
            dataset_path location and their respective meta data.

        data_subset : pandas dataframe, a subset of rows of the all_data
            attribute filtered on meta data such as exposure time, etc.

    Methods
    -------

        select_fits_dataset : initialise a dataset by a root directory
            location containing fits files.

        select_fits_subset : select a subset of fits data by binning,
            exposure, reduction, observatory, jd.

    """

    def __init__(self):
        self._dataset_path = None
        self.datastream = None
        self._data_path = None
        self.data_release = None
        self.target = None
        self.dimension = None
        self.pre_processing = None
        self._columns = [
            'Datastream', 'DataRelease', 'Target', 'Dimension',
            'PreProcessing', 'Binning', 'Exposure', 'Reduction',
            'Observatory', 'JD', 'Date', 'FilePath']
        self.all_data = pd.DataFrame()
        self.data_subset = pd.DataFrame()

    def __repr__(self):
        return 'Spectral data handler for {}.'.format(self.target)

    def select_fits_dataset(self, dataset_path, datastream,
                            data_release, target, dimension):
        """ Initialise fits files dataset. """
        self.data_release = data_release
        self.target = target
        self.dimension = dimension

        print('Looking through fits files on disk.')

        # Location of dataset.
        self._dataset_path = dataset_path
        self.datastream = datastream
        if datastream == 'GlobalJetWatch':
            self._data_path = os.path.join(
                self._dataset_path, self.datastream, 'Results',
                self.data_release, self.target, self.dimension)
        else:
            self._data_path = os.path.join(
                self._dataset_path, self.datastream, 'Results',
                self.target)

        # Load in data files.
        found_data = pd.Series()
        for path, directories, files in os.walk(self._data_path):
            for name in files:
                if name.endswith('.fits'):
                    found_data = found_data.append(
                        pd.Series(os.path.join(path, name)),
                        ignore_index=True)

        # Reformat into df.
        if datastream == 'GlobalJetWatch':
            self.all_data[self._columns] = found_data.apply(
                lambda f_path: unpack_gjw_fits_file_name_to_series(f_path, self))
        elif datastream == 'GMOS':
            self.all_data[self._columns] = found_data.apply(
                lambda f_path: unpack_gmos_fits_file_to_series(f_path, self))
        elif datastream == 'STIS':
            self.all_data[self._columns] = found_data.apply(
                lambda f_path: unpack_stis_fits_file_to_series(f_path, self))
        elif datastream == 'UVES':
            self.all_data[self._columns] = found_data.apply(
                lambda f_path: unpack_uves_fits_file_to_series(f_path, self))

        print('Found {} individual fits files.'.format(
            len(self.all_data)))

    def select_fits_subset(self, pre_processing=None, binning=None,
                           exposure=(None, 'max'), reduction=None,
                           observatory=None, jd=(None, 'exact')):
        """ Select a subset of the loaded fits dataset. """
        self.data_subset = self.all_data.copy()

        # Select on pre-processing if provided.
        if pre_processing:
            self.data_subset = self.data_subset.loc[
                self.data_subset['PreProcessing'] == pre_processing]

        # Select on binning if provided.
        if binning:
            self.data_subset = self.data_subset.loc[
                self.data_subset['Binning'] == binning]

        # Select on exposure if provided.
        if exposure[0]:
            if exposure[1] == 'exact':
                self.data_subset = self.data_subset.loc[
                    self.data_subset['Exposure'] == exposure[0]]
            elif exposure[1] == 'min':
                self.data_subset = self.data_subset.loc[
                    self.data_subset['Exposure'] > exposure[0]]
            elif exposure[1] == 'max':
                self.data_subset = self.data_subset.loc[
                    self.data_subset['Exposure'] < exposure[0]]
            elif exposure[1] == 'between':
                self.data_subset = self.data_subset.loc[
                    (self.data_subset['Exposure'] > exposure[0][0]) &
                    (self.data_subset['Exposure'] < exposure[0][1])]

        # Select on reduction if provided.
        if reduction:
            self.data_subset = self.data_subset.loc[
                self.data_subset['Reduction'] == reduction]

        # Select on observatory if provided.
        if observatory:
            self.data_subset = self.data_subset.loc[
                self.data_subset['Observatory'] == observatory]

        # Select on jd if provided.
        if jd[0]:
            if jd[1] == 'exact':
                if isinstance(jd[0], list):
                    self.data_subset = self.data_subset.loc[
                        self.data_subset['JD'].isin(jd[0])]
                else:
                    self.data_subset = self.data_subset.loc[
                        self.data_subset['JD'] == jd[0]]
            elif jd[1] == 'except':
                if isinstance(jd[0], list):
                    self.data_subset = self.data_subset.drop(
                        self.data_subset[self.data_subset[
                            'JD'].isin(jd[0])].index)
                else:
                    self.data_subset = self.data_subset.drop(
                        self.data_subset[self.data_subset['JD']
                                         == jd[0]].index)
            elif jd[1] == 'after':
                self.data_subset = self.data_subset.loc[
                    self.data_subset['JD'] > jd[0]]
            elif jd[1] == 'before':
                self.data_subset = self.data_subset.loc[
                    self.data_subset['JD'] < jd[0]]
            elif jd[1] == 'between':
                self.data_subset = self.data_subset.loc[
                    (self.data_subset['JD'] > jd[0][0]) &
                    (self.data_subset['JD'] < jd[0][1])]

        # Sort by JD and reset index.
        self.data_subset = self.data_subset.sort_values(
            by=['JD'], ascending=True).reset_index(drop=True)

        # Check spectra found:
        if len(self.data_subset) == 0:
            raise FileExistsError(
                'No spectra found for that location or subset.')

        print('Selecting {} spectra from fits dataset.'.format(
            len(self.data_subset)))
