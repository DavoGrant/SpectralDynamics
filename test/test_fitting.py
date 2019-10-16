import unittest
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from extract.helpers import construct_gaussian


class TestSpectralFitting(unittest.TestCase):
    """ Test spectral fitting. """

    def __init__(self, *args, **kwargs):
        super(TestSpectralFitting, self).__init__(*args, **kwargs)
        self.draw = False

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

    def test_curve_fit(self):
        """ Test curve fit. """
        # Synthetic two emission Gaussian profile.
        bins = 200
        wavelength = np.linspace(6450, 6550, bins)
        comp_1 = construct_gaussian(wavelength, 6500, 2, 100)
        comp_2 = construct_gaussian(wavelength, 6510, 3, 70)
        flux = comp_1 + comp_2
        flux_errors = np.sqrt(abs(flux))
        flux += np.random.normal(0, flux_errors, bins)

        # Fit.
        popt, pcov = curve_fit(
            self._multi_component_model,
            wavelength, flux, sigma=flux_errors,
            p0=(6503, 1.5, 90, 6508, 4, 75),
            bounds=((6490, 0, 0, 6500, 0, 0),
                    (6510, 5, 200, 6520, 6, 150)))

        total_components = []
        comp_params = np.split(np.array(popt), len(popt) / 3)
        for model_params in comp_params:
            comp = construct_gaussian(
                wavelength, model_params[0],
                model_params[1], model_params[2])
            total_components.append(comp.tolist())
        total_fit = np.asarray(total_components).sum(axis=0)

        if self.draw:
            # Observations
            plt.scatter(wavelength, flux)

            # Components.
            for c in total_components:
                plt.plot(wavelength, c)

            # Total fit.
            plt.plot(wavelength, total_fit)

            # Show.
            plt.tight_layout()
            plt.show()

        # Goodness of fit.
        chisq = (((flux - total_fit) / flux_errors) ** 2).sum()
        reduced_chisq = chisq / (len(wavelength) - len(popt))

        # Sanity checks.
        self.assertLess(reduced_chisq, 1.5)


if __name__ == '__main__':
    unittest.main()
