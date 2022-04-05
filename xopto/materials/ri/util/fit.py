# -*- coding: utf-8 -*-
################################ Begin license #################################
# Copyright (C) Laboratory of Imaging technologies,
#               Faculty of Electrical Engineering,
#               University of Ljubljana.
#
# This file is part of PyXOpto.
#
# PyXOpto is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PyXOpto is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PyXOpto. If not, see <https://www.gnu.org/licenses/>.
################################# End license ##################################

from typing import Tuple, List
import sys

from scipy.optimize import minimize
import numpy as np

from xopto.materials.ri.util import model


class Fit:
    def __init__(self, kind: model.Model,
                 wavelengths: np.ndarray, n: np.ndarray,
                 x0: Tuple or List = None, verbose: bool = False):
        '''
        Fit one of the models to the measured values of the refractive index
        as a function of wavelength.

        Parameters
        ----------
        kind: model.Model
            Model of the refractive index as a function of wavelength.
            Use one of the predefined models:

            - "Conrady_1': Conrady 1 model :math:`n^{2} = A_{1} + A_{2}/\\lambda + A_{3}/\\lambda^{3.5}`
            - "Conrady_2': Conrady 2 model :math:`n^{2} = A_{1} + A_{2}/\\lambda^{2} + A_{3}/\\lambda^{3.5}`
            - "ConradyEx_1': Extended Conrady 1 model :math:`n^{2} = A_{1} + A_{2}/\\lambda + A_{3}/\\lambda^{A_{4}}`
            - "ConradyEx_2': Extended Conrady 2 model :math:`n^{2} = A_{1} + A_{2}/\\lambda^{2} + A_{3}/\\lambda^{A_{4}}`
            - "Exponential": Exponential model :math:`n = a + b*e^{-\\lambda/c}`
            - "Cauchy': Cauchy’s model :math:`n^{2} = A_{1} + A_{2}\\lambda^{2} + A_{3}/\\lambda^{2} + A_{4}/\\lambda^{4} + A_{5}/\\lambda^{6} + A_6/\\lambda^{8}`
            - "Sellmeier_1": 1 term Sellmeier formula :math:`n^{2} = 1 + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1})`
            - "SellmeierEx_1": 1 term extended Sellmeier formula :math:`n^{2} = A + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1})`
            - "Sellmeier_2": 2 term Sellmeier formula :math:`n^{2} = 1 + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2})`
            - "SellmeierEx_2": 2 term extended Sellmeier formula :math:`n^{2} = A + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2})`
            - "Sellmeier_3": 3 term Sellmeier formula :math:`n^{2} = 1 + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2}) + B_{3}\\lambda^{2}/(\\lambda^{2} - C_{3})`
            - "SellmeierEx_3": 3 term extended Sellmeier formula :math:`n^{2} = A + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2}) + B_{3}\\lambda^{2}/(\\lambda^{2} - C_{3})`
            - "Sellmeier_4": 4 term Sellmeier formula :math:`n^{2} = 1 + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2}) + B_{3}\\lambda^{2}/(\\lambda^{2} - C_{3} + B_{4}\\lambda^{2}/(\\lambda^{2} - C_{4})`
            - "SellmeierEx_4": 4 term extended Sellmeier formula :math:`n^{2} = A + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2}) + B_{3}\\lambda^{2}/(\\lambda^{2} - C_{3} + B_{4}\\lambda^{2}/(\\lambda^{2} - C_{4})`
            - "Sellmeier_5": 5 term Sellmeier formula :math:`n^{2} = 1 + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2}) + B_{3}\\lambda^{2}/(\\lambda^{2} - C_{3} + B_{4}\\lambda^{2}/(\\lambda^{2} - C_{4} + B_{5}\\lambda^{2}/(\\lambda^{2} - C_{5})`
            - "SellmeierEx_5": 5 term extended Sellmeier formula :math:`n^{2} = A + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2}) + B_{3}\\lambda^{2}/(\\lambda^{2} - C_{3} + B_{4}\\lambda^{2}/(\\lambda^{2} - C_{4} + B_{5}\\lambda^{2}/(\\lambda^{2} - C_{5})`
            - "Herzberger_3_2": Herzberger 4 x 2 formula :math:`n = A + B\\lambda^{2} + C\\lambda^{4} + D/(\\lambda^{2} - 0.028) + E/(\\lambda^{2} - 0.028)^{2}`
            - "HerzbergerEx_3_2": Extended Herzberger 4 x 2 formula :math:`n = A + B\\lambda^{2} + C\\lambda^{4} + D/(\\lambda^{2} - E) + F/(\\lambda^{2} - G)^{2}`
            - "Herzberger_4_2": Herzberger 4 x 2 formula :math:`n = A + B\\lambda^{2} + C\\lambda^{4} + D\\lambda^{6} + E/(\\lambda^{2} - 0.028) + F/(\\lambda^{2} - 0.028)^{2}`
            - "HerzbergerEx_4_2": Extended Herzberger 4 x 2 formula :math:`n = A + B\\lambda^{2} + C\\lambda^{4} + D\\lambda^{6} + E/(\\lambda^{2} - F) + G/(\\lambda^{2} - H)^{2}`

        wavelengths: np.ndarray vector
            Wavelengths at which the refractive index was measured. The
            wavelengths must be ordered in ascending order.
        n: np.ndarray vector
            The measured values of refractive index at the specified
            wavelengths.
        x0: list, tuple or np.ndarray vector
            Initial guess of the model parameters.
        verbose: bool
            Verbose mode.

        Note
        ----
        For better convergence of the fit process use an appropriate
        preprocessor with the supplied refractive index model.
        A preprocessor that scales the wavelengths by the inverse of the first
        (shortest) wavelength usually produces good results, e.g.
        :code:`Sellmeier_1(pp=Scale(1.0/wavelengths[0]))`
        '''
        if not isinstance(kind, model.Model):
            raise TypeError(
                'Unsupported refractive index model type: "{}"! ' \
                'Must be an instance of model.Model!'.format(
                    kind.__class__.__name__)
            )

        self._verbose = bool(verbose)

        wavelengths = np.asarray(wavelengths, dtype=np.float64)
        n = np.asarray(n, dtype=np.float64)

        if x0 is None:
            x0 = kind.guess(wavelengths, n)

        self._fit_model = kind

        if self._verbose:
            print('Initial condition x0:', x0)

        self._wavelengths = np.asarray(wavelengths, dtype=np.float64)
        self._n = np.asarray(n, dtype=np.float64)

        results = minimize(
            lambda x: (
                (self._fit_model(self._wavelengths, params=x) - self._n)**2
            ).sum(),
            x0, method='L-BFGS-B')
        self._x = results['x']

        if self._verbose:
            print('Results of the model fit:', self._x)

        self._model = type(self._fit_model)(self._x, self._fit_model.pp)

    def _get_wavelengths(self) -> np.ndarray:
        return self._wavelengths

    wavelengths = property(_get_wavelengths, None, None,
        'Wavelengths at which the refractive index was measured '
        'Returns values as passed to the constructor.'
    )

    def _get_n(self) -> np.ndarray:
        return self._n

    n = property(
        _get_n, None, None,
        'The measured values of refractive index at the specified wavelengths. '
        'The values are returned as passed to the constructor.'
    )

    def error(self) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Fit error at the wavelengths passed to the constructor.

        Returns
        -------
        err: np.ndarray of length wavelengths.size
            Fit errors defined as :math:`(n_{fit} - n_{measured})`.
        rel_err: np.ndarray of length wavelengths.size
            Relative fit error defined as
            :math:`(n_{fit} - n_{measured})/n_{measured}`.
        '''
        n_fit = self(self._wavelengths)
        err = n_fit - self._n
        rerr = err/self._n

        return err, rerr

    def _get_model(self) -> model.Model:
        return self._model

    model = property(
        _get_model, None, None,
        'Model instance obtained by the fit to the measured data.'
    )
        

    def visualize(self, wavelengths: np.ndarray = None, show=False):
        '''
        Plots the quality of fit.

        Parameters
        ----------
        wavelengths: np.ndarray vector
            Wavelengths (m) at which to visualize the fit. If None, the
            wavelengths passed to the constructor are used.
        show: bool
            If nonzero, calls pp.show().
        '''
        from matplotlib import pyplot as pp

        if wavelengths is None:
            wavelengths = self._wavelengths

        err, rerr = self.error()

        f, (ax1, ax2) = pp.subplots(2, 1)
        ax1.plot(self._wavelengths*1e9, self._n, 'ok')
        n_visualize = self(wavelengths)
        ax1.plot(wavelengths*1e9, n_visualize)

        ax2.plot(self._wavelengths*1e9, 100.0*rerr, '-o')

        ax1.set_title('Model fit - "{}"'.format(self.model.name))
        ax1.grid()
        ax1.set_xlabel('Wavelength (nm)')

        ax2.set_title('Relative fit error as (fit - measured)/measured (%)')
        ax2.set_xlabel('Wavelength (nm)')
        ax2.grid()

        pp.tight_layout()

        if show:
            pp.show()

    def __call__(self, wavelengths: float or np.ndarray) -> float or np.ndarray:
        '''
        Evaluate the refractive index model at the given wavelengths of light.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength of light (m).
        '''
        return self._model(wavelengths)

    def __str__(self):
        return 'Fit of {:s} model:\n  \"{:s}\"\nwhere wn is:\n  \"{:}\"'.format(
            self._model.name, *self._model.render())

    def __repr__(self):
        return '{} # id {}'.format(self.__str__(), id(self))


class ThermalFit:
    def __init__(self, kind: model.Model,
                 wavelengths: np.ndarray, n: np.ndarray,
                 temperatures: np.ndarray, pk_order: int = 2,
                 x0: Tuple or List[float] = None):
        '''
        Fit a temperature dependent refractive index model. Each parameter of
        the selected refractive index model is modelled as a polynomial
        function of the temperature.

        Parameters
        ----------
        kind: model.Model
            Refractive index model instance.
        wavelengths: np.ndarray vector
            Wavelengths of light (m). The first wavelength in the vector
            is used by the model preprocessor to normalize/divide all
            the wavelengths.
        n: numpy_ndarray of shape (wavelengths.size, temperatures.size)
            The measured refractive indices.
        temperatures: np.ndarray vector
            Temperatures (should be in K, but accepts any units).
        pk_order: int
            Order of the polynomial that models the temperature dependence of
            the refractive index model parameters.
        x0: list, tuple ndarray vector
            Initial guess of the refractive index model parameters.
        '''
        self._fit_model = kind 

        self._wavelengths = wavelengths = np.asarray(
            wavelengths, dtype=np.float64)
        self._temperatures = temperatures = np.asarray(
            temperatures, dtype=np.float64)
        self._n = n = np.asarray(n, dtype=np.float64)

        fits = []
        params = []
        for index, temperature in enumerate(temperatures):
            fits.append(Fit(kind, wavelengths, n[:, index], x0=x0))
            params.append(fits[-1].model.params)

        params = np.vstack(params)
        self._fits = fits

        self._ri_model = self._fits[0].model

        # polyfit for temperatures
        # polynomial coefficients of each model parmeter are in columns of pk
        self._pk = np.polyfit(temperatures, params, pk_order)

    def error(self) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Fit error at the wavelengths and temperatures passed to the constructor.

        Returns
        -------
        err: np.ndarray of shape (wavelengths.size, temperatures.size)
            Fit errors defined as :math:`(n_{fit} - n_{measured})`.
        rel_err: np.ndarray of shape (wavelengths.size, temperatures.size)
            Relative fit error defined as
            :math:`(n_{fit} - n_{measured})/n_{measured}`.
        '''
        n_fit = np.zeros((self._wavelengths.size, self._temperatures.size))
        for index, temperature in enumerate(self._temperatures):
            n_fit[:, index] = self(self._wavelengths, temperature)
        err = n_fit - self._n
        rel_err = err/self._n

        return err, rel_err

    def visualize(self, wavelengths: np.ndarray = None, show: bool = False):
        '''
        Plots the quality of fit.

        Parameters
        ----------
        wavelengths: np.ndarray vector
            Wavelengths (m) at which to visualize the fit. If None, the
            wavelengths passed to the constructor are used.
        show: bool
            If nonzero, calls pp.show().
        '''
        from matplotlib import pyplot as pp

        if wavelengths is None:
            wavelengths = self._wavelengths

        err, rerr = self.error()

        f, (ax1, ax2) = pp.subplots(2, 1)
        for index, temperature in enumerate(self._temperatures):
            ax1.plot(self._wavelengths*1e9, self._n, 'ok')
            n_visualize = self(wavelengths, temperature)
            ax1.plot(wavelengths*1e9, n_visualize,
                     label='T {:.1f}'.format(temperature))

            ax2.plot(self._wavelengths*1e9, 100.0*rerr[:, index], '-o')

        ax1.set_title('Model fit - "{}"'.format(self._ri_model.name))
        ax1.grid()
        ax1.set_xlabel('Wavelength (nm)')
        ax1.legend()
        ax2.set_title('Relative fit error as (fit - measured)/measued (%)')
        ax2.set_xlabel('Wavelength (nm)')
        ax2.grid()

        pp.tight_layout()

        if show:
            pp.show()

    def polynomial_coefficients(self) -> np.ndarray:
        '''
        Returns the polynomial coefficients that model the temperature
        dependence of the refractive index model parameters.
        Use np.polyval to determine the values of polynomial coefficients
        at a given temperature (use the same temperture units as for
        the temperatures that were passed to the constructor).

        Returns
        -------
        pk: np.ndarray of shape (pk_order + 1, num_terms*2)
            The polynomial coefficients of each model parameter
            are stored in columns.
        '''
        return self._pk

    def _get_wavelengths(self) -> np.ndarray:
        return self._wavelengths

    property(_get_wavelengths, None, None,
        'Wavelengths of light at which the refractive index was measured. '\
        'The values are returned as passed to the constructor.'
    )

    def _get_n(self) -> np.ndarray:
        return self._n

    n = property(
        _get_n, None, None,
        'The measured values of refractive index as passed to the constructor.'
    )

    def _get__temperatures(self) -> np.ndarray:
        return self._temperatures

    temperatures = property(
        _get_n, None, None,
        'The temperatures at which the refractive index was measured. '
        'Returns the values as passed to the constructor.'
    )

    def __call__(self, wavelength: float, temperature: float = None):
        '''
        Computes the refractive index at the given wavelength and temperature.

        Parameters
        ----------
        wavelength: float
            Wavelength of light (m).
        temperature: float
            Temperature in same units as used for the constructor parameter
            temperatures. If None, the value of the first element of the
            parameter temperatures passed to the constructor is used.
        '''
        if temperature is None:
            temperature = self._temperatures[0]

        return self._ri_model(
            wavelength, params=np.polyval(self._pk, temperature)
        )

if __name__ == '__main__':
    from xopto.materials import ri

    w = np.linspace(400e-9, 800e-9, 5)
    n = ri.water.default(w)

    # Fit the wavelength dependence at a single temperature
    f = Fit(model.Conrady_2(None, model.Scale(1/w[0])), w, n)
    f.visualize(show=True)
    print(f.model.render())

    # Fit the temperature dependence - polynomial modelling of the 
    # refractive index model parameters.
    T = [293, 298, 303]
    n = np.array([ri.siliglass.default(w, temperature=T[0]),
                  ri.siliglass.default(w, temperature=T[1]),
                  ri.siliglass.default(w, temperature=T[2])]).T

    ft = ThermalFit(model.Conrady_2(None, model.Scale(1/w[0])), w, n, T)
    ft.visualize(show=True)
