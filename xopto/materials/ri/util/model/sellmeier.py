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

import numpy as np

from .base import Scale, Normalize, Model


def _ri(params: np.ndarray, pp: Scale or Normalize,
        wavelengths: np.ndarray, num_terms: int) -> np.ndarray:
    '''
    Compute the refractive index for the given model parameters
    and wavelengths of light using the Sellmeier equation.

    Parameters
    ----------
    params: np.ndarray
        Model parameters.
    pp: Scale or Normalize
        Wavelength preprocessor.
    wavelengths: np.ndarray
        Wavelengths of light (m).
    num_terms: int
        Number of terms.

    Returns
    -------
    n: np.ndarray
        Refractive index estimated for the given model parameters and
        at the given wavelengths of light. 
    '''
    wn = pp(wavelengths)
    inv_wn2 = 1.0/(wn)**2

    result = None
    for term in range(num_terms):
        if term == 0:
            result = params[2*term]/(1.0 - params[2*term + 1]*inv_wn2)
        else:
            result += params[2*term]/(1.0 - params[2*term + 1]*inv_wn2)
    result += 1.0

    return np.sqrt(result)

def _ri_ex(params: np.ndarray, pp: Scale or Normalize,
        wavelengths: np.ndarray, num_terms: int) -> np.ndarray:
    '''
    Compute the refractive index for the given model parameters
    and wavelengths of light using the extended Sellmeier equation.

    Parameters
    ----------
    params: np.ndarray
        Model parameters.
    pp: Scale or Normalize
        Wavelength preprocessor.
    wavelengths: np.ndarray
        Wavelengths of light (m).
    num_terms: int
        Number of terms.

    Returns
    -------
    n: np.ndarray
        Refractive index estimated for the given model parameters and
        at the given wavelengths of light. 
    '''
    wn = pp(wavelengths)
    inv_wn2 = 1.0/(wn)**2

    result = None
    for term in range(num_terms):
        if term == 0:
            result = params[2*term + 1]/(1.0 - params[2*term + 2]*inv_wn2)
        else:
            result += params[2*term + 1]/(1.0 - params[2*term + 2]*inv_wn2)
    result += params[0]

    return np.sqrt(result)

class Sellmeier_1(Model):
    def __init__(self, params, pp: Scale or Normalize = None, **kwargs):
        '''
        One term Sellmeier model of the refractive index.

        :math:`n^{2} = 1 + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1})`

        Parameters
        ----------
        params: np.ndarray
            Default model parameters.
        pp: Scale or Normalize
            Wavelength preprocessor instance.
        kwargs: dict
            Parameters passed to the baseclass.
        '''
        formatstr = 'n*n = 1.0 + {:.8e}/(1.0 - {:.8e}/wn**2'
        super().__init__(formatstr, params, pp, **kwargs)

    def guess(self, wavelengths: np.ndarray, n: np.ndarray) -> list:
        '''
        Returns an initial guess of the model parameters for optimization/fit.

        Parameters
        ----------
        wavelengths: np.ndarray
            The wavelengths of light at which the values of refractive index
            are defined.
        n: np.ndarray
            The values of refractive index at the given wavelengths of light.

        Returns
        -------
        params0: np.ndarray
            Initial guess for the values of the model parameters
        '''
        return [0.0]*2

    def ri(self, params: np.ndarray, wavelengths: np.ndarray):
        '''
        Compute the refractive index for the given model parameters
        and wavelengths of light.

        Parameters
        ----------
        params: np.ndarray
            Model parameters.
        wavelengths: np.ndarray
            Wavelengths of light (m).

        Returns
        -------
        n: np.ndarray
            Refractive index estimated for the given model parameters and
            at the given wavelengths of light. 
        '''
        return _ri(params, self.pp, wavelengths, 1)


class SellmeierEx_1(Model):
    def __init__(self, params, pp: Scale or Normalize = None, **kwargs):
        '''
        Extended one term Sellmeier model of the refractive index.

        :math:`n^{2} = A + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1})`

        Parameters
        ----------
        params: np.ndarray
            Default model parameters.
        pp: Scale or Normalize
            Wavelength preprocessor instance.
        kwargs: dict
            Parameters passed to the baseclass.
        '''
        formatstr = 'n*n = {:.8e} + {:.8e}/(1.0 - {:.8e}/wn**2)'
        super().__init__(formatstr, params, pp, **kwargs)

    def guess(self, wavelengths: np.ndarray, n: np.ndarray) -> list:
        '''
        Returns an initial guess of the model parameters for optimization/fit.

        Parameters
        ----------
        wavelengths: np.ndarray
            The wavelengths of light at which the values of refractive index
            are defined.
        n: np.ndarray
            The values of refractive index at the given wavelengths of light.

        Returns
        -------
        params0: np.ndarray
            Initial guess for the values of the model parameters
        '''
        return [1.0] + [0.0]*2

    def ri(self, params: np.ndarray, wavelengths: np.ndarray):
        '''
        Compute the refractive index for the given model parameters
        and wavelengths of light.

        Parameters
        ----------
        params: np.ndarray
            Model parameters.
        wavelengths: np.ndarray
            Wavelengths of light (m).

        Returns
        -------
        n: np.ndarray
            Refractive index estimated for the given model parameters and
            at the given wavelengths of light. 
        '''
        return _ri_ex(params, self.pp, wavelengths, 1)


class Sellmeier_2(Model):
    def __init__(self, params, pp: Scale or Normalize = None, **kwargs):
        '''
        Two term Sellmeier model of the refractive index.

        :math:`n^{2} = 1 + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2})`

        Parameters
        ----------
        params: np.ndarray
            Default model parameters.
        pp: Scale or Normalize
            Wavelength preprocessor instance.
        kwargs: dict
            Parameters passed to the baseclass.
        '''
        formatstr = 'n*n = 1.0 + ' \
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2)'
        super().__init__(formatstr, params, pp, **kwargs)

    def guess(self, wavelengths: np.ndarray, n: np.ndarray) -> list:
        '''
        Returns an initial guess of the model parameters for optimization/fit.

        Parameters
        ----------
        wavelengths: np.ndarray
            The wavelengths of light at which the values of refractive index
            are defined.
        n: np.ndarray
            The values of refractive index at the given wavelengths of light.

        Returns
        -------
        params0: np.ndarray
            Initial guess for the values of the model parameters
        '''
        return [0.0]*4

    def ri(self, params: np.ndarray, wavelengths: np.ndarray):
        '''
        Compute the refractive index for the given model parameters
        and wavelengths of light.

        Parameters
        ----------
        params: np.ndarray
            Model parameters.
        wavelengths: np.ndarray
            Wavelengths of light (m).

        Returns
        -------
        n: np.ndarray
            Refractive index estimated for the given model parameters and
            at the given wavelengths of light. 
        '''
        return _ri(params, self.pp, wavelengths, 2)


class SellmeierEx_2(Model):
    def __init__(self, params, pp: Scale or Normalize = None, **kwargs):
        '''
        Extended two term Sellmeier model of the refractive index.

        :math:`n^{2} = A + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2})`

        Parameters
        ----------
        params: np.ndarray
            Default model parameters.
        pp: Scale or Normalize
            Wavelength preprocessor instance.
        kwargs: dict
            Parameters passed to the baseclass.
        '''
        formatstr = 'n*n = {:.8e} + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2)'
        super().__init__(formatstr, params, pp, **kwargs)

    def guess(self, wavelengths: np.ndarray, n: np.ndarray) -> list:
        '''
        Returns an initial guess of the model parameters for optimization/fit.

        Parameters
        ----------
        wavelengths: np.ndarray
            The wavelengths of light at which the values of refractive index
            are defined.
        n: np.ndarray
            The values of refractive index at the given wavelengths of light.

        Returns
        -------
        params0: np.ndarray
            Initial guess for the values of the model parameters
        '''
        return [1.0] + [0.0]*4

    def ri(self, params: np.ndarray, wavelengths: np.ndarray):
        '''
        Compute the refractive index for the given model parameters
        and wavelengths of light.

        Parameters
        ----------
        params: np.ndarray
            Model parameters.
        wavelengths: np.ndarray
            Wavelengths of light (m).

        Returns
        -------
        n: np.ndarray
            Refractive index estimated for the given model parameters and
            at the given wavelengths of light. 
        '''
        return _ri_ex(params, self.pp, wavelengths, 2)


class Sellmeier_3(Model):
    def __init__(self, params, pp: Scale or Normalize = None, **kwargs):
        '''
        Three term Sellmeier model of the refractive index.

        :math:`n^{2} = 1 + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2}) + B_{3}\\lambda^{2}/(\\lambda^{2} - C_{3})`

        Parameters
        ----------
        params: np.ndarray
            Default model parameters.
        pp: Scale or Normalize
            Wavelength preprocessor instance.
        kwargs: dict
            Parameters passed to the baseclass.
        '''
        formatstr = 'n*n = 1.0 + ' \
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2)'
        super().__init__(formatstr, params, pp, **kwargs)

    def guess(self, wavelengths: np.ndarray, n: np.ndarray) -> list:
        '''
        Returns an initial guess of the model parameters for optimization/fit.

        Parameters
        ----------
        wavelengths: np.ndarray
            The wavelengths of light at which the values of refractive index
            are defined.
        n: np.ndarray
            The values of refractive index at the given wavelengths of light.

        Returns
        -------
        params0: np.ndarray
            Initial guess for the values of the model parameters
        '''
        return [0.0]*6

    def ri(self, params: np.ndarray, wavelengths: np.ndarray):
        '''
        Compute the refractive index for the given model parameters
        and wavelengths of light.

        Parameters
        ----------
        params: np.ndarray
            Model parameters.
        wavelengths: np.ndarray
            Wavelengths of light (m).

        Returns
        -------
        n: np.ndarray
            Refractive index estimated for the given model parameters and
            at the given wavelengths of light. 
        '''
        return _ri(params, self.pp, wavelengths, 3)


class SellmeierEx_3(Model):
    def __init__(self, params, pp: Scale or Normalize = None, **kwargs):
        '''
        Extended three term Sellmeier model of the refractive index.

        :math:`n^{2} = A + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2}) + B_{3}\\lambda^{2}/(\\lambda^{2} - C_{3})`

        Parameters
        ----------
        params: np.ndarray
            Default model parameters.
        pp: Scale or Normalize
            Wavelength preprocessor instance.
        kwargs: dict
            Parameters passed to the baseclass.
        '''
        formatstr = 'n*n = {:.8e} + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2)'
        super().__init__(formatstr, params, pp, **kwargs)

    def guess(self, wavelengths: np.ndarray, n: np.ndarray) -> list:
        '''
        Returns an initial guess of the model parameters for optimization/fit.

        Parameters
        ----------
        wavelengths: np.ndarray
            The wavelengths of light at which the values of refractive index
            are defined.
        n: np.ndarray
            The values of refractive index at the given wavelengths of light.

        Returns
        -------
        params0: np.ndarray
            Initial guess for the values of the model parameters
        '''
        return [1.0] + [0.0]*6

    def ri(self, params: np.ndarray, wavelengths: np.ndarray):
        '''
        Compute the refractive index for the given model parameters
        and wavelengths of light.

        Parameters
        ----------
        params: np.ndarray
            Model parameters.
        wavelengths: np.ndarray
            Wavelengths of light (m).

        Returns
        -------
        n: np.ndarray
            Refractive index estimated for the given model parameters and
            at the given wavelengths of light. 
        '''
        return _ri_ex(params, self.pp, wavelengths, 3)


class Sellmeier_4(Model):
    def __init__(self, params, pp: Scale or Normalize = None, **kwargs):
        '''
        Four term Sellmeier model of the refractive index.

        :math:`n^{2} = 1 + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2}) + B_{3}\\lambda^{2}/(\\lambda^{2} - C_{3} + B_{4}\\lambda^{2}/(\\lambda^{2} - C_{4})`

        Parameters
        ----------
        params: np.ndarray
            Default model parameters.
        pp: Scale or Normalize
            Wavelength preprocessor instance.
        kwargs: dict
            Parameters passed to the baseclass.
        '''
        formatstr = 'n*n = 1.0 + ' \
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2)'
        super().__init__(formatstr, params, pp, **kwargs)

    def guess(self, wavelengths: np.ndarray, n: np.ndarray) -> list:
        '''
        Returns an initial guess of the model parameters for optimization/fit.

        Parameters
        ----------
        wavelengths: np.ndarray
            The wavelengths of light at which the values of refractive index
            are defined.
        n: np.ndarray
            The values of refractive index at the given wavelengths of light.

        Returns
        -------
        params0: np.ndarray
            Initial guess for the values of the model parameters
        '''
        return [0.0]*8

    def ri(self, params: np.ndarray, wavelengths: np.ndarray):
        '''
        Compute the refractive index for the given model parameters
        and wavelengths of light.

        Parameters
        ----------
        params: np.ndarray
            Model parameters.
        wavelengths: np.ndarray
            Wavelengths of light (m).

        Returns
        -------
        n: np.ndarray
            Refractive index estimated for the given model parameters and
            at the given wavelengths of light. 
        '''
        return _ri(params, self.pp, wavelengths, 4)


class SellmeierEx_4(Model):
    def __init__(self, params, pp: Scale or Normalize = None, **kwargs):
        '''
        Four term Sellmeier model of the refractive index.

        :math:`n^{2} = A + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2}) + B_{3}\\lambda^{2}/(\\lambda^{2} - C_{3} + B_{4}\\lambda^{2}/(\\lambda^{2} - C_{4})`

        Parameters
        ----------
        params: np.ndarray
            Default model parameters.
        pp: Scale or Normalize
            Wavelength preprocessor instance.
        kwargs: dict
            Parameters passed to the baseclass.
        '''
        formatstr = 'n*n = {:.8e} + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2)'
        super().__init__(formatstr, params, pp, **kwargs)

    def guess(self, wavelengths: np.ndarray, n: np.ndarray) -> list:
        '''
        Returns an initial guess of the model parameters for optimization/fit.

        Parameters
        ----------
        wavelengths: np.ndarray
            The wavelengths of light at which the values of refractive index
            are defined.
        n: np.ndarray
            The values of refractive index at the given wavelengths of light.

        Returns
        -------
        params0: np.ndarray
            Initial guess for the values of the model parameters
        '''
        return [1.0] + [0.0]*8

    def ri(self, params: np.ndarray, wavelengths: np.ndarray):
        '''
        Compute the refractive index for the given model parameters
        and wavelengths of light.

        Parameters
        ----------
        params: np.ndarray
            Model parameters.
        wavelengths: np.ndarray
            Wavelengths of light (m).

        Returns
        -------
        n: np.ndarray
            Refractive index estimated for the given model parameters and
            at the given wavelengths of light. 
        '''
        return _ri_ex(params, self.pp, wavelengths, 4)


class Sellmeier_5(Model):
    def __init__(self, params, pp: Scale or Normalize = None, **kwargs):
        '''
        Five term Sellmeier model of the refractive index.

        :math:`n^{2} = 1 + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2}) + B_{3}\\lambda^{2}/(\\lambda^{2} - C_{3} + B_{4}\\lambda^{2}/(\\lambda^{2} - C_{4} + B_{5}\\lambda^{2}/(\\lambda^{2} - C_{5})`

        Parameters
        ----------
        params: np.ndarray
            Default model parameters.
        pp: Scale or Normalize
            Wavelength preprocessor instance.
        kwargs: dict
            Parameters passed to the baseclass.
        '''
        formatstr = 'n*n = 1.0 + ' \
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2)'
        super().__init__(formatstr, params, pp, **kwargs)

    def guess(self, wavelengths: np.ndarray, n: np.ndarray) -> list:
        '''
        Returns an initial guess of the model parameters for optimization/fit.

        Parameters
        ----------
        wavelengths: np.ndarray
            The wavelengths of light at which the values of refractive index
            are defined.
        n: np.ndarray
            The values of refractive index at the given wavelengths of light.

        Returns
        -------
        params0: np.ndarray
            Initial guess for the values of the model parameters
        '''
        return [0.0]*10

    def ri(self, params: np.ndarray, wavelengths: np.ndarray):
        '''
        Compute the refractive index for the given model parameters
        and wavelengths of light.

        Parameters
        ----------
        params: np.ndarray
            Model parameters.
        wavelengths: np.ndarray
            Wavelengths of light (m).

        Returns
        -------
        n: np.ndarray
            Refractive index estimated for the given model parameters and
            at the given wavelengths of light. 
        '''
        return _ri(params, self.pp, wavelengths, 5)


class SellmeierEx_5(Model):
    def __init__(self, params, pp: Scale or Normalize = None, **kwargs):
        '''
        Five term Sellmeier model of the refractive index.

        :math:`n^{2} = A + B_{1}\\lambda^{2}/(\\lambda^2 - C_{1}) + B_{2}\\lambda^{2}/(\\lambda^{2} - C_{2}) + B_{3}\\lambda^{2}/(\\lambda^{2} - C_{3} + B_{4}\\lambda^{2}/(\\lambda^{2} - C_{4} + B_{5}\\lambda^{2}/(\\lambda^{2} - C_{5})`

        Parameters
        ----------
        params: np.ndarray
            Default model parameters.
        pp: Scale or Normalize
            Wavelength preprocessor instance.
        kwargs: dict
            Parameters passed to the baseclass.
        '''
        formatstr = 'n*n = {:.8e} + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2) + '\
                    '{:.8e}/(1.0 - {:.8e}/wn**2)'
        super().__init__(formatstr, params, pp, **kwargs)

    def guess(self, wavelengths: np.ndarray, n: np.ndarray) -> list:
        '''
        Returns an initial guess of the model parameters for optimization/fit.

        Parameters
        ----------
        wavelengths: np.ndarray
            The wavelengths of light at which the values of refractive index
            are defined.
        n: np.ndarray
            The values of refractive index at the given wavelengths of light.

        Returns
        -------
        params0: np.ndarray
            Initial guess for the values of the model parameters
        '''
        return [1.0] + [0.0]*10

    def ri(self, params: np.ndarray, wavelengths: np.ndarray):
        '''
        Compute the refractive index for the given model parameters
        and wavelengths of light.

        Parameters
        ----------
        params: np.ndarray
            Model parameters.
        wavelengths: np.ndarray
            Wavelengths of light (m).

        Returns
        -------
        n: np.ndarray
            Refractive index estimated for the given model parameters and
            at the given wavelengths of light. 
        '''
        return _ri_ex(params, self.pp, wavelengths, 5)
