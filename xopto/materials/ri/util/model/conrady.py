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


class Conrady_1(Model):
    def __init__(self, params, pp: Scale or Normalize = None, **kwargs):
        '''
        Conrady 1 model of the refractive index.

        :math:`n^{2} = A_{1} + A_{2}/\\lambda + A_{3}/\\lambda^{3.5}`

        Parameters
        ----------
        params: np.ndarray
            Default model parameters.
        pp: xopto.materials.ri.util.model.base.Scale or xopto.materials.ri.util.model.base.Normalize
            Wavelength preprocessor instance.
        kwargs: dict
            Parameters passed to the baseclass.
        '''
        
        formatstr = 'n = {:.8e} + {:.8e}/wn +{:.8e}/wn**3.5'
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
        return [n[0], 0.0, 0.0]

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
        wn = self.pp(wavelengths)
        return params[0] + params[1]/wn + params[2]/wn**3.5


class ConradyEx_1(Model):
    def __init__(self, params, pp: Scale or Normalize = None, **kwargs):
        '''
        Extended Conrady 1 model of the refractive index.

        :math:`n^{2} = A_{1} + A_{2}/\\lambda + A_{3}/\\lambda^{A_{4}}`

        Parameters
        ----------
        params: np.ndarray
            Default model parameters.
        pp: xopto.materials.ri.util.model.base.Scale or xopto.materials.ri.util.model.base.Normalize
            Wavelength preprocessor instance.
        kwargs: dict
            Parameters passed to the baseclass.
        '''
        
        formatstr = 'n = {:.8e} + {:.8e}/wn +{:.8e}/wn**{:.8e}'
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
        return [n[0], 0.0, 0.0, 3.5]

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
        wn = self.pp(wavelengths)
        return params[0] + params[1]/wn + params[2]/wn**params[3]


class Conrady_2(Model):
    def __init__(self, params, pp: Scale or Normalize = None, **kwargs):
        '''
        Conrady 2 model of the refractive index.

        :math:`n^{2} = A_{1} + A_{2}/\\lambda^{2} + A_{3}/\\lambda^{3.5}`

        Parameters
        ----------
        params: np.ndarray
            Default model parameters.
        pp: xopto.materials.ri.util.model.base.Scale or xopto.materials.ri.util.model.base.Normalize
            Wavelength preprocessor instance.
        kwargs: dict
            Parameters passed to the baseclass.
        '''
        
        formatstr = 'n = {:.8e} + {:.8e}/wn**2 +{:.8e}/wn**3.5'
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
        return [n[0], 0.0, 0.0]

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
        wn = self.pp(wavelengths)
        return params[0] + params[1]/wn**2 + params[2]/wn**3.5


class ConradyEx_2(Model):
    def __init__(self, params, pp: Scale or Normalize = None, **kwargs):
        '''
        Extended Conrady 2 model of the refractive index.

        :math:`n^{2} = A_{1} + A_{2}/\\lambda^{2} + A_{3}/\\lambda^{A_{4}}`

        Parameters
        ----------
        params: np.ndarray
            Default model parameters.
        pp: xopto.materials.ri.util.model.base.Scale or xopto.materials.ri.util.model.base.Normalize
            Wavelength preprocessor instance.
        kwargs: dict
            Parameters passed to the baseclass.
        '''
    
        formatstr = 'n = {:.8e} + {:.8e}/wn**2 +{:.8e}/wn**{:.8e}'
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
        return [n[0], 0.0, 0.0, 3.5]

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
        wn = self.pp(wavelengths)
        return params[0] + params[1]/wn**2 + params[2]/wn**params[3]
