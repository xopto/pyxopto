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

from typing import Tuple

import numpy as np

class Normalize:
    def __init__(self, src_range: Tuple[float, float],
                 dest_range: Tuple[float, float]=(-1.0, 1.0)):
        '''
        Linearly transforms the input data range so that they tightly fit the
        specified output range.

        Parameters
        ----------
        src_range: numpy.ndarray, list or tuple of 2 values
            The range of input values as [min, max]. Can also be a full
            array of values in which case the minimum and maximum values
            are taken as the range.
        dest_range: numpy.ndarray, list or tuple of 2 values
            The range of normalized values as [min, max].
        '''
        self._dest_range = (float(dest_range[0]), float(dest_range[1]))
        src_range = np.asarray(src_range, dtype=np.float64)
        self._src_range = (float(src_range.min()), float(src_range.max()))

        out_span = self._dest_range[1] - self._dest_range[0]
        in_span = self._src_range[1] - self._src_range[0]
        self._k = out_span/in_span

    def __call__(self, data: np.ndarray) -> np.ndarray:
        '''
        Applies the normalization/transformation to the data.

        Parameters
        ----------
        data: numpy.ndarray, list or tuple
            Data to transfor.

        Returns
        -------
        normalized: numpy.ndarray, list or tuple
            The input data transformed to the output range specified
            in the constructor call.
        '''
        data = np.asarray(data, dtype=np.float64)
        return self._dest_range[0] + (data - self._low)*self._k

    def undo(self, data: np.ndarray) -> np.ndarray:
        '''
        Apply inverse of the normalization/transformation to the data.

        Parameters
        ----------
        data: numpy.ndarray, list or tuple
            Data to transfor.

        Returns
        -------
        denormalized: numpy.ndarray, list or tuple
            The input data inversely transformed to the input range specified
            in the constructor call.
        '''
        data = np.asarray(data, dtype=np.float64)
        return (data - self._dest_range[0])/self._k + self._low

    def render(self, input='wavelength', output='wn'):
        '''
        Render the preprocessor equation to a string.

        Parameters
        ----------
        input: str
            Input symbol as a string.
        output: str
            Output symbol as a string.
        '''
        return '{:s} = {:.8e} + ({:s} - {:.8e})*{:.8e}'.format(
            output, self._dest_range[0], input, self._low, self._k)

    def __str__(self):
        return 'Normalize(values={}, interval={}) '\
               '# Transformation: "{} + (x - {})*{}"'.format(
                   self._values, self._interval,
                   self._dest_range[0], self._low, self._k)

    def __repr__(self):
        return '{} # id {}'.format(self.__str__(), id(self))


class Scale:
    def __init__(self, factor: float, offset: float = 0.0):
        '''
        Constructs an object that transforms the input data as
        (data - offset)*factor.

        Parameters
        ----------
        factor: float
            The multiplicative term.
        offset: float
            The additive term.
        '''
        self._factor = float(factor)
        self._offset = float(offset)

    def __call__(self, data: np.ndarray) -> np.ndarray:
        '''
        Applies the normalization/transformation to the data.

        Parameters
        ----------
        data: numpy.ndarray, list or tuple
            Data to transfor.

        Returns
        -------
        scaled: numpy.ndarray, list or tuple
            The input data shifted and scaled by the values specified
        '''
        return (np.asarray(data, dtype=np.float64) - self._offset)*self._factor

    def undo(self, data: np.ndarray) -> np.ndarray:
        '''
        Apply inverse of the transformation to the data.

        Parameters
        ----------
        data: numpy.ndarray, list or tuple
            Data to transfor.

        Returns
        -------
        descaled: numpy.ndarray, list or tuple
            The input data inversely transformed.
        '''
        return np.asarray(data, dtype=np.float64)/self._factor + self._offset

    def render(self, input='wavelength', output='wn'):
        '''
        Render the preprocessor equation to a string.

        Parameters
        ----------
        input: str
            Input symbol as a string.
        output: str
            Output symbol as a string.
        '''
        if self._offset == 0.0:
            res = '{:s} = {:s}*{:.8e}'.format(output, input, self._factor)
        else:
            res = '{:s} = ({:s} - {:.8e})*{:.8e}'.format(
                output, input, self._offset, self._factor)
        return res

    def __str__(self):
        return 'Scale(factor={offset}, offset={factor}) '\
               '# Transformation:  "(x - {offset})*{factor}"'.format(
                   offset=self._offset, factor=self._factor)

    def __repr__(self):
        return '{} # id {}'.format(self.__str__(), id(self))


class Model:
    def __init__(self, formatstr:str, params: np.ndarray or None = None,
                 pp: Scale or Normalize = None):
        '''
        Creates a refractive index model with the give wavelength preprocessor.

        Parameters
        ----------
        formatstr: str
            Format string that can be used to print the model equation.
        params: np.ndarray
            Default values of the model parameters.
        pp: Scale or Normalize
            Wavelength preprocesor.
        '''
        if params is not None:
            params = np.asarray(params)
        if pp is None:
            pp = Scale(1.0)

        self._pp = pp
        self._params = params
        self._formatstr = formatstr

    def guess(self, wavelengths: np.ndarray, n: np.ndarray) -> list:
        '''
        Returns an initial guess for optimization/fit.

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
        if self.params is None:
            raise RuntimeError('Cannot derive the number of model parameters!')

        return [0.0]*len(self.params)

    def _get_pp(self) -> Scale or Normalize:
        return self._pp

    pp = property(_get_pp, None, None,
                  'Preprocessor that scales/normalizes the wavelengths')

    def _get_params(self) -> Scale or Normalize:
        return self._params

    params = property(_get_params, None, None,
                     'Default parameter values of the model.')

    def _get_formatstr(self) -> str:
        return self._formatstr

    formatstr = property(_get_formatstr, None, None,
                         'Format string for rendering the model equation.')

    def __call__(self, wavelengths:np.ndarray, params=None) -> np.ndarray:
        '''
        Evaluate the model for the given wavelengths. Use the given parameters
        or default values from constructor if params is None.

        Parameters
        ----------
        wavelength: np.ndarray
            Wavelengths of light (m) at which to evaluate the model.
        params: np.ndarray or None
            Parameters of the model. If none, uses default values that were
            passed to the constructor.

        Returns
        -------
        n: np.ndarray
            Refractive index at the given wavelengths as estimated by the model.
        '''
        if params is None:
            params = self.params
        return self.ri(params, wavelengths)

    def _get_name(self):
        return self.__class__.__name__

    name = property(_get_name, None, None, 'Model name.')

    def render(self, params: np.ndarray = None) -> str:
        '''
        Render the model equation with the given set of parameters
        (use default parameter values passed to the constructor if None).

        Parameters
        ----------
        params: np.ndarray or None
            Model parameters. If None, use the default parameter values as
            passed to the constructor.

        Returns
        -------
        eq: str
            Model equation rendered as a string.
        pp: str
            Preprocessor equation rendered as a string.
        '''
        if params is None:
            if self.params is None:
                raise ValueError(
                    'Cannot derive the values of model parameters!')
            params = self.params

        return self._formatstr.format(*np.asarray(params).tolist()), \
            self.pp.render()
