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
import builtins

import numpy as np

from xopto.mcbase import mcobject
from xopto.mcbase.mcworker import ClWorker


class PreprocessorFunction(mcobject.McObject):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def apply(self, x: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Base class of preprocessor functions that can be used to scale or
        normalize datasets.

        Parameters
        ----------
        x: np.ndarray
            Input data array.
        out: np.ndarray
            Data array for the results. Must be of the same
            size and type as the x data array.
        
        Returns
        -------
        y: np.ndarray
            Preprocessed input data.
        '''
        if out is None:
            out = np.array(x)
        else:
            out = np.copyto(out, x)

        return out

    def undo(self, y: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Undo preprocessing from the output data.

        Parameters
        ----------
        y: np.ndarray
            Preprocessed output data.
        out: np.ndarray
            Data array for the results. Must be of the same
            size and type as the y data array.

        Returns
        -------
        x: np.ndarray
            Input data obtained by applying inverse preprocessing to the
            output data.
        '''
        if out is None:
            out = np.array(y)
        else:
            np.copyto(out, y)

        return out

    def train(self, x: np.ndarray):
        '''
        Determine the required statistics of input values.

        Parameters
        ----------
        x: np.ndarray
            A representative training dataset.
        '''
        pass

    def todict(self) -> dict:
        '''
        Export instance to a dict object.

        Returns
        -------
        data: dict
            Instance exported to a dict object.
        '''
        return {'type': self.__class__.__name__}

    @classmethod
    def fromdict(cls, data: dict) -> 'Preprocessor':
        '''
        Create a new instance from data that were exported to a dict
        object with the
        :py:meth:`xopto.util.regression.pp.PreprocessorFunction.todict` method.
        '''
        data = dict(data)
        type_name = data.pop('type')
        T = globals().get(type_name, cls)
        return T(**data)

    def __repr__(self):
        return '{} # {}'.format(self.__str__(), id(self))


class LogScale(PreprocessorFunction):
    @staticmethod
    def cl_implementation(worker:  ClWorker) -> str:
        '''
        Implementation of the logarithmic scaling preprocessor in OpenCL.
        '''
        return '\n'.join((
            'inline void pp_logscale_apply(',
            '		mc_size_t start, mc_size_t stop, mc_size_t step, ',
            '		mc_fp_t *inout)',
            '{',
            '	for(mc_size_t pos = start; pos < stop; pos += step){',
            '		inout[pos] = mc_log(inout[pos]);',
            '	}',
            '};',
            '',
            'inline void pp_logscale_undo(',
            '		mc_size_t start, mc_size_t stop, mc_size_t step, ',
            '		mc_fp_t *inout)',
            '{',
            '	for(mc_size_t pos = start; pos < stop; pos += step){',
            '		inout[pos] = mc_exp(inout[pos]);',
            '	}',
            '};',
        ))

    def cl_render_call(self, worker: ClWorker, selector: slice,
                       inout_var: str, indent: str = None,
                       undo: bool = False) -> List[str]:
        '''
        Render preprocessor call to a list of strings, each holding one line
        of OpenCL code.

        Parameters
        ----------
        worker: ClWorker
            OpeCL worker instance.
        selector: slice
            Input data slice that selects the input items to preprocess.
        inout_var: str
            Name of the RW variable that holds the data.
        indent: str
            Indentation string. Defaults to 4 spaces.
        undo: bool
            If True render the undo call, else the apply call.

        Returns
        -------
        call: List[str]
            The rendered preprocessor call as a string.
        '''
        T = worker.types
        indent = '    ' if indent is None else indent
        return [
            'pp_logscale_{:s}('.format('undo' if undo else 'apply'),
            '{indent:s}{start:s}, {stop:s}, {step:s},'.format(
                start=T.size_t_str(selector.start),
                stop=T.size_t_str(selector.stop),
                step=T.size_t_str(selector.step), indent=indent),
            '{indent:s}{inout_var:s}'.format(
                inout_var=inout_var, indent=indent),
            ');'
        ]

    def apply(self, x: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Apply natural logarithm to the input data.

        Parameters
        ----------
        x: np.ndarray
            Input data array.
        out: np.ndarray
            Data array for the results. Must be of the same
            size and type as the x data array.
        
        Returns
        -------
        y: np.ndarray
            Logarithmically scaled input data.
        '''
        result = np.log(x, out=out)
        return result

    def undo(self, x: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Undo natural logarithm from the output data.

        Parameters
        ----------
        y: np.ndarray
            Preprocessed output data.
        out: np.ndarray
            Data array for the results. Must be of the same
            size and type as the y data array.

        Returns
        -------
        x: np.ndarray
            Input data obtained by applying inverse preprocessing to the
            output data.
        '''
        result = np.exp(x, out=out)
        return result

    def __str__(self):
        return 'LogScale()'


class Log10Scale(PreprocessorFunction):
    @staticmethod
    def cl_implementation(worker:  ClWorker) -> str:
        '''
        Implementation of the log10 scaling preprocessor in OpenCL.
        '''
        return '\n'.join((
            'inline void pp_log10scale_apply(',
            '		mc_size_t start, mc_size_t stop, mc_size_t step, ',
            '		mc_fp_t const *inout)',
            '{',
            '	for(mc_size_t pos = start; pos < stop; pos += step){',
            '		inout[pos] = mc_log10(inout[pos]);',
            '	}',
            '};',
            '',
            'inline void pp_log10scale_undo(',
            '		mc_size_t start, mc_size_t stop, mc_size_t step, ',
            '		mc_fp_t *inout)',
            '{',
            '	for(mc_size_t pos = start; pos < stop; pos += step){',
            '		inout[pos] = mc_pow(FP_10, inout[pos]);',
            '	}',
            '};',
        ))

    def cl_render_call(self, worker: ClWorker, selector: slice,
                       inout_var: str, indent: str = None,
                       undo: bool = False) -> List[str]:
        '''
        Render preprocessor call to a list of strings, each holding one line
        of OpenCL code.

        Parameters
        ----------
        worker: ClWorker
            OpeCL worker instance.
        selector: slice
            Input data slice that selects the input items to preprocess.
        inout_var: str
            Name of the RW variable that holds the data.
        indent: str
            Indentation string. Defaults to 4 spaces.
        undo: bool
            If True render the undo call, else the apply call.

        Returns
        -------
        call: List[str]
            The rendered preprocessor call as a string.
        '''
        T = worker.types
        indent = '    ' if indent is None else indent
        return [
            'pp_log10scale_{:s}('.format('undo' if undo else 'apply'),
            '{indent:s}{start:s}, {stop:s}, {step:s},'.format(
                start=T.size_t_str(selector.start),
                stop=T.size_t_str(selector.stop),
                step=T.size_t_str(selector.step), indent=indent),
            '{indent:s}{inout_var:s}'.format(
                inout_var=inout_var, indent=indent),
            ');'
        ]

    def apply(self, x: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Apply log10 function to the input data.

        Parameters
        ----------
        x: np.ndarray
            Input data array.
        out: np.ndarray
            Data array for the results. Must be of the same
            size and type as the x data array.
        
        Returns
        -------
        y: np.ndarray
            Logarithmically scaled input data.
        '''
        result = np.log10(x, out=out)
        return result

    def undo(self, x: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Undo log10 from the output data.

        Parameters
        ----------
        y: np.ndarray
            Preprocessed output data.
        out: np.ndarray
            Data array for the results. Must be of the same
            size and type as the y data array.

        Returns
        -------
        x: np.ndarray
            Input data obtained by applying inverse preprocessing to the
            output data.
        '''
        result = np.power(x, 10, out=out)
        return result

    def __str__(self):
        return 'Log10Scale()'


class ExpScale(PreprocessorFunction):
    @staticmethod
    def cl_implementation(worker:  ClWorker) -> str:
        '''
        Implementation of the log10 scaling preprocessor in OpenCL.
        '''
        return '\n'.join((
            'inline void pp_expscale_apply(',
            '		mc_size_t start, mc_size_t stop, mc_size_t step, ',
            '		mc_fp_t *inout){',
            '	for(mc_size_t pos = start; pos < stop; pos += step){',
            '		output[pos] = mc_exp(inout[pos]);',
            '	}',
            '};',
            '',
            'inline void pp_expscale_undo(',
            '		mc_size_t start, mc_size_t stop, mc_size_t step, ',
            '		mc_fp_t *inout){',
            '	for(mc_size_t pos = start; pos < stop; pos += step){',
            '		inout[pos] = mc_log(inout[pos]);',
            '	}',
            '};',
        ))

    def cl_render_call(self, worker: ClWorker, selector: slice,
                       inout_var: str, indent: str = None,
                       undo: bool = False) -> List[str]:
        '''
        Render preprocessor call to a list of strings, each holding one
        line of OpenCL code.

        Parameters
        ----------
        worker: ClWorker
            OpeCL worker instance.
        selector: slice
            Input data slice that selects the input items to preprocess.
        inout_var: str
            Name of the RW variable that holds the data.
            Indentation string. Defaults to 4 spaces.
        undo: bool
            If True render the undo call, else the apply call.

        Returns
        -------
        call: List[str]
            The rendered preprocessor call as a string.
        '''
        T = worker.types
        indent = '    ' if indent is None else indent
        return [
            'pp_expscale_{:s}('.format('undo' if undo else 'apply'),
            '{indent:s}{start:s}, {stop:s}, {step:s},'.format(
                start=T.size_t_str(selector.start),
                stop=T.size_t_str(selector.stop),
                step=T.size_t_str(selector.step), indent=indent),
            '{indent:s}{inout_var:s}'.format(
                inout_var=inout_var, indent=indent),
            ');'
        ]

    def apply(self, x: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Apply exponential function to the input data.

        Parameters
        ----------
        x: np.ndarray
            Input data array.
        out: np.ndarray
            Data array for the results. Must be of the same
            size and type as the x data array.
        
        Returns
        -------
        y: np.ndarray
            Exponentially scaled input data.
        '''
        result = np.exp(x, out=out)
        return result

    def undo(self, x: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Undo exponential scale from the output data.

        Parameters
        ----------
        y: np.ndarray
            Preprocessed output data.
        out: np.ndarray
            Data array for the results. Must be of the same
            size and type as the y data array.

        Returns
        -------
        x: np.ndarray
            Input data obtained by applying inverse preprocessing to the
            output data.
        '''
        result = np.log(x, out=out)
        return result

    def __str__(self):
        return 'ExpScale()'


class Normalize(PreprocessorFunction):
    @staticmethod
    def cl_implementation(worker:  ClWorker) -> str:
        '''
        Implementation of the normalization preprocessor in OpenCL.
        '''
        return '\n'.join((
            'inline void pp_normalize_apply(',
            '		mc_size_t start, mc_size_t stop, mc_size_t step, ',
            '		mc_fp_t scale, mc_fp_t offset, ',
            '		mc_fp_t *inout)',
            '{',
            '	for(mc_size_t pos = start; pos < stop; pos += step){',
            '		inout[pos] = inout[pos]*scale + offset;',
            '	}',
            '};',
            '',
            'inline void pp_normalize_undo(',
            '		mc_size_t start, mc_size_t stop, mc_size_t step, ',
            '		mc_fp_t inv_scale, mc_fp_t offset, ',
            '		mc_fp_t *inout)',
            '{',
            '	for(mc_size_t pos = start; pos < stop; pos += step){',
            '		inout[pos] = (inout[pos] - offset)*inv_scale;',
            '	}',
            '};',
        ))

    def cl_render_call(self, worker: ClWorker, selector: slice,
                       inout_var: str, indent: str = None,
                       undo: bool = False) -> List[str]:
        '''
        Render preprocessor call to a list of strings, each holding one line
        of OpenCL code.

        Parameters
        ----------
        worker: ClWorker
            OpeCL worker instance.
        selector: slice
            Input data slice that selects the input items to preprocess.
        inout_var: str
            Name of the RW variable that holds the data.
        indent: str
            Indentation string. Defaults to 4 spaces.
        undo: bool
            If True render the undo call, else the apply call.

        Returns
        -------
        call: List[str]
            The rendered preprocessor call as a string.
        '''
        T = worker.types
        indent = '    ' if indent is None else indent
        return [
            'pp_normalize_{:s}('.format('undo' if undo else 'apply'),
            '{indent:s}{start:s}, {stop:s}, {step:s},'.format(
                start=T.size_t_str(selector.start),
                stop=T.size_t_str(selector.stop),
                step=T.size_t_str(selector.step), indent=indent),
            '{indent:s}{scale:s}, {offset:s},'.format(
                scale=T.fp_str(1.0/self._scale if undo else self._scale),
                offset=T.fp_str(self._offset),
                indent=indent),
            '{indent:s}{inout_var:s}'.format(
                inout_var=inout_var, indent=indent),
            ');'
        ]

    def __init__(self, output_range: Tuple[float, float] = (0.0, 1.0),
                       input_range: Tuple[float, float] = None):
        '''
        Preprocessor that normalizes the values to the specified range.

        Parameters
        ----------
        output_range: Tuple(float, float)
            Range of the output values.
        input_range: Tuple(float, float)
            Range of the input values. If None, the range must be determined by
            calling the
            :py:meth:`~xopto.util.regression.sequential.Normalize.train`
            method.
        '''
        self._scale = 1.0
        self._offset = 0.0

        self._output_range = (float(output_range[0]), float(output_range[1]))

        self._input_range = None
        if input_range is not None:
            self._input_range = (float(input_range[0]), float(input_range[1]))
            self._prepare()

    def _get_input_range(self) -> Tuple[float, float]:
        return self._input_range
    input_range = property(_get_input_range, None, None,
                           'Range of value in the input data.')

    def _get_output_range(self) -> Tuple[float, float]:
        return self._output_range
    output_range = property(_get_output_range, None, None,
                           'Range of values in the output data.')

    def train(self, x: np.ndarray):
        '''
        Determine the range of input values.

        Parameters
        ----------
        x: np.ndarray
            A representative training dataset.
        '''
        x = np.asarray(x)
        self._input_range = (float(x.min()), float(x.max()))
        self._prepare()

    def apply(self, x: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Apply normalization to the input data.

        Parameters
        ----------
        x: np.ndarray
            Array of input data to be normalized.
        out: np.ndarray
            Data array for the results. Must be of the same
            size and type as the x data array.

        Returns
        -------
        y: np.ndarray
            Input data normalized to the output range.
        '''
        if self._scale is None or self._offset is None:
            raise ValueError('The range of input values is undefined!')

        if out is None:
            out = np.array(x)
        else:
            np.copyto(out, x)
        out *= self._scale
        out += self._offset

        return out

    def undo(self, y: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Apply denormalization to the output data.

        Parameters
        ----------
        y: np.ndarray
            Array of normalized output data to be denormalized.

        Returns
        -------
        x: np.ndarray
            Denormalized output data.
        '''
        if out is None:
            out = np.array(y)
        else:
            np.copyto(out, y)
        out -= self._offset
        out /= self._scale

        return out

    def _prepare(self):
        '''
        Internal method that precalculates the scaling factor and offset.
        '''
        delta_input = self._input_range[1] - self._input_range[0]
        delta_output = self._output_range[1] - self._output_range[0]
        if delta_input == 0.0:
            self._scale = 1.0
        else:
            self._scale = delta_output/delta_input

        self._offset = -self._input_range[0]*self._scale + self._output_range[0]

    def todict(self) -> dict:
        '''
        Export instance to a dict object.

        Returns
        -------
        data: dict
            Instance data exported to a dict.
        '''
        return {
            'type': self.__class__.__name__,
            'input_range': self._input_range,
            'output_range': self._output_range
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'Normalize':
        '''
        Create a new instance from data exported to dict.

        Parameters
        ----------
        data: dict
            Data that were exported to a dict by the
            :py:meth:`~xopto.util.regression.sequential.Normalize.todict`
            method.
        '''
        return PreprocessorFunction.fromdict(data)

    def __str__(self):
        return 'Normalize(output_range={}, input_range={})'.format(
            self.output_range, self.input_range)


class SNV(PreprocessorFunction):
    @staticmethod
    def cl_implementation(worker:  ClWorker) -> str:
        '''
        Implementation of the normalization preprocessor in OpenCL.
        '''
        return '\n'.join((
            'inline void pp_snv_apply(',
            '		mc_size_t start, mc_size_t stop, mc_size_t step, ',
            '		mc_fp_t scale, mc_fp_t offset, ',
            '		mc_fp_t *inout)',
            '{',
            '	for(mc_size_t pos = start; pos < stop; pos += step){',
            '		inout[pos] = inout[pos]*scale + offset;',
            '	}',
            '};',
            '',
            'inline void pp_snv_undo(',
            '		mc_size_t start, mc_size_t stop, mc_size_t step, ',
            '		mc_fp_t inv_scale, mc_fp_t offset, ',
            '		mc_fp_t *inout)',
            '{',
            '	for(mc_size_t pos = start; pos < stop; pos += step){',
            '		inout[pos] = (inout[pos] - offset)*inv_scale;',
            '	}',
            '};',
        ))

    def cl_render_call(self, worker: ClWorker, selector: slice,
                       inout_var: str, indent: str = None,
                       undo: bool = False) -> List[str]:
        '''
        Render preprocessor call to a list of strings, each holding one line
        of OpenCL code.

        Parameters
        ----------
        worker: ClWorker
            OpeCL worker instance.
        selector: slice
            Input data slice that selects the input items to preprocess.
        inout_var: str
            Name of the RW variable that holds the data.
        indent: str
            Indentation string. Defaults to 4 spaces.
        undo: bool
            If True render the undo call, else the apply call.

        Returns
        -------
        call: List[str]
            The rendered preprocessor call as a string.
        '''
        T = worker.types
        indent = '    ' if indent is None else indent
        return [
            'pp_snv_{:s}('.format('undo' if undo else 'apply'),
            '{indent:s}{start:s}, {stop:s}, {step:s},'.format(
                start=T.size_t_str(selector.start),
                stop=T.size_t_str(selector.stop),
                step=T.size_t_str(selector.step), indent=indent),
            '{indent:s}{scale:s}, {offset:s},'.format(
                T.fp_str(1.0/self._scale if undo else self._scale),
                T.fp_str(self._offset), indent=indent),
            '{indent:s}{inout_var:s}'.format(
                inout_var=inout_var, indent=indent),
            ');'
        ]

    def __init__(self, output_mean: float = 0.0, output_std: float = 1.0,
                 input_mean: float = None, input_std: float = None):
        '''
        Preprocessor that normalizes data
        to the specified mean and standard deviation (std).

        Parameters
        -----------
        output_mean: float
            The mean value of variables after preprocessing.
        output_std: float
            The standard deviation of variables after preprocessing.
        input_mean: float
            The expected mean value of the input variables.
        input_std: float
            The expected standard deviation of the input variables.
        '''
        self._output_mean = float(output_mean)
        self._output_std = float(output_std)

        self._input_mean = self._input_std = None

        self._scale = 1.0
        self._offset = 0.0

        if input_mean is not None:
            self._input_mean = float(input_mean)

        if input_std is not None:
            self._input_std = float(input_std)

        if self._input_mean is not None and self._input_std is not None:
            self._prepare()

    def _get_input_mean(self) -> float:
        return self._input_mean
    input_mean = property(_get_input_mean, None, None,
                          'Mean of the input values.')

    def _get_input_std(self) -> float:
        return self._input_std
    input_std = property(_get_input_std, None, None,
                         'Standard deviation of the input values.')

    def _get_output_mean(self) -> float:
        return self._output_mean
    output_mean = property(_get_output_mean, None, None,
                           'Mean of the output values.')

    def _get_output_std(self) -> float:
        return self._output_std
    output_std = property(_get_output_std, None, None,
                          'Standard deviation of the output values.')

    def train(self, x: np.ndarray):
        '''
        Determine the mean and std of the input values.

        Parameters
        ----------
        x: np.ndarray
            A representative training dataset.
        '''
        x = np.asarray(x)
        self._input_mean = x.mean()
        self._input_std = x.std()
        self._prepare()

    def apply(self, x: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Apply normalization to the input data.

        Parameters
        ----------
        x: np.ndarray
            Array of input data to be normalized.
        out: np.ndarray
            Data array for the results. Must be of the same
            size and type as the x data array.

        Returns
        -------
        y: np.ndarray
            Input data normalized to the output range.
        '''
        if self._input_mean is None or self._input_std is None:
            raise ValueError(
                'The mean and/or std of input values are undefined!')

        if out is None:
            out = np.array(x)
        else:
            np.copyto(out, x)
        out *= self._scale
        out += self._offset

        return out

    def undo(self, y: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Apply denormalization to the output data.

        Parameters
        ----------
        y: np.ndarray
            Array of normalized output data to be denormalized.
        out: np.ndarray
            Data array for the results. Must be of the same
            size and type as the y data array.

        Returns
        -------
        x: np.ndarray
            Denormalized output data.
        '''
        if out is None:
            out = np.array(y)
        else:
            np.copyto(out, y)
        out -+ self._offset
        out /= self._scale

        return out

    def _prepare(self):
        '''
        Internal method that precalculates the scaling factor and offset.
        '''
        if self._input_std == 0.0:
            self._scale = 1.0
        else:
            self._scale = self._output_std/self._input_std
        self._offset = -self._input_mean*self._scale + self._output_mean

    def todict(self) -> dict:
        '''
        Export instance to a dict object.

        Returns
        -------
        data: dict
            Instance data exported to a dict.
        '''
        return {
            'type': self.__class__.__name__,
            'input_mean': self._input_mean,
            'input_std': self._input_std,
            'output_mean': self._output_mean,
            'output_std': self._output_std,
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'Normalize':
        '''
        Create a new instance from data exported to dict.

        Parameters
        ----------
        data: dict
            Data that were exported to a dict by the
            :py:meth:`~xopto.util.regression.sequential.SNV.todict`
            method.
        '''
        return PreprocessorFunction.fromdict(data)

    def __str__(self):
        return 'SNV(output_mean={}, output_std={}, ' \
                    'input_mean={}, input_stf={})'.format(
                        self.output_mean, self.output_std,
                        self.input_mean, self.input_std)


class PreprocessorItem:
    def cl_render_call(self, worker: ClWorker,
                       inout_var: str, n: int, undo: bool = False,
                       indent: str = None) -> List[str]:
        '''
        Render preprocessor item call to a list of strings, each representing
        one line of OpenCL code.

        Parameters
        ----------
        worker: ClWorker
            OpeCL worker instance.
        inout_var: str
            Name of the RW variable that holds the data.
        n: int
            The length of the input/output data vector.
        indent: str
            Indentation string. Defaults to 4 spaces.
        undo: bool
            If True render the undo call, else the apply call.

        Returns
        -------
        call: List[str]
            The rendered preprocessor call as a list of code lines.
        '''
        sequence = self.sequence
        if undo:
            sequence = sequence[::-1]
        selector = self.selector[-1]
        start, stop, step = selector.start, selector.stop, selector.step
        start = 0 if start is None else start
        stop = n if stop is None else stop
        step = 1 if step is None else step
        selector = slice(start, stop, step)

        lines = []
        for func_index, func in enumerate(sequence):
            func_lines = func.cl_render_call(
                worker, selector, inout_var,
                undo=undo, indent=indent)
            lines.extend(func_lines)
        
        return lines

    def __init__(self, sequence: PreprocessorFunction or \
                                 Tuple[PreprocessorFunction],
                 selector: slice or Tuple[slice] = None):
        '''
        Preprocessor that applies a sequence of preprocessor functions to
        the selected data.

        Parameters
        ----------
        sequence: PreprocessorFunction or Tuple[PreprocessorFunction]
            A sequence of preprocessor functions to apply to the data.
        selector: slice or Tuple[slice]
            A slice selector for the input data. The preprocessor is only
            applied to the selected data.
        '''
        if isinstance(sequence, PreprocessorFunction):
            sequence = (sequence, )
 
        if selector is None:
            selector = slice(None)

        if isinstance(selector, slice):
            selector = (selector, )

        self._sequence = sequence
        self._selector = selector

    def _get_sequence(self) -> Tuple[PreprocessorFunction]:
        return self._sequence
    sequence = property(_get_sequence, None, None,
                        'Sequence of preprocessor functions.')

    def _get_selector(self) -> Tuple[slice]:
        return self._selector
    selector = property(_get_selector, None, None, 'Data selector.')

    def apply(self, x: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Apply the preprocessor to the input data.

        Parameters
        ----------
        x: np.ndarray
            Input data.
        out: np.ndarray
            Data array for the results. Must be of the same
            size and type as the x data array.

        Returns
        -------
        y: np.ndarray
            Selected input data that were preprocessed by the sequence
            of preprocessor functions.
        '''
        selected = x[..., self._selector]
        for func in self._sequence:
            selected = func.apply(selected, out=out)

        return selected

    def undo(self, y: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Apply inverse of the preprocessor to the selected output data.

        Parameters
        ----------
        y: np.ndarray
            Output data.
        out: np.ndarray
            Data array for the results. Must be of the same
            size and type as the y data array.

        Returns
        -------
        x: np.ndarray
            Selected output data that were preprocessed by the reversed
            sequence of the inverse preprocessor functions.
        '''
        selected = y[..., self._selector]
        for func in self._sequence[::-1]:
            selected = func.undo(selected, out= out)

        return selected

    def train(self, x: np.ndarray):
        '''
        Train the sequence of the preprocessor functions.

        Parameters
        ----------
        x: np.ndarray
            A representative dataset for the training.
        '''
        selected = x[..., self._selector]
        for func in self._sequence:
            func.train(selected)
            selected = func.apply(selected)

    def todict(self) -> dict:
        '''
        Export preprocessor to a dict object.

        Returns
        -------
        preprocessor: dict
            Preprocessor exported to a dict.
        '''
        s = self._selector
        return {
            'type': self.__class__.__name__,
            'sequence': [item.todict() for item in self._sequence],
            'selector': [{'type': 'slice', 'start': s.start, 'stop': s.stop,
                          'step':s.step} for s in self._selector],
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'PreprocessorItem':
        '''
        Lod a preprocessor item from a dict.

        Parameters
        ----------
        data: dict
            :py:class:`PreprocessorItem` instance exported to a dict.

        Returns
        -------
        item: PreprocessorItem
            An instance of :py:class:`PreprocessorItem` created from the
            data in the dict.
        '''
        data = dict(data)
        sequence_list = data.pop('sequence')
        sequence = [PreprocessorFunction.fromdict(item)
                    for item in sequence_list]
        selector_list = data.pop('selector')
        selector = []
        for item in selector_list:
            item_data = dict(item)
            type_name = item_data.pop('type')
            T = getattr(builtins, type_name, globals().get(type))
            if T is slice:
                slice(item_data.get('start'), item_data.get('stop'),
                      item_data.get('step'))
            else:
                selector.append(T(**item_data))
        data['sequence'] = sequence
        data['selector'] = selector
        type_name = data.pop('type')
        T = globals().get(type_name, cls)
        return T(**data)

    def __str__(self):
        return 'PreprocessorItem(sequence={}, selector={})'.format(
            self.sequence, self.selector)

    def __repr__(self):
        return '{} # {}'.format(self.__str__(), id(self))


class Preprocessor:
    def cl_render_call(self, worker: ClWorker,
                       prefix: str, inout_var: str,
                       n: int, undo: bool = False,
                       indent: str = None, width: int = 80) -> List[str]:
        '''
        Render preprocessor item call to a list of strings, each representing
        one line of OpenCL source code.

        Parameters
        ----------
        worker: ClWorker
            OpeCL worker instance.
        prefix: str
            Preprocessor function prefix.
        inout_var: str
            Name of the RW variable that holds the data.
        n: int
            The length of the input/output data vector.
        undo: bool
            If True render the undo call, else the apply call.
        indent: str
            Indentation string. Defaults to 4 spaces.
        width: str
           Maximum line width.

        Returns
        -------
        call: List[str]
            The rendered preprocessor call as a list of code lines.
        '''
        lines = []
        lines.extend([
            'void {prefix}pp_{apply_or_undo}('.format(
                prefix=prefix,
                apply_or_undo = 'undo' if undo else 'apply'),
            '{indent}mc_fp_t *{inout_var})'.format(
                indent=indent, inout_var=inout_var),
            '{'
        ])

        body = []
        for item_index, item in enumerate(self.items):
            body.append(
                '/* preprocessor item {} of {} */'.format(
                    item_index + 1, len(self.items)
                )
            )
            body.extend(
                item.cl_render_call(
                    worker, inout_var, n, undo=undo, indent=indent)
            )
        for line in body:
            lines.append(indent + line)

        lines.append('};')

        return lines

    def __init__(self, items: PreprocessorItem or
                              Tuple[PreprocessorItem] = None):
        '''
        Preprocessor that is applied to the input data.

        Parameters
        ----------
        items: PreprocessorItem or Tuple[PreprocessorItem]
            Preprocessor items.
        '''
        if items is None:
            items = tuple()
        if isinstance(items, PreprocessorItem):
            items = (items, )

        self._items = tuple(items)

    def _get_items(self) -> Tuple[PreprocessorItem]:
        return self._items
    items = property(_get_items, None, None, 'Preprocessor items.')

    def apply(self, x: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Apply preprocessor to the input data.

        Parameters
        ----------
        x: np.ndarray
            The input data.
        out: np.ndarray
            Data array for the results. Must be of the same
            size and type as the x data array.

        Returns
        -------
        y: np.ndarray
            Preprocessed input data.
        '''
        if out is None:
            out = np.array(x)
        else:
            np.copyto(out, x)

        for item in self._items:
            out[..., item.selector] = item.apply(x)
        
        return out

    def undo(self, y: np.ndarray, out: np.ndarray = None) -> np.ndarray:
        '''
        Undo preprocessing from the output data.

        Parameters
        ----------
        y: np.ndarray
            Preprocessed output data.
        out: np.ndarray
            Data array for the results. Must be of the same
            size and type as the y data array.

        Returns
        -------
        x: np.ndarray
            Input data obtained by applying inverse preprocessing to the
            output data.
        '''
        if out is None:
            out = np.array(y)
        else:
            np.copyto(out, y)

        for item in self._items:
            out[..., item.selector] = item.undo(y)

        return out

    def train(self, x: np.ndarray):
        '''
        Train the preprocessor items.

        Parameters
        ----------
        x: np.ndarray
            A representative input training dataset.
        '''
        for item in self._items:
            item.train(x)

    def todict(self) -> dict:
        '''
        Export preprocessor to a dict object.

        Returns
        -------
        data: dict
            Preprocessor exported to a dict object.
        '''
        return {
            'type': self.__class__.__name__,
            'items': [item.todict()  for item in self._items]
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'Preprocessor':
        '''
        Lod a preprocessor item from a dict.

        Parameters
        ----------
        data: dict
            :py:class:`Preprocessor` instance exported to a dict.

        Returns
        -------
        item: PreprocessorItem
            An instance of :py:class:`Preprocessor` created from the
            data in the dict.
        '''
        data = dict(data)
        type_name = data.pop('type')
        itmes_list = data.pop('items')
        items = [PreprocessorItem.fromdict(item) for item in itmes_list]
        data['items'] = items
        T = globals().get(type_name, cls)
        return T(**data)

    def __str__(self):
        return 'Preprocessor(items={})'.format(self.items)

    def __repr__(self):
        return '{} # {}'.format(self.__str__(), id(self))


if __name__ == '__main__':
    data = np.array([[1.0, 2, 3, 4, 5],[10, 10**2, 10**3, 10**4, 10**5]])

    pp = Preprocessor(
        (
            PreprocessorItem((Log10Scale(), Normalize()), slice(0, 1)),
            PreprocessorItem((Log10Scale(), ), slice(1, 2))
        )
    )

    d = pp.todict()
    Preprocessor.fromdict(d)

    pp.train(data)
    pp_data = pp.apply(data)
    data_ = pp.undo(pp_data)

    def fp_str(value: float) -> str:
        return '{:.8e}'.format(float(value))

    weights = np.random.rand(20, 40)
    bias = np.random.rand(weights.shape[0])

    def generate_transformation(weights, bias, indent=0, width=80):
        lines = []
        
        for line_index, line in enumerate(weights):
            items = []
            for weight_index, weight in enumerate(line):
                items.append('{:s}*input[{:d}]'.format(fp_str(weight), weight_index))
            if bias is not None:
                items.append(fp_str(bias[line_index]))
            header = 'output[{:d}] = '.format(line_index)
            offset = len(header)
            active_line = indent + header
            joiner = ' + '
            for item_index, item in enumerate(items):
                last = item_index >= len(items) - 1
                if item_index > 0:
                    active_line += joiner
                if len(active_line) + len(item) >= width - 1:
                    lines.append(active_line)
                    active_line = indent + ' '*offset

                active_line += item
            lines.append(active_line + ';')
        print('inline void mua_sequential_layer_0(fp_t const *input, fp_t *output){')
        print('\n'.join(lines))
        print('}')

    generate_transformation(weights, bias, indent=' '*4, width=80)
