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

from xopto.mcbase.mcrun import RunMinWeightBase, RunMinPacketsBase
from xopto.mcvox import mc



class _RunMinWeightDetector(RunMinWeightBase):
    def __init__(self, location: str, min_weight: float, batch_size: int,
                 selection: slice or int):
        '''
        Base class for running MC simulations until the desired weight
        is collected by the selected detector.

        Parameters
        ----------
        location: str
            Detector location ("top", "bottom" or "specular")
        min_weight: float
            Minimum weight that needs to be collected by the detector.
        batch_size: int
            Number of packets to launch in a single simulation.
        selection: slice or int
            A slice / region or index of the detector that is used to compute
            the collected weight. If None (default), the entire detector
            contributes to the total weight.
        '''
        self._location = location
        super().__init__(min_weight, batch_size, selection)

    def process_batch(self, result: mc.MC_RESULT_TYPE) \
            -> Tuple[float, int, mc.MC_RESULT_TYPE]:
        '''
        Process one MC simulation.

        Note
        ----
        Do not call this method directly. The method is called by
        :py:meth:`run`.

        Parameters
        ----------
        result: mc.MC_RESULT_TYPE
            A tuple as returned by the :py:meth:`mcvox.mc.Mc.run` method.
    
        Returns
        -------
        total_weight: float
            Total weight collected by the selected detector.
        n: int
            The number of packets collected by the detector or None.
        result: mc.MC_RESULT_TYPE
            Input simulation result ready for the next Monte Carlo simulation.
        '''
        detectors = result[2]
        if detectors is None or type(getattr(detectors, self._location)) == \
                mc.mcdetector.DetectorDefault:
            raise RuntimeError(
                'This MC simulator instance does not use the configured'
                ' ("{}") detector!'.format(self._location))

        detector = getattr(detectors, self._location)

        raw_data = detector.raw
        if self.selection is not None:
            raw_data = raw_data[self.selection]
        total_weight = np.sum(raw_data)

        return total_weight, None, result

    def _get_location(self) -> str:
        return self._location
    location = property(_get_location, None, None,
                        'Detector location.')


class RunMinWeightTop(_RunMinWeightDetector):
    def __init__(self, min_weight: float, batch_size: int,
                 selection: slice or int = None):
        '''
        Object for running MC simulations until the desired weight is collected
        by the top detector.

        Parameters
        ----------
        min_weight: float
            Minimum weight that needs to be collected by the detector.
        batch_size: int
            Number of packets to launch in a single simulation.
        selection: slice or int
            A slice / region or index of the detector that is used to compute
            the collected weight. If None (default), the entire detector
            contributes to the total weight.
        '''
        super().__init__('top', min_weight, batch_size, selection)

class RunMinWeightBottom(_RunMinWeightDetector):
    def __init__(self, min_weight: float, batch_size: int,
                 selection: slice or int = None):
        '''
        Object for running MC simulations until the desired weight is collected
        by the bottom detector.

        Parameters
        ----------
        min_weight: float
            Minimum weight that needs to be collected by the detector.
        batch_size: int
            Number of packets to launch in a single simulation.
        selection: slice or int
            A slice / region or index of the detector that is used to compute
            the collected weight. If None (default), the entire detector
            contributes to the total weight.
        '''
        super().__init__('bottom', min_weight, batch_size, selection)

class RunMinWeightSpecular(_RunMinWeightDetector):
    def __init__(self, min_weight: float, batch_size: int,
                 selection: slice or int = None):
        '''
        Object for running MC simulations until the desired weight is collected
        by the specular detector.

        Parameters
        ----------
        min_weight: float
            Minimum weight that needs to be collected by the detector.
        batch_size: int
            Number of packets to launch in a single simulation.
        selection: slice or int
            A slice / region or index of the detector that is used to compute
            the collected weight. If None (default), the entire detector
            contributes to the total weight.
        '''
        super().__init__('specular', min_weight, batch_size, selection)


class RunMinWeightTrace(RunMinWeightBase):
    def __init__(self, min_weight: float, batch_size: int):
        '''
        Object for running MC simulations until the desired weight is collected
        by the trace, i.e. the sum of terminal weights of the traced packets.

        Parameters
        ----------
        min_weight: float
            Minimum total terminal weight of the packets that needs to be
            collected by the trace.
        batch_size: int
            Number of packets to launch in a single simulation.
        '''
        super().__init__(min_weight, batch_size, None)

    def process_batch(self, result: mc.MC_RESULT_TYPE) \
            -> Tuple[float, int, mc.MC_RESULT_TYPE]:
        '''
        Process one MC simulation.

        Note
        ----
        Do not call this method directly. The method is called by
        :py:meth:`run`.

        Parameters
        ----------
        result: mc.MC_RESULT_TYPE
            A tuple as returned by the :py:meth:`mcvox.mc.Mc.run` method.
    
        Returns
        -------
        total_weight: float
            Total weight collected by the selected detector.
        n: int
            The number of packets collected by the detector or None.
        result: mc.MC_RESULT_TYPE
            Input simulation result ready for the next Monte Carlo simulation.
        '''
        trace = result[0]
        if trace is None:
            raise RuntimeError('This MC simulator instance does not use trace!')

        total_weight = np.sum(trace.terminal['w'])
        num_packets = len(trace)

        return total_weight, num_packets, result


class RunMinPacketsTrace(RunMinPacketsBase):
    def __init__(self, min_packets: float, batch_size: int):
        '''
        Object for running MC simulations until the desired number of packets
        is collected in the (filtered) trace.

        Parameters
        ----------
        min_packets: float
            Minimum number of packet traces that needs to be collected.
        batch_size: int
            Number of packets to launch in a single simulation.
        '''
        super().__init__(min_packets, batch_size, None)

    def process_batch(self, result: mc.MC_RESULT_TYPE) \
            -> Tuple[float, int, mc.MC_RESULT_TYPE]:
        '''
        Process one MC simulation.

        Note
        ----
        Do not call this method directly. The method is called by
        :py:meth:`run`.

        Parameters
        ----------
        result: mc.MC_RESULT_TYPE
            A tuple as returned by the :py:meth:`mcvox.mc.Mc.run` method.
    
        Returns
        -------
        total_weight: float
            Total weight collected by the selected detector.
        n: int
            The number of packets collected by the detector or None.
        result: mc.MC_RESULT_TYPE
            Input simulation result ready for the next Monte Carlo simulation.
        '''
        trace = result[0]
        if trace is None:
            raise RuntimeError('This MC simulator instance does not use trace!')

        total_weight = np.sum(trace.terminal['w'])
        num_packets = len(trace)

        return total_weight, num_packets, result
