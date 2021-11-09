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

from xopto.mcbase import mcobject


class RunMinWeightBase:
    def __init__(self, min_weight: float, batch_size: int,
                 selection: slice or int = None):
        '''
        Base class for running MC simulations until the desired weight
        is collected by the selected detector.

        Parameters
        ----------
        min_weight: float
            Minimum weight that needs to be collected by the detector.
        batch_size: int
            Number of packets to launch in a single simulation.
        '''
        if selection is None:
            selection = slice(None)
        self._min_weight = float(min_weight)
        self._batch_size = int(batch_size)
        self._num_packets = 0
        self._selection = selection

    def run(self, mc_obj: mcobject.McObject, min_packets: int = None,
            out: tuple = None, *args, **kwargs) -> tuple:
        '''
        Repeat simulations in batches until the required minimum weight
        is collected by the selected detector and optionally the minimum number
        of packets is launched.

        Parameters
        ----------
        mc_obj: mcobject.McObject
            Simulator instance.
        min_packets: int
            Optional minimum number of packets to launch and simulate. 
            If None (default), there is no restriction on te minimum number of
            packets that need to be simulated.
        out: tuple
            Output from the previous run of the simulator.
        *args: tuple
            Positional arguments passed to the simulator.
        **kwargs: dict
            Keyword arguments passed to the simulator.

        Returns
        -------
        result: tuple
            Accumulated results as returned by the Mont Carlo simulator
            instance.
        '''
        if min_packets is None:
            min_packets = 0

        total_num_packets = 0
        total_weight = 0.0
        while out is None or self.min_weight > total_weight or \
                total_num_packets < min_packets:
            out = mc_obj.run(self.batch_size, out=out, *args, **kwargs)
            total_weight, num_packets, out = self.process_batch(out)
            total_num_packets += self.batch_size

        self._num_packets = total_num_packets
        self._total_weight = total_weight

        return out

    def process_batch(self, result: tuple) -> Tuple[float, int, tuple]:
        '''
        Implement this method in a subclass.
        Parameters
        ----------
        result: tuple
            Result as returned by the Monte Carlo simulation.

        Returns
        -------
        total_weight: float
            The total accumulated weight.
        num_packets: int
            The total number of packets accumulated after processing this batch. 
        result: tuple
            Simulation result ready for the next Monte Carlo simulation. 
        '''
        return 0.0, 0, result

    def _get_batch_selection(self) -> int:
        return self._selection
    selection = property(_get_batch_selection, None, None,
                         'A slice / region or index of the detector that is '
                         'used to compute the collected weight. If None '
                         '(default), the entire detector contributes to the '
                         'total weight.')

    def _get_batch_size(self) -> int:
        return self._batch_size
    batch_size = property(_get_batch_size, None, None,
                          'Number of packets processed / launched in one '
                          'simulation run / batch.')

    def _get_min_weight(self) -> int:
        return self._min_weight
    min_weight = property(_get_min_weight, None, None,
                          'The minimum total weight that needs to be collected '
                          'by the selected detector.')

    def _get_num_packets(self) -> int:
        return self._num_packets
    n = property(_get_num_packets, None, None,
                 'The number of packets that were launched during '
                 'the last call to the :py:meth:`run` method.')

    def _get_total_weight(self) -> float:
        return self._total_weight
    weight = property(_get_total_weight, None, None,
                      'The weight collected by the detector during '
                      'the last call to the :py:meth:`run` method.')


class RunMinPacketsBase:
    def __init__(self, min_packets: float, batch_size: int,
                 selection: slice or int):
        '''
        Object for running MC simulations until the desired number of packets
        is collected by the selected detector.

        Parameters
        ----------
        min_packets: int
            Minimum required number of packets.
        min_weight: float
            Minimum weight that needs to be collected by the detector.
        batch_size: int
            Number of packets to launch in a single simulation.
        '''
        if selection is None:
            selection = slice(None)
        self._min_packets = int(min_packets)
        self._batch_size = int(batch_size)
        self._num_packets = 0
        self._selection = selection

    def run(self, mc_obj: mcobject.McObject, min_packets: int = None,
            out: tuple = None, *args, **kwargs) -> tuple:
        '''
        Repeat simulations in batches until the required minimum weight
        is collected by the selected detector and optionally the minimum number
        of packets is launched.

        Parameters
        ----------
        mc_obj: mcobject.McObject
            Simulator instance.
        min_packets: int
            Optional minimum number of packets to launch and simulate. 
            If None (default), there is no restriction on te minimum number of
            packets that need to be simulated.
        out: tuple
            Output from the previous run of the simulator.
        *args: tuple
            Positional arguments passed to the simulator.
        **kwargs: dict
            Keyword arguments passed to the simulator.

        Returns
        -------
        result: tuple
            Accumulated results as returned by the Mont Carlo simulator
            instance.
        '''
        if min_packets is None:
            min_packets = 0

        total_num_packets = 0
        total_weight = 0.0
        while out is None or self.min_weight > total_weight or \
                total_num_packets < min_packets:
            out = mc_obj.run(self.batch_size, out=out, *args, **kwargs)
            total_weight, num_packets, out = self.process_batch(out)
            total_num_packets += self.batch_size

        self._num_packets = total_num_packets
        self._total_weight = total_weight

        return out

    def process_batch(self, result: tuple) -> Tuple[float, int, tuple]:
        '''
        Implement this method in a subclass. Return True if the batch is
        complete.

        Parameters
        ----------
        result: tuple
            Result as returned by the Monte Carlo simulation.

        Returns
        -------
        total_weight: float
            The total accumulated weight.
        num_packets: int
            The total number of packets accumulated after processing this batch. 
        result: tuple
            Simulation result ready for the next Monte Carlo simulation. 
        '''
        return 0.0, 0, result

    def _get_batch_selection(self) -> int:
        return self._batch_size
    selection = property(_get_batch_selection, None, None,
                         'A slice / region or index of the detector that is '
                         'used to compute the collected weight. If None '
                         '(default), the entire detector contributes to the '
                         'total weight.')

    def _get_batch_size(self) -> int:
        return self._batch_size
    batch_size = property(_get_batch_size, None, None,
                          'Number of packets processed / launched in one '
                          'simulation run / batch.')

    def _get_min_packets(self) -> int:
        return self._min_weight
    min_packets = property(_get_min_packets, None, None,
                          'The minimum number of packets that needs to be '
                          'collected by the selected detector.')

    def _get_num_packets(self) -> int:
        return self._num_packets
    n = property(_get_num_packets, None, None,
                 'The number of packets that were launched during '
                 'the last call to the :py:meth:`run` method.')

    def _get_total_weight(self) -> float:
        return self._total_weight
    weight = property(_get_total_weight, None, None,
                      'The weight collected by the detector during '
                      'the last call to the :py:meth:`run` method.')
