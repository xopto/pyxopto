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

import time
import threading

import numpy as np

from xopto.mcbase import mcobject
import pyopencl as cl

class ProgressMonitor:
    def __init__(self, mcsim: mcobject.McObject, interval: float = 0.5):
        '''
        Initialize Monte Carlo progress monitor.

        Parameters
        ----------
        mcsim: mcobject.McObject
            Monte Carlo simulator instance.
        target: int
            Optional target/final value for the monitored quantity, e.g. the
            number of launched packets.
        interval: float
            Polling interval in seconds.
        '''
        self._mcsim = mcsim
        self._interval = max(0.1, float(interval))
        self._processed = 0
        self._threads = 0
        self._target = 0
        self._terminate_on_stop = False

        self._track = False
        self._stop = False
        self._condition = threading.Semaphore(0)

        self._thread = threading.Thread(target=self._proc, args=(mcsim,))
        self._thread.start()

    def start(self, target: int, terminate: bool = True) -> 'ProgressMonitor':
        '''
        Starts monitoring the progress of the related Monte Carlo
        simulator instance.

        Parameters
        ----------
        target: int
            Optional target/final value for the monitored quantity, e.g. the
            number of launched packets.
        terminate: bool
            Flag that terminates the monitor if set to True.
            Setting this flag to False will allow to reuse the monitor.
            When using the monitor in a for loop, it is much more efficient to
            reuse the monitor, since each instance of a monitor is using
            a background thread that needs to be created and initialized
            before the monitoring can star.
            
        Returns
        -------
        self: ProgressMonitor
            Returns self.

        Notes
        -----
        Note that a terminated monitor cannot be restarted. A RuntimeError
        will be raised if an attept is made to start a terminated monitor.
        '''
        if self._stop:
            raise RuntimeError(
                'A terminated progress monitor can not be started!')

        self._target = int(target)
        self._processed = 0
        self._threads = 0
        self._track = True
        self._condition.release()
        self._terminate_on_stop = bool(terminate)
        return self

    def resume(self, target: int = None):
        '''
        Resume the monitor with the last state.

        Parameters
        ----------
        target: int
            Optional new target/final value for the monitored quantity, e.g.
            the number of launched packets, that will replace the existing
            target value.
        '''
        if target is not None:
            self._target = int(target)

        self._track = True

    def progress(self) -> float:
        '''
        Returns progress as a floating point number from 0.0 to 1.0, where
        1.0 is returned when the monitored task is complete.

        Returns
        -------
        progress: float
            Returns a value between 0.0 and 1.0 (complete).
        '''
        return min(self._processed/self._target, 1.0)

    def target(self) -> int:
        '''
        Returns the current target value of the monitored quantity.

        Returns
        -------
        target: int
            Target/final value for the monitored quantity, e.g the number of
            launched packets.
        '''
        return self._target

    def processed(self) -> int:
        '''
        Returns the number of processed items.
        
        Returns
        -------
        processed: int
            The number of processed items.
        '''
        return min(self._target, self._processed)

    def threads(self) -> int:
        '''
        Returns the number of OpenCL threads that are working
        on the job. Returns None, if this functionality is not
        supported.
        
        Returns
        -------
        n: int
            the number of OpenCL threads that are working on the job or
            None if not supported.
        '''
        return self._threads

    def stop(self):
        '''
        Stop monitoring the progres of the related Monte Carlo simulator
        instance. Note that the monitoring can be restarted.
        '''
        self._track = False
        self._processed = self._target
        if self._terminate_on_stop:
            self.terminate()

    def terminate(self):
        '''
        Terminate the progress monitor. Note that a terminated progress monitor
        cannot be restarted.
        '''
        self._stop = True
        self._track = False
        self._condition.release()

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._terminate_on_stop:
            self.terminate()

    def _proc(self, mcsim: mcobject.McObject):
        num_processed = np.zeros([1], dtype=mcsim.types.np_cnt)
        num_threads = np.zeros([1], dtype=np.uint32)
        queue = cl.CommandQueue(mcsim.cl_context)

        while not self._stop:
            #print('\nloop 1\n')
            self._condition.acquire()
            while self._track and self._target > self._processed:
                #print('\nloop 2\n')
                cl_num_packets = mcsim.cl_buffers.get('num_processed_packets')
                cl_num_kernels = mcsim.cl_buffers.get('num_kernels')
                if cl_num_packets is not None and cl_num_kernels is not None:
                    cl.enqueue_copy(queue, num_processed, cl_num_packets)
                    cl.enqueue_copy(queue, num_threads, cl_num_kernels)
                    self._threads = num_threads[0]
                    if self._processed != num_processed[0]:
                        self._processed = num_processed[0]
                        self.report()
                    if self._target <= self._processed:
                        print()
                        
                time.sleep(self._interval)

    def report(self):
        '''
        Overload this method for custom handling of the progress report.
        This function is called each time the value of the progress
        changes. Note that the polling interval is set with the
        constructor  
        :py:class:`~xopto.mcbase.mcprogress.ProgressMonitor` parameter
        :code:`interval`.
        '''
        N = 42
        n = int(self.progress()*N)
        print('|{}>{}| {:d}%'.format(
            '-'*n, ' '*(N - n), int(100.0*self.progress())), end='\r')
