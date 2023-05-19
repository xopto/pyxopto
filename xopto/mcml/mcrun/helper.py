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

from typing import Tuple, List, Dict, Callable
import os
import time
import argparse

from xopto.mcml import mc

import numpy as np

class McRunHelper:
    @staticmethod
    def cli_input(first: int = 0, n: int = 1000, batch: int = 1000,
                  packets: int = 100e6,
                  device: str = None, device_index: int = 0,
                  root_dir: str = None, mc_dir: str = None,
                  processed_dir: str = None,
                  verbose: bool = False, label: str = None,):
        '''
        Customize MC simulations with command line interface (CLI) arguments.
        Command line parameters are initialized with the default values
        that are passed to the function.

        Parameters
        ----------
        first: int
            Index of the first sample that will be simulated.
        n: int
            Number of samples that will be simulated from and including the
            sample at the specified offset.
        batch: int
            Number of samples in batch that will be saved to a single file.
        packets: int
            Number of packets to launch in each simulation.
        device: str
            Unique part of tht OpenCL device name that will run the simulations.
        device_index: int
            Zero-based index of the OpenCL device.
        root_dir: str
            Root directory of the dataset. Defaults to current working directory.
        mc_dir: str
            Output subdirectory for the MC simulations. Defaults to "mc".
        processed_dir: str
            Output subdirectory for the postprocessed data. Defaults to "processed".
        verbose: bool
            Turn on/off verbose mode.
        label: str
            Label displayed by the command line help.

        Returns
        -------
        cli_args: dict
            A dictionary of input arguments (keys) and their values. The
            values are initialized with the input arguments and updated if
            a corresponding CLI parameter is passed in the call.
        '''
        if device is None:
            device = ''
        if label is None:
            label = 'Command line arguments for MC simulations'
        if mc_dir is None:
            mc_dir = 'mc'
        if processed_dir is None:
            processed_dir = 'processed'
        if root_dir is None:
            root_dir = os.getcwd()

        parser = argparse.ArgumentParser(description=label)
        parser.add_argument('-d', '--device', dest='device',
                            type=str, default=str(device),
                            help='Unique part of the OpenCL device/platform name.')
        parser.add_argument('-i', '--index', dest='device_index',
                            type=int, default=int(device_index),
                            help='OpenCL device index.')

        parser.add_argument('-f', '--first', dest='first',
                            type=int, default=int(first),
                            help='Zero-based index of the sample that '
                                'will be simulated first.')
        parser.add_argument('-n', '--number', dest='number',
                            type=int, default=int(n),
                            help='Number of samples that will be simulated '
                                'simulated.')
        parser.add_argument('-b', '--batch', dest='batch',
                            type=int, default=int(batch),
                            help='Number of samples in a batch that will be saved'
                                'into a single file.')

        parser.add_argument('-p', '--packets', dest='packets',
                            type=int, default=int(packets),
                            help='Number of packets to launch in each simulation.')

        parser.add_argument('--root-dir', dest='root_dir',
                            type=str, default=str(root_dir),
                            help='Root directory of the dataset.')

        parser.add_argument('--output-dir', dest='mc_dir',
                            type=str, default=str(mc_dir),
                            help='Output subdirectory for the MC simulations.')

        parser.add_argument('--processed-dir', dest='processed_dir',
                            type=str, default=str(processed_dir),
                            help='Output subdirectory for the postprocessed data.')

        parser.add_argument('-v', '--verbose', dest='verbose',
                            action='store_true', default=bool(verbose),
                            help='Turn on/off verbose mode.')

        args = parser.parse_args()

        return {'first': args.first, 'n': args.number, 'batch': args.batch,
                'packets': args.packets, 'root_dir': args.root_dir,
                'mc_dir': args.mc_dir, 'processed_dir': args.processed_dir,
                'device': args.device if args.device else None,
                'device_index': args.device_index, 'verbose': args.verbose}

    def __init__(self, *args, **kwargs):
        '''
        Creates a new instance of simulator helper that can be used
        to customize and run simulations for various configurations.

        args: tuple
            Positional and keyword arguments passed to the MC
            constructor :py:meth:`xopto.mcml.mc.Mc`.
        kwargs: dict
            Keyword arguments passed to the MC
            constructor :py:meth:`xopto.mcml.mc.Mc`.

        Note
        ----
        The base class implementation calls the
        :py:meth:`McRunHelper.create_mc` method to create and initialize
        a new managed simulator instance.
        '''
        self._mc_obj = self.create_mc(*args, **kwargs)

    def create_layers(self) -> mc.mclayer.Layers:
        '''
        Creates and returns the layer stack.

        Returns
        -------
        layers: xopto.mcml.mclayer.layer.Layers
            The created layer stack.

        Note
        ----
        Overload this method to create a custom layer stack. The default
        implementation creates a basic semi-infinite slab.
        '''
        return mc.mclyer.Layers((
            mc.mclayer.Layer(d=float('inf'),
                             mua=0.0e2, mus=0.0e2,
                             n=1.0, pf=mc.mcpf.Hg(0.8)),
            mc.mclayer.Layer(d=float('inf'),
                             mua=0.0e2, mus=0.0e2,
                             n=1.0, pf=mc.mcpf.Hg(0.8)),
            mc.mclayer.Layer(d=float('inf'),
                             mua=0.0e2, mus=0.0e2,
                             n=1.0, pf=mc.mcpf.Hg(0.8))
        ))

    def create_source(self) -> mc.mcsource.Source:
        '''
        Creates and returns a photon packet source.

        Returns
        -------
        source: xopto.mcml.mcsource.base.Source
            The created packet source.

        Note
        ----
        Overload this method to create custom sources. The default
        implementation creates an instance of
        :py:class:`xopto.mcml.mcsource.line.Line` source.
        '''
        return mc.mcsource.Line()

    def create_surface(self) -> mc.mcsurface.SurfaceLayouts or None:
        '''
        Creates and returns a surface layout or None if surface layout is
        not used.

        Returns
        -------
        layouts: xopto.mcml.mcsurface.base.SurfaceLayouts or None
            The created surface layouts or None.

        Note
        ----
        Overload this method to create custom surface layouts. Default
        implementation returns None, i.e. surface layouts are not used.
        '''
        return None

    def create_detectors(self) -> mc.mcdetector.Detectors or None:
        '''
        Creates and returns a surface detectors or None if surface detectors
        are not used.

        Returns
        -------
        detectors: xopto.mcml.mcdetector.base.Detectors or None
            The created surface detectors or None.

        Note
        ----
        Overload this method to create custom surface detectors. Default
        implementation returns a radial detector
        :py:class:`~xopto.mcml.mcdetector.radial.Radial`
        on the top and bottom 
        sample surfaces with a range of 5 mm and concentric ring size of 10 μm
        (500 rings). Specular detector uses a total reflectance
        detector :py:class:`~xopto.mcml.mcdetector.total.Total`
        '''
        return mc.mcdetector.Detectors(
            top=mc.mcdetector.Radial(0.0, 5.0, 500),
            bottom=mc.mcdetector.Radial(0.0, 5.0, 500),
            specular=mc.mcdetector.Total()
        )

    def create_fluence(self) -> mc.mcfluence.Fluence or None:
        '''
        Creates and returns a flunce detector or None if fluence detector
        is not used.

        Returns
        -------
        detectors: xopto.mcbase.mcfluence.fluence.Fluence or None
            The created fluence detector or None.

        Note
        ----
        Overload this method to create a custom fluence detector. Default
        implementation returns None, i.e. fluence detector is not used.
        '''
        return None

    def create_trace(self) -> mc.mctrace.Trace or None:
        '''
        Creates and returns a trace detector or None if trace detector
        is not used.

        Returns
        -------
        detectors: xopto.mcbase.mctrace.Trace or None
            The created trace detector or None.

        Note
        ----
        Overload this method to create a custom trace detector. Default
        implementation returns None, i.e. trace detector is not used.
        '''
        return None

    def create_mc(self, *args, **kwargs) -> mc.Mc:
        '''
        Creates an instance of MC simulator.

        args: tuple
            Positional and keyword arguments passed to the MC constructor
            :py:meth:`xopto.mcml.mc.Mc`.
        kwargs: dict
            Keyword arguments passed to the MC constructor
            :py:meth:`xopto.mcml.mc.Mc`.

        Returns
        -------
        mc_obj: xopto.mcml.mc.Mc
            The created simulator instance.

        Note
        ----
        The MC components are created in the following order:

        1. The packet source (:py:meth:`McRunHelper.create_source`)
        2. The surface layouts (:py:meth:`McRunHelper.create_surface`)
        3. The layer stack (:py:meth:`McRunHelper.create_layers`)
        4. The surface detectors (:py:meth:`McRunHelper.create_detectors`)
        5. The fluence detector (:py:meth:`McRunHelper.create_fluence`)
        6. The trace detector (:py:meth:`McRunHelper.create_trace`)

        '''
        source = self.create_source()
        surface = self.create_surface()
        layers = self.create_layers()
        detectors = self.create_detectors()
        fluence = self.create_fluence()
        trace = self.create_trace()
        
        mc_obj = mc.Mc(
            layers, source, detectors,
            trace=trace, fluence=fluence, surface=surface,
            *args, **kwargs)

        return mc_obj

    def update_layers(self):
        '''
        Updates the layer stack for the next MC simulation.

        Note
        ----
        Overload this method to customize the update process of the layer
        stack. Default implementation does nothing.
        '''
        pass

    def update_source(self):
        '''
        Updates the packet source for the next MC simulation.

        Note
        ----
        Overload this method to customize the update process of the packet
        source. Default implementation does nothing.
        '''

        pass

    def update_surface(self):
        '''
        Updates the surface layouts for the next MC simulation.

        Note
        ----
        Overload this method to customize the update process of the surface
        layouts. Default implementation does nothing.
        '''
        pass

    def update_detectors(self):
        '''
        Updates the surface detectors for the next MC simulation.

        Note
        ----
        Overload this method to customize the update process of the surface
        detectors. Default implementation does nothing.
        '''
        pass

    def update_fluence(self):
        '''
        Updates the fluence detector the next MC simulation.

        Note
        ----
        Overload this method to customize the update process of the fluence
        detector. Default implementation does nothing.
        '''
        pass

    def update_trace(self):
        '''
        Updates the trace detector for the next MC simulation.

        Note
        ----
        Overload this method to customize the update process of the trace
        detector. Default implementation does nothing.
        '''
        pass

    def update_rmax(self):
        '''
        Update the maximum simulation radius.

        Note
        ----
        Overload this method to customize the update process of the simulation
        radius. Default implementation does nothing.
        '''
        pass

    def verbose_print(self):
        '''
        Called when verbose progress mode is enabled.

        Note
        ----
        Overload this method to customize the update process of the trace
        detector. Default implementation does nothing.
        '''
        pass

    def update(self):
        '''
        Updates the MC simulator for the next simulation.
        The default update order is:

        1. The layer stack (:py:meth:`McRunHelper.update_layers`).
        2. The packet source (:py:meth:`McRunHelper.update_source`).
        3. The packet surface layout (:py:meth:`McRunHelper.update_surface`).
        4. The surface detectors. (:py:meth:`McRunHelper.update_detectors`).
        5. The fluence detector (:py:meth:`McRunHelper.update_fluence`).
        6. The trace detector (:py:meth:`McRunHelper.update_trace`).
        7. The simulation radius (:py:meth:`McRunHelper.update_rmax`).
        '''
        self.update_layers()
        self.update_source()
        self.update_surface()
        self.update_detectors()
        self.update_fluence()
        self.update_trace()
        self.update_rmax()

    def _get_mc_obj(self) -> mc.Mc:
        return self._mc_obj
    mc_obj = property(_get_mc_obj, None, None, 'Monte Carlo simulator instance')

    def collect_mc_config(self) -> dict:
        '''
        Collects the entire MC configuration and creates a dict with the
        following keys:

        - "source" - the packet source configuration
        - "surface" - the surface layout configuration (optional component)
        - "layers" - the layer stack configuration
        - "detectors" - the surface detectors (optional component)
        - "fluence" - the fluence detector configuration (optional component)
        - "trace" - the trace detector configuration (optional component)
        - "rmax" - the simulation radius
        - "run_report" - the simulation run report dict as returned by
                         the :py:attr:`xopto.mcml.mc.Mc.run_report` property.

        The values of unused optional components are set to None.

        Returns
        -------
        cfg: dict
            Complete configuration of the MC simulator as a dict with the
            listed keys. Depending on the configuration of the MC simulator,
            some or all of the values under the keys of optional components
            might be None.
        '''
        mc_obj = self.mc_obj
        return {
                'source': mc_obj.source.todict(),
                'surface': mc_obj.surface.todict()
                           if mc_obj.surface is not None else None,
                'layers': mc_obj.layers.todict(),
                'detectors': mc_obj.detectors.todict()
                             if mc_obj.detectors is not None else None,
                'fluence': mc_obj.fluence.todict()
                           if mc_obj.fluence is not None else None,
                'trace': mc_obj.trace.todict()
                           if mc_obj.trace is not None else None,
                'rmax': mc_obj.rmax,
                'run_report': dict(mc_obj.run_report),
        }
    
    def collect_detectors(self, result: mc.mcdetector.Detectors) -> dict:
        '''
        Collects the simulation results from surface detectors and
        returns a dict with keys "reflectance" and "transmittance".
        Depending on the configuration of the MC simulator, the value of one
        or both keys might be None

        Parameters
        ----------
        result: xopto.mcml.mcdetector.base.Detectors
            Results of simulations or None if surface detectors were
            not configured/used.

        Returns
        -------
        data: dict
            Collected data of the surface detectors as a dict with keys
            "reflectance", "transmittance" nad "specular". 
        '''
        data = {'reflectance': None, 'transmittance': None, 'specular': None}
        if result is not None:
            if type(result.top) != mc.mcdetector.DetectorDefault:
                data['reflectance'] = result.top.reflectance
            if type(result.bottom) != mc.mcdetector.DetectorDefault:
                data['transmittance'] = result.bottom.reflectance
            if type(result.specular) != mc.mcdetector.DetectorDefault:
                data['specular'] = result.specular.reflectance
        return data

    def collect_fluence(self, result: mc.mcfluence.Fluence) -> dict or None:
        '''
        Collects the simulator results from fluence detector.

        Parameters
        ----------
        result: xopto.mcbase.mcfluence.fluence.Fluence
            Results of fluence simulations or None if fluence detector was
            not configured/used.

        Returns
        -------
        data: dict or None
            Collected data of the fluence detector as a dict with key
            "data". None if fluence detector was not configured/used.
        '''
        data = {'data': None}
        if result is not None:
            data['data'] = result.data
        return data

    def collect_trace(self, result: mc.mctrace.Trace) -> dict or None:
        '''
        Collects the simulator results from trace detector.

        Parameters
        ----------
        result: xopto.mcbase.mctrace.Trace
            Results of trace simulations or None if trace detector was
            not configured/used.

        Returns
        -------
        data: dict or None
            Collected data of the trace detector as a dict with key
            "data" and "n". None if trace detector was not configured/used.
        '''
        if result is not None:
            return {'data': result.data, 'n': result.n}
        return None

    def run_one(self, nphotons: int or Callable[['McRunHelper'], int],
                verbose: bool = False, *args, **kwargs) -> dict:
        '''
        Run one simulation.

        Parameters
        ----------
        nphotons: int
            Run the MC simulation with the specified number of packets.
            Can be a callable that takes the instance of :py:class:`McRunHelper`
            and returns the number of packets.
        verbose: bool
            Enables verbose mode.
        *args, **kwargs: 
            Positional and keyword arguments passed to the
            :py:meth:`~xopto.mcml.mc.Mc.run` method.

        Returns
        -------
        data: dict
            Simulation results as dict with keys:
            - "epoch_timestamp" - epoch timestamp
            - "mc" - simulator configuration
            - "detectors" - simulation results for surface detectors
            - "fluence" - simulation results for fluence detector
            - "trace" - trace results for trace detector

        '''
        if not isinstance(nphotons, (int, float, np.integer)):
            nphotons = nphotons(self)

        self.update()

        if verbose:
            self.verbose_print()

        trace_res, fluence_res, detectors_res = \
            self.mc_obj.run(nphotons, *args, **kwargs)

        return self.collect_one(
            nphotons, trace_res, fluence_res, detectors_res)

    def collect_one(self, nphotons: int,
                    trace_res: mc.mctrace.Trace = None,
                    fluence_res: mc.mcfluence.Fluence = None,
                    detectors_res: mc.mcdetector.Detectors = None) -> dict:
        '''
        Collect simulation results and configuration from on simulation
        run.

        Parameters
        ----------
        nphotons: int
            Number of packets that were used in the simulation.
        trace_res: xopto.mcbase.mctrace.Trace or None
            Simulation results for trace as returned
            by the :py:meth:`~xopto.mcml.mc.Mc.run` method.
        fluence_res: xopto.mcbase.mcfluence.fluence.Fluence or None
            Simulation results for fluence as returned
            by the :py:meth:`~xopto.mcml.mc.Mc.run` method.
        detectors_res: xopto.mcml.mcdetector.base.Detectors or None
            Simulation results for surface detectors as returned
            by the :py:meth:`~xopto.mcml.mc.Mc.run` method.

        Returns
        -------
        data: dict
            Collected simulation configuration and results as a dict wit
            keys:

            - 'epoch_time' - epoch timestamp as returned by the
              :py:meth:`time.time` method.
            - 'mc' - simulator configuration as returned by the
              :py:meth:`McRunHelper.collect_mc_config` method
            - 'num_packets' - number of launched packets
            - 'detectors' - collected surface detector result (optional)
              returned by the :py:meth:`McRunHelper.collect_detectors` method.
            - 'fluence' - collected fluence detector results (optional) as
              returned by the :py:meth:`McRunHelper.collect_fluence` method.
            - 'trace' - collected trace detector results (optional) as
              returned by the :py:meth:`McRunHelper.collect_trace` method.

            Unused optional items are set to None.
        '''
        data = {
            'epoch_timestamp': time.time(),
            'mc': self.collect_mc_config(),
            'num_packets': nphotons,
        }
        data['detectors'] = self.collect_detectors(detectors_res)
        data['fluence'] = self.collect_fluence(fluence_res)
        data['trace'] = self.collect_trace(trace_res)

        return data
        
    def run_batch(self, size: int, nphotons: int or Callable[['McRunHelper'], int],
                  first: int = 0, verbose: bool = False,
                  *args, **kwargs) -> List[Dict]:
        '''
        Run a batch of simulations. The MC configuration is updated before
        each simulation by a call to the :py:meth:`McRunHelper.update`.

        Parameters
        ----------
        size: int
            Batch size, i.e. the number of simulations to run.
        nphotons: int
            Run the MC simulation with the specified number of packets.
            Can be a callable that takes the instance of :py:class:`McRunHelper`
            and returns the number of packets.
        first: int
            Index of the first sample in the batch.
        verbose: bool
            Print verbose information using the
            :py:meth:`McRunHelper.verbose_print`.
        *args, **kwargs: 
            Positional and keyword arguments passed to the
            :py:meth:`~xopto.mcml.mc.Mc.run` method. 

        Returns
        -------
        results: List[Dict]
            A list of simulation results as returns by the
            :py:meth`McRunHelper.run_one` method. An extra key "index"
            is added to each item. The value of this key is the zero-based
            index of the sample computed from the input argument "first + i",
            where "i" is the zero-based sample index within the batch. 
        '''
        batch_data = []
        batch_start_t = time.perf_counter()
        rate = float('nan')
        for i in range(size):
            if verbose:
                print('{}-{}: {}/{} - {:.1f} samples/h'.format(
                    first, first + size, i + 1, size, rate))
            data = self.run_one(nphotons,verbose=verbose, *args, **kwargs)
            data['index'] = first + i
            batch_data.append(data)
            dt_start = time.perf_counter() - batch_start_t
            rate = (i + 1)/dt_start*3600

        return batch_data
