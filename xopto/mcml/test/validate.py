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

import os
import os.path
import pickle
import time

import numpy as np
from scipy import io
from scipy.integrate import simps

from xopto.mcbase.mctest import McTest

from xopto.mcml import mc
from xopto.mcml import mctypes
from xopto.mcml import mcpf
from xopto.mcml import mclayer
from xopto.mcml import mcdetector
from xopto.mcml import mctrace
from xopto.mcml import mcsource
from xopto.mcml import mcoptions
from xopto.mcml import mcsurface
from xopto.mcml.mcutil import fiber as fiberutil
from xopto.mcml.mcutil import fourier

from xopto.pf import Hg, Gk
from xopto.pf.util import GkMap

from xopto.cl import clinfo

from xopto.util import convolve
from xopto.pf import Hg


os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
DATA_DIR = datadir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'reference'))

def mc_options(cfg: dict) -> dict:
    '''
    Extract options that can be passed to the Monte Carlo simulator
    constructor :py:meth:`xopto.mcml.mc.Mc`.

    Parameters
    ----------
    cfg: dict
        Group of options in a dict.
    Returns
    -------
    options: dict
        The options from cfg that can be passed to the simulator constructor.
    '''
    return {
        'options': cfg.get('options', []),
        'cl_devices': cfg.get('cl_devices', []),
        'cl_build_options': cfg.get('cl_build_options', []),
        'types': cfg.get('types', mctypes.McDataTypesSingle)
    }

def run_options(cfg: dict) -> dict:
    '''
    Extract options that can be passed to the Monte Carlo run
    method :py:meth:`xopto.mcml.mc.Mc.run`.

    Parameters
    ----------
    cfg: dict
        Group of options in a dict.
    Returns
    -------
    options: dict
        The options from cfg that can be passed to the simulator run method.
    '''
    return {
        'nphotons': cfg.get('nphotons', 10e6),
        'wgsize': cfg.get('wgsize'),
        'exportsrc': cfg.get('exportsrc'),
        'verbose': cfg.get('verbose', False)
    }

def _show_mua_musr_pf_param_grid(
        pf_param_name: str, grid: tuple or list, data: np.ndarray,
        sds: list or tuple = None, show_mean: bool = True,
        units: str = None, title: str = None, maximize: bool = False):
    '''
    A helper function for displaying 3D lut data with a singe
    scattering phase function parameter.

    Parameters
    ----------
    pf_param_name: str
        Scattering phase function parameter name.
    grid: list or tuple
        List or tuple of vectors forming the domain, e.g. (mua, musr, gamma).
        The full grid is formed by np.meshgrid(*grid, indexing='ij')
    data: np.ndarray
        4D numpy array of spatially resolved reflectance or corresponding
        errors organized as [mua, musr, gamma, sds].
    sds: list of tuple
        List of source-detector separations (m).
    show_mean: bool
        If set to True, show mean value of the slices.
    units: str
        Units of the mean value.
    title: str
        Figure title.
    maximize: bool
        Maximize the figure on screen.
    '''
    import matplotlib.pyplot as pp

    if units is None:
        units = ''
    if pf_param_name is None:
        pf_param_name = ''

    data_min = np.min(data, axis=(0, 1, 2))
    data_max = np.max(data, axis=(0, 1, 2))

    axes = []
    images = []
    titles = []
    fig = pp.figure()
    num_sds = data.shape[-1]
    mua, musr, pf_param = grid
    mua_mur_extent = (musr[0]*1e-2, musr[-1]*1e-2, mua[0]*1e-2, mua[-1]*1e-2)
    
    if title is not None:
        pp.suptitle(title)

    pf_param_min = pf_param.min()
    pf_param_max = pf_param.max()

    if pf_param.size > 1:
        pf_param_step = (pf_param_max - pf_param_min)/(pf_param.size - 1)
    else:
        pf_param_step = 0.0

    for sds_index in range(num_sds):
        data_slice = data[:, :, 0, sds_index]
        axes.append(pp.subplot(2, 3, sds_index + 1))
        pp.xlabel('musr (1/cm)')
        pp.ylabel('mua (1/cm)')
        if sds is not None:
            if show_mean:
                titles.append(pp.title('SDS: {:.1f} μm, mean:{:.3f}{}'.format(
                    sds[sds_index]*1e6, data_slice.mean(), units)))
            else:
                titles.append(
                    pp.title('SDS: {:.1f} μm'.format(sds[sds_index]*1e6)))

        images.append(
            pp.imshow(data_slice, extent=mua_mur_extent, aspect='auto'))
        cb = pp.colorbar()
        cb.mappable.set_clim(data_min[sds_index], data_max[sds_index])

    axes.append(pp.subplot(2, 3, 6))
    pp.plot(pf_param, '.g')
    pp.xlabel('N/A')
    pp.ylabel(pf_param_name)
    pf_param_pos_plot, = pp.plot(0, pf_param[0],
                                 'o', fillstyle='none', color='r')
    pp.title('Click on plot to select a mua-musr slice')

    if maximize:
        _pp_maximize()

    #pp.tight_layout()

    def on_pf_param_pos_plot_click(event):
        if event.inaxes == axes[-1]:
            # print(event.xdata, event.ydata)

            pf_param_index = np.round((event.ydata - pf_param_min)/pf_param_step)
            pf_param_index = int(min(max(pf_param_index, 0), pf_param.size - 1))

            pf_param_pos_plot.set_data(
                pf_param_index, pf_param[pf_param_index])

            for sds_index in range(num_sds):
                data_slice = data[:, :, pf_param_index, sds_index]
                images[sds_index].set_array(data_slice)
                if sds is not None:
                    if show_mean:
                        titles[sds_index].set_text(
                            'SDS: {:.1f} μm, mean:{:.3f}{}'.format(
                                sds[sds_index]*1e6, data_slice.mean(), units)
                        )
                    else:
                        titles[sds_index].set_text(
                            'SDS: {:.1f} μm'.format(sds[sds_index]*1e6)
                        )
            pp.draw()

    fig.canvas.mpl_connect('button_press_event', on_pf_param_pos_plot_click)

    return fig


class SingleLayerLineSourceRadialProfile(McTest):
    def __init__(self, usepflut: bool = False, verbose: bool = False,
                 **kwargs):
        '''
        Test of the Monte Carlo simulator for a singel layer sample and
        a normally incident thin line source. The radial reflectance
        distribution is computed for a single set of optical properties and
        compared to the reference data.

        Parameters
        ----------
        usepflut: bool
            Run the simulator with lookup table-based sampling of the
            scattering phase function.
        verbose: bool
            Turn on/off verbose output.
        kwargs: dict
            Parameters passed to the Monte Carlo simulator constructor
            :py:meth:`xopto.mcml.mc.Mc`.
        '''

        self._verbose = verbose
        self._use_specular = kwargs.get('specular', False)

        datafile = os.path.join(DATA_DIR, 'cuda_singlelayer_linesource_radial.pkl')
        with open(datafile, 'rb') as fid:
            self._data = data = pickle.load(fid)

        if usepflut:
            pf = [mcpf.LutEx(Hg, data['layers'][0]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][1]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][2]['pf']['pfparams'])]
        else:
            pf = [mcpf.Hg(*data['layers'][0]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][1]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][2]['pf']['pfparams'])]

        layers = mclayer.Layers([
            mclayer.Layer(pf=pf[0], **data['layers'][0]['basic']),
            mclayer.Layer(pf=pf[1], **data['layers'][1]['basic']),
            mclayer.Layer(pf=pf[2], **data['layers'][2]['basic'])
        ])

        source = mcsource.Line(**data['source'])

        detector = mcdetector.Radial(
            mcdetector.RadialAxis(**data['detector']['axis']),
            cosmin=data['detector']['cosmin']
        )
        specular = None
        if self._use_specular:
            specular = mcdetector.Total()
        detectors = mcdetector.Detectors(top=detector, specular=specular)

        mc_obj = mc.Mc(layers, source, detectors, **mc_options(kwargs))
        mc_obj.rmax = data['mc']['rmax']
        self._reflectance = None

        super().__init__(mc_obj, data['reflectance'],
                         'Single layer line source radial reflectance profile')

    def run(self, nphotons: int or float = 10e6, **kwargs) -> float:
        '''
        Run the test.

        Parameters
        ----------
        nphotons: int
            The number of photon packets to use for each simulation run.
        kwargs: dict
            Parameters passed to the Monte carlo run method
            :py:meth:`xopto.mcml.mc.Mc.run`.

        Returns
        -------
        tmc: float
            Time (s) consumed by the simulator core.
        '''
        t_start = time.perf_counter()
        _, _, detectors_res = self.mc.run(nphotons, **kwargs)

        if self._verbose:
            dt = time.perf_counter() - t_start
            print(self.progress_str(1, 1, dt))

        t_mc = time.perf_counter() - t_start

        simulated = detectors_res.top.reflectance
        simulated.shape = self.reference.shape
        self.simulated = simulated
        self._r = detectors_res.top.r

        self._rel_error = (simulated - self.reference)/self.reference*100.0
        self.error = np.mean(self._rel_error)

        return t_mc

    def visualize(self):
        '''
        Visualize the results of the test.
        '''
        if self.simulated is None:
            return


        from matplotlib import pyplot as pp
        fig, axis = pp.subplots(1, 2)
        pp.suptitle('Test of: {}'.format(self.description))
        ax = axis.flat[0]
        ax.semilogy(self._r, self.simulated, '-r', label='Simulated')
        ax.semilogy(self._r, self.reference, '-g', label='Reference')
        ax.set_xlabel('Radius (m)')
        ax.set_ylabel('Reflectance')
        ax.legend(loc='upper right')

        ax = axis.flat[1]
        ax.plot(self._r, self._rel_error)
        ax.set_ylabel('Relative error (%)')
        ax.set_xlabel('Radius (m)')

    def passed(self) -> bool:
        '''
        Returns True if the simulator passed the test, else False.
        '''
        if self.error is not None:
            return np.abs(self.error).max() < 0.5
        return False


class SingleLayerUniformFiberRadial(McTest):
    def __init__(self, usepflut: bool = False, verbose: bool = False,
                 **kwargs):
        '''
        Test of the Monte Carlo simulator for a singel layer sample and
        a normally incident optical fiber source. The reflectance is computed
        for a wide range of optical properties and compared to reference data.
        A Radial detector is used to collect the reflectance.
        The comparison to reference data is made for linearly placed fibers.

        Parameters
        ----------
        usepflut: bool
            Run the simulator with lookup table-based sampling of the
            scattering phase function.
        verbose: bool
            Turn on/off verbose output.
        kwargs: dict
            Parameters passed to the Monte Carlo simulator constructor
            :py:meth:`xopto.mcml.mc.Mc`.
        '''

        self._verbose = verbose
        self._use_specular = kwargs.get('specular', False)

        datafile = os.path.join(DATA_DIR, 'cuda_singlelayer_uniformfiber.pkl')
        with open(datafile, 'rb') as fid:
            self._data = data = pickle.load(fid)

        if usepflut:
            pf = [mcpf.LutEx(Hg, data['layers'][0]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][1]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][2]['pf']['pfparams'])]
        else:
            pf = [mcpf.Hg(*data['layers'][0]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][1]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][2]['pf']['pfparams'])]

        layers = mclayer.Layers([
            mclayer.Layer(pf=pf[0], **data['layers'][0]['basic']),
            mclayer.Layer(pf=pf[1], **data['layers'][1]['basic']),
            mclayer.Layer(pf=pf[2], **data['layers'][2]['basic'])
        ])

        source = mcsource.UniformFiberNI(
            fiberutil.MultimodeFiber(**data['source'])
        )

        detector = mcdetector.Radial(
            mcdetector.RadialAxis(**data['detector']['axis']),
            cosmin=data['detector']['cosmin']
        )
        specular = None
        if self._use_specular:
            specular = mcdetector.Total()
        detectors = mcdetector.Detectors(top=detector, specular=specular)

        mc_obj = mc.Mc(layers, source, detectors, **mc_options(kwargs))
        mc_obj.rmax = data['mc']['rmax']
        self._fiber_reflectance = None

        super().__init__(mc_obj, data['reflectance'],
                         'Single layer uniform fiber radial detector '
                         'lookup table')

    def run(self, nphotons: int or float = 10e6, **kwargs) -> float:
        '''
        Run the test.

        Parameters
        ----------
        nphotons: int
            The number of photon packets to use for each simulation run.
        kwargs: dict
            Parameters passed to the Monte carlo run method
            :py:meth:`xopto.mcml.mc.Mc.run`.

        Returns
        -------
        tmc: float
            Time (s) consumed by the simulator core.
        '''
        mua_vector = self._data['lut']['mua']
        musr_vector = self._data['lut']['musr']
        mua, musr = np.meshgrid(mua_vector, musr_vector, indexing='ij')

        g = self._data['layers'][1]['pf']['pfparams'][0]
        N = mua.size

        t_start = time.perf_counter()

        reflectance =[]
        for i in range(N):
            self.mc.layers[1].mua = mua.flat[i]
            self.mc.layers[1].mus = musr.flat[i]/(1.0 - g)
            _, _, detectors_res = self.mc.run(nphotons, **kwargs)
            reflectance.append(detectors_res.top.reflectance)

            if self._verbose:
                dt = time.perf_counter() - t_start
                print(self.progress_str(i + 1, N, dt), end='\r')
                if i == N - 1:
                    print()
        t_mc = time.perf_counter() - t_start

        sds = self._data['detector']['fibers']['sds']
        dcore = self._data['detector']['fibers']['dcore']
        simulated = convolve.fiber_reflectance(
            detectors_res.top.r, reflectance, sds=sds, dcore=dcore)
        simulated.shape = self.reference.shape
        self.simulated = simulated

        self._rel_error = (simulated - self.reference)/self.reference*100.0
        self.error = self._rel_error

        return t_mc

    def visualize(self):
        '''
        Visualize the results of the test.
        '''
        if self.simulated is None:
            return
        
        sds = self._data['detector']['fibers']['sds']
        nfib = len(sds)

        from matplotlib import pyplot as pp
        musr_vector = self._data['lut']['musr']
        mua_vector = self._data['lut']['mua']
        extent = (musr_vector[0]*1e-2, musr_vector[-1]*1e-2,
                  mua_vector[0]*1e-2, mua_vector[-1]*1e-2)
        
        nrows, ncols = int(np.ceil(nfib/3)), 3
        fig, axis = pp.subplots(nrows, ncols)
        pp.suptitle('Test of: {} - Rel. error(SDS) (mean={:.2f})%'.format(
            self.description, self._rel_error.mean()))

        for i in range(nfib):
            ax = axis.flat[i]
            pos = ax.imshow(self._rel_error[..., i], extent=extent,
                            aspect='auto', origin='lower')
            ax.set_xlabel('Red. scattering (1/cm)')
            ax.set_ylabel('Absorption (1/cm)')
            ax.set_title('SDS={:.1f} um, error={:.2f}%'.format(
                sds[i]*1e6, self._rel_error[...,i].mean()))
            fig.colorbar(pos, ax=ax)
        
        ax = axis.flat[-1]
        labels = []
        for d in sds:
            labels.append('{:.1f} um'.format(d*1e6))
        ax.set_title('Rel. error (%) as f(SDS)')
        new_shape = (int(self._rel_error.size//nfib), nfib)
        ax.boxplot(np.reshape(self._rel_error, new_shape) , labels=labels)
        ax.set_ylabel('Rel. error (%)')

    def passed(self) -> bool:
        '''
        Returns True if the simulator passed the test, else False.
        '''
        if self.error is not None:
            return np.abs(self.error.mean()) < 0.5
        return False


class SingleLayerUniformFiberCartesian(McTest):
    def __init__(self, usepflut: bool = False, verbose: bool = False,
                 **kwargs):
        '''
        Test of the Monte Carlo simulator for a singel layer sample and
        a normally incident optical fiber source. The reflectance is computed
        for a wide range of optical properties and compared to reference data.
        A Cartesian detector is used to collect the reflectance.
        The comparison to reference data is made for linearly placed fibers.

        Parameters
        ----------
        usepflut: bool
            Run the simulator with lookup table-based sampling of the
            scattering phase function.
        verbose: bool
            Turn on/off verbose output.
        kwargs: dict
            Parameters passed to the Monte Carlo simulator constructor
            :py:meth:`xopto.mcml.mc.Mc`.
        '''

        self._verbose = verbose
        self._use_specular = kwargs.get('specular', False)

        datafile = os.path.join(DATA_DIR, 'cuda_singlelayer_uniformfiber.pkl')
        with open(datafile, 'rb') as fid:
            self._data = data = pickle.load(fid)

        if usepflut:
            pf = [mcpf.LutEx(Hg, data['layers'][0]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][1]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][2]['pf']['pfparams'])]
        else:
            pf = [mcpf.Hg(*data['layers'][0]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][1]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][2]['pf']['pfparams'])]

        layers = mclayer.Layers([
            mclayer.Layer(pf=pf[0], **data['layers'][0]['basic']),
            mclayer.Layer(pf=pf[1], **data['layers'][1]['basic']),
            mclayer.Layer(pf=pf[2], **data['layers'][2]['basic'])
        ])

        source = mcsource.UniformFiberNI(
            fiberutil.MultimodeFiber(**data['source'])
        )

        nr = data['detector']['axis']['n']
        rmin = data['detector']['axis']['start']
        rmax = data['detector']['axis']['stop']
        dr = rmax/nr
        cosmin = data['detector']['cosmin']
        nx = np.ceil(nr*2.2/2.0)*2.0
        ny = np.ceil(nr*2.1/2.0)*2.0
        detector = mcdetector.Cartesian(
            mcdetector.Axis(-dr*nx*0.5, dr*nx*0.5, nx),
            mcdetector.Axis(-dr*ny*0.5, dr*ny*0.5, ny),
            cosmin=cosmin
        )
        specular = None
        if self._use_specular:
            specular = mcdetector.Total()
        detectors = mcdetector.Detectors(top=detector, specular=specular)

        self._raxis = mcdetector.RadialAxis(rmin, rmax, nr)

        mc_obj = mc.Mc(layers, source, detectors, **mc_options(kwargs))
        mc_obj.rmax = data['mc']['rmax']
        self._fiber_reflectance = None

        super().__init__(mc_obj, data['reflectance'],
                         'Single layer uniform fiber cartesian detector '
                         'lookup table')

    def run(self, nphotons: int or float = 10e6, **kwargs) -> float:
        '''
        Run the test.

        Parameters
        ----------
        nphotons: int
            The number of photon packets to use for each simulation run.
        kwargs: dict
            Parameters passed to the Monte carlo run method
            :py:meth:`xopto.mcml.mc.Mc.run`.

        Returns
        -------
        tmc: float
            Time (s) consumed by the simulator core.
        '''
        mua_vector = self._data['lut']['mua']
        musr_vector = self._data['lut']['musr']
        mua, musr = np.meshgrid(mua_vector, musr_vector, indexing='ij')

        g = self._data['layers'][1]['pf']['pfparams'][0]
        N = mua.size

        t_start = time.perf_counter()
        t_mc = 0.0
        reflectance =[]
        for i in range(N):
            t1 = time.perf_counter()
            self.mc.layers[1].mua = mua.flat[i]
            self.mc.layers[1].mus = musr.flat[i]/(1.0 - g)
            _, _, detectors_res = self.mc.run(nphotons, **kwargs)
            t_mc = t_mc + time.perf_counter() - t1

            reflectance.append(detectors_res.top.radial(self._raxis).reflectance)

            if self._verbose:
                dt = time.perf_counter() - t_start
                print(self.progress_str(i + 1, N, dt), end='\r')
                if i == N - 1:
                    print()

        sds = self._data['detector']['fibers']['sds']
        dcore = self._data['detector']['fibers']['dcore']
        simulated = convolve.fiber_reflectance(
            self._raxis.centers, reflectance, sds=sds, dcore=dcore)
        simulated.shape = self.reference.shape
        self.simulated = simulated

        self._rel_error = (simulated - self.reference)/self.reference*100.0
        self.error = self._rel_error

        return t_mc

    def visualize(self):
        '''
        Visualize the results of the test.
        '''
        if self.simulated is None:
            return
        
        sds = self._data['detector']['fibers']['sds']
        nfib = len(sds)

        from matplotlib import pyplot as pp
        musr_vector = self._data['lut']['musr']
        mua_vector = self._data['lut']['mua']
        extent = (musr_vector[0]*1e-2, musr_vector[-1]*1e-2,
                  mua_vector[0]*1e-2, mua_vector[-1]*1e-2)
        
        nrows, ncols = int(np.ceil(nfib/3)), 3
        fig, axis = pp.subplots(nrows, ncols)
        pp.suptitle('Test of: {} - Rel. error(SDS) (mean={:.2f})%'.format(
            self.description, self._rel_error.mean()))

        for i in range(nfib):
            ax = axis.flat[i]
            pos = ax.imshow(self._rel_error[..., i], extent=extent,
                            aspect='auto', origin='lower')
            ax.set_xlabel('Red. scattering (1/cm)')
            ax.set_ylabel('Absorption (1/cm)')
            ax.set_title('SDS={:.1f} um, error={:.2f}%'.format(
                sds[i]*1e6, self._rel_error[..., i].mean()))
            fig.colorbar(pos, ax=ax)
        
        ax = axis.flat[-1]
        labels = []
        for d in sds:
            labels.append('{:.1f} um'.format(d*1e6))
        ax.set_title('Rel. error (%) as f(SDS)')
        new_shape = (int(self._rel_error.size//nfib), nfib)
        ax.boxplot(np.reshape(self._rel_error, new_shape) , labels=labels)
        ax.set_ylabel('Rel. error (%)')

    def passed(self) -> bool:
        '''
        Returns True if the simulator passed the test, else False.
        '''
        if self.error is not None:
            return np.abs(self.error.mean()) < 0.5
        return False


class DoubleLayerUniformFiberRadial(McTest):
    def __init__(self, usepflut: bool = False, verbose: bool = False,
                 **kwargs):
        '''
        Test of the Monte Carlo simulator for a double layer sample and
        a normally incident optical fiber source. The reflectance is computed
        for a wide range of optical properties and compared to reference data.
        A Cartesian detector is used to collect the reflectance.
        The comparison to reference data is made for linearly placed fibers.

        Parameters
        ----------
        usepflut: bool
            Run the simulator with lookup table-based sampling of the
            scattering phase function.
        verbose: bool
            Turn on/off verbose output.
        kwargs: dict
            Parameters passed to the Monte Carlo simulator constructor
            :py:meth:`xopto.mcml.mc.Mc`.
        '''
        self._verbose = verbose
        self._use_specular = kwargs.get('specular', False)

        datafile = os.path.join(DATA_DIR, 'cuda_doublelayer_uniformfiber.pkl')
        with open(datafile, 'rb') as fid:
            self._data = data = pickle.load(fid)

        if usepflut:
            pf = [mcpf.LutEx(Hg, data['layers'][0]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][1]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][2]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][3]['pf']['pfparams'])]
        else:
            pf = [mcpf.Hg(*data['layers'][0]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][1]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][2]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][3]['pf']['pfparams'])]

        layers = mclayer.Layers([
            mclayer.Layer(pf=pf[0], **data['layers'][0]['basic']),
            mclayer.Layer(pf=pf[1], **data['layers'][1]['basic']),
            mclayer.Layer(pf=pf[2], **data['layers'][2]['basic']),
            mclayer.Layer(pf=pf[3], **data['layers'][3]['basic'])
        ])

        source = mcsource.UniformFiberNI(
            fiberutil.MultimodeFiber(**data['source'])
        )

        detector = mcdetector.Radial(
            mcdetector.RadialAxis(**data['detector']['axis']),
            cosmin=data['detector']['cosmin']
        )
        specular = None
        if self._use_specular:
            specular = mcdetector.Total()
        detectors = mcdetector.Detectors(top=detector, specular=specular)

        mc_obj = mc.Mc(layers, source, detectors, **mc_options(kwargs))
        mc_obj.rmax = data['mc']['rmax']
        self._fiber_reflectance = None

        super().__init__(mc_obj, data['reflectance'],
                         'Double layer uniform fiber radial detector '
                         'lookup table')

    def run(self, nphotons: int or float = 10e6, **kwargs) -> float:
        '''
        Run the test.

        Parameters
        ----------
        nphotons: int
            The number of photon packets to use for each simulation run.
        kwargs: dict
            Parameters passed to the Monte carlo run method
            :py:meth:`xopto.mcml.mc.Mc.run`.

        Returns
        -------
        tmc: float
            Time (s) consumed by the simulator core.
        '''
        mua_1 = self._data['lut'][0]['mua']
        musr_1 = self._data['lut'][0]['musr']

        mua_2_vector = self._data['lut'][1]['mua']
        musr_2_vector = self._data['lut'][1]['musr']
        mua_2, musr_2 = np.meshgrid(mua_2_vector, musr_2_vector, indexing='ij')

        g_1 = self._data['layers'][1]['pf']['pfparams'][0]
        g_2 = self._data['layers'][2]['pf']['pfparams'][0]
        N = mua_2.size

        t_start = time.perf_counter()

        self.mc.layers[1].mua = mua_1
        self.mc.layers[1].mus = musr_1/(1.0 - g_1)

        reflectance =[]
        for i in range(N):
            self.mc.layers[2].mua = mua_2.flat[i]
            self.mc.layers[2].mus = musr_2.flat[i]/(1.0 - g_2)
            _, _, detectors_res = self.mc.run(nphotons, **kwargs)
            reflectance.append(detectors_res.top.reflectance)

            if self._verbose:
                dt = time.perf_counter() - t_start
                print(self.progress_str(i + 1, N, dt), end='\r')
                if i == N - 1:
                    print()
        t_mc = time.perf_counter() - t_start

        sds = self._data['detector']['fibers']['sds']
        dcore = self._data['detector']['fibers']['dcore']
        simulated = convolve.fiber_reflectance(
            detectors_res.top.r, reflectance, sds=sds, dcore=dcore)
        simulated.shape = self.reference.shape
        self.simulated = simulated

        self._rel_error = (simulated - self.reference)/self.reference*100.0
        self.error = np.mean(self._rel_error, (0,1))

        return t_mc

    def visualize(self):
        '''
        Visualize the results of the test.
        '''
        if self.simulated is None:
            return
        
        sds = self._data['detector']['fibers']['sds']
        nfib = len(sds)

        from matplotlib import pyplot as pp
        musr_vector = self._data['lut'][1]['musr']
        mua_vector = self._data['lut'][1]['mua']
        extent = (musr_vector[0]*1e-2, musr_vector[-1]*1e-2,
                  mua_vector[0]*1e-2, mua_vector[-1]*1e-2)
        
        nrows, ncols = int(np.ceil(nfib/3)), 3
        fig, axis = pp.subplots(nrows, ncols)
        pp.suptitle('Test of: {} - Rel. error(SDS) (mean={:.2f})%'.format(
            self.description, self._rel_error.mean()))

        for i in range(nfib):
            ax = axis.flat[i]
            pos = ax.imshow(self._rel_error[..., i], extent=extent,
                            aspect='auto', origin='lower')
            ax.set_xlabel('Red. scattering (1/cm)')
            ax.set_ylabel('Absorption (1/cm)')
            ax.set_title('SDS={:.1f} um, error={:.2f}%'.format(
                sds[i]*1e6, self._rel_error[..., i].mean()))
            fig.colorbar(pos, ax=ax)
        
        ax = axis.flat[-1]
        labels = []
        for d in sds:
            labels.append('{:.1f} um'.format(d*1e6))
        ax.set_title('Rel. error (%) as f(SDS)')
        new_shape = (int(self._rel_error.size//nfib), nfib)
        ax.boxplot(np.reshape(self._rel_error, new_shape) , labels=labels)
        ax.set_ylabel('Rel. error (%)')

    def passed(self) -> bool:
        '''
        Returns True if the simulator passed the test, else False.
        '''
        if self.error is not None:
            return np.abs(self.error).max() < 0.5
        return False


class SingleLayerSfdiIncidence30VariableNa(McTest):
    def __init__(self, usepflut: bool = False, verbose: bool = False,
                 **kwargs):
        '''
        Test of the Monte Carlo simulator for a single layer sample and
        a 30 deg. incidence and variable detector NA.
        A :py:class:`detector.SymmetricX` detector is used to collect the
        backscattered light.

        Parameters
        ----------
        usepflut: bool
            Run the simulator with lookup table-based sampling of the
            scattering phase function.
        verbose: bool
            Turn on/off verbose output.
        kwargs: dict
            Parameters passed to the Monte Carlo simulator constructor
            :py:meth:`xopto.mcml.mc.Mc`.
        '''
        self._verbose = verbose
        self._use_specular = kwargs.get('specular', False)

        self._usepflut = usepflut

        filename =os.path.join(
            DATA_DIR,
            'single_layer_sfdi_30deg_incidence_variable_na.pkl'
        )
        with open(filename, 'rb') as fid:
            self._data = pickle.load(fid)

        reference = np.array([self._data['amplitude'], self._data['phase']])

        self._mc_options = mc_options(kwargs)

        super().__init__(None, reference,
                        'SFDI single layer 30 deg incidence variable detector '
                        'NA SymmetriX detector')

    def run(self, nphotons: int or float = 10e6, **kwargs) -> float:
        '''
        Run the test.

        Parameters
        ----------
        nphotons: int
            The number of photon packets to use for each simulation run.
        kwargs: dict
            Parameters passed to the Monte carlo run method
            :py:meth:`xopto.mcml.mc.Mc.run`.

        Returns
        -------
        tmc: float
            Time (s) consumed by the simulator core.
        '''
        data = self._data

        if nphotons is None:
            nphotons = int(4e9)

        rmax = 150e-3

        if self._usepflut:
            pf = [mcpf.LutEx(Hg, data['layers'][0]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][1]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][2]['pf']['pfparams'])]
        else:
            pf = [mcpf.Hg(*data['layers'][0]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][1]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][2]['pf']['pfparams'])]

        layers = mclayer.Layers([
            mclayer.Layer(pf=pf[0], **data['layers'][0]['basic']),
            mclayer.Layer(pf=pf[1], **data['layers'][1]['basic']),
            mclayer.Layer(pf=pf[2], **data['layers'][2]['basic'])
        ])

        # the photon packet source
        source = mcsource.Line(**data['source'])

        frequencies = data['frequencies']
        NA = data['lut']['na']
        rmax = data['mc']['rmax']

        sf_cplx = np.zeros((NA.size, frequencies.size), dtype=np.cdouble)
        t_mc = 0.0
        t_start = time.perf_counter()
        N = NA.size
        for i in range(N):
            t1 = time.perf_counter()

            cosmin = np.sqrt(1.0 - NA[i]**2)
            detector = mcdetector.SymmetricX(
                mcdetector.SymmetricAxis(**data['detector']['axis']),
                cosmin=cosmin
            )
            specular = None
            if self._use_specular:
                specular = mcdetector.Total()
            detectors = mcdetector.Detectors(top=detector, specular=specular)

            mc_obj = mc.Mc(layers, source, detectors, **self._mc_options)
            mc_obj.rmax = rmax

            _, _, detectors_res = mc_obj.run(nphotons, **kwargs)
    
            sf_cplx[i] = fourier.discreteSimpson(
                frequencies, detectors_res.top.x, detectors_res.top.reflectance,
                uneven=True)

            dt = time.perf_counter() - t1
            t_mc += dt
            if self._verbose:
                dt = time.perf_counter() - t_start
                print(self.progress_str(i + 1, N, dt), end='\r')
                if i == N - 1:
                    print()

        self._amplitude = np.abs(sf_cplx)
        self._phase = np.angle(sf_cplx)
        self._amplitude_error = \
            (self._amplitude - data['amplitude'])/data['amplitude']
        self._phase_error = (self._phase - data['phase'])

        self.simulated = np.array([self._amplitude, self._phase])
        self.error = np.array([self._amplitude_error, self._phase_error])

        if self._verbose:
            amplitude_error = \
                np.mean(np.abs(self._amplitude_error[4:]), 1).max()
            phase_error = np.mean(np.abs(self._phase_error[4:]), 1).max()
            print('SFDI errors for NA >= 10°:')
            print('   relative amplitude error {:.2f}%.'.format(
                amplitude_error*100.0))
            print('   phase error {:.2f}°.'.format(
                np.rad2deg(phase_error)))

        return t_mc

    def visualize(self):
        import matplotlib
        from matplotlib import pyplot as pp
        from cycler import cycler

        frequencies = self._data['frequencies']
        NA = self._data['lut']['na']
        n = NA.size

        colors = np.vstack([
            np.interp(np.arange(n), [0, (n - 1)/2, n - 1], [1.0, 0.0, 0.0]),
            np.interp(np.arange(n), [0, (n - 1)/2, n - 1], [0.0, 1.0, 0.0]),
            np.interp(np.arange(n), [0, (n - 1)/2, n - 1], [0.0, 0.0, 1.0]),
        ])
        matplotlib.rcParams['axes.prop_cycle'] = cycler(color=colors.T.tolist())

        fig, axes = pp.subplots(3, 1)
        ax = axes.flat[0]
        ax.set_title('Relative amplitude error')
        ax.plot(frequencies*1e-3, self._amplitude_error.T)

        ax = axes.flat[1]
        ax.set_title('Phase error (rad)')
        ax.plot(frequencies*1e-3, self._phase_error.T)
        ax.set_xlabel('Frequency (1/mm)')

        ax = axes.flat[2]
        for index in range(n):
            acceptance_deg = np.rad2deg(np.arcsin(NA[index]))
            ax.plot([], [], label='NA={:.1f}°'.format(acceptance_deg))
        ax.legend(ncol=8)
        ax.set_axis_off()

    def passed(self) -> bool:
        '''
        Returns True if the simulator passed the test, else False.
        '''
        if self.error is not None:
            # SFDI errors for NA >= 10° must be less than 2.5 %
            amplitude_error = self._amplitude_error
            amplitude_rerror = np.mean(np.abs(amplitude_error[4:]), 1).max()

            phase_error = self._phase_error
            phase_rerror = np.mean(np.abs(phase_error[4:]), 1).max()

            print(amplitude_rerror, phase_rerror)
            return  amplitude_rerror < 0.025 and phase_rerror < 0.025

        return False


class SingleLayerUniformFiberTrace(McTest):
    def __init__(self, usepflut: bool = False, verbose: bool = False,
                 **kwargs: dict):
        '''
        Test of the Monte Carlo simulator for a single layer sample and
        trace analysis.
        The results of trace analysis are compared to the data collectred by
        the :py:class`mcdetector.Radial` detector.

        Parameters
        ----------
        usepflut: bool
            Run the simulator with lookup table-based sampling of the
            scattering phase function.
        verbose: bool
            Turn on/off verbose output.
        kwargs: dict
            Parameters passed to the Monte Carlo simulator constructor
            :py:meth:`xopto.mcml.mc.Mc`.
        '''
        self._verbose = verbose
        self._use_specular = kwargs.get('specular', False)

        filename = os.path.join(DATA_DIR, 'single_layer_uniformfiber_trace.pkl')
        with open(filename, 'rb') as fid:
            self._data = data = pickle.load(fid)

        if usepflut:
            pf = [mcpf.LutEx(Hg, data['layers'][0]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][1]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][2]['pf']['pfparams'])]
        else:
            pf = [mcpf.Hg(*data['layers'][0]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][1]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][2]['pf']['pfparams'])]

        layers = mclayer.Layers([
            mclayer.Layer(pf=pf[0], **data['layers'][0]['basic']),
            mclayer.Layer(pf=pf[1], **data['layers'][1]['basic']),
            mclayer.Layer(pf=pf[2], **data['layers'][2]['basic'])
        ])

        source = mcsource.UniformFiber(
            fiberutil.MultimodeFiber(**data['source'])
        )

        detector = mcdetector.Radial(
            mcdetector.RadialAxis(**data['detector']['axis']),
            cosmin=data['detector']['cosmin']
        )
        specular = None
        if self._use_specular:
            specular = mcdetector.Total()
        detectors = mcdetector.Detectors(top=detector, specular=specular)

        trace = mctrace.Trace(**data['trace'])

        mc_obj = mc.Mc(layers, source, detectors, trace, **mc_options(kwargs))
        mc_obj.rmax = data['mc']['rmax']

        self._data = data

        super().__init__(mc_obj, None, 'Single layer uniform fiber trace')
        self._reflectance = None

    def run(self, nphotons: int or float = 1e6, **kwargs) -> float:
        '''
        Run the test.

        Parameters
        ----------
        nphotons: int
            The number of photon packets to use for each simulation run.
        kwargs: dict
            Parameters passed to the Monte carlo run method
            :py:meth:`xopto.mcml.mc.Mc.run`.

        Returns
        -------
        tmc: float
            Time (s) consumed by the simulator core.
        '''
        mua = self._data['lut']['mua']
        musr = self._data['lut']['musr']
        g1 = self._data['layers'][1]['pf']['pfparams'][0]
        N = mua.size
        nr = self.mc.detectors.top.n
    
        r = self.mc.detectors.top.edges
        dr = r[1] - r[0]
        cosmin = self.mc.detectors.top.cosmin
        arings = np.pi*(r[1:]**2 - r[:-1]**2)

        reflectance_trace = np.zeros((N, nr))
        reflectance_radial = np.zeros_like(reflectance_trace)
        pkt_inds = np.arange(nphotons)
        err = np.zeros([nr])
        t_mc = 0.0


        t_start = tstart = time.perf_counter()

        for i in range(N):
            self.mc.layers[1].mua = mua[i]
            self.mc.layers[1].mus = musr[i]/(1.0 - g1)

            t1 = time.perf_counter()
            trace_res, _, detectors_res = self.mc.run(nphotons, **kwargs)
            t_mc += time.perf_counter() - t1

            x, y, z = trace_res.data['x'], trace_res.data['y'], trace_res.data['z']
            pz, w, n = trace_res.data['pz'], trace_res.data['w'], trace_res.n
            termind = np.minimum(n - 1, trace_res.maxlen - 1)

            xterm, yterm, zterm = x[pkt_inds, termind], y[pkt_inds, termind], \
                z[pkt_inds, termind]
            pzterm, wterm = pz[pkt_inds, termind], w[pkt_inds, termind]
            rind = np.minimum(np.floor(np.sqrt(xterm**2 + yterm**2)/dr), nr - 1)
            mask = np.logical_and(zterm <= 0.0, np.abs(pzterm) >= cosmin)
            rindok = rind[mask]
            wok = wterm[mask]
            for ind in range(nr):
                reflectance_trace[i, ind] = np.sum(wok[rindok == ind])
            reflectance_trace[i] = reflectance_trace[i]/(nphotons*arings)
            reflectance_radial[i] = detectors_res.top.reflectance

            if self._verbose:
                dt = time.perf_counter() - t_start
                print(self.progress_str(i + 1, N, dt), end='\r')
                if i == N - 1:
                    print()

        self.simulated = reflectance_trace
        self.reference = reflectance_radial

        error = np.zeros_like(reflectance_radial)
        mask = reflectance_radial != 0
        error[mask] = (reflectance_trace[mask] - reflectance_radial[mask])/ \
            reflectance_radial[mask]
        error[np.logical_and(reflectance_radial == 0.0,
                             reflectance_trace != 0.0)] = np.inf
        self.error = error

        return t_mc

    def visualize(self):
        from matplotlib import pyplot as pp

        errstr = 'Maximum relative trace reflectance error {:.4f}%'.format(
            np.abs(self.error).max()*100.0)

        fig, ax = pp.subplots()
        ax.set_title(errstr)
        ax.plot(self.mc.detectors.top.r, self.simulated.T, '-xr', label='Trace')
        ax.plot(self.mc.detectors.top.r, self.reference.T, '-+g', label='Radial')
        ax.set_xlabel('Radius (m)')
        ax.set_ylabel('Reflectance (1/m^2)')
        ax.legend()

    def passed(self):
        '''
        Returns True if the simulator passed the test, else False.
        '''
        if self.error is not None:
            return np.abs(self.error).max()*100.0 < 0.1
        return False


class SingleLayerSixLinearProbe(McTest):
    def __init__(self, usepflut: bool = False, verbose: bool = False,
                 **kwargs):
        '''
        Test of the Monte Carlo simulator for a singel layer sample and
        a normally incident six-linear optical fiber probe.

        Parameters
        ----------
        usepflut: bool
            Run the simulator with lookup table-based sampling of the
            scattering phase function.
        verbose: bool
            Turn on/off verbose output.
        kwargs: dict
            Parameters passed to the Monte Carlo simulator constructor
            :py:meth:`xopto.mcml.mc.Mc`.
        '''

        self._verbose = verbose
        self._use_specular = kwargs.get('specular', False)

        datafile = os.path.join(DATA_DIR, 'single_layer_hg_lut_six_linear_probe.pkl')
        with open(datafile, 'rb') as fid:
            self._data = data = pickle.load(fid)

        if usepflut:
            pf = [mcpf.LutEx(Hg, data['layers'][0]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][1]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][2]['pf']['pfparams'])]
        else:
            pf = [mcpf.Hg(*data['layers'][0]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][1]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][2]['pf']['pfparams'])]

        layers = mclayer.Layers([
            mclayer.Layer(pf=pf[0], **data['layers'][0]['basic']),
            mclayer.Layer(pf=pf[1], **data['layers'][1]['basic']),
            mclayer.Layer(pf=pf[2], **data['layers'][2]['basic'])
        ])
        source_data = data['source']
        fiber = fiberutil.MultimodeFiber(**source_data.pop('fiber'))
        source = mcsource.UniformFiber(fiber, **source_data)

        surface_data = data['surface']
        fiber = fiberutil.MultimodeFiber(**surface_data.pop('fiber'))
        surface = mcsurface.SurfaceLayouts(
            top=mcsurface.LinearArray(fiber, **surface_data)
        )

        detector_data = data['detector']
        fiber = fiberutil.MultimodeFiber(**detector_data.pop('fiber'))
        detector = mcdetector.Detectors(
            top=mcdetector.LinearArray(fiber, **detector_data)
        )

        specular = None
        if self._use_specular:
            specular = mcdetector.Total()
        detectors = mcdetector.Detectors(top=detector, specular=specular)

        mc_obj = mc.Mc(layers, source, detectors, surface=surface,
                       **mc_options(kwargs))
        mc_obj.rmax = data['mc']['rmax']
        self._reflectance = None

        super().__init__(mc_obj, data['reflectance'],
                         'Single layer six-linear probe HG 2D mua musr lut')

    def run(self, nphotons: int or float = 10e6, **kwargs) -> float:
        '''
        Run the test.

        Parameters
        ----------
        nphotons: int
            The number of photon packets to use for each simulation run.
        kwargs: dict
            Parameters passed to the Monte carlo run method
            :py:meth:`xopto.mcml.mc.Mc.run`.

        Returns
        -------
        tmc: float
            Time (s) consumed by the simulator core.
        '''
        mua_vector = self._data['lut']['mua']
        musr_vector = self._data['lut']['musr']
        mua, musr = np.meshgrid(mua_vector, musr_vector, indexing='ij')
        g = self._data['layers'][1]['pf']['pfparams'][0]
        N = mua.size

        t_start = time.perf_counter()

        reflectance =[]
        for i in range(N):
            self.mc.layers[1].mua = mua.flat[i]
            self.mc.layers[1].mus = musr.flat[i]/(1.0 - g)
            _, _, detectors_res = self.mc.run(nphotons, **kwargs)
            reflectance.append(detectors_res.top.reflectance[1:])

            if self._verbose:
                dt = time.perf_counter() - t_start
                print(self.progress_str(i + 1, N, dt), end='\r')
                if i == N - 1:
                    print()
        t_mc = time.perf_counter() - t_start

        simulated = np.array(reflectance)
        simulated.shape = self.reference.shape
        self.simulated = simulated

        self._rel_error = (simulated - self.reference)/self.reference*100.0
        self.error = self._rel_error

        return t_mc

    def visualize(self):
        '''
        Visualize the results of the test.
        '''
        if self.simulated is None:
            return
        
        nfib = self._data['detector']['n'] - 1
        sds = np.arange(1, nfib + 1)*self._data['detector']['spacing']

        from matplotlib import pyplot as pp
        musr_vector = self._data['lut']['musr']
        mua_vector= self._data['lut']['mua']
        extent = (musr_vector[0]*1e-2, musr_vector[-1]*1e-2,
                  mua_vector[0]*1e-2, mua_vector[-1]*1e-2)
        
        nrows, ncols = int(np.ceil(nfib/3)), 3
        fig, axis = pp.subplots(nrows, ncols)
        pp.suptitle('Test of: {} - Rel. error(SDS) (mean={:.2f})%'.format(
            self.description, self._rel_error.mean()))

        for i in range(nfib):
            ax = axis.flat[i]
            pos = ax.imshow(self._rel_error[..., i], extent=extent,
                            aspect='auto', origin='lower')
            ax.set_xlabel('Red. scattering (1/cm)')
            ax.set_ylabel('Absorption (1/cm)')
            ax.set_title('SDS={:.1f} um, error={:.2f}%'.format(
                sds[i]*1e6, self._rel_error[...,i].mean()))
            fig.colorbar(pos, ax=ax)
        
        ax = axis.flat[-1]
        labels = []
        for d in sds:
            labels.append('{:.1f}'.format(d*1e6))
        ax.set_title('Rel. error (%) as f(SDS)')
        new_shape = (int(self._rel_error.size//nfib), nfib)
        ax.boxplot(np.reshape(self._rel_error, new_shape) , labels=labels, )
        ax.set_ylabel('Rel. error (%)')
        ax.set_xlabel('Source detector separation (um)')

    def passed(self) -> bool:
        ''''
        Returns True if the simulator passed the test, else False.
        '''
        if self.error is not None:
            return np.abs(self.error.mean()) < 0.5
        return False


class SingleLayerNineLinearProbe(McTest):
    def __init__(self, usepflut: bool = False, verbose: bool = False,
                 **kwargs):
        '''
        Test of the Monte Carlo simulator for a singel layer sample and
        a normally incident nine-linear optical fiber probe.

        Parameters
        ----------
        usepflut: bool
            Run the simulator with lookup table-based sampling of the
            scattering phase function.
        verbose: bool
            Turn on/off verbose output.
        kwargs: dict
            Parameters passed to the Monte Carlo simulator constructor
            :py:meth:`xopto.mcml.mc.Mc`.
        '''

        self._verbose = verbose
        self._use_specular = kwargs.get('specular', False)

        datafile = os.path.join(DATA_DIR, 'single_layer_hg_lut_nine_linear_probe.pkl')
        with open(datafile, 'rb') as fid:
            self._data = data = pickle.load(fid)

        if usepflut:
            pf = [mcpf.LutEx(Hg, data['layers'][0]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][1]['pf']['pfparams']),
                  mcpf.LutEx(Hg, data['layers'][2]['pf']['pfparams'])]
        else:
            pf = [mcpf.Hg(*data['layers'][0]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][1]['pf']['pfparams']),
                  mcpf.Hg(*data['layers'][2]['pf']['pfparams'])]

        layers = mclayer.Layers([
            mclayer.Layer(pf=pf[0], **data['layers'][0]['basic']),
            mclayer.Layer(pf=pf[1], **data['layers'][1]['basic']),
            mclayer.Layer(pf=pf[2], **data['layers'][2]['basic'])
        ])
        source_data = data['source']
        fiber = fiberutil.MultimodeFiber(**source_data.pop('fiber'))
        source = mcsource.UniformFiber(fiber, **source_data)

        surface_data = data['surface']
        fiber = fiberutil.MultimodeFiber(**surface_data.pop('fiber'))
        surface = mcsurface.SurfaceLayouts(
            top=mcsurface.LinearArray(fiber, **surface_data)
        )

        detector_data = data['detector']
        fiber = fiberutil.MultimodeFiber(**detector_data.pop('fiber'))
        detector = mcdetector.Detectors(
            top=mcdetector.LinearArray(fiber, **detector_data)
        )

        specular = None
        if self._use_specular:
            specular = mcdetector.Total()
        detectors = mcdetector.Detectors(top=detector, specular=specular)

        mc_obj = mc.Mc(layers, source, detectors, surface=surface,
                       **mc_options(kwargs))
        mc_obj.rmax = data['mc']['rmax']
        self._reflectance = None

        super().__init__(mc_obj, data['reflectance'],
                         'Single layer nine-linear probe HG 2D mua musr lut')

    def run(self, nphotons: int or float = 10e6, **kwargs) -> float:
        '''
        Run the test.

        Parameters
        ----------
        nphotons: int
            The number of photon packets to use for each simulation run.
        kwargs: dict
            Parameters passed to the Monte carlo run method
            :py:meth:`xopto.mcml.mc.Mc.run`.

        Returns
        -------
        tmc: float
            Time (s) consumed by the simulator core.
        '''
        mua_vector = self._data['lut']['mua']
        musr_vector = self._data['lut']['musr']
        mua, musr = np.meshgrid(mua_vector, musr_vector, indexing='ij')

        g = self._data['layers'][1]['pf']['pfparams'][0]
        N = mua.size

        t_start = time.perf_counter()

        reflectance =[]
        for i in range(N):
            self.mc.layers[1].mua = mua.flat[i]
            self.mc.layers[1].mus = musr.flat[i]/(1.0 - g)
            _, _, detectors_res = self.mc.run(nphotons, **kwargs)
            reflectance.append(detectors_res.top.reflectance[1:])

            if self._verbose:
                dt = time.perf_counter() - t_start
                print(self.progress_str(i + 1, N, dt), end='\r')
                if i == N - 1:
                    print()
        t_mc = time.perf_counter() - t_start

        simulated = np.array(reflectance)
        simulated.shape = self.reference.shape
        self.simulated = simulated

        self._rel_error = (simulated - self.reference)/self.reference*100.0
        self.error = self._rel_error

        return t_mc

    def visualize(self):
        '''
        Visualize the results of the test.
        '''
        if self.simulated is None:
            return
        
        nfib = self._data['detector']['n'] - 1
        sds = np.arange(1, nfib + 1)*self._data['detector']['spacing']

        from matplotlib import pyplot as pp
        musr_vector = self._data['lut']['musr']
        mua_vector = self._data['lut']['mua']
        extent = (musr_vector[0]*1e-2, musr_vector[-1]*1e-2,
                  mua_vector[0]*1e-2, mua_vector[-1]*1e-2)
        
        nrows, ncols = int(np.ceil(nfib/3)), 3
        fig, axis = pp.subplots(nrows, ncols)
        pp.suptitle('Test of: {} - Rel. error(SDS) (mean={:.2f})%'.format(
            self.description, self._rel_error.mean()))

        for i in range(nfib):
            ax = axis.flat[i]
            pos = ax.imshow(self._rel_error[..., i], extent=extent,
                            aspect='auto', origin='lower')
            ax.set_xlabel('Red. scattering (1/cm)')
            ax.set_ylabel('Absorption (1/cm)')
            ax.set_title('SDS={:.1f} um, error={:.2f}%'.format(
                sds[i]*1e6, self._rel_error[...,i].mean()))
            fig.colorbar(pos, ax=ax)
        
        ax = axis.flat[-1]
        labels = []
        for d in sds:
            labels.append('{:.1f}'.format(d*1e6))
        ax.set_title('Rel. error (%) as f(SDS)')
        new_shape = (int(self._rel_error.size//nfib), nfib)
        ax.boxplot(np.reshape(self._rel_error, new_shape) , labels=labels, )
        ax.set_ylabel('Rel. error (%)')
        ax.set_xlabel('Source detector separation (um)')

    def passed(self) -> bool:
        ''''
        Returns True if the simulator passed the test, else False.
        '''
        if self.error is not None:
            return np.abs(self.error.mean()) < 0.5
        return False


class SingleLayerGkLutSixLinearProbe(McTest):
    def __init__(self, usepflut: bool = False, verbose: bool = False,
                 **kwargs):
        '''
        Test of the Monte Carlo simulator for a singel layer sample,
        3D GK gamma lookup table and a normally incident six-linear
        optical fiber probe.

        Parameters
        ----------
        usepflut: bool
            Run the simulator with lookup table-based sampling of the
            scattering phase function.
        verbose: bool
            Turn on/off verbose output.
        kwargs: dict
            Parameters passed to the Monte Carlo simulator constructor
            :py:meth:`xopto.mcml.mc.Mc`.
        '''

        self._verbose = verbose
        self._use_specular = kwargs.get('specular', False)

        datafile = os.path.join(DATA_DIR, 'single_layer_gk_lut_six_linear_probe.pkl')
        with open(datafile, 'rb') as fid:
            self._data = data = pickle.load(fid)

        if usepflut:
            pf = [mcpf.LutEx(Gk, data['layers'][0]['pf']['pfparams']),
                  mcpf.LutEx(Gk, data['layers'][1]['pf']['pfparams']),
                  mcpf.LutEx(Gk, data['layers'][2]['pf']['pfparams'])]
        else:
            pf = [mcpf.Gk(*data['layers'][0]['pf']['pfparams']),
                  mcpf.Gk(*data['layers'][1]['pf']['pfparams']),
                  mcpf.Gk(*data['layers'][2]['pf']['pfparams'])]

        layers = mclayer.Layers([
            mclayer.Layer(pf=pf[0], **data['layers'][0]['basic']),
            mclayer.Layer(pf=pf[1], **data['layers'][1]['basic']),
            mclayer.Layer(pf=pf[2], **data['layers'][2]['basic'])
        ])
        source_data = data['source']
        fiber = fiberutil.MultimodeFiber(**source_data.pop('fiber'))
        source = mcsource.UniformFiber(fiber, **source_data)

        surface_data = data['surface']
        fiber = fiberutil.MultimodeFiber(**surface_data.pop('fiber'))
        surface = mcsurface.SurfaceLayouts(
            top=mcsurface.LinearArray(fiber, **surface_data)
        )

        detector_data = data['detector']
        fiber = fiberutil.MultimodeFiber(**detector_data.pop('fiber'))
        detector = mcdetector.Detectors(
            top=mcdetector.LinearArray(fiber, **detector_data)
        )

        specular = None
        if self._use_specular:
            specular = mcdetector.Total()
        detectors = mcdetector.Detectors(top=detector, specular=specular)

        mc_obj = mc.Mc(layers, source, detectors, surface=surface,
                       **mc_options(kwargs))
        mc_obj.rmax = data['mc']['rmax']
        self._reflectance = None

        super().__init__(mc_obj, data['reflectance'],
                         'Single layer six-linear probe GG 3D mua musr lut')

    def run(self, nphotons: int or float = 10e6, **kwargs) -> float:
        '''
        Run the test.

        Parameters
        ----------
        nphotons: int
            The number of photon packets to use for each simulation run.
        kwargs: dict
            Parameters passed to the Monte carlo run method
            :py:meth:`xopto.mcml.mc.Mc.run`.

        Returns
        -------
        tmc: float
            Time (s) consumed by the simulator core.
        '''
        mua_vector = self._data['lut']['mua']
        musr_vector = self._data['lut']['musr']
        gamma_vector = self._data['lut']['gamma']
        g = self._data['g']
        pf_params_lut = self._data['lut'].get('pf_params_lut')
        gkmap = GkMap.fromfile()

        if pf_params_lut is None:
            for index, gamma in enumerate(gamma_vector):
                if self._verbose:
                    print('Computing GK parameters for (g={}, gamma={}), '
                        '{}/{}'.format(g, gamma, index + 1, gamma_vector.size),
                        end='\r')
                    if index == len(gamma_vector):
                        print()
                pf_params_lut[(g, gamma)] = gkmap.invgamma(g, gamma)

        mua, musr, gamma = np.meshgrid(
            mua_vector, musr_vector, gamma_vector, indexing='ij')

        N = mua.size

        t_start = time.perf_counter()

        reflectance =[]
        for i in range(N):
            pf_params = pf_params_lut.get((g, gamma.flat[i]))

            self.mc.layers[1].mua = mua.flat[i]
            self.mc.layers[1].mus = musr.flat[i]/(1.0 - g)
            self.mc.layers[1].pf.g = pf_params[0]
            self.mc.layers[1].pf.a = pf_params[1]
            _, _, detectors_res = self.mc.run(nphotons, **kwargs)
            reflectance.append(detectors_res.top.reflectance[1:])

            if self._verbose:
                dt = time.perf_counter() - t_start
                print(self.progress_str(i + 1, N, dt), end='\r')
                if i == N - 1:
                    print()
        t_mc = time.perf_counter() - t_start

        simulated = np.array(reflectance)
        simulated.shape = self.reference.shape
        self.simulated = simulated

        self._rel_error = (simulated - self.reference)/self.reference*100.0
        self.error = self._rel_error

        return t_mc

    def visualize(self):
        '''
        Visualize the results of the test.
        '''
        if self.simulated is None:
            return
        
        nfib = self._data['detector']['n'] - 1
        sds = np.arange(1, nfib + 1)*self._data['detector']['spacing']

        from matplotlib import pyplot as pp
        musr_vector = self._data['lut']['musr']
        mua_vector = self._data['lut']['mua']
        gamma_vector = self._data['lut']['gamma']

        title = 'Test of: {} - Rel. error(SDS) (mean={:.2f})%'.format(
            self.description, self._rel_error.mean()
        )

        _show_mua_musr_pf_param_grid(
            pf_param_name='gamma', grid=(mua_vector, musr_vector, gamma_vector),
            data=self._rel_error, sds=sds, show_mean=True, units='%',
            title='')

        return
        extent = (musr_vector[0]*1e-2, musr_vector[-1]*1e-2,
                  mua_vector[0]*1e-2, mua_vector[-1]*1e-2)
        
        nrows, ncols = int(np.ceil(nfib/3)), 3
        fig, axis = pp.subplots(nrows, ncols)
        pp.suptitle('Test of: {} - Rel. error(SDS) (mean={:.2f})%'.format(
            self.description, self._rel_error.mean()))

        for i in range(nfib):
            ax = axis.flat[i]
            pos = ax.imshow(self._rel_error[..., i].mean(-1), extent=extent,
                            aspect='auto', origin='lower')
            ax.set_xlabel('Red. scattering (1/cm)')
            ax.set_ylabel('Absorption (1/cm)')
            ax.set_title('SDS={:.1f} um, error={:.2f}%'.format(
                sds[i]*1e6, self._rel_error[...,i].mean()))
            fig.colorbar(pos, ax=ax)
        
        ax = axis.flat[-1]
        labels = []
        for d in sds:
            labels.append('{:.1f}'.format(d*1e6))
        ax.set_title('Rel. error (%) as f(SDS)')
        new_shape = (int(self._rel_error.size//nfib), nfib)
        ax.boxplot(np.reshape(self._rel_error, new_shape) , labels=labels, )
        ax.set_ylabel('Rel. error (%)')
        ax.set_xlabel('Source detector separation (um)')

    def passed(self) -> bool:
        ''''
        Returns True if the simulator passed the test, else False.
        '''
        if self.error is not None:
            return np.abs(self.error.mean()) < 0.5
        return False

def _append_report(test: McTest, reports: list = None) -> list:
    '''
    A private support function for collectig test reports.
    '''
    if reports is None:
        reports = []

    if test.passed():
        fmt_str = '{:s} PASSED!'
    else:
        fmt_str = '{:s} FAILED!'
    reports.append(fmt_str.format(test.__class__.__name__))

    return reports


if __name__ == '__main__':
    import argparse
    import os.path

    from xopto import USER_TMP_PATH    

    nphotons = 10e6
    cldevices = ['nvidia', 'amd', 'cpu', 'hd']

    parser = argparse.ArgumentParser(
        description='Input parameters',
        usage='python -m xopto.mcml.test.validate [<args>]\n')
    
    parser.add_argument(
        '-l', '--lut', action='store_true',
        required=False,
        help='Use lookup table-based sampling of the scattering angle.')

    parser.add_argument(
        '-n', '--nphotons', type=float, default=10e6,
        metavar='NUM_PHOTONS',
        help='Number of photons to use in MC simulations.')

    parser.add_argument(
        '-d', '--device', type=str, default='', metavar='OPENCL_DEVICE_NAME',
        required=False,
        help='String that uniquely identifies one OpenCL device. '\
             'Existing value in the configuration file is overwritten.')

    parser.add_argument(
        '-p', '--precision', type=str, default='single', metavar='FLOAT_PRECISION',
        required=False, choices=('single', 'double'),
        help='String that selects single or double floating-point precision.')

    parser.add_argument(
        '-c', '--counter_bits', type=int, default=32,
        metavar='PACKET_COPUNTER_BITS', required=False, choices=(32, 64),
        help='String that selects single or double floating-point precision.')

    parser.add_argument(
        '-b', '--ballistic', action='store_true', default=False, required=False,
        help='Enable ballistic implementation of the Monte Carlo kernel.')


    args = parser.parse_args()
    nphotons = max(0, int(args.nphotons))
    usepflut = bool(args.lut)
    custom_types = []
    ballistic = bool(args.ballistic)
    if args.device:
        cldevices = [args.device]
    if args.precision in ('double', 'd'):
        custom_types.append(mctypes.McDouble)
        cl_build_options = ['-cl-mad-enable']
    else:
        custom_types.append(mctypes.McSingle)
        cl_build_options = ['-cl-mad-enable', '-cl-fast-relaxed-math']
    if args.counter_bits == 64:
        custom_types.append(mctypes.McCnt64)
    else:
        custom_types.append(mctypes.McCnt32)

    options = {
        'nphotons': nphotons,
        'usepflut':usepflut,
        'wgsize':None,
        'verbose':True,
        'visualize':True,
        'cl_devices':clinfo.device(cldevices),
        'cl_build_options': cl_build_options,
        'types': mctypes.McDataTypes.custom(*custom_types),
        'options':[
            mcoptions.McUseNativeMath.on,
            mcoptions.McFloatLutMemory('constant'),
            mcoptions.McDebugMode.off,
            mcoptions.McUseBallisticKernel(ballistic)
        ],
        'specular': True, 
        'exportsrc': os.path.join(USER_TMP_PATH, 'mcml_validate.h')
    }

    reports = []

    test = SingleLayerLineSourceRadialProfile(**options)
    test.run(**run_options(options))
    if options['visualize']: test.visualize()
    if not test.passed(): print('SingleLayerLineSourceRadialProfile test FAILED!')
    reports = _append_report(test, reports)

    test = SingleLayerUniformFiberRadial(**options)
    test.run(**run_options(options))
    if options['visualize']: test.visualize()
    if not test.passed(): print('SingleLayerUniformFiberRadial test FAILED!')
    reports = _append_report(test, reports)

    test = DoubleLayerUniformFiberRadial(**options)
    test.run(**run_options(options))
    if options['visualize']: test.visualize()
    if not test.passed(): print('DoubleLayerUniformFiberRadial test FAILED!')
    reports = _append_report(test, reports)

    test = SingleLayerUniformFiberCartesian(**options)
    test.run(**run_options(options))
    if options['visualize']: test.visualize()
    if not test.passed(): print('SingleLayerUniformFiberCartesian test FAILED!')
    reports = _append_report(test, reports)

    test = SingleLayerSfdiIncidence30VariableNa(**options)
    test.run(**run_options(options))
    if options['visualize']: test.visualize()
    if not test.passed(): print('SingleLayerSfdiIncidence30VariableNa test FAILED!')
    reports = _append_report(test, reports)

    test = SingleLayerUniformFiberTrace(**options)
    trace_options = run_options(options)
    # reduce the number of photon packets for trace validation
    trace_options['nphotons'] = min(100000, trace_options.get('nphotons', 10000000))
    test.run(**trace_options)
    if options['visualize']: test.visualize()
    if not test.passed(): print('SingleLayerUniformFiberTrace test FAILED!')
    reports = _append_report(test, reports)

    test = SingleLayerSixLinearProbe(**options)
    test.run(**run_options(options))
    if options['visualize']: test.visualize()
    if not test.passed(): print('SingleLayerSixLinearProbe test FAILED!')
    reports = _append_report(test, reports)

    test = SingleLayerNineLinearProbe(**options)
    test.run(**run_options(options))
    if options['visualize']: test.visualize()
    if not test.passed(): print('SingleLayerNineLinearProbe test FAILED!')
    reports = _append_report(test, reports)

    test = SingleLayerGkLutSixLinearProbe(**options)
    test.run(**run_options(options))
    if options['visualize']: test.visualize()
    if not test.passed(): print('SingleLayerGkLutSixLinearProbe test FAILED!')
    reports = _append_report(test, reports)
    
    print('#'*60)
    print('Validation test summary:')
    print('\n    '.join(reports))

    if options['visualize']:
        from matplotlib import pyplot as pp
        pp.show()
