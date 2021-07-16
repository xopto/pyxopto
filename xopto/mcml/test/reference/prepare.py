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

import os.path
import pickle

import numpy as np
import scipy.io

from xopto import PICKLE_PROTOCOL

ROOT_DIR = os.path.dirname(__file__)
from pyopto.pf.util import GkMap

def single_layer_line_source_radial_profile():
    out = {}
    data = scipy.io.loadmat(
        os.path.join(ROOT_DIR, 'old_format', 'Rr_CUDAMCML.mat')
    )

    # sample layer
    n1 = 1.33
    mua1 = 0.01*1e2
    musr1 = 60*1e2
    g1 = 0.85
    d1 = 8e-3

    # surrounding medium
    ntop = nbottom = 1.33

    # source

    # detector
    nr = 340
    dr = 5e-6
    cos_min = 0.0

    # reflectance
    reflectance_profile = np.squeeze(data['Rr'])*1e4 # units are 1/cm^2

    # mc simulation radius
    rmax = 10e-3


    out['layers'] = [
        {'basic': {'n': ntop, 'd': float('inf'), 'mua':0.0, 'mus': 0.0 },
         'pf': {'type': 'Hg', 'pfparams': (g1,)}},
        {'basic': {'n': n1, 'd': d1, 'mua':mua1, 'mus':musr1/(1.0 - g1)},
         'pf': {'type': 'Hg', 'pfparams': (g1,)}},
        {'basic': {'n': nbottom, 'd': float('inf'), 'mua': 0.0, 'mus': 0.0},
         'pf': {'type': 'Hg', 'pfparams': (g1,)}}
    ]

    out['source'] = {}
    out['detector'] = {
        'axis':{'start': 0.0, 'stop': nr*dr, 'n':nr}, 
        'cosmin': cos_min
    }
    out['reflectance'] = reflectance_profile
    out['mc'] = {'rmax': rmax}

    filename = os.path.join(ROOT_DIR, 'cuda_singlelayer_linesource_radial')
    with open(filename + '.pkl', 'wb') as fid:
        pickle.dump(out, fid, protocol=PICKLE_PROTOCOL)


def single_layer_uniform_fiber_lut():
    out = {}
    data = scipy.io.loadmat(
        os.path.join(ROOT_DIR, 'old_format', 'CUDA_reference.mat')
    )

    # sample layer
    n1 = float(data['CUDA_lut'][0, 0]['layer1_n'])
    d1 = float(data['CUDA_lut'][0, 0]['layer1_d'])*1e-2
    g1 = float(data['CUDA_lut'][0, 0]['layer1_g'])
    muA1 = data['CUDA_lut'][0, 0]['layer1_mua']*1e2
    muSr1 = data['CUDA_lut'][0, 0]['layer1_musr']*1e2

    # surrounding medium
    ntop = float(data['CUDA_lut'][0, 0]['layer_top'])
    nbottom = float(data['CUDA_lut'][0, 0]['layer_bottom'])

    # source
    fib_na = float(data['CUDA_lut'][0, 0]['fiber_NA'])
    fib_sds = data['CUDA_lut'][0, 0]['fiber_SDS'].flatten()*1e-2
    fib_d = 2*float(data['CUDA_lut'][0, 0]['fiber_radius'])*1e-2

    # detector
    nr = 340
    dr = 5e-6
    cos_min = float(np.sqrt(1.0 - (fib_na/ntop)**2))

    # reflectance
    R1 = data['CUDA_lut'][0, 0]['reflectance_fiber1']
    R2 = data['CUDA_lut'][0, 0]['reflectance_fiber2']
    R3 = data['CUDA_lut'][0, 0]['reflectance_fiber3']
    R4 = data['CUDA_lut'][0, 0]['reflectance_fiber4']
    R5 = data['CUDA_lut'][0, 0]['reflectance_fiber5']
    Rref = np.dstack([R1, R2, R3, R4, R5])

    # mc simulation radius
    rmax = float(data['CUDA_lut'][0, 0]['termination_radius']*1e-2)

    out['layers'] = [
        {'basic': {'n': ntop, 'd': float('inf'), 'mua':0.0, 'mus': 0.0 },
         'pf': {'type': 'Hg', 'pfparams': (g1,)}},
        {'basic': {'n': n1, 'd': d1, 'mua':0.0, 'mus':0.0},
         'pf': {'type': 'Hg', 'pfparams': (g1,)}},
        {'basic': {'n': nbottom, 'd': float('inf'), 'mua': 0.0, 'mus': 0.0},
         'pf': {'type': 'Hg', 'pfparams': (g1,)}}
    ]
    out['source'] = {'na': fib_na, 'dcore': fib_d, 'dcladding': fib_d, 'ncore': ntop}
    out['detector'] = {
        'axis':{'start': 0.0, 'stop': nr*dr, 'n':nr}, 
        'cosmin': cos_min,
        'fibers':{'sds': fib_sds.tolist(), 'dcore': fib_d}
    }
    mua = np.sort(np.unique(muA1))
    musr = np.sort(np.unique(muSr1))
    out['lut'] = {'mua': np.asarray(mua), 'musr': np.asarray(musr)}
    out['reflectance'] = np.asarray(Rref)
    out['mc'] = {'rmax': rmax}

    filename = os.path.join(ROOT_DIR, 'cuda_singlelayer_uniformfiber')
    with open(filename + '.pkl', 'wb') as fid:
        pickle.dump(out, fid, protocol=PICKLE_PROTOCOL)


def double_layer_uniform_fiber_lut():
    out = {}
    data = scipy.io.loadmat(
        os.path.join(ROOT_DIR, 'old_format', 'CUDA_reference2.mat')
    )

    # sample layer 1
    n1 = float(data['CUDA_lut'][0, 0]['layer1_n'])
    d1 = float(data['CUDA_lut'][0, 0]['layer1_d'])*1e-2
    g1 = float(data['CUDA_lut'][0, 0]['layer1_g'])
    muA1 = data['CUDA_lut'][0, 0]['layer1_mua']*1e2
    muSr1 = data['CUDA_lut'][0, 0]['layer1_musr']*1e2

    # sample layer 2
    n2 = float(data['CUDA_lut'][0, 0]['layer2_n'])
    d2 = float(data['CUDA_lut'][0, 0]['layer2_d'])*1e-2
    g2 = float(data['CUDA_lut'][0, 0]['layer2_g'])
    muA2 = data['CUDA_lut'][0, 0]['layer2_mua']*1e2
    muSr2 = data['CUDA_lut'][0, 0]['layer2_musr']*1e2

    # sourrounding medium
    ntop = float(data['CUDA_lut'][0, 0]['layer_top'])
    nbottom = float(data['CUDA_lut'][0, 0]['layer_bottom'])

    # optical fiber source
    fib_na = float(data['CUDA_lut'][0, 0]['fiber_NA'])
    fib_sds = data['CUDA_lut'][0, 0]['fiber_SDS'].flatten()*1e-2
    fib_d = 2*float(data['CUDA_lut'][0, 0]['fiber_radius'])*1e-2

    # detector
    nr = 340
    dr = 5e-6
    cos_min = float(np.sqrt(1.0 - (fib_na/ntop)**2))

    # reflectance
    R1 = data['CUDA_lut'][0, 0]['reflectance_fiber1']
    R2 = data['CUDA_lut'][0, 0]['reflectance_fiber2']
    R3 = data['CUDA_lut'][0, 0]['reflectance_fiber3']
    R4 = data['CUDA_lut'][0, 0]['reflectance_fiber4']
    R5 = data['CUDA_lut'][0, 0]['reflectance_fiber5']
    Rref = np.dstack([R1, R2, R3, R4, R5])

    # mc simulation radius
    rmax = float(data['CUDA_lut'][0, 0]['termination_radius']*1e-2)

    out['layers'] = [
        {'basic': {'n': ntop, 'd': float('inf'), 'mua':0.0, 'mus': 0.0 },
         'pf': {'type': 'Hg', 'pfparams': (g1,)}},
        {'basic': {'n': n1, 'd': d1, 'mua':0.0, 'mus':0.0},
         'pf': {'type': 'Hg', 'pfparams': (g1,)}},
        {'basic': {'n': n2, 'd': d2, 'mua':0.0, 'mus':0.0},
         'pf': {'type': 'Hg', 'pfparams': (g2,)}},
        {'basic': {'n': nbottom, 'd': float('inf'), 'mua': 0.0, 'mus': 0.0},
         'pf': {'type': 'Hg', 'pfparams': (g1,)}}
    ]
    out['source'] = {'na': fib_na, 'dcore': fib_d, 'dcladding': fib_d, 'ncore': ntop}
    out['detector'] = {
        'axis':{'start': 0.0, 'stop': nr*dr, 'n':nr}, 
        'cosmin': cos_min,
        'fibers':{'sds': fib_sds.tolist(), 'dcore': fib_d}
    }
    mua1 = np.sort(np.unique(muA1))
    musr1 = np.sort(np.unique(muSr1))
    mua2 = np.sort(np.unique(muA2))
    musr2 = np.sort(np.unique(muSr2))
    out['lut'] = [{'mua': mua1, 'musr': musr1},
                  {'mua': mua2, 'musr': musr2}]
    out['reflectance'] = np.asarray(Rref)
    out['mc'] = {'rmax': rmax}

    filename = os.path.join(ROOT_DIR, 'cuda_doublelayer_uniformfiber')
    with open(filename + '.pkl', 'wb') as fid:
        pickle.dump(out, fid, protocol=PICKLE_PROTOCOL)


def single_layer_sfdi_30deg_incidence_variable_na():
    data = data = np.load(
        os.path.join(
            ROOT_DIR, 'old_format',
            'single_layer_sfdi_30deg_incidence_multi_na_ilm.npz'
        )
    )

    # sample layer
    n1 = float(data['n'])
    g1 = float(data['g'])
    mua1 = float(data['mua'])
    musr1 = float(data['musr'])
    d1 = float('inf')

    # source
    incidence_angle = float(data['incidence_angle']) # radians
    
    #detector
    NA = np.sin(data['acceptance_angle']) # this is a vector

    # sourrounding medium
    ntop = 1.0
    nbottom = 1.0

    # spatial frequencies 1/m
    frequencies = np.array(data['frequencies'])

    # amplitude and phase
    amplitude = np.array(data['amplitude'].T)
    phase = np.array(data['phase'].T)

    # simulation radius
    rmax = 150e-3

    out = {}

    out['layers'] = [
        {'basic': {'n': ntop, 'd': float('inf'), 'mua':0.0, 'mus': 0.0 },
         'pf': {'type': 'Hg', 'pfparams': (g1,)}},
        {'basic': {'n': n1, 'd': d1, 'mua':mua1, 'mus':musr1/(1.0 - g1)},
         'pf': {'type': 'Hg', 'pfparams': (g1,)}},
        {'basic': {'n': nbottom, 'd': float('inf'), 'mua': 0.0, 'mus': 0.0},
         'pf': {'type': 'Hg', 'pfparams': (g1,)}}
    ]
    out['source'] = {
        'direction': (np.sin(incidence_angle), 0.0, np.cos(incidence_angle))
    }
    out['detector'] = {
        'axis':{'center': 0.0, 'range': 150e-3, 'n_half':4000, 'logscale': True}
    }
    out['frequencies'] = frequencies
    
    out['lut'] = {'na': NA}

    out['amplitude'] = amplitude
    out['phase'] = phase

    out['mc'] = {'rmax': rmax}

    filename = os.path.join(
        ROOT_DIR,
        'single_layer_sfdi_30deg_incidence_variable_na'
    )
    with open(filename + '.pkl', 'wb') as fid:
        pickle.dump(out, fid, protocol=PICKLE_PROTOCOL)


def single_layer_uniform_fiber_trace():
    # the sample layer
    g1 = 0.85
    n1 = 1.33

    # the source fiber
    na = 0.22
    dcore = 200e-6
    dcladding = 220e-6
    ncore = 1.452
    d1 = 8e-3

    # detector
    nr = 170
    dr = 10e-6
    cosmin = np.sqrt(1.0 - (na/ncore)**2)

    # the surrounding medium
    ntop = nbottom = ncore

    # simulation radius
    rmax = 10e-3

    out = {}
    out['layers'] = [
        {'basic': {'n': ntop, 'd': float('inf'), 'mua':0.0, 'mus': 0.0 },
        'pf': {'type': 'Hg', 'pfparams': (g1,)}},
        {'basic': {'n': n1, 'd': d1, 'mua':0.0, 'mus':0.0},
        'pf': {'type': 'Hg', 'pfparams': (g1,)}},
        {'basic': {'n': nbottom, 'd': float('inf'), 'mua': 0.0, 'mus': 0.0},
        'pf': {'type': 'Hg', 'pfparams': (g1,)}}
    ]
    out['source'] = {
        'na': na, 'dcore': dcore, 'dcladding': dcladding, 'ncore': ncore
    }
    out['detector'] = {
        'axis': {'start': 0.0, 'stop': nr*dr, 'n': nr},
        'cosmin': cosmin
    }
    out['trace'] = {'maxlen': 100}
    out['lut'] = {'mua': np.array([0.0, 0.0, 20.0, 20.0, 10.0])*1e2,
                  'musr': np.array([0.0, 60.0, 60.0, 0.0, 30.0])*1e2}
    out['mc'] = {'rmax': rmax}

    filename = os.path.join(
        ROOT_DIR,
        'single_layer_uniformfiber_trace'
    )
    with open(filename + '.pkl', 'wb') as fid:
        pickle.dump(out, fid, protocol=PICKLE_PROTOCOL)

def single_layer_hg_lut_6_linear_probe():
    out = {}

    out['layers'] = [
        {'basic': {'n': 1.0, 'd': float('inf'), 'mua':0.0, 'mus': 0.0},
         'pf': {'type': 'Hg', 'pfparams': (0.8,)}},
        {'basic': {'n': 1.33, 'd': 8e-3, 'mua': 0.0, 'mus': 0.0},
         'pf': {'type': 'Hg', 'pfparams': (0.8,)}},
        {'basic': {'n': 1.33, 'd': float('inf'), 'mua': 0.0, 'mus': 0.0},
         'pf': {'type': 'Hg', 'pfparams': (0.8,)}}
    ]

    fiber = {'dcore': 200e-6, 'dcladding': 220e-6, 'ncore': 1.452, 'na': 0.22}

    out['source'] = {'fiber': fiber, 'position': (-220e-6*2.5, 0.0, 0.0)}
    out['surface'] = {
        'fiber': fiber, 'n': 6, 'spacing': 220e-6,
        'cutout': (220e-6*6, 220e-6), 'cutoutn':1.6,
        'diameter': 6.1e-3, 'reflectivity': 0.43
    }
    out['detector'] = {'fiber': fiber, 'n':6, 'spacing': 220e-6}

    mua = np.linspace(0.0, 25e2, 35)
    musr = np.linspace(10e2, 100e2, 35)
    out['lut'] = {'mua': mua, 'musr': musr}

    data = np.load(
        os.path.join(
            ROOT_DIR, 'old_format',
            'single_layer_hg_lut_six_linear_probe_ex_data.npz'
        )
    )
    out['reflectance'] = data['reflectance']
    out['mc'] = {'rmax': 5e-3}

    filename = os.path.join(
        ROOT_DIR, 'single_layer_hg_lut_six_linear_probe')
    with open(filename + '.pkl', 'wb') as fid:
        pickle.dump(out, fid, protocol=PICKLE_PROTOCOL)

def single_layer_hg_lut_9_linear_probe():
    out = {}

    out['layers'] = [
        {'basic': {'n': 1.0, 'd': float('inf'), 'mua':0.0, 'mus': 0.0},
         'pf': {'type': 'Hg', 'pfparams': (0.8,)}},
        {'basic': {'n': 1.33, 'd': 8e-3, 'mua': 0.0, 'mus': 0.0},
         'pf': {'type': 'Hg', 'pfparams': (0.8,)}},
        {'basic': {'n': 1.33, 'd': float('inf'), 'mua': 0.0, 'mus': 0.0},
         'pf': {'type': 'Hg', 'pfparams': (0.8,)}}
    ]

    fiber = {'dcore': 200e-6, 'dcladding': 220e-6, 'ncore': 1.452, 'na': 0.22}

    out['source'] = {'fiber': fiber, 'position': (-220e-6*4, 0.0, 0.0)}
    out['surface'] = {
        'fiber': fiber, 'n': 9, 'spacing': 220e-6,
        'cutout': (220e-6*9, 220e-6), 'cutoutn':1.6,
        'diameter': 6.1e-3, 'reflectivity': 0.43
    }
    out['detector'] = {'fiber': fiber, 'n':9, 'spacing': 220e-6}

    mua = np.linspace(0.0, 25e2, 35)
    musr = np.linspace(10e2, 100e2, 35)
    out['lut'] = {'mua': mua, 'musr': musr}

    data = np.load(
        os.path.join(
            ROOT_DIR, 'old_format',
            'single_layer_hg_lut_nine_linear_probe_ex_data.npz'
        )
    )
    out['reflectance'] = data['reflectance']
    out['mc'] = {'rmax': 5e-3}

    filename = os.path.join(
        ROOT_DIR, 'single_layer_hg_lut_nine_linear_probe')
    with open(filename + '.pkl', 'wb') as fid:
        pickle.dump(out, fid, protocol=PICKLE_PROTOCOL)


def single_layer_gk_lut_6_linear_probe():
    out = {}

    out['layers'] = [
        {'basic': {'n': 1.0, 'd': float('inf'), 'mua':0.0, 'mus': 0.0},
         'pf': {'type': 'Gk', 'pfparams': (0.8, 0.5)}},
        {'basic': {'n': 1.33, 'd': 8e-3, 'mua': 0.0, 'mus': 0.0},
         'pf': {'type': 'Gk', 'pfparams': (0.8, 0.5)}},
        {'basic': {'n': 1.33, 'd': float('inf'), 'mua': 0.0, 'mus': 0.0},
         'pf': {'type': 'Gk', 'pfparams': (0.8, 0.5)}}
    ]
    out['g'] = 0.8
    source_fiber = {
        'dcore': 200e-6, 'dcladding': 220e-6, 'ncore': 1.452, 'na': 0.18}
    detector_fiber = {
        'dcore': 200e-6, 'dcladding': 220e-6, 'ncore': 1.452, 'na': 0.22}

    out['source'] = {'fiber': source_fiber, 'position': (-220e-6*2.5, 0.0, 0.0)}
    out['surface'] = {
        'fiber': detector_fiber, 'n': 6, 'spacing': 220e-6,
        'cutout': (220e-6*6, 220e-6), 'cutoutn':1.6,
        'diameter': 6.0e-3, 'reflectivity': 0.57
    }
    out['detector'] = {'fiber': detector_fiber, 'n':6, 'spacing': 220e-6}

    mua = np.linspace(0.0, 25e2, 30)
    musr = np.linspace(5e2, 70e2, 30)
    gamma = np.linspace(1.6, 2.3, 20)
    out['lut'] = {'mua': mua, 'musr': musr, 'gamma': gamma}

    gkmap = GkMap.fromfile()

    pf_params_lut = {}
    g = float(out['g'])
    print()
    for index, gamma_value in enumerate(gamma):
        gamma_value = float(gamma_value)
        print('Computing GK parameters for (g={}, gamma={}), '
                '{}/{}'.format(g, gamma_value, index + 1, gamma.size),
                end='\r')
        if index == len(gamma) - 1:
            print()
        pf_params_lut[(g, gamma_value)] = gkmap.invgamma(out['g'], gamma_value)

    out['lut']['pf_params_lut'] = pf_params_lut

    data = np.load(
        os.path.join(
            ROOT_DIR, 'old_format',
            'single_layer_gk_gamma_lut_six_linear_probe_ex_data.npz'
        )
    )
    reference = data['reflectance']
    reference.shape = (mua.size, musr.size, gamma.size, 5)
    out['reflectance'] = reference
    out['mc'] = {'rmax': 5e-3}

    filename = os.path.join(
        ROOT_DIR, 'single_layer_gk_lut_six_linear_probe')
    with open(filename + '.pkl', 'wb') as fid:
        pickle.dump(out, fid, protocol=PICKLE_PROTOCOL)

if __name__ == '__main__':
    single_layer_uniform_fiber_lut()
    double_layer_uniform_fiber_lut()
    single_layer_line_source_radial_profile()
    single_layer_sfdi_30deg_incidence_variable_na()
    single_layer_uniform_fiber_trace()
    single_layer_hg_lut_6_linear_probe()
    single_layer_hg_lut_9_linear_probe()
    single_layer_gk_lut_6_linear_probe()
