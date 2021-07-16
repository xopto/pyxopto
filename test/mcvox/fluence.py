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

import numpy as np

from xopto.mcvox import mc
from xopto.cl import clinfo
from xopto import make_user_dirs, USER_TMP_PATH
from xopto.mcvox.mcutil.fiber import MultimodeFiber


debug = False

exportsrc = os.path.join(USER_TMP_PATH, 'mcvox.h')
make_user_dirs()

nphotons = 1e6
r_term = 100e-3
dz = 0.02e-3
dr = 0.1e-3
nz = 500
nr = 200
n = 1.33
mua = 10e2
mus = 125e2
d = 1.5e-3
g = 0.96

######################
r_vessels = 0.200e-3
z0 = 0.25e-3
y0 = 0.0e-3


cl_device = clinfo.gpus()[0]

# define the scattering phase function Henyey-Greenstein
# "mcpf" is a submodule assigned for handling OpenCL structures
# that are passed to the Monte Carlo simulator
# a variety of phase functions are currently possible:
# Henyey-Greenstein, Gegenbauer kernel or Reynolds-McCromick,
# modified Henyey-Greenstein, modified product of cosines, product of
# cosines, Rayleigh scattering and lut (lookup-table based). The latter
# offers arbitrary implementation of a scattering phase function
# over the polar angle theta, while the azimuth angle is always
# sampled over the 2pi interval
pf = mc.mcpf.Hg(0.7)

#source = mc.mcsource.Line()
source = mc.mcsource.UniformBeam(
    50e-6,
    position=(-0.5e-3, 0.0, 0.5e-3),
    direction=(1.0, 0.0, 0.0)
)
source = mc.mcsource.GaussianBeam(
    [25e-6, 5e-6],
    position=(-0.5e-3, 0.0, 0.5e-3),
    direction=(1.0, 0.0, 0.0)
)

source = mc.mcsource.UniformFiber(
    MultimodeFiber(dcore=200e-6, dcladding=220e-6, ncore=1.452, na=0.22),
    position=(-0.5e-3, 0.0, 0.5e-3),
    direction=(1.0, 0.0, 0.0)
)

xaxis = mc.mcfluence.Axis(-0.5e-3, 0.5e-3, 200)
yaxis = mc.mcfluence.Axis(-0.5e-3, 0.5e-3, 200)
zaxis = mc.mcfluence.Axis( 0.0e-3, 1.0e-3, 200)

# define vessel

voxels = mc.mcgeometry.Voxels(
    xaxis=xaxis, yaxis=yaxis, zaxis=zaxis
)

z, y, x = voxels.meshgrid()
logical_inds = (z - 0.5e-3)**2 + (y - 0.0)**2 <= 0.25e-3**2

materials = mc.mcmaterial.Materials([
    mc.mcmaterial.Material(n=1.33, mua=0.0e2, mus=50e2, pf=pf),
    mc.mcmaterial.Material(n=1.33, mua=0.5e2, mus=50e2, pf=pf)
])

fluence = mc.mcfluence.Fluence(
    xaxis=xaxis, yaxis=yaxis, zaxis=zaxis, mode='deposition'
)

mc_obj = mc.Mc(voxels=voxels, materials=materials, fluence=fluence,
    source=source, cl_devices=cl_device,
    options=[mc.mcoptions.McDebugMode(debug)],
    #cl_build_options=['-cl-opt-disable', '-cl-strict-aliasing']
)

mc_obj.rmax = r_term
mc_obj.voxels.material[logical_inds] = 1

nphotons = 10e6 if not debug else 1
fluence, detector = mc_obj.run(nphotons, verbose=True, exportsrc=exportsrc)[1:]

#import matplotlib.pyplot as pp
#pp.imshow(fluence.data.mean(axis=-1))
#pp.show()
fluence.plot(axis='x', show=False)
fluence.plot(axis='xproj')
