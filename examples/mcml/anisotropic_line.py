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

import matplotlib.pyplot as pp
import numpy as np

from xopto.mcml import mc
from xopto import make_user_dirs, USER_TMP_PATH

# used for exporting the MCML kernel code for inspection/debugging
exportsrc = os.path.join(USER_TMP_PATH, 'mcml_anisotropic.h')
make_user_dirs()

# Parameter g of Hga is internally a 3x3 tensor but can be initialized or set by:
# (a) a scalar value (applied to all the 3 diagonal elements of the tensor,
#                     all the non diagonal elements are set to 0)
# (b) a vector of length 3 (applied to the 3 diagonal elements of the tensor,
#                           all the non diagonal elements are set to 0))
# (c) a full tensor of shape 3x3
pf = mc.mcpf.Hga([0.9, 0.01, 0.8])
cl_device = mc.clinfo.gpus()[0]

# Example of an isotropic layer stack (regular pyxopto). Note that the
# isotropic layers can also take the new anisotropic scattering phase function
# Hga but the values of mua and mus are scalar.
isotropic_layers = mc.mclayer.Layers([
    mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=pf),
    mc.mclayer.Layer(d=np.inf, n=1.33, mua=0.1e2, mus=100e2, pf=pf),
    mc.mclayer.Layer(d=0.0, n=1.0, mua=1e2, mus=0.0, pf=pf),
])

# Example of an anisotropic layer stack.
# Parameters mua and mus are internally 3x3 tensors but can be initializer or set by:
# (a) a scalar value (applied to all the 3 diagonal elements of the tensor,
#                     all the non diagonal elements are set to 0)
# (b) a vector of length 3 (applied to the 3 diagonal elements of the tensor,
#                           all the non diagonal elements are set to 0))
# (c) a full tensor of shape 3x3
anisotropic_layers = mc.mclayer.Layers([
    mc.mclayer.AnisotropicLayer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=pf),
    mc.mclayer.AnisotropicLayer(d=np.inf, n=1.33, mua=0.1e2, mus=50e2, pf=pf),
    mc.mclayer.AnisotropicLayer(d=0.0, n=1.0, mua=1e2, mus=0.0, pf=pf),
])

# picks one of the two layer stacks for simulations 
layers = anisotropic_layers

# The anisotropi MCML supports all the sources of the isotropic MCML.
theta = np.deg2rad(0)
source = mc.mcsource.Line(
    (0.0, 0.0, 0.0),
    (np.sin(theta), 0.0, np.cos(theta))
)

# The anisotropic MCML supports the same fluence objects as the isotropic MCML.
fluence = mc.mcfluence.Fluence(
    mc.mcfluence.Axis(-2.5e-3, 2.5e-3, 101),
    mc.mcfluence.Axis(-2.5e-3, 2.5e-3, 101),
    mc.mcfluence.Axis( 0.0e-3, 5.0e-3, 101)
)

# The anisotropic MCML supports the same detector objects as the isotropic MCML.
detectors = mc.mcdetector.Detectors(
    top=mc.mcdetector.Cartesian(
        xaxis=mc.mcdetector.Axis(-2.5e-3, 2.5e-3, 101),
        yaxis=mc.mcdetector.Axis(-2.5e-3, 2.5e-3, 101)
    ),
    specular=mc.mcdetector.Total()
)
# The simulator interface of the anisotropic MCML is the same as of the isotropic MCML.
mc_obj = mc.Mc(layers, source, detectors, fluence=fluence,
               cl_devices=cl_device)
mc_obj.rmax = 25e-3

nphotons = 1e6
_, fluence_res, detectors_res = mc_obj.run(
    nphotons, exportsrc=exportsrc, verbose=True)

specular = float(detectors_res.specular.reflectance)

detectors_res.top.plot(show=False)
pp.title('Reflectance')

fluence_res.plot(axis='xy', show=False)

fluence_res.plot()
