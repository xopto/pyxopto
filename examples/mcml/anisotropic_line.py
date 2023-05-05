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


debug = False

exportsrc = os.path.join(USER_TMP_PATH, 'mcml_anisotropic.h')
make_user_dirs()

# g can be a scalar value, a vector of length 3 or a full tensor of size 3x3
pf = mc.mcpf.Hga([0.9, 0.01, 0.8])
cl_device = mc.clinfo.gpus()[0]

isotropic_layers = mc.mclayer.Layers([
    mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=pf),
    mc.mclayer.Layer(d=np.inf, n=1.33, mua=0.1e2, mus=100e2, pf=pf),
    mc.mclayer.Layer(d=0.0, n=1.0, mua=1e2, mus=0.0, pf=pf),
])

# mus, mua can be a vector of length 3 or a full tensor matrix of size 3x3
anisotropic_layers = mc.mclayer.Layers([
    mc.mclayer.AnisotropicLayer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=pf),
    mc.mclayer.AnisotropicLayer(d=np.inf, n=1.33, mua=0.1e2, mus=50e2, pf=pf),
    mc.mclayer.AnisotropicLayer(d=0.0, n=1.0, mua=1e2, mus=0.0, pf=pf),
])

layers = anisotropic_layers

theta = np.deg2rad(0)
source = mc.mcsource.Line(
    (0.0, 0.0, 0.0),
    (np.sin(theta), 0.0, np.cos(theta))
)

fluence = mc.mcfluence.Fluence(
    mc.mcfluence.Axis(-2.5e-3, 2.5e-3, 101),
    mc.mcfluence.Axis(-2.5e-3, 2.5e-3, 101),
    mc.mcfluence.Axis( 0.0e-3, 5.0e-3, 101)
)

detectors = mc.mcdetector.Detectors(
    top=mc.mcdetector.Cartesian(
        xaxis=mc.mcdetector.Axis(-2.5e-3, 2.5e-3, 101),
        yaxis=mc.mcdetector.Axis(-2.5e-3, 2.5e-3, 101)
    ),
    specular=mc.mcdetector.Total()
)

mc_obj = mc.Mc(layers, source, detectors, fluence=fluence, cl_devices=cl_device,
           options=[mc.mcoptions.McDebugMode(debug),
                    ('MC_USE_EVENTS', 1)])
mc_obj.rmax = 25e-3

nphotons = 1e6 if not debug else 1
_, fluence_res, detectors_res = mc_obj.run(
    nphotons, exportsrc=exportsrc, verbose=True)

specular = float(detectors_res.specular.reflectance)

#pp.figure()
#pp.title('Top surface reflectance (specular {:.1f}%)'.format(specular*100.0))
#pp.semilogy(detectors_res.top.r*1e3, detectors_res.top.reflectance)
#pp.xlabel('r (mm)')
detectors_res.top.plot(show=False)
pp.title('Reflectance')

fluence_res.plot(axis='xy', show=False)

fluence_res.plot()
