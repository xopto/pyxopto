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

from xopto.mcvox import mc
from xopto import make_user_dirs, USER_TMP_PATH


debug = False

exportsrc = os.path.join(USER_TMP_PATH, 'mcvox_isotropicpoint.h')
make_user_dirs()

pf = mc.mcpf.Hg(0.8)
cl_device = mc.clinfo.gpus()[0]

materials = mc.mcmaterial.Materials([
    mc.mcmaterial.Material(n=1.0, mua=0.0, mus=0.0, pf=pf),
    mc.mcmaterial.Material(n=1.33, mua=1e2, mus=100e2, pf=pf),
])

voxels = mc.mcgeometry.Voxels(
    mc.mcfluence.Axis(-2.5e-3, 2.5e-3, 100),
    mc.mcfluence.Axis(-2.5e-3, 2.5e-3, 100),
    mc.mcfluence.Axis( 0.0e-3, 5.0e-3, 100)
)

theta = np.deg2rad(45)
source = mc.mcsource.Line(
    (0.0, 0.0, 0.0),
    (np.sin(theta), 0.0, np.cos(theta))
)

fluence = mc.mcfluence.Fluence(
    mc.mcfluence.Axis(-2.5e-3, 2.5e-3, 100),
    mc.mcfluence.Axis(-2.5e-3, 2.5e-3, 100),
    mc.mcfluence.Axis( 0.0e-3, 5.0e-3, 100)
)

detectors = mc.mcdetector.Detectors(
    top=mc.mcdetector.Radial(
        mc.mcdetector.Axis(0e-3, 5e-3, 1000)
    ),
    specular=mc.mcdetector.Total()
)

mc_obj = mc.Mc(voxels, materials, source, detectors, fluence=fluence,
           cl_devices=cl_device,
           options=[mc.mcoptions.McDebugMode(debug)])
mc_obj.rmax = 25e-3

mc_obj.voxels.material[:] = 1

nphotons = 1e6 if not debug else 1
_, fluence_res, detectors_res = mc_obj.run(
    nphotons, exportsrc=exportsrc, verbose=True)

specular = float(detectors_res.specular.reflectance)

pp.figure()
pp.title('Top surface reflectance (specular {:.1f}%)'.format(specular*100.0))
pp.semilogy(detectors_res.top.r*1e3, detectors_res.top.reflectance)
pp.xlabel('r (mm)')

fluence_res.plot(axis='xy', show=False)

fluence_res.plot()
