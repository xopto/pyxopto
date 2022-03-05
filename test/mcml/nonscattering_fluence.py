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
from xopto.mcbase import mcoptions

from xopto.mcml import mc
from xopto import make_user_dirs, USER_TMP_PATH


debug = False

exportsrc = os.path.join(USER_TMP_PATH, 'mcml_nonscattering_fluence.h')
make_user_dirs()

pf = mc.mcpf.Hg(0.8)
cl_device = mc.clinfo.gpus()[0]
layers = mc.mclayer.Layers([
    mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=pf),
    mc.mclayer.Layer(d=1.0, n=1.0, mua=4.0, mus=0.0, pf=pf),
    mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=pf),
])

source = mc.mcsource.Line()

detectors = mc.mcdetector.Detectors(
    bottom=mc.mcdetector.Total()
)
fluence = mc.mcfluence.Fluence(
    mc.mcfluence.Axis(-0.5e-3, 0.5e-3, 1),
    mc.mcfluence.Axis(-0.5e-3, 0.5e-3, 1),
    mc.mcfluence.Axis( 0.0e-3, 1.0,    1000)
)

nphotons = 500e6

meth = mc.mcoptions.McMethod
for method, name in ((meth.ar, 'Albedo Rejection'), (meth.aw, 'Albedo Weight'),
                     (meth.mbl, 'Microscopic Beer-Lambert')):
    print('\nProcessing with MC method:', name)
    print('-'*40)
    mc_obj = mc.Mc(layers, source, detectors, fluence=fluence,
                cl_devices=cl_device, options=[method])
    mc_obj.rmax = 2.0

    _, fluence_res, detectors_res = mc_obj.run(
        nphotons, exportsrc=exportsrc, verbose=True)
    mc_transmittance = detectors_res.bottom.transmittance

    analytical_transmittance = np.exp(-layers[1].mua*layers[1].d)
    print('-'*40)
    print('Transmittance (MC, analytical): ({:.6f}, {:.6f})'.format(
        float(mc_transmittance), float(analytical_transmittance)))

    analytical_deposition = -np.diff(np.exp(-layers[1].mua*fluence_res.zaxis.edges))
    mc_deposition = (fluence_res.raw/nphotons).flat

    fig, (ax1, ax2) = pp.subplots(2, 1)
    pp.suptitle(name)
    ax1.set_title('Energy deposition in nonscattering medium')
    ax1.plot(fluence_res.z, mc_deposition, 'r-', label="MC")
    ax1.plot(fluence_res.z, analytical_deposition, 'g--', label='Analytical')
    ax1.legend()
    ax2.set_title('Relative error (%)')
    ax2.plot(fluence_res.z,
            (mc_deposition - analytical_deposition)/analytical_deposition*100.0)
    ax2.set_xlabel('z (mm)')
    pp.tight_layout()

pp.show()
