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
import tempfile

import matplotlib as mpl
import matplotlib.pyplot as pp
from matplotlib.animation import FuncAnimation
import numpy as np

from xopto.mcml import mc
from xopto import make_user_dirs, USER_TMP_PATH


exportsrc = os.path.join(USER_TMP_PATH, 'mcml_fluencecyl_tr.h')
make_user_dirs()

cl_device = mc.clinfo.gpus()[0]

pf = mc.mcpf.Hg(0.8)

layers = mc.mclayer.Layers([
    mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=pf),
    mc.mclayer.Layer(d=1.0, n=1.0, mua=1.0e2, mus=50.0e2, pf=pf),
    mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=pf),
])

source = mc.mcsource.Line()

detectors = mc.mcdetector.Detectors(
    bottom=mc.mcdetector.Total()
)
fluence = mc.mcfluence.FluenceCylt(
    raxis=mc.mcfluence.Axis(0.0e-3, 1.0e-3, 101),
    fiaxis=mc.mcfluence.Axis(0.0e-3, 2.0*np.pi, 101),
    zaxis=mc.mcfluence.Axis( 0.0e-3, 1.0e-3, 101),
    taxis=mc.mcfluence.Axis( 0.0e-9, 0.01e-9, 200)
)

nphotons = 10e6

mc_obj = mc.Mc(layers, source, detectors, fluence=fluence,
               cl_devices=cl_device)
mc_obj.rmax = 2.0e-3

_, fluence_res, detectors_res = mc_obj.run(
    nphotons, exportsrc=exportsrc, verbose=True)

fluence_res.plot(axis='fi')

# VISUALIZE THE RESULTS
# plot the deposited energy in 1/m^3.
rz_proj = fluence_res.data.sum(axis=1)
rz_proj[rz_proj <= 0.0] = rz_proj[rz_proj > 0.0].min()
extent = [*fluence_res.raxis.span, *fluence_res.zaxis.span]
rz_proj_log = np.log10(rz_proj)
vmin, vmax =  np.log10(rz_proj[rz_proj > 0.0].min()), np.log10(rz_proj.max())

fig, ax = pp.subplots()
im = ax.imshow(rz_proj_log[..., 0], extent=extent, origin='lower',
               vmin=vmin, vmax=vmax)
colorbar = pp.colorbar(im)
colorbar.set_label(r'$\frac{1}{m^3s}$', rotation='horizontal')

title = ax.set_title('Deposited energy @ 0 ps')
ax.set_xlabel('r (mm)')
ax.set_ylabel('z (mm)')

pp.savefig(tempfile.TemporaryFile(), format='png')

ticks = colorbar.ax.get_yticklabels()
new_ticks = []
for tick in ticks:
    value = None
    if isinstance(tick, mpl.text.Text):
        value = tick.get_text()
        if value:
            value = float(value)
    elif isinstance(tick, (int, float)):
        value = float(tick)
    if isinstance(value, (int, float)):
        new_ticks.append('$10^{' + '{:.1f}'.format(value) + '}$')
colorbar.ax.set_yticklabels(new_ticks)

def init():
    return im, title

def update(frame):
    title.set_text('Deposited energy @ {:.1f} ps'.format(
        fluence_res.t[frame]*1e12))
    im.set_array(rz_proj_log[..., frame])
    return im, title

ani = FuncAnimation(fig, update, frames=np.arange(fluence_res.taxis.n),
                    init_func=init, blit=False)
pp.show()
