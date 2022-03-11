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

# This example shows how to simulate deposited energy using the 
# voxelized Monte Carlo method for a simple case of a two-layered 
# skin model with an embedded blood vessel.

import time
import tempfile

from xopto.mcvox import mc
from xopto.cl import clinfo

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pp
from matplotlib.animation import FuncAnimation

# DEFINE A COMPUTATIONAL DEVICE
# select the first available GPU device
cl_device = clinfo.gpu()

num_packets = 5e6

# DEFINE VOXELS AND TISSUE STRUCTURE

# define each axis
xaxis = mc.mcgeometry.Axis(start=-0.5e-3, stop=0.5e-3, n=201)
yaxis = mc.mcgeometry.Axis(start=-0.5e-3, stop=0.5e-3, n=201)
zaxis = mc.mcgeometry.Axis(start=0.0, stop=1.0e-3, n=201)
taxis = mc.mcgeometry.Axis(start=0.0, stop=0.005e-9, n=50)

# define voxels
voxels = mc.mcgeometry.Voxels(
    xaxis=xaxis, 
    yaxis=yaxis, 
    zaxis=zaxis
)

# get mesh
z, y, x = voxels.meshgrid()

# set logical indices for material assignment
r_vessels = 100e-6
x0 = 0.0
z0 = zaxis.n*zaxis.step*0.5

logical_epi = z <= 60e-6
logical_der = z > 60e-6
logical_vessel = (x-x0)**2 + (z-z0)**2 <= r_vessels**2

# DEFINE MATERIALS
# surrounding medium
material_sm = mc.mcmaterial.Material(
    n=1.33,
    mua=0.0001e2,
    mus=1.0e2,
    pf=mc.mcpf.Hg(1.0)
)

# epidermis
material_epi = mc.mcmaterial.Material(
    n=1.33,
    mua=16.5724e2,
    mus=375.9398e2,
    pf=mc.mcpf.Hg(0.9)
)

# dermis
material_der = mc.mcmaterial.Material(
    n=1.33,
    mua=0.4585e2,
    mus=356.5406e2,
    pf=mc.mcpf.Hg(0.9)
)

# blood
material_bl = mc.mcmaterial.Material(
    n=1.0, # 1.33,
    mua=230.5427e2,
    mus=93.9850e2,
    pf=mc.mcpf.Hg(0.9)
)

materials = mc.mcmaterial.Materials([
    material_sm,
    material_epi,
    material_der,
    material_bl
])

# DEFINE PENCIL BEAM SOURCE
source = mc.mcsource.Line()
#source = mc.mcsource.GaussianBeam(sigma=10e-6)

# DEFINE FLUENCE OBJECT
fluence = mc.mcfluence.Fluencet(
    xaxis=xaxis, 
    yaxis=yaxis, 
    zaxis=zaxis,
    taxis=taxis,
    mode='deposition'
)

# DEFINE MC OBJECT FOR MONTE CARLO SIMULATIONS AND ASSIGN MATERIALS
mc_obj = mc.Mc(
    voxels=voxels,
    materials=materials,
    fluence=fluence,
    source=source,
    cl_devices=cl_device,
)

# assign the materials to voxels
mc_obj.voxels.material[logical_epi] = 1
mc_obj.voxels.material[logical_der] = 2
mc_obj.voxels.material[logical_vessel] = 3

# RUN THE MONTE CARLO SIMULATIONS
_, fluence_res, _ = mc_obj.run(num_packets, verbose=True)

# VISUALIZE THE RESULTS
# plot the deposited energy in 1/m^3.
xz_proj = fluence_res.data.sum(axis=1)
xz_proj[xz_proj <= 0.0] = xz_proj[xz_proj > 0.0].min()
extent = [*fluence_res.xaxis.span, *fluence_res.zaxis.span]
xz_proj_log = np.log10(xz_proj)

fig, ax = pp.subplots()
im = ax.imshow(xz_proj_log[..., 0], extent=extent, origin='upper')
colorbar = pp.colorbar(im)
colorbar.set_label(r'$\frac{1}{m^3s}$', rotation='horizontal')

title = ax.set_title('Deposited energy @ 0 ps')
ax.set_xlabel('x (mm)')
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
    im.set_array(xz_proj_log[..., frame])
    return im, title

ani = FuncAnimation(fig, update, frames=np.arange(fluence_res.taxis.n),
                    init_func=init, blit=False)
pp.show()
