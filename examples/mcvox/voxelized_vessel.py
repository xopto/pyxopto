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

from xopto.mcvox import mc
from xopto.mcvox.mcutil import axis
from xopto.cl import clinfo

import numpy as np
import matplotlib.pyplot as pp

# DEFINE A COMPUTATIONAL DEVICE
# select the first device from the provided platform
cl_device = clinfo.gpu(platform='nvidia')

# DEFINE VOXELS AND TISSUE STRUCTURE
# number of bins and size
nx = 201
ny = 201
nz = 200
binsize = 5e-6

# define each axis
xaxis = axis.Axis(
    start=-nx/2 * binsize, 
    stop=nx/2 * binsize, 
    n=nx
)

yaxis = axis.Axis(
    start=-ny/2 * binsize, 
    stop=ny/2 * binsize, 
    n=ny
)

zaxis = axis.Axis(
    start=0.0, 
    stop=nz * binsize,
    n=nz
)

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
z0 = nz/2 * binsize

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
    n=1.33,
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

# DEFINE FLUENCE OBJECT
fluence = mc.mcfluence.Fluence(
    xaxis=xaxis, 
    yaxis=yaxis, 
    zaxis=zaxis,
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
nphotons = 1e8
deposit = mc_obj.run(nphotons, verbose=True)[1]

# VISUALIZE THE RESULTS
# plot the deposited energy in 1/m^3.
pp.figure()
pp.title('Deposited energy')
extent = [1e3*(xaxis.centers[0] - binsize/2), 1e3*(xaxis.centers[-1] - binsize/2),
    1e3*(zaxis.centers[-1] + binsize/2), 1e3*(zaxis.centers[0] - binsize/2)
]
pp.imshow(np.log10(deposit.data[:,100,:]/nphotons), 
    extent=extent,
    origin='upper'
)
pp.xlabel('x (mm)')
pp.ylabel('z (mm)')
cbar = pp.colorbar()
cbar.set_label('Log10(A)')
pp.show()

# use the slicer view
deposit.plot(axis='y')