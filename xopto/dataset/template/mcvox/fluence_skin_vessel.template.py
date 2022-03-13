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

from xopto.mcvox import mc
from xopto.cl import clinfo

import numpy as np

# This example shows how to simulate deposited energy using the 
# voxelized Monte Carlo method for a simple case of a two-layered 
# skin model with an embedded blood vessel.

root_dataset_dir = '{{ root_dataset_dir }}'

cl_index = int(os.environ.get('CL_INDEX', 0))
cl_device = clinfo.device(
    os.environ.get('CL_DEVICE', {{ cl_name }}),
    index=cl_index
)

cl_build_options = {{ cl_build_options or None }}

overwrite = int(os.environ.get('MC_OVERWRITE', False))

num_packets = int(os.environ.get('MC_NUM_PACKETS', {{ num_packets }}))
batch_packets = int(os.environ.get('MC_BATCH_PACKETS', {{ batch_packets }}))

# DEFINE VOXELS AND TISSUE STRUCTURE
# number of bins and size
nx = {{ voxelization.nx }}
ny = {{ voxelization.ny }}
nz = {{ voxelization.nz }}
voxel_size = {{ voxelization.voxel_size }}

# define each axis
xaxis = mc.mcgeometry.Axis(
    start=-nx/2*voxel_size, 
    stop=nx/2*voxel_size, 
    n=nx
)

yaxis = mc.mcgeometry.Axis(
    start=-ny/2*voxel_size, 
    stop=ny/2*voxel_size, 
    n=ny
)

zaxis = mc.mcgeometry.Axis(
    start=0.0, 
    stop=nz*voxel_size,
    n=nz
)

# define voxels
voxels = mc.mcgeometry.Voxels(
    xaxis=xaxis, 
    yaxis=yaxis, 
    zaxis=zaxis
)

# DEFINE MATERIALS
# surrounding medium
material_surrounding = mc.mcmaterial.Material(
    n={{ materials.surrounding.n }},
    mua={{ materials.surrounding.mua }},
    mus={{ materials.surrounding.mus }},
    pf=mc.mcpf.Hg({{ materials.surrounding.g }})
)

# epidermis
material_epidermis = mc.mcmaterial.Material(
    n={{ materials.epidermis.n }},
    mua={{ materials.epidermis.mua }},
    mus={{ materials.epidermis.mus }},
    pf=mc.mcpf.Hg({{ materials.epidermis.g }})
)

# dermis
material_dermis = mc.mcmaterial.Material(
    n={{ materials.dermis.n }},
    mua={{ materials.dermis.mua }},
    mus={{ materials.dermis.mus }},
    pf=mc.mcpf.Hg({{ materials.dermis.g }})
)

# blood
material_blood = mc.mcmaterial.Material(
    n={{ materials.blood.n }},
    mua={{ materials.blood.mua }},
    mus={{ materials.blood.mus }},
    pf=mc.mcpf.Hg({{ materials.blood.g }})
)

materials = mc.mcmaterial.Materials([
    material_surrounding,
    material_epidermis,
    material_dermis,
    material_blood
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

detectors = mc.mcdetector.Detectors(
    top=mc.mcdetector.Cartesian(xaxis, yaxis),
    bottom=mc.mcdetector.Cartesian(xaxis, yaxis)
)

mcoptions = [getattr(mc.mcoptions.McMethod, '{{ method }}'),
             mc.mcoptions.McUseFluenceCache({{ cache }})]

# DEFINE MC OBJECT FOR MONTE CARLO SIMULATIONS AND ASSIGN MATERIALS
mc_obj = mc.Mc(
    voxels=voxels,
    materials=materials,
    detectors=detectors,
    fluence=fluence,
    source=source,
    options = mcoptions,
    cl_devices=cl_device,
    cl_build_options=cl_build_options
)
mc_obj.rmax = {{ rmax }}

# get mesh
z, y, x = voxels.meshgrid()

# compute logical indices for material assignment
r_vessel = {{ vessel_diameter }}*0.5
epidermis_mask = z <= {{ epidermis_thickness }}
dermis_mask = z > {{ epidermis_thickness }}

if batch_packets <= 0:
    batch_packets = num_packets

vessel_depth = {{ vessel_depth }}
for index, depth in enumerate(vessel_depth):
    print('Processing vessel position z={:.1f} um {}/{}'.format(
        depth*1e6, index + 1, vessel_depth.size ))

    vessel_mask = (x - 0.0)**2 + (z - depth)**2 <= r_vessel**2

    # assign the materials to voxels
    mc_obj.voxels.material.fill(0)
    mc_obj.voxels.material[epidermis_mask] = 1
    mc_obj.voxels.material[dermis_mask] = 2
    mc_obj.voxels.material[vessel_mask] = 3

    # voxels changed ... update of kernel voxel data required
    mc_obj.voxels.update()

    filename = '2-layer-skin-{:.0f}um-vessel-{:.0f}um-depth-deposition'.format(
        r_vessel*1e6*2.0, depth*1e6).replace('.', '_') + '.npz'
    output_dir = os.path.join(root_dataset_dir, 'mcvox', 'fluence')
    full_filename = os.path.join(output_dir, filename)

    if os.path.exists(full_filename) and not overwrite:
        print('Output file "{}" already exists and will not be '
              'overwritten!'.format(full_filename))
        continue

    # RUN THE MONTE CARLO SIMULATIONS
    packets_done = 0
    out = None
    while packets_done < num_packets:
        out = mc_obj.run(batch_packets, out)
        packets_done += batch_packets
        print('Processed {}/{}'.format(
            int(packets_done), int(num_packets)), end='\r')
    print()

    _, fluence_res, detectors_res = out

    data = {
        'num_packets': num_packets,
        'batch_packets': batch_packets,
        'voxels': voxels.todict(),
        'voxel_data': voxels.material,
        'materials': materials.todict(),
        'source': source.todict(),
        'detectors': detectors.todict(),
        'fluence': fluence.todict(),
        'fluence_data': fluence_res.data,
        'reflectance': detectors_res.top.reflectance,
        'transmittance': detectors_res.bottom.transmittance,
        'method': '{{ method }}'
    }

    os.makedirs(output_dir, exist_ok=True)
    np.savez_compressed(full_filename, **data)
