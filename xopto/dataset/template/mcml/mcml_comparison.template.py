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
import time

import numpy as np
from xopto.mcml import mc
from xopto.cl import clinfo

STORAGE_ROOT = '{{ root_dataset_dir }}'

cl_index = int(os.environ.get('CL_INDEX', {{ cl_index }}))
cl_device = clinfo.device(
    os.environ.get('CL_DEVICE', {{ cl_name }}),
    index=cl_index
)

cl_build_options = {{ cl_build_options or None }}

overwrite = os.environ.get('MC_OVERWRITE', False) 

num_packets = int(os.environ.get('MC_NUM_PACKETS', {{ num_packets }}))

PF_DIR_FORMAT = lambda g: 'g-{:.2f}'.format(g).replace('.', '_')
FILENAME_FORMAT = lambda mua, musr: 'mua-{:.2f}-musr-{:.2f}-invcm'.format(
    mua*1e-2, musr*1e-2).replace('.', '_') + '.npz'

mua_values = {{ sample.mua }}
musr_values = {{ sample.musr }}
g_values = {{ sample.g }}

layers = mc.mclayer.Layers([
    {%- for layer in sample.layers %}
    mc.mclayer.Layer(d={{ layer.d }}, n={{ layer.n }},
                     mua={{ layer.mua }}, mus={{ layer.mus }}, pf=mc.mcpf.Hg({{ layer.g }})),
    {%- endfor %}
])

source = mc.mcsource.Line()

detectors = mc.mcdetector.Detectors(
    top=mc.mcdetector.Radial(
        mc.mcdetector.Axis(0.0, 5.0e-3, 500)
    ),
    specular = mc.mcdetector.Total(),
    bottom=mc.mcdetector.Radial(
        mc.mcdetector.Axis(0.0, 5.0e-3, 500)
    )
)

fluence = mc.mcfluence.FluenceRz(
    mc.mcfluence.Axis(0.0, 5.0e-3, 500),
    mc.mcfluence.Axis(0.0, 5.0e-3, 500),
    mode='deposition',
)
{% if fluence_k -%}fluence.k = {{ fluence_k }}{%- endif %}

mcoptions = [getattr(mc.mcoptions.McMethod, '{{ method }}')]

mc_obj = mc.Mc(layers, source, detectors, fluence=fluence,
               cl_build_options=cl_build_options, options=mcoptions,
               cl_devices=cl_device{% if double_precision %}, types=mc.mctypes.McDataTypesDouble{% endif %})
mc_obj.rmax = 1000e-3
pos = 0
N = len(mua_values)*len(musr_values)*len(g_values)
progress_template = \
    '{:>4d}/{:<4d} mua={:5.2f} 1/cm; musr={:5.2f} 1/cm; pf=Hg({:.2f})'

for g in g_values:
    mc_obj.layers[1].pf.g = g
    pf_param_dir = PF_DIR_FORMAT(g)
    save_dir = os.path.join(STORAGE_ROOT,
                            'mcml_comparison', '{{ sample.dir }}',
                            'line', 'radial', 'hg', pf_param_dir)
    os.makedirs(save_dir, exist_ok=True)

    print(mc_obj.layers)
    print(mc_obj.source)
    print(mc_obj.detectors)
    print(mc_obj.fluence)
    print()

    for mua in mua_values:
        for musr in musr_values:
            pos += 1
            report = progress_template.format(pos, N, mua*1e-2, musr*1e-2, g)
            print(report, end='\r')

            filename = FILENAME_FORMAT(mua, musr)
            full_filename = os.path.join(save_dir, filename)
            if os.path.exists(filename) and not overwrite:
                print('Skipping simulation for existing file "{}"!'.format(
                    full_filename))
                continue

            mc_obj.layers[1].mua = mua
            mc_obj.layers[1].mus = musr/(1.0 - g)
            if pos == 1:
                # build the OpenCL source here - don't want to include in timing
                _ = mc_obj.run(1)
            t_start = time.perf_counter()
            _, res_fluence, res_detectors = mc_obj.run(num_packets)
            run_time = time.perf_counter() - t_start
            data = {
                'rmax': float(mc_obj.rmax),
                'num_packets': int(num_packets),
                'layers': mc_obj.layers.todict(),
                'source': mc_obj.source.todict(),
                'detectors': mc_obj.detectors.todict(),
                'fluence': res_fluence.todict(),
                'fluence_data': res_fluence.data,
                'reflectance': res_detectors.top.reflectance,
                'total_reflectance': \
                    res_detectors.top.raw.sum()/res_detectors.top.nphotons,
                'specular': res_detectors.specular.reflectance,
                'transmittance': res_detectors.bottom.reflectance,
                'total_transmittance': \
                    res_detectors.bottom.raw.sum()/res_detectors.bottom.nphotons,
                'run_time': run_time,
                'cl_device': mc_obj.cl_device[0].name,
                'cl_index': {{ cl_index }},
                'method': '{{ method }}'
            }
            np.savez_compressed(full_filename, **data)

print()
