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

{% macro render_args_kwargs(args=[], kwargs={}) -%}
{%- for item in args %}{{ item }}{% if not loop.last or kwargs %}, {% endif %}{% endfor -%}
{%- for key, value in kwargs.items() %}{{ key }}={{ value }}{% if not loop.last %}, {% endif %}{% endfor -%}
{%- endmacro -%}

import os
import os.path

import numpy as np

from xopto.cl import clinfo
from xopto.mcml import mc
from xopto.mcml.mcutil.fiber import MultimodeFiber

root_dataset_dir = '{{ root_dataset_dir }}'

cl_index = int(os.environ.get('CL_INDEX', {{ cl_index }}))
cl_device = clinfo.device(
    os.environ.get('CL_DEVICE', {{ cl_device }}),
    index=cl_index
)

cl_build_options = {{ cl_build_options or None }}

overwrite = int(os.environ.get('MC_OVERWRITE', False))

sv_num_packets = int(os.environ.get('MC_NUM_PACKETS', {{ sv_num_packets }}))
batch_packets = int(os.environ.get('MC_BATCH_PACKETS', {{ batch_packets }}))

layers = mc.mclayer.Layers([
    {%- for layer in sample.layers %}
    mc.mclayer.Layer(d={{ layer.d }}, mua={{ layer.mua }}, mus={{ layer.mus }}, n={{ layer.n }},
                     pf=mc.mcpf.Hg({{ layer.g }})),
    {%- endfor %}
])
layers[0].n = {{ source.n_above }}
layers[-1].n = {{ source.n_bellow }}

fiber = MultimodeFiber(
    {{ render_args_kwargs(source.fiber.args, source.fiber.kwargs) }})

top_layout = mc.mcsurface.LinearArray(
    fiber, n=2, spacing={{ source.kwargs.spacing }},
    reflectivity={{ source.kwargs.reflectivity}},
    diameter={{ source.kwargs.diameter }})
source = mc.mcsource.UniformFiber(
    fiber, position=(*top_layout.fiber_position(0), 0.0))

fp_eps = np.finfo(np.float32).eps
fiber_cosmin = (1.0 - (source.fiber.na/source.fiber.ncore)**2)**0.5
trace = mc.mctrace.Trace(
    maxlen={{ trace.kwargs.maxlen }},
    filter=mc.mctrace.Filter(
        z=(-fp_eps, fp_eps),
        pz=(-1.0, -fiber_cosmin),
        r=(0.0, {{ source.fiber.kwargs.dcore }}*0.5, (top_layout.fiber_position(1)))
    )
)

surface = mc.mcsurface.SurfaceLayouts(
    top=top_layout
)

mcoptions = [getattr(mc.mcoptions.McMethod, {{ method }})]

mc_obj = mc.Mc(layers, source, surface=surface, trace=trace,
               cl_devices=cl_device, cl_build_options=cl_build_options,
               options=mcoptions)
mc_obj.rmax = {{ rmax }}

sv = mc.mcsv.SamplingVolume(
    xaxis=mc.mcsv.Axis({{ render_args_kwargs(sv.xaxis.args, sv.xaxis.kwargs) }}),
    yaxis=mc.mcsv.Axis({{ render_args_kwargs(sv.yaxis.args, sv.yaxis.kwargs) }}),
    zaxis=mc.mcsv.Axis({{ render_args_kwargs(sv.zaxis.args, sv.zaxis.kwargs) }})
)

sample_name = '{{ sample.dir }}'
source_name = '{{ source.dir }}'

sds = 'sds-{:.0f}um'.format({{ source.kwargs.spacing }}*1e6).replace('.', '_')
hg_g = 'g-{:.2f}'.format(layers[1].pf.g).replace('.', '_')
filename = 'mua-{:.2f}-musr-{:.2f}-invcm'.format(
    layers[1].mua*1e-2,
    layers[1].mus*(1.0 - layers[1].pf.g)*1e-2
).replace('.', '_')

output_dir = os.path.join(
    root_dataset_dir, 'mcml', sample_name, 'sv','reflectance',
    source_name, sds, 'hg', hg_g)
filename_full = os.path.join(output_dir, filename + '.npz')

if os.path.exists(filename_full) and not overwrite:
    print('Output file "{}" already exists and will '
          'not be overwritten!'.format(filename_full))
else:
    print('Sampling volume for reflectance configuration of two optical fibers:')
    print('    Fiber: {:.0f} um | SDS: {:.0f} um'.format(
        {{source.fiber.kwargs.dcore}}*1e6, {{source.kwargs.spacing}}*1e6))
    print('    Hg({:.2f}) | mua={:.2f} 1/cm | musr={:.2f} 1/cm'.format(
        layers[1].pf.g, layers[1].mua*1e-2,
        layers[1].mus*1e-2*(1.0 - layers[1].pf.g)))

    reached_detector = 0
    sv_res = sv
    while reached_detector <= sv_num_packets:
        trace_res, _, _ = mc_obj.run(batch_packets)
        sv_res = mc_obj.sampling_volume(trace_res, sv_res)
        reached_detector += len(trace_res)
        print('    Collected paths: {}/{}'.format(
            reached_detector, sv_num_packets), end='\r')
    print()

    data = {
        'rmax': {{ rmax }},
        'batch_packets': batch_packets,
        'sv_num_packets': sv_num_packets,
        'layers': layers.todict(),
        'source': source.todict(),
        'sv': sv_res.todict(),
        'sv_data': sv_res.data,
        'method': {{ method }}
    }

    os.makedirs(output_dir, exist_ok=True)
    np.savez_compressed(filename_full, **data)
