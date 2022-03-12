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

from xopto.mcml import mc
from xopto.mcml.mcutil.fiber import MultimodeFiber
from xopto import pf
from xopto.cl import clinfo
from xopto.util import hankel, fourier

import numpy as np

{% macro render_args_kwargs(args=[], kwargs={}) -%}
{%- for item in args %}{{ item }}{% if not loop.last or kwargs %}, {% endif %}{% endfor -%}
{%- for key, value in kwargs.items() %}{{ key }}={{ value }}{% if not loop.last %}, {% endif %}{% endfor -%}
{%- endmacro -%}

{%- macro render_detector_axis(args=[], kwargs={}) -%}
mc.mcdetector.Axis({{ render_args_kwargs(args, kwargs) }})
{%- endmacro -%}

{%- macro render_detector_symmetricaxis(args=[], kwargs={}) -%}
mc.mcdetector.SymmetricAxis({{ render_args_kwargs(args, kwargs) }})
{%- endmacro -%}

{%- macro render_fluence_axis(args=[], kwargs={}) -%}
mc.mcfluence.Axis({{ render_args_kwargs(args, kwargs) }})
{%- endmacro -%}

{%- macro render_fluence(fluence) -%}
mc.mcfluence.Fluence(
    xaxis={{ render_fluence_axis(flence.xaxis.args, flence.xaxis.kwargs) }},
    yaxis={{ render_fluence_axis(flence.yaxis.args, flence.yaxis.kwargs) }},
    zaxis={{ render_fluence_axis(flence.zaxis.args, flence.zaxis.kwargs) }},
    {{ render_args_kwargs(fluence.args, fluence.kwargs) }}
)
{%- endmacro -%}

{%- macro render_multimode_fiber(args=[], kwargs={}) -%}
mcutil.fiber.MultimodeFiber({{ render_args_kwargs(args, kwargs) }})
{%- endmacro -%}

storage_dir_common = os.path.join(
    '{{ root_dataset_dir }}',
    'mcml', '{{ sample.dir }}', '{{ source.dir }}'
)

cl_index = int(os.environ.get('CL_INDEX', {{ cl_index }}))
cl_device = clinfo.device(
    os.environ.get('CL_DEVICE', {{ cl_device }}),
    index=cl_index
)

cl_build_options = {{ cl_build_options or None }}

overwrite = int(os.environ.get('MC_OVERWRITE', False))

num_packets = int(os.environ.get('MC_NUM_PACKETS', {{ num_packets }}))

layers = mc.mclayer.Layers([
    {%- for layer in sample.layers %}
    mc.mclayer.Layer(
        d={{ layer.d }}, n={% if loop.first %}{{source.n_above }}{% elif loop.last %}{{ source.n_bellow }}{% else %}{{ layer.n }}{% endif %}, mua={{ layer.mua }}, mus={{ layer.mus }}, pf=mc.mcpf.LutEx(pf.Mie, {{ pf.default }}, lutsize={{ pf.lutsize }})),
    {%- endfor %}
])

surface = soure = top_layout = top_detector = None

{% if source.type.startswith('Line') and not source.type.startswith('LinearArray') -%}
{%- if source.is_sfdi -%}
source = mc.mcsource.Line({{ render_args_kwargs(source.args, source.kwargs) }})
{%- if source.detector == 'radial_sfdi' %}
top_detector = mc.mcdetector.Radial(
    {{ render_detector_axis(detector.radial_sfdi.top.raxis.args, detector.radial_sfdi.top.raxis.kwargs) }},
    {{ render_args_kwargs(detector.radial_sfdi.top.args, detector.radial_sfdi.top.kwargs) }}
)
detector_dir = '{{ detector.radial_sfdi.dir }}'
{%- elif source.detector == 'symmetricx_sfdi' %}
top_detector = mc.mcdetector.SymmetricX(
    {{ render_detector_symmetricaxis(detector.symmetricx_sfdi.top.xaxis.args, detector.symmetricx_sfdi.top.xaxis.kwargs) }},
    {{ render_args_kwargs(detector.symmetricx_sfdi.top.args, detector.symmetricx_sfdi.top.kwargs) }}
)
detector_dir = '{{ detector.symmetricx_sfdi.dir }}'
{%- else %}
raise ValueError('Unsupported SFDI detector!')
{%- endif %}
{%- else %}
source = mc.mcsource.Line({{ render_args_kwargs(source.args, source.kwargs) }})
top_detector = mc.mcdetector.Radial(
    {{ render_detector_axis(detector.radial.top.raxis.args, detector.radial.top.raxis.kwargs) }},
    {{ render_args_kwargs(detector.radial.top.args, detector.radial.top.kwargs) }}
)
detector_dir = 'radial'
{%- endif %}
{%- elif source.type.startswith('UniformBeam') -%}
source = mc.mcsource.UniformBeam({{ render_args_kwargs(source.args, source.kwargs) }})
top_detector = mc.mcdetector.Radial(
    {{ render_detector_axis(detector.radial.top.raxis.args, detector.radial.top.raxis.kwargs) }},
    {{ render_args_kwargs(detector.radial.top.args, detector.radial.top.kwargs) }}
)
detector_dir = 'radial'
{%- elif source.type.startswith('GaussianBeam') -%}
source = mc.mcsource.GaussianBeam({{ render_args_kwargs(source.args, source.kwargs) }})
top_detector = mc.mcdetector.Radial(
    {{ render_detector_axis(detector.radial.top.raxis.args, detector.radial.top.raxis.kwargs) }},
    {{ render_args_kwargs(detector.radial.top.args, detector.radial.top.kwargs) }}
)
detector_dir = 'radial'
{%- elif source.type.startswith('UniformFiber') -%}
fiber = MultimodeFiber({{ render_args_kwargs(source.fiber.args, source.fiber.kwargs) }})
source = mc.mcsource.UniformFiber(fiber, {{ render_args_kwargs(source.args, source.kwargs) }})
top_detector = mc.mcdetector.Radial(
    {{ render_detector_axis(detector.radial.top.raxis.args, detector.radial.top.raxis.kwargs) }},
    {% if detector.radial.top.args or detector.radial.top.kwargs %}
    {{ render_args_kwargs(detector.radial.top.args, detector.radial.top.kwargs) }},
    {%- endif -%}
    cosmin=(1.0 - (fiber.na/fiber.ncore)**2)**0.5
)
detector_dir = 'radial'
{%- elif source.type.startswith('SixAroundOne') -%}
fiber = MultimodeFiber({{ render_args_kwargs(source.fiber.args, source.fiber.kwargs) }})
top_detector = mc.mcdetector.SixAroundOne(fiber, spacing={{ source.kwargs.spacing }})
top_layout = mc.mcsurface.SixAroundOne(
    fiber, spacing={{ source.kwargs.spacing }}, reflectivity={{ source.kwargs.reflectivity }}, diameter={{ source.kwargs.diameter }})
source = mc.mcsource.UniformFiber(
    fiber, position=[*top_layout.fiber_position(0), 0.0])
detector_dir = 'probe'
{%- elif source.type.startswith('LinearArray') -%}
fiber = MultimodeFiber({{ render_args_kwargs(source.fiber.args, source.fiber.kwargs) }})
top_layout = mc.mcsurface.LinearArray(
    fiber, n={{ source.kwargs.n }}, spacing={{ source.kwargs.spacing }}, reflectivity={{ source.kwargs.reflectivity}}, diameter={{ source.kwargs.diameter }})
top_detector = mc.mcdetector.LinearArray(fiber, n={{ source.kwargs.n }}, spacing={{ source.kwargs.spacing }})
source = mc.mcsource.UniformFiber(
    fiber, position=[*top_layout.fiber_position(0), 0.0])
detector_dir = 'probe'
{%- else -%}
raise TypeError('Unsupported source!')
{%- endif %}

{% if source.is_sfdi -%}
{%- if source.detector == 'radial_sfdi' -%}
bottom_detector = mc.mcdetector.Radial(
    {{ render_detector_axis(detector.radial_sfdi.bottom.raxis.args, detector.radial_sfdi.bottom.raxis.kwargs) }},
    {{ render_args_kwargs(detector.radial_sfdi.bottom.args, detector.radial_sfdi.bottom.kwargs) }}
)
{%- elif source.detector == 'symmetricx_sfdi' -%}
bottom_detector = mc.mcdetector.SymmetricX(
    {{ render_detector_symmetricaxis(detector.symmetricx_sfdi.bottom.raxis.args, detector.symmetricx_sfdi.bottom.raxis.kwargs) }},
    {{ render_args_kwargs(detector.symmetricx_sfdi.bottom.args, detector.symmetricx_sfdi.bottom.kwargs) }}
)
{%- else -%}
raise('Unsupported SFDI detector!')
{%- endif %}
{%- else %}
bottom_detector = mc.mcdetector.Radial(
    {{ render_detector_axis(detector.radial.bottom.raxis.args, detector.radial.bottom.raxis.kwargs) }},
    {{ render_args_kwargs(detector.radial.bottom.args, detector.radial.bottom.kwargs) }}
)
{%- endif %}
detectors = mc.mcdetector.Detectors(
    top=top_detector, bottom=bottom_detector,
    specular=mc.mcdetector.Total() 
)

if top_layout is not None:
    surface = mc.mcsurface.SurfaceLayouts(top = top_layout)

mcoptions = [mc.mcoptions.McFloatLutMemory.constant_mem,
             getattr(mc.mcoptions.McMethod, '{{ method }}')]

mc_obj = mc.Mc(layers, source, detectors, surface=surface,
               cl_devices=cl_device, cl_build_options=cl_build_options,
               options=mcoptions)
mc_obj.rmax = {{ rmax }}

mua_values = {{ sample.mua }}
musr_values = {{ sample.musr }}

{% if source.is_sfdi -%}
spatial_frequencies = {{ spatial_frequencies }}

def sfd_transform(frequencies: np.ndarray, r_x: np.ndarray,
                  reflectance: np.ndarray) -> np.ndarray:
    {%- if source.detector == 'symmetricx_sfdi' %}
    return fourier.discrete_simpson(frequencies, r_x, reflectance, uneven=True)
    {%- else %}
    return hankel.discrete_simpson(frequencies, r_x, reflectance, uneven=True)
    {%- endif %}
{%- endif %}

def process_mua_musr(storage_dir : str, g: float, pf_info: str,
                     num_packets: int, pos: int = 0,
                     overwrite: bool = False) -> int:
    for mua in {{ sample.mua }}:
        for musr in {{ sample.musr }}:
            pos += 1
            print('{:>5d}/{:<5d}: {:s} | '
                  'g1={:.4f} | mua={:.2f} 1/cm | musr={:.2f} 1/cm '.format(
                pos, N, pf_info, g1, mua*1e-2, musr*1e-2), end='\r')

            filename = 'mua-{:.2f}-musr-{:.2f}-invcm'.format(mua*1e-2, musr*1e-2)
            full_filename = os.path.join(
                storage_dir, filename.replace('.', '_') + '.npz')
            if not overwrite and os.path.exists(full_filename):
                print('Skipping simulation for existing file "{}"\n'.format(
                    full_filename))
                continue

            layers[1].mua = mua
            layers[1].mus = musr/(1.0 - g)
            _, _, detectors_res = mc_obj.run(num_packets)
            data = {
                'layers': mc_obj.layers.todict(),
                'source': mc_obj.source.todict(),
                'surface': mc_obj.surface.todict() \
                    if mc_obj.surface is not None else None,
                'detectors': mc_obj.detectors.todict(),
                'reflectance': detectors_res.top.reflectance,
                'transmittance': detectors_res.bottom.transmittance,
                'total_reflectance': detectors_res.top.raw.sum()/ \
                    max(detectors_res.top.nphotons, 1),
                'total_transmittance': detectors_res.bottom.raw.sum()/ \
                    max(detectors_res.top.nphotons, 1),
                'method': '{{ method }}'
            }
            {%- if source.is_sfdi %}
            data['sfd_reflectance'] = sfd_transform(
                spatial_frequencies,
                {% if source.detector == 'symmetricx_sfdi' %}detectors_res.top.x{% else %}detectors_res.top.r{% endif %},
                detectors_res.top.reflectance)
            data['sfd_frequency'] = spatial_frequencies
            {%- endif %}
            np.savez_compressed(full_filename, **data)
    return pos

wavelengths = {{ pf.wavelengths }}
diameters = {{ pf.diameters }}
nparticle = {{ pf.nparticle }}
nmedium = {{ pf.nmedium }}
nfibercore = {{ pf.nfibercore }}
probe_reflectivity = {{ pf.probe_reflectivity }}

N = {{ pf.wavelengths | length }}*{{ pf.diameters | length }}*{{ sample.n_mua }}*{{ sample.n_musr  }}
pos = 0
for diameter in diameters:
    for wavelength_index, wavelength in enumerate(wavelengths):
        pf_args = (float(nparticle[wavelength_index]),
                   float(nmedium[wavelength_index]),
                   float(diameter), float(wavelength))
        g1 = pf.Mie(*pf_args).g(1)
        layers[1].pf = mc.mcpf.LutEx(pf.Mie, pf_args, lutsize={{ pf.lutsize }})
        layers[1].n = pf_args[1]
        {% if source.type.startswith('UniformFiber') %}
        mc_obj.source.fiber.ncore = nfibercore[wavelength_index]
        mc_obj.detectors.top.cosmin = \
            (1.0 - (mc_obj.source.fiber.na/mc_obj.source.fiber.ncore)**2)**0.5
        {% elif source.type.startswith('LinearArray') or source.type.startswith('SixAroundOne') %}
        mc_obj.source.fiber.ncore = nfibercore[wavelength_index]
        mc_obj.surface.top.fiber.ncore = nfibercore[wavelength_index]
        mc_obj.surface.top.reflectivity = probe_reflectivity[wavelength_index]
        mc_obj.detectors.top.fiber.n = nfibercore[wavelength_index]
        {% endif %}
        storage_dir = os.path.join(
            storage_dir_common, detector_dir, '{{ pf.dir }}',
            '{{ pf.param_dir[0] }}'.format(diameter*1e6).replace('.', '_'),
            '{{ pf.param_dir[1] }}'.format(wavelength*1e9).replace('.', '_')
        )
        os.makedirs(storage_dir, exist_ok=True)
        pf_info = 'Mie({:.3f}, {:.3f}, {:.2e}, {:.2e})'.format(*pf_args)
        pos = process_mua_musr(
            storage_dir, g1, pf_info, num_packets, pos, overwrite)

print()
