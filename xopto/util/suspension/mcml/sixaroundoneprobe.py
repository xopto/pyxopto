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

from typing import Callable, List
import time

from xopto.mcml import mc
from xopto.materials import ri
from xopto.mcml.mcutil import fiber

from xopto.util.suspension import Suspension


class SixAroundOneProbe:
    def __init__(
            self,
            suspension: Suspension,
            surrounding_ri: None or float or Callable([float, float], float) = None,
            fiber_core_ri: None or float or Callable[[float, float], float] = None,
            fiber_core_d: float = 400e-6,
            fiber_cladding_d: float = 420e-6,
            fiber_na: float or Callable[[float, float], float] = 0.22,
            fiber_spacing: float or None = None,
            src_fiber: int  = 0,
            probe_diameter: float = 6.1e-3,
            probe_reflectivity: float or Callable[[float], float] = 0.6,
            probe_cutout: bool = True,
            epoxy_fill_ri: None or float or Callable[[float, float], float] = None,
            cl_build_options: None or List[str or mc.cloptions.ClBuildOption] = None,
            options: None or List[mc.mcoptions.McOption] = None,
            cl_device: str or mc.cl.Device or List[mc.cl.Device] or
                       mc.cl.Context or mc.cl.CommandQueue = None):
        '''
        Parameters
        ----------
        suspension: Suspension
            Suspension instance use with the Monte Carlo simulations.
        surrounding_ri:  float or Callable[[float, float], float]
            Refractive index of the medium above the probe tip surface.
            A callable that takes wavelength and temperature or a fixed
            floating-point value that is independent of wavelength and
            temperature.
            Default implementation uses the refractive index of the suspension
            medium (method medium_ri).
        fiber_core_ri: float or Callable[[float, float], float]
            Refractive index of the fiber core.
            A callable that takes wavelength and temperature or a fixed
            floating-point value that is independent of wavelength and
            temperature.
            Default implementation of the fused silica refractive index
            is used by default.
        fiber_core_d: float
            Diameter of the fiber core.
        fiber_cladding_d: float
            Outer diameter of the fiber cladding.
        fiber_na: float or Callable[[float, float], float]
            Numerical aperture of the fibers as a callable that takes
            wavelength and temperature or a fixed floating-point value that
            is independent of wavelength and temperature.
        fiber_spacing: float or None
            Spacing of optical fibers in the probe. If None, the fibers are
            tightly packed.
        src_fiber: int
            Index of the source fiber. Defaults to the central fiber
            with index 0.
        probe_diameter: float
            diameter of the optical fiber probe.
        probe_reflectivity: float or Callable[[float, float], float]]
            Reflectivity of the stainless steel fiber probe tip as a callable
            that takes wavelength and temperature or a fixed floating-point
            value that is independent of wavelength and temperature.
        probe_cutout: bool
            If True, the fibers of the probe are placed in a rectangular cutout
            filled with epoxy.
        epoxy_fill_ri: float or Callable[[float, float], float]
            Refractive index of the epoxy fill that surrounds the fibers in the
            cutout.
            This value is only used if the probe_cutout is set to True.
            A callable that takes wavelength and temperature or is a fixed
            floating-point value that is independent of wavelength and
            temperature.
            A constant value of 1.6 is used by default.
        cl_device: None or str or cl.Device or List[cl.Device] or cl.Context or cl.CommandQueue
            OpenCL device that will run the simulations. The first available
            GPU is used if None.
        cl_build_options: None or List[str or cloptions.ClBuildOption]
            OpenCL build options. See the Mont Carlo constructor
            :py:meth:`~xopto.mcml.mc.Mc.__init__` for more details.
            Default build options are set to:

            .. code-block:: python

                :py:code:`[mc.cloptions.FastRelaxedMath]`.

        options: None or List[mcoptions.McOption]
            Monte Carlo simulator options. See the Mont Carlo constructor
            :py:meth:`~xopto.mcml.mc.Mc.__init__` for more details.
            Default options are set to:

            .. code-block:: python

                :py:code:`[mc.mcoptions.McFloatLutMemory.constant_mem]`.

        '''
        self._suspension = suspension

        default_temperature = 293.15
        init_wavelength, init_temperature = 550e-9, default_temperature

        if cl_build_options is None:
            cl_build_options=[mc.cloptions.FastRelaxedMath]
        if options is None:
            options = [mc.mcoptions.McFloatLutMemory.constant_mem]

        if isinstance(fiber_core_ri, (float, int)):
            fiber_core_ri_value = float(fiber_core_ri)
            fiber_core_ri = \
                lambda w, t=default_temperature: fiber_core_ri_value
        if isinstance(fiber_na, (float, int)):
            fiber_na_value = float(fiber_na)
            fiber_na = lambda w, t=default_temperature: fiber_na_value
        if isinstance(epoxy_fill_ri, (float, int)):
            epoxy_fill_ri_value = float(epoxy_fill_ri)
            epoxy_fill_ri = \
                lambda w, t=default_temperature: epoxy_fill_ri_value
        if isinstance(probe_reflectivity, (float, int)):
            probe_reflectivity_value = float(probe_reflectivity)
            probe_reflectivity = \
                lambda w, t=default_temperature: probe_reflectivity_value

        if surrounding_ri is None:
            surrounding_ri = suspension.medium_ri
        if fiber_core_ri is None:
            fiber_core_ri = ri.glass.fusedsilica.default
        if epoxy_fill_ri is None:
            epoxy_fill_ri = lambda w, t=default_temperature: 1.6

        if cl_device is None:
            cl_device = mc.clinfo.gpu()

        fiber_core_d = float(fiber_core_d)
        fiber_cladding_d = float(fiber_cladding_d)
        src_fiber = int(src_fiber)
        probe_diameter = float(probe_diameter)
    
        if fiber_spacing is None:
            fiber_spacing = 3.0*fiber_cladding_d

        self._epoxy_fill_ri = epoxy_fill_ri
        self._surrounding_ri = epoxy_fill_ri

        self._fiber_core_ri = fiber_core_ri
        self._fiber_na = fiber_na

        self._probe_reflectivity = probe_reflectivity

        multimode_fiber = fiber.MultimodeFiber(
            dcore=fiber_core_d, dcladding=fiber_cladding_d,
            ncore=fiber_core_ri(init_wavelength, init_temperature),
            na=fiber_na(init_wavelength, init_temperature)
        )

        if probe_cutout:
            cutout = 2*fiber_spacing + fiber_cladding_d
        else:
            cutout = 0.0

        top_surface_layout = mc.mcsurface.SixAroundOne(
            fiber=multimode_fiber,
            spacing=fiber_spacing,
            cutout=cutout,
            cutoutn=epoxy_fill_ri(init_wavelength, init_temperature),
            reflectivity=probe_reflectivity(init_wavelength, init_temperature),
            diameter=probe_diameter
        )

        surface_layouts = mc.mcsurface.SurfaceLayouts(
            top=top_surface_layout
        )

        source = mc.mcsource.UniformFiber(
            fiber=multimode_fiber,
            position=top_surface_layout.fiber_position(src_fiber) + (0.0,)
        )

        top_detector = mc.mcdetector.SixAroundOne(
            fiber=multimode_fiber,
            spacing=fiber_spacing,
            position=top_surface_layout.position
        )

        detectors = mc.mcdetector.Detectors(
            top=top_detector
        )

        mcpf_obj = self._suspension.mcpf(init_wavelength, init_temperature)

        layers = mc.mclayer.Layers([
            mc.mclayer.Layer(d=float('inf'), mua=0.0, mus=0.0,
                             n=surrounding_ri(init_wavelength),
                             pf=mcpf_obj),
            mc.mclayer.Layer(d=float('inf'), mua=0.0, mus=0.0,
                             n=suspension.medium_ri(init_wavelength),
                             pf=mcpf_obj),
            mc.mclayer.Layer(d=float('inf'), mua=0.0, mus=0.0,
                             n=surrounding_ri(init_wavelength),
                             pf=mcpf_obj)
        ])

        self._mc_obj = mc.Mc(
            layers, source, detectors, surface=surface_layouts,
            cl_build_options=cl_build_options, options=options,
            cl_devices=cl_device)
        self._mc_obj.rmax = 25e-3

    def run(self, wavelength: float, temperature=293.15, *args, **kwargs):
        '''
        Updates the optical with values computed the provided wavelength and
        temperature, and runs the MC simulation.

        Parameters
        ----------
        wavelength: float
            Wavelength of light (m) that is used to update the optical
            properties.
        temperature: float
            Temperature of the suspension in (K).
        *args, **kwargs:
            Positional and keyword arguments passed to the MC object.

        Returns
        -------
        data: dict
            Simulation results and configuration as a dict with keys
            "reflectance", "wavelength", "temperature" and "mc".
        '''
        wavelength = float(wavelength)
        temperature = float(temperature)

        # recompute wavelength and temperature dependent values
        medium_ri = float(self._suspension.medium_ri(wavelength, temperature))
        surrounding_ri = float(self._surrounding_ri(wavelength, temperature))
        fiber_na = float(self._fiber_na(wavelength, temperature))
        fiber_core_ri = float(self._fiber_core_ri(wavelength, temperature))
        probe_reflectivity = float(
            self._probe_reflectivity(wavelength, temperature))
        epoxy_fill_ri = float(self._epoxy_fill_ri(wavelength, temperature))

        # compute scattering phase function LUT using cache
        mcpf_obj = self._suspension.mcpf(wavelength, temperature)

        mc_obj = self.mc_obj

        # update the layers
        mc_obj.layers[0].pf = mcpf_obj
        mc_obj.layers[0].n = surrounding_ri
        #
        mc_obj.layers[1].pf = mcpf_obj
        mc_obj.layers[1].n = medium_ri
        mc_obj.layers[1].mus = self._suspension.mus(wavelength, temperature)
        mc_obj.layers[1].mua = self._suspension.medium_mua(
            wavelength, temperature)
        #
        mc_obj.layers[2].pf = mcpf_obj
        mc_obj.layers[2].n = surrounding_ri

        # update the source
        mc_obj.source.fiber.na = fiber_na
        mc_obj.source.fiber.ncore = fiber_core_ri

        # update the top surface layout
        mc_obj.surface.top.fiber.na = fiber_na
        mc_obj.surface.top.fiber.ncore = fiber_core_ri
        mc_obj.surface.top.ncutout = epoxy_fill_ri
        mc_obj.surface.top.reflectivity = probe_reflectivity

        # update the reflectance detector
        mc_obj.detectors.top.fiber.na = fiber_na
        mc_obj.detectors.top.fiber.ncore = fiber_core_ri

        _, _, detectors_res = mc_obj.run(*args, **kwargs)

        data = {
            'epoch_timestamp': time.time(),
            'mc': {
                'source': mc_obj.source.todict(),
                'surface': mc_obj.surface.todict(),
                'layers': mc_obj.layers.todict(),
                'detectors': mc_obj.detectors.todict(),
                'num_packets': int(detectors_res.top.nphotons),
                'rmax': mc_obj.rmax,
                'run_report': mc_obj.run_report
            },
            'reflectance': detectors_res.top.reflectance,
            'wavelength': wavelength,
            'temperature': temperature,
        }

        return data

    def _get_mc_obj(self) -> mc.Mc:
        return self._mc_obj
    mc_obj = property(_get_mc_obj, None, None, 'MC simulator instance.')

    def _get_suspension(self) -> Suspension:
        return self._suspension
    suspension = property(_get_suspension, None, None,
                          'The underlying suspension.')

