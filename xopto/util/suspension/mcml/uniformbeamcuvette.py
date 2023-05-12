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

from typing import Callable, Tuple, List
import time

from xopto.mcml import mc
from xopto.materials import ri
from xopto.mcml.mcutil import fiber

from xopto.util.suspension import Suspension

import numpy as np


class UniformBeamCuvette:
    def __init__(
            self,
            suspension: Suspension,
            cuvette_ri: None or float or Callable([float, float], float) = None,
            cuvette_mua: float or Callable[[float, float], float] = 0.0,
            cuvette_wall: float = 1.25e-3,
            cuvette_path: float = 10.0e-3,
            surrounding_ri: None or float or Callable([float, float], float) = None,
            surrounding_mua: float or Callable([float, float], float) = 0.0,
            beam_d: float = 4.0e-3,
            beam_na: float = 0.01,
            beam_fi_theta: Tuple[float, float] = (0.0, 0.0),
            detector_d: float = 4.0e-3,
            detector_na: float = 0.01,
            detector_fi_theta: Tuple[float, float] = (0.0, 0.0),
            detector_lateral_offset: Tuple[float, float]= (0.0, 0.0),
            detector_axial_offset: float = 0.0,
            cl_build_options: None or List[str or mc.cloptions.ClBuildOption] = None,
            options: None or List[mc.mcoptions.McOption] = None,
            cl_device: str or mc.cl.Device or List[mc.cl.Device] or
                       mc.cl.Context or mc.cl.CommandQueue = None):
        '''
        Parameters
        ----------
        suspension: Suspension
            Suspension instance use with the Monte Carlo simulations.
        cuvette_ri:  float or Callable[[float, float], float]
            Refractive index of the cuvette.
            A callable that takes wavelength and temperature or a fixed
            floating-point value that is independent of wavelength and
            temperature.
            Default implementation uses fused silica.
        cuvette_mua: float or Callable[[float, float], float]
            Absorption coefficient (1/m) of the cuvette.
            A callable that takes wavelength and temperature or a fixed
            floating-point value that is independent of wavelength and
            temperature.
            Default implementation uses 0.0 (no absorption).
        cuvette_wall: float
            Cuvette wall thickness (m).
        cuvette_path: float
            Clear path length of the cuvette (m).
        surrounding_ri:  float or Callable[[float, float], float]
            Refractive index of the medium surrounding the cuvette.
            A callable that takes wavelength and temperature or a fixed
            floating-point value that is independent of wavelength and
            temperature.
            Default implementation uses 1 (air).
        surrounding_mua: float or Callable[[float, float], float]
            Absorption coefficient (1/m) of the surrounding medium.
            A callable that takes wavelength and temperature or a fixed
            floating-point value that is independent of wavelength and
            temperature.
            Default implementation uses 0.0 (no absorption).
        beam_d: float
            Beam diameter (m) at the position of the source.
        beam_na: float or Callable[[float, float], float]
            Numerical aperture of the beam as a callable that takes
            wavelength and temperature or a fixed floating-point value that
            is independent of wavelength and temperature.
            Divergence of the beam be computed as :math:`\\asin(NA)`.
        beam_fi_theta: Tuple[float, float]
            Beam azimuth and incidence angles.
        detector_d: float
            Detector diameter.
        detector_na: float or Callable[[float, float], float]
            Numerical aperture of the detector as a callable that takes
            wavelength and temperature or a fixed floating-point value that
            is independent of wavelength and temperature.
            Detector acceptance angle can be computed as :math:`\\asin(NA)`.
        detector_fi_theta: float
            Detector azimuth and incidence angles (rad).
        detector_lateral_offset: Tuple[float, float]
            Later displacement of the detector from the axis of the
            source as a tuple (x, y).
        detector_axial_offset: float
            Axial distance (m) of the detector from the transmittance side
            of the cuvette (must be zero or a positive value). 
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
        if detector_axial_offset < 0.0:
            raise ValueError(
                'Distance of the detector from the cuvette surface '
                '(detector_axial_offset) must be a positive value!')

        self._suspension = suspension

        default_temperature = 293.15
        init_wavelength, init_temperature = 550e-9, default_temperature

        if cl_build_options is None:
            cl_build_options=[mc.cloptions.FastRelaxedMath]
        if options is None:
            options = [mc.mcoptions.McFloatLutMemory.constant_mem]

        if surrounding_ri is None:
            surrounding_ri = 1.0

        if surrounding_mua is None:
            surrounding_mua = 0.0

        if cuvette_ri is None:
            cuvette_ri = ri.glass.fusedsilica.default

        if cuvette_mua is None:
            cuvette_mua = 0.0

        if isinstance(cuvette_ri, (float, int)):
            cuvette_ri_value = float(cuvette_ri)
            cuvette_ri = \
                lambda w, t=default_temperature: cuvette_ri_value

        if isinstance(cuvette_mua, (float, int)):
            cuvette_mua_value = float(cuvette_mua)
            cuvette_mua = \
                lambda w, t=default_temperature: cuvette_mua_value

        if isinstance(surrounding_ri, (float, int)):
            surrounding_ri_value = float(surrounding_ri)
            surrounding_ri = \
                lambda w, t=default_temperature: surrounding_ri_value

        if isinstance(surrounding_mua, (float, int)):
            surrounding_mua_value = float(surrounding_mua)
            surrounding_mua = \
                lambda w, t=default_temperature: surrounding_mua_value

        if isinstance(beam_na, (float, int)):
            beam_na_value = float(beam_na)
            beam_na = lambda w, t=default_temperature: beam_na_value

        if isinstance(detector_na, (float, int)):
            detector_na_value = float(detector_na)
            detector_na = lambda w, t=default_temperature: detector_na_value

        if cl_device is None:
            cl_device = mc.clinfo.gpu()

        self._surrounding_ri = surrounding_ri
        self._surrounding_mua = surrounding_mua

        self._cuvette_ri = cuvette_ri
        self._cuvette_mua = cuvette_mua
        self._cuvette_wall = float(cuvette_wall)
        self._cuvette_path = float(cuvette_path)

        self._beam_d = float(beam_d)
        self._beam_fi_theta = (float(beam_fi_theta[0]),
                               float(beam_fi_theta[1]))
        self._beam_na = beam_na

        self._detector_d = float(detector_d)
        self._detector_fi_theta = (float(detector_fi_theta[0]),
                                   float(detector_fi_theta[1]))
        self._detector_na = detector_na
        self._detector_lateral_offset = (float(detector_lateral_offset[0]),
                                         float(detector_lateral_offset[1]))
        self._detector_axial_offset = float(detector_axial_offset)

        # create initialize the simulator objects
        surrounding_ri_value = \
            self._surrounding_ri(init_wavelength, init_temperature)
        surrounding_mua_value = \
            self._surrounding_mua(init_wavelength, init_temperature)

        cuvette_ri_value = self._cuvette_ri(init_wavelength, init_temperature)
        cuvette_mua_value = self._cuvette_mua(init_wavelength, init_temperature)

        suspension_ri_value = \
            self._suspension.medium_ri(init_wavelength, init_temperature)
        suspension_mua_value = \
            self._suspension.mua(init_wavelength, init_temperature)
        suspension_mus_value = \
            self._suspension.mus(init_wavelength, init_temperature)

        beam_na_value = self._beam_na(init_wavelength, init_temperature)
        detector_na_value = self._detector_na(init_wavelength, init_temperature)

        source_fiber = fiber.MultimodeFiber(
            dcore=self._beam_d, dcladding=self._beam_d,
            ncore=surrounding_ri_value,
            na=beam_na_value
        )

        detector_fiber = fiber.MultimodeFiber(
            dcore=self._detector_d, dcladding=self._detector_d,
            ncore=surrounding_ri_value,
            na=detector_na_value
        )

        fi, theta = self._beam_fi_theta
        source = mc.mcsource.UniformFiber(
            fiber=source_fiber,
            direction=(
                np.sin(theta)*np.cos(fi), 
                np.sin(theta)*np.sin(fi), 
                np.cos(theta))
        )

        fi, theta = self._detector_fi_theta
        bottom_detector = mc.mcdetector.LinearArray(
            fiber=detector_fiber,
            n=1,
            spacing=0.0,
            position=self._detector_lateral_offset,
            direction=(
                np.sin(theta)*np.cos(fi), 
                np.sin(theta)*np.sin(fi), 
                np.cos(theta)
            )
        )

        detectors = mc.mcdetector.Detectors(
            bottom=bottom_detector
        )

        mcpf_obj = suspension.mcpf(init_wavelength, init_temperature)

        layers = [
            mc.mclayer.Layer(d=float('inf'),
                             mua=surrounding_mua_value, mus=0.0,
                             n=surrounding_ri_value, pf=mcpf_obj),

            mc.mclayer.Layer(d=self._cuvette_wall,
                             mua=cuvette_mua_value, mus=0.0,
                             n=cuvette_ri_value, pf=mcpf_obj),
            mc.mclayer.Layer(d=self._cuvette_path,
                             mua=suspension_mua_value, mus=suspension_mus_value,
                             n=suspension_ri_value, pf=mcpf_obj),
            mc.mclayer.Layer(d=self._cuvette_wall,
                             mua=cuvette_mua_value, mus=0.0,
                             n=cuvette_ri_value, pf=mcpf_obj),

            mc.mclayer.Layer(d=float('inf'),
                             mua=surrounding_mua_value, mus=0.0,
                             n=surrounding_ri_value, pf=mcpf_obj)
        ]
        if self._detector_axial_offset > 0.0:
            layers.insert(
                4,
                mc.mclayer.Layer(
                    d=self._detector_axial_offset,
                    mua=surrounding_mua_value, mus=0.0,
                    n=surrounding_ri_value, pf=mcpf_obj),
            )
        layers = mc.mclayer.Layers(layers)

        self._mc_obj = mc.Mc(
            layers, source, detectors,
            cl_build_options=cl_build_options, options=options,
            cl_devices=cl_device)
        self._mc_obj.rmax = max(25e-3, 2.0*layers.thickness())

    def run(self, wavelength: float, temperature=293.15, *args, **kwargs) \
            -> dict:
        '''
        Updates the the optical with values computed at the provided wavelength
        and temperature, and runs the MC simulation.

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

        mc_obj = self.mc_obj

        # recompute wavelength and temperature dependent values
        surrounding_ri = float(self._surrounding_ri(wavelength, temperature))
        surrounding_mua = float(self._surrounding_mua(wavelength, temperature))

        beam_na = float(self._beam_na(wavelength, temperature))
        detector_na = float(self._detector_na(wavelength, temperature))

        cuvette_ri = float(self._cuvette_ri(wavelength, temperature))
        cuvette_mua = float(self._cuvette_mua(wavelength, temperature))

        suspension_mus = self._suspension.mus(wavelength, temperature)
        suspension_mua = self._suspension.mua(wavelength, temperature)
        suspension_ri = self._suspension.medium_ri(wavelength, temperature)

        # compute scattering phase function LUT using cache
        mcpf_obj = self._suspension.mcpf(wavelength, temperature)

        # update the layer stack
        # surrounding medium on the incident side
        mc_obj.layers[0].pf = mcpf_obj
        mc_obj.layers[0].n = surrounding_ri
        mc_obj.layers[0].mua = surrounding_mua
        # cuvette wall
        mc_obj.layers[1].pf = mcpf_obj
        mc_obj.layers[1].n = cuvette_ri
        mc_obj.layers[1].mua = cuvette_mua
        # sample inside the cuvette
        mc_obj.layers[2].pf = mcpf_obj
        mc_obj.layers[2].n = suspension_ri
        mc_obj.layers[2].mus = suspension_mus
        mc_obj.layers[2].mua = suspension_mua
        # cuvette wall
        mc_obj.layers[3].pf = mcpf_obj
        mc_obj.layers[3].n = cuvette_ri
        mc_obj.layers[3].mua = cuvette_mua
        # axial offset layer (between the cuvette and detector)
        if self._detector_axial_offset > 0.0:
            mc_obj.layers[-2].pf = mcpf_obj
            mc_obj.layers[-2].n = surrounding_ri
            mc_obj.layers[-2].mua = surrounding_mua
        # surrounding medium on the transmittance side
        mc_obj.layers[-1].pf = mcpf_obj
        mc_obj.layers[-1].n = surrounding_ri
        mc_obj.layers[-1].mua = surrounding_mua

        # update the source
        mc_obj.source.fiber.na = beam_na
        mc_obj.source.fiber.ncore = surrounding_ri

        # update the transmittance detector
        mc_obj.detectors.bottom.fiber.na = detector_na
        mc_obj.detectors.bottom.fiber.ncore = surrounding_ri

        _, _, detectors_res = mc_obj.run(*args, **kwargs)

        data = {
            'epoch_timestamp': time.time(),
            'mc': {
                'source': mc_obj.source.todict(),
                'layers': mc_obj.layers.todict(),
                'detectors': mc_obj.detectors.todict(),
                'num_packets': int(detectors_res.bottom.nphotons),
                'rmax': mc_obj.rmax,
                'run_report': mc_obj.run_report
            },
            'transmittance': detectors_res.bottom.reflectance,
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

if __name__ == '__main__':
    from xopto.pf.distribution import Normal

    cuvette = UniformBeamCuvette(
        Suspension(Normal(1.0e-6, 0.01e-6))
    )

    r = cuvette.run(550e-9, nphotons=100)