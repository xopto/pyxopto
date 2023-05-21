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

from xopto.mcml import mc
from xopto.materials import skin
from xopto.util import color
from xopto.materials import ri

import numpy as np
from matplotlib import pyplot as pp


# creating a default 3-layer skin model
model = skin.Skin3()

# using the skin model to create the layer stack (5 MC layers, i.e. the first
# and the last layer in the stack represent the surrounding medium)
layers = model.create_mc_layers(550e-9)

# customization of the skin model
model[0].n = ri.skin.default
model[0].d = 0.100e-3
model[0].melanin = 0.02
model[0].baseline_absorption = 0.0
model[1].d = 1.5e-3
model[1].blood = 0.005
model[1].spo2 = 0.97
model[1].baseline_absorption = 0.0
model[2].blood = 0.01
model[2].spo2 = 0.97
model[2].baseline_absorption = 0.0

# creating the source of photon packets
incidence_angle = np.deg2rad(45.0) # 0 deg for perpendicular incidence
source = mc.mcsource.Line(
    direction=[np.sin(incidence_angle), 0.0, np.cos(incidence_angle)])

# RGB color space and observer selection
RGB = color.cie.SRGB
#observer = color.cie.CIE1931 # 2 deg field of view
observer = color.cie.CIE1964 # 10 deg field of view
observer_cosmin = np.cos(observer.fov/2)

# creating total reflectance surface detector
detectors = mc.mcdetector.Detectors(
    top=mc.mcdetector.Total(cosmin=observer_cosmin),
    specular=mc.mcdetector.Total(cosmin=observer_cosmin)
)
# include specular component in color computation
with_specular = False

# selecting the first available OpenCL GPU device
gpu = mc.clinfo.gpu()

# creating a Monte Carlo simulator instance
mc_obj = mc.Mc(layers, source, detectors, cl_devices=gpu)
mc_obj.rmax = 50.0e-3

# standard wavelength range with a 10 nm spectral resolution
wavelengths = color.cie.STANDARD_WAVELENGTHS_10_NM

# running MC simulations for each wavelength
reflectance = np.zeros_like(wavelengths)
for index, wavelength in enumerate(wavelengths):
    print('Processing wavelength {:d}/{:<10d}'.format(
        index + 1, wavelengths.size), end='\r')
    # update the MC layers with the skin model
    model.update_mc_layers(mc_obj.layers, wavelength)
    # running the MC simulation with 1,000,000 photon packets
    _, _, detectors_res = mc_obj.run(1e6)
    # extract reflectance
    r = detectors_res.top.reflectance
    if with_specular: # include specular component if required
        r += detectors_res.specular.reflectance
    # using RGB color space illuminant spectral power density
    reflectance[index] = r*RGB.illuminant.spd(wavelength)

print()
# compute color, correct illuminant luminosity for Lambertian reflector and
# solid acceptance angle of the CIE observer
rgb = RGB.from_spectrum(wavelengths, reflectance,
                        observer=observer,
                        normalize=(1.0 - observer_cosmin))
rgb_uint8 = np.round(np.clip(rgb, 0.0, 1.0)*255).astype(np.uint8)
print('Raw normalized RGB components:', rgb)

fig, ax = pp.subplots(1, 2)
pp.suptitle('Skin color - {} color space'.format(RGB.name))
ax[0].plot(wavelengths*1e9, reflectance, label='Skin')
ax[0].plot(wavelengths*1e9,
           RGB.illuminant.spd(wavelengths)*(1.0 - observer_cosmin),
           label='Lambertian-{}'.format(RGB.illuminant.name))
ax[0].set_xlabel('Wavelength (nm)')
ax[0].set_ylabel('Reflectance (a.u.)')
ax[0].set_title('Illuminant {:s} & observer {:s}'.format(
    RGB.illuminant.name, observer.name))
ax[0].legend()

ax[1].imshow(np.reshape(rgb_uint8, [1, 1, 3]))
ax[1].set_title('Skin color: RGB=({:d}, {:d}, {:d})'.format(*rgb_uint8))
ax[1].set_axis_off()
pp.tight_layout()
pp.show()
