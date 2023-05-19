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

import numpy as np
from matplotlib import pyplot as pp


# creating a default 3-layer skin model
model = skin.Skin3()

# using the skin model to create the layer stack (5 MC layers, i.e. the first
# and the last layer in the stack represent the surrounding medium)
layers = model.create_mc_layers(550e-9)

# use low melanin and blood content
model[0].melanin = 0.001
model[1].blood = 0.02
model[2].blood = 0.01

# creating the source of photon packets
source = mc.mcsource.Line(direction=[np.sin(np.pi/4), 0.0, np.cos(np.pi/4)])

# using CIE1964 observer, field of view is 10 deg
cosmin = np.cos(np.deg2rad(10.0/2))

# creating total reflectance surface detector
detectors = mc.mcdetector.Detectors(
    top=mc.mcdetector.Total(cosmin=cosmin),
    specular=mc.mcdetector.Total(cosmin=cosmin)
)
# include specular component in color computation
with_specular = False

# selecting the first available OpenCL GPU device
gpu = mc.clinfo.gpu()

# creating a Monte Carlo simulator
mc_obj = mc.Mc(layers, source, detectors, cl_devices=gpu)
mc_obj.rmax = 25.0e-3

# standard wavelength range
wavelengths = color.cie.STANDARD_WAVELENGTHS_10_NM

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
    if with_specular:
        r += detectors_res.specular.reflectance
    # using SRGB color space illuminant spectral power density
    reflectance[index] = r*color.cie.SRGB.illuminant.spd(wavelength)

print()
# compute color, correct illuminant luminosity for Lambertian reflector and
# solid acceptance angle of the CIE observer
rgb = color.cie.SRGB.from_spectrum(wavelengths, reflectance,
                                   observer=color.cie.CIE1964,
                                   normalize=1.0/(1.0 - cosmin))
print('Raw normalized SRGB components:', rgb)
fig, ax = pp.subplots(1, 2)
ax[0].plot(wavelengths*1e9, reflectance)
ax[0].set_xlabel('Wavelength (nm)')
ax[0].set_ylabel('Reflectance (a.u.)')
ax[0].set_title('Reflectance for illuminant {:s}'.format(
    color.cie.SRGB.illuminant.name))
rgb_uint8 = np.round(np.clip(rgb, 0.0, 1.0)*255).astype(np.uint8)
r, g, b = rgb_uint8 
rgb_uint8.shape = (1, 1, 3)
ax[1].imshow(rgb_uint8)
ax[1].set_title('Skin color: RGB=({:d}, {:d}, {:d})'.format(r, g, b))
pp.show()