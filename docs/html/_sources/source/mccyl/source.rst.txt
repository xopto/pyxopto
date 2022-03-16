.. ****************************** Begin license ********************************
.. Copyright (C) Laboratory of Imaging technologies,
..               Faculty of Electrical Engineering,
..               University of Ljubljana.
..
.. This file is part of PyXOpto.
..
.. PyXOpto is free software: you can redistribute it and/or modify
.. it under the terms of the GNU General Public License as published by
.. the Free Software Foundation, either version 3 of the License, or
.. (at your option) any later version.
..
.. PyXOpto is distributed in the hope that it will be useful,
.. but WITHOUT ANY WARRANTY; without even the implied warranty of
.. MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
.. GNU General Public License for more details.
..
.. You should have received a copy of the GNU General Public License
.. along with PyXOpto. If not, see <https://www.gnu.org/licenses/>.
.. ******************************* End license *********************************

.. include:: ../common.rst

.. _mccyl-source-label:

Sources
=======

All the photon packet sources that can be used in the Monte Carlo
simulator are implemented by subclassing
:py:class:`xopto.mccyl.mcsource.base.Source`. The 
:py:mod:`xopto.mccyl.mcsource` module includes a number of commonly
used sources:

* :py:class:`xopto.mccyl.mcsource.point.IsotropicPoint` implements an isotropic
  point source.

* :py:class:`xopto.mccyl.mcsource.line.Line` implements an infinitely narrow and
  collimated source.

* :py:class:`xopto.mccyl.mcsource.uniformbeam.UniformBeam` implements a 
  uniform collimated source of elliptical cross section.

* :py:class:`xopto.mccyl.mcsource.gaussianbeam.GaussianBeam` implements
  a collimated Gaussian source of elliptical cross section.

All sources have a position attribute that defines the position of the
geometrical center of the source. This position is used as a reference point
for terminating photon packets that leave the simulation radius set
through the :py:attr:`xopto.mccyl.mc.Mc.rmax` property. Note that by default
all the sources point along the positive direction of the :math:`x` axis and
enter the outer layer of the sample on the negative :math:`x` axis at
:math:`x=d/2` .

.. note::

    The refractive index of the source (not all the sources implement
    this property) is used to compute the initial photon packet weight (the
    specular reflectance at the source-sample boundary is subtracted from 1.0).
    The sources that do not implement a refractive index, inherit the value from
    the surrounding medium.

The :py:class:`~xopto.mccyl.mcsource.uniformbeam.UniformBeam`,
:py:class:`xopto.mccyl.mcsource.gaussianbeam.GaussianBeam` and all the
fiber sources

Isotropic point
---------------

This source can be placed inside and outside the sample. The following example
shows how to create an isotropic point source in the center of the sample:

.. code-block:: python

    from xopto.mccyl import mc

    src = mc.mcsource.IsotropicPoint((0.0, 0.0, 0.0e-3))

If the source is positioned outside of the sample, the entry point
into the sample is determined by propagating the packet from the source
position along the launch direction. The MC simulation will start after
refracting the packet into the sample and subtracting the specular reflectance
at the sample boundary from the initial weight of the packet.
If the photon packet does not intersect the sample, the
initial weight will be set to 0 and the packet
will be launched from the center of the sample :math:`(x,y,z)=(0, 0, 0)`. Such
zero-weight packets are immediately terminated and have no contribution to the
fluence and surface detector, however will be included in the trace (have no 
effect on the sampling volume or other trace-based analysis due to the
zero-weight).
Note that in case the position lies within the sample, it will be used as the
launch point and the packets will retain the full initial weight.

Line
----

The following example shows how to create a Line source that enters the sample
1 |nbsp| mm of the :math:`x` axis and propagates in the direction of the
:math:`x` axis:

.. code-block:: python

    from xopto.mccyl import mc
    import numpy as np

    ang = np.deg2rad(10.0)
    src = mc.mcsource.Line(position=(0.0, 1.0e-3, 0.0), direction=(1.0, 0.0, 0.0))

The packets are always launched from the outer surface of the sample. The
source position and direction are used to determine the launch point at the
sample surface. From there, the packet is first refracted into the sample
and the surface reflectance is subtracted from the initial packet weight.
If a specular surface :ref:`detector <mccyl-detector-label>` is used, the
reflectance is deposited into that detector. 

Uniform beam
---------------
A uniform beam of elliptical cross section can be created with
:py:class:`xopto.mccyl.mcsource.uniformbeam.UniformBeam` source. The beam
diameter along the y and z axis is controlled by the :code:`diameter` parameter.
In case the diameter is give as a scalar :code:`float` value, the cross
section of the beam becomes circular. Note that the diameters are applied in
the coordinate system of the beam. Optionally, the beam can be repositioned
and tilted through the :code:`position` and :code:`direction` parameters.
The following example creates a uniform beam with incidence along the
:math:`axis` and a circular cross section diameter of 1 |nbsp| mm.

.. code-block:: python

    from xopto.mccyl import mc

    src = mc.mcsource.UniformBeam(1.0e-3)

The packets are always launched from the top surface of the sample. The
source position and direction are used to determine the launch point at the
sample surface. From there, the packet is first refracted into the sample
and the surface reflectance is subtracted from the initial packet weight.
If a specular surface :ref:`detector <mccyl-detector-label>` is used, the
reflectance is deposited into that detector. 

Gaussian beam
-------------

A Gaussian beam can be created with
:py:class:`xopto.mccyl.mcsource.gaussianbeam.GaussianBeam` source. The standard
deviation (:math:`\sigma`) of the beam can be independently set along
the x and y axis.
In case :math:`\sigma` is give as a scalar :code:`float` value, the value is
applied along the y and z axis. The beam is clipped at a distance of
:math:`5\sigma` from the beam axis. The clip distance can be customized through
the :code:`clip` parameter. Note that :math:`\sigma` is applied in
the coordinate system of the beam. Optionally, the beam can be repositioned
and tilted through the :code:`position` and :code:`direction` parameters.
The following example creates a
:py:class:`~xopto.mccyl.mcsource.gaussianbeam.GaussianBeam` source with a
perpendicular incidence and with :math:`\sigma_y` 1 |nbsp| mm and
:math:`\sigma_z` 2 |nbsp| mm.

.. code-block:: python

    from xopto.mccyl import mc

    src = mc.mcsource.GaussianBeam([1.0e-3, 2.0e-3])

The packets are always launched from the surface of the sample. The
source position and direction are used to determine the launch point at the
sample surface. From there, the packet is first refracted into the sample
and the reflectance is subtracted from the initial packet weight.
If a specular surface :ref:`detector <mccyl-detector-label>` is used, the
reflectance is deposited into that detector. 

