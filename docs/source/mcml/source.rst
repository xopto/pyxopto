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

.. _mcml-source-label:

Sources
=======

All the photon packet sources that can be used in the Monte Carlo
simulator are implemented by subclassing
:py:class:`xopto.mcml.mcsource.base.Source`. The 
:py:mod:`xopto.mcml.mcsource` module includes a number of commonly
used sources:

* :py:class:`xopto.mcml.mcsource.point.IsotropicPoint` implements an isotropic
  point source.

* :py:class:`xopto.mcml.mcsource.line.Line` implements an infinitely narrow and
  collimated source.

* :py:class:`xopto.mcml.mcsource.uniformbeam.UniformBeam` implements a 
  uniform collimated source of elliptical cross section.

* :py:class:`xopto.mcml.mcsource.gaussianbeam.GaussianBeam` implements
  a collimated Gaussian source of elliptical cross section.

* :py:mod:`xopto.mcml.mcsource.fiber` module implements several optical fiber sources:
  
    * :py:mod:`xopto.mcml.mcsource.fiber.UniformFiber` implements a fiber
      source that emits uniformly within the NA from each point of the fiber
      core surface.
    * :py:mod:`xopto.mcml.mcsource.fiber.UniformFiberLut` implements a fiber
       source that follows a measured or nonparametric angular emission.
    * :py:mod:`xopto.mcml.mcsource.fiber.LambertianFiber` implements a fiber
      source that emits like a Lambertian surface within the NA of the core.

* :py:mod:`xopto.mcml.mcsource.fiberni` module implements the same sources as
  the :py:mod:`xopto.mcml.mcsource.fiber` module but for a normal source
  incidence.

* :py:mod:`xopto.mcml.mcsource.rectangular` module implements the following sources:

    * :py:mod:`xopto.mcml.mcsource.rectangular.UniformRectangular` implements
      a rectangular source that emits uniformly within the NA from each point
      of the source surface.
    * :py:mod:`xopto.mcml.mcsource.rectangular.UniformRectangularLut` implements
      a rectangular source that follows a measured or nonparametric angular
      emission characteristics.
    * :py:mod:`xopto.mcml.mcsource.rectangular.LambertianRectangular` implements
      a rectangular source that emits like a Lambertian surface within the NA
      of the source.

All sources have a position attribute that defines the position of the
geometrical center of the source. This position is used as a reference point
for terminating photon packets that leave the simulation radius set
through the :py:attr:`xopto.mcml.mc.Mc.rmax` property.

Note that the refractive index of the source (not all the sources implement
this property) is used to compute the initial photon packet weight (the
specular reflectance at the source-sample boundary is subtracted from 1.0).
The sources that do not implement a refractive index, inherit the value from
the surrounding medium.

The :py:class:`~xopto.mcml.mcsource.uniformbeam.UniformBeam`,
:py:class:`xopto.mcml.mcsource.gaussianbeam.GaussianBeam` and all the
fiber sources

Isotropic point
---------------

This source can be placed inside and above the sample. The following example
shows how to create an isotropic point source 5 |nbsp| mm bellow the top
surface of the sample:

.. code-block:: python

    from xopto.mcml import mc

    src = mc.mcsource.IsotropicPoint((0.0, 0.0, 5.0e-3))

If the source is positioned above the sample (z <= 0), the entry point
into the sample is determined by propagating the packet from the source
position along the launch direction. The MC simulation will start after
refracting the packet into the sample and subtracting the specular reflectance
at the sample boundary from the initial weight of the packet.
If the photon packet does not intersect the sample, the
initial weight will be set to 0 (reflectance to 1) and the packet
will be launched with the z coordinate set to 0. If a specular surface
:ref:`detector <mcml-detector-label>` is used, the reflectance is deposited into 
that detector. Such zero-weight packets are
immediately terminated and have no contribution to the fluence and surface
detectors, however will be included in the trace (have no effect on
the sampling volume or other trace-based analysis due to the zero-weight).
Note that in case the position lies within the sample (z > 0),
it will be used as the launch point and the packets will retain the
full initial weight.

Line
----

The following example shows how to create a Line source with a
10 :superscript:`o` incidence (tilted along the x axis) that is located at
the center of the top sample surface :math:`(0, 0, 0)`:

.. code-block:: python

    from xopto.mcml import mc
    import numpy as np

    ang = np.deg2rad(10.0)
    src = mc.mcsource.Line(direction=(np.sin(ang), 0.0, np.cos(ang)))

The packets are always launched from the top surface of the sample. The
source position and direction are used to determine the launch point at the
sample surface. From there, the packet is first refracted into the sample
and the surface reflectance is subtracted from the initial packet weight.
If a specular surface :ref:`detector <mcml-detector-label>` is used, the
reflectance is deposited into that detector. 

Uniform beam
---------------
A uniform beam of elliptical cross section can be created with
:py:class:`xopto.mcml.mcsource.uniformbeam.UniformBeam` source. The beam
diameter along the x and y axis is controlled by the :code:`diameter` parameter.
In case the diameter is give as a scalar :code:`float` value, the cross
section of the beam becomes circular. Note that the diameters are applied in
the coordinate system of the beam. Optionally, the beam can be repositioned
and tilted through the :code:`position` and :code:`direction` parameters.
The following example creates a uniform beam with a perpendicular incidence
and a circular cross section diameter of 1 |nbsp| mm.

.. code-block:: python

    from xopto.mcml import mc

    src = mc.mcsource.UniformBeam(1.0e-3)

The packets are always launched from the top surface of the sample. The
source position and direction are used to determine the launch point at the
sample surface. From there, the packet is first refracted into the sample
and the surface reflectance is subtracted from the initial packet weight.
If a specular surface :ref:`detector <mcml-detector-label>` is used, the
reflectance is deposited into that detector. 

Gaussian beam
-------------

A Gaussian beam can be created with
:py:class:`xopto.mcml.mcsource.gaussianbeam.GaussianBeam` source. The standard
deviation (:math:`\sigma`) of the beam can be independently set along
the x and y axis.
In case :math:`\sigma` is give as a scalar :code:`float` value, the value is
applied along the x and y axis. The beam is clipped at a distance of
:math:`5\sigma` from the beam axis. The clip distance can be customized through
the :code:`clip` parameter. Note that :math:`\sigma` is applied in
the coordinate system of the beam. Optionally, the beam can be repositioned
and tilted through the :code:`position` and :code:`direction` parameters.
The following example creates a
:py:class:`~xopto.mcml.mcsource.gaussianbeam.GaussianBeam` source with a
perpendicular incidence and with :math:`\sigma_x` 1 |nbsp| mm and
:math:`\sigma_y` 2 |nbsp| mm.

.. code-block:: python

    from xopto.mcml import mc

    src = mc.mcsource.GaussianBeam([1.0e-3, 2.0e-3])

The packets are always launched from the top surface of the sample. The
source position and direction are used to determine the launch point at the
sample surface. From there, the packet is first refracted into the sample
and the reflectance is subtracted from the initial packet weight.
If a specular surface :ref:`detector <mcml-detector-label>` is used, the
reflectance is deposited into that detector. 

Fiber
-----
This example shows how to create an optical fiber source with an :math:`NA=0.22`,
a fused silica core of diameter 200 |nbsp| μm and outer diameter of the
cladding 220 |nbsp| μm. Note that in this example we also use the
:py:mod:`xopto.materials.ri` module to calculate the refractive index of
fused silica at 550 |nbsp| nm and the :py:mod:`xopto.mcbase.mcutil.fiber` module for
creating a multimode fiber :py:class:`~xopto.mcbase.mcutil.fiber.MultimodeFiber`.

.. code-block:: python

    from xopto.materials.ri
    from xopto.mcml import mc
    from xopto.mcml.mcutil import fiber

    fib = fiber.MultimodeFiber(
        dcore=200e-6,
        dcladding=220e-6,
        ncore=ri.glass.fusedsilica.default(550e-9), 
        na=0.22
    )
    src = mc.mcsource.UniformFiber(fib)

Any of the source parameters can be accessed and changed through the
class properties. In the following example we change the :math:`NA` and
the diameters of the fiber core and cladding.

.. code-block:: python

    src.fiber.na = 0.23
    src.fiber.dcladding = 420e-6
    src.fiber.dcore = 400e-6

The packets are always launched from the top surface of the sample. The
fiber position and direction are used to determine the launch point at the
sample surface. If the incidence is not perpendicular, the cross-section of
the fiber still forms a tight contact with the sample surface (the fiber 
is cut at an angle and has an elliptic cross-section). From there, the packet
is launched according to the NA of the fiber core and the refractive index
of the sample. The reflectance at the fiber core-sample boundary is subtracted
from the initial packet weight. If a specular surface
:ref:`detector <mcml-detector-label>` is used, the reflectance is deposited
into that detector. 

Rectangular
-----------
Rectangular sources are similar to optical fiber sources, however with a 
rectangular emission surface:

* :py:class:`xopto.mcml.mcsource.rectangular.UniformRectangular`
  emits uniformly within the NA of the source from each point of the source surface.
* :py:class:`xopto.mcml.mcsource.rectangular.LambertianRectangular`
  emits like a Lambertian surface within the NA of the source.
* :py:class:`xopto.mcml.mcsource.rectangular.UniformRectangularLut`
  follows a nonparametric/measured angular emission characteristics.

The following example creates a 
:py:class:`~xopto.mcml.mcsource.rectangular.UniformRectangular` source of
size :math:`(x, y) = (1, 2)` |nbsp| mm, NA 0.22 and refractive index 1.452.

.. code-block:: python

    from xopto.mcml import mc

    mc.UniformRectangular(1.0e-3, 2.0e-3, n=1.452, na=0.22)
