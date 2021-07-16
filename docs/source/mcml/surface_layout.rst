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

.. _mcml-surface-label:

.. include:: ../common.rst

Surface layouts
===============

Surface layouts can be used to describe complex interfaces at the top and
bottom sample surfaces, such as the surface of the optical fiber probe tip.
The surface layouts are subclasses of 
:py:class:`xopto.mcml.mcsurface.base.SurfaceLayoutBase`. Note that the top
and bottom sample surface can be assigned a different surface layout.

Readily available surface layots are conveniently imported into the
:py:mod:`xopto.mcml.mc` and :py:mod:`xopto.mcml.mcsurface` modules and include:

* :py:mod:`xopto.mcml.mcsurface.lambertian.LambertianReflector`
  implements a uniform surface with optical properties varying
  from an ideal first surface mirror to a an ideal Lambertian surface with
  an optional absorption.

* The surface layouts of optical fiber probes are implemented in:

    * 'py:mod:`xopto.mcml.mcsurface.probe.sixaroundone.SixAroundOne` implements a
      parameterized layout of a six-around-one optical fiber probe.
    * 'py:mod:`mcml.mcsurface.probe.lineararray.LinearArray` implements a
      parameterized layout of probe with linear placement of optical fibers.
    * 'py:mod:`mcml.mcsurface.probe.fiberarray.FiberArray` implements a general
      parameterized layout of differnt optical fibers in an optical fiber probe.

In this example we create a layout of a tightly packed six-around-one optical
fiber probe with 7 multimode fused silica fibers. The multimode optical fibers
have a core diameter of 400 |nbsp| μm, a cladding diameter of 420 |nbsp| μm
and a numerical aperture 0.22. The stainless
steel tip of the optical fiber probe has a 6 |nbsp| mm diameter and a
reflectivity of 65%. Note that we utilize the
:py:mod:`xopto.mcbase.mcutil.fiber` module to create an instance of
:py:class:`~xopto.mcbase.mcutil.fiber.MultimodeFiber` that represents a
multimode fiber.

.. code:: python

    from xopto.mcml.mcutil import fiber
    from xopto.mcml import mc

    fiber.MultimodeFiber(
        400e-6,
        420e-6,
        ncore=ri.glass.fusedsilica.default(550e-9)
    )

    mc.mcsurface.SixAroundOne(fiber, diameter=6.0e-3, reflectivity=0.65)

Use in Monte Carlo simulator
----------------------------

The complex layouts at the top and / or bottom sample surface are passed to the
Monte Carlo simulator through
:py:class:`xopto.mcml.mcsurface.base.SurfaceLayouts`. In the following example
we set the layout at the top sample surface to a six-around-one probe
:py:class:`~xopto.mcml.mcsurface.probe.sixaroundone.SixAroundOne` and the
layout at the bottom sample surface to a Lambertian reflector
:py:class:`~xopto.mcml.mcsurface.lambertian.LambertianReflector` that absorbs
30% of the incident light and reflects the remaining 70% as an ideal Lambertian
reflector. Note that it is not required to populate  both the top and bottom
sample surface layouts as illustrated in this example.

.. code-block:: python

    from xopto.mcml.mcutil import fiber
    from xopto.mcml import mc

    fiber.MultimodeFiber(
        400e-6,
        420e-6,
        ncore=ri.glass.fusedsilica.default(550e-9)
    )

    layouts = mc.mcsurface.SurfaceLayouts(
        top=mc.mcsurface.SixAroundOne(fiber, diameter=6.0e-3, reflectivity=0.65),
        bottom=mc.mcsurface.LambertianReflector(reflectance=0.7) 
    )

