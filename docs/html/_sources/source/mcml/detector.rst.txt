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

.. _mcml-detector-label:

Detectors
=========

Detectors are used to collect the weights of photon packets that escape the
sample at the top or bottom sample surface or to collect the weight of photon
packets that is specularly reflected when launching with the selected source.
The three detector locations can be populated and passed to the Monte Carlo
simulator through the :py:class:`xopto.mcml.mcdetector.base.Detectors`.
All the detector types that can be used in the Monte Carlo simulator
are implemented by subclassing :py:class:`xopto.mcml.mcdetector.base.Detector`.
The :py:mod:`xopto.mcml.mcdetector` module includes a number different detectors:

* :py:class:`xopto.mcml.mcdetector.total.Total` implements a 
  detector that collects the weight of photon packets into a single accumulator.

* :py:class:`xopto.mcml.mcsdetector.radial.Radial` implements a radial detector
  that collects the photon packets through concentric annular rings.

* :py:class:`xopto.mcml.mcdetector.cartesian.Cartesian` implements a 
  Cartesian detector that collects the photon packets through a grid of
  rectangular detectors.

* :py:class:`xopto.mcml.mcdetector.symmetric.SymmetricX` implements a 
  Cartesian detector that is symmetric across the x axis (uses the absolute
  value of y coordinate).

* :py:mod:`xopto.mcml.mcdetector.probe` implements several several standard
  and custom parameterized layouts of optical fibers in optical fiber probes.

    * :py:class:`xopto.mcml.mcdetector.probe.sixaroundone.SixAroundOne` implements
      the common six-around-one layout of optical fibers.
    * :py:class:`xopto.mcml.mcdetector.probe.lineararray.LinearArray` implements
      a linear layout of optical fibers, assuming that all the fibers are equal.
    * :py:class:`xopto.mcml.mcdetector.probe.fiberarray.FiberArray` implements
      a general layout of optical fibers, where each optical fiber can have
      different properties.

The individual detectors are conveniently imported into the
:py:mod:`xopto.mcml.mc` and :py:mod:`xopto.mcml.mcdetector` modules.

Radial
------

The following example shows how to create a
:py:class:`~xopto.mcml.mcdetector.radial.Radial` detector centered at
:math:`(x,y)=(0,0)` that collects photon packets from :math:`r=0` |nbsp| mm
to :math:`r=1` |nbsp| mm through concentric annular rings of thickness
10 |nbsp| μm. Note that the distribution of detectors along the radial
direction is defined by :py:class:`~xopto.mcbase.mcutil.axis.Axis`.

.. code-block:: python

    from xopto.mcml import mc

    det = mc.mcdetector.Radial(mc.mcdetector.Axis(0.0, 1.0e-3, 100))

To use a logarithmic spacing of annular rings, set the input parameter
:code:`logscale=True`.

.. code-block:: python

    from xopto.mcml import mc

    det = mc.mcdetector.Radial(mc.mcdetector.Axis(0.0, 1.0e-3, 100, logscale=True))

The acceptance cone of the detector can be limited by specifying the minimum
acceptance angle cosine. The angle is by default computed relative to the
top or bottom sample surface normal. This behavior can be customized by
setting the reference direction through the :code:`direction` input argument
to a custom unit-length vector. Note that the acceptance cone is applied after
the photon packet leaves the sample and is refracted into the surrounding
medium.

The following example creates a radial detector that collects photon packets
within ±10 |nbsp| :superscript:`o` of sample surface normal.

To use a logarithmic spacing of annular rings, set the input parameter
:code:`logscale=True`.

.. code-block:: python

    from xopto.mcml import mc
    import numpy as np

    cosmin = np.cos(np.deg2rad(10.0))
    det = mc.mcdetector.Radial(
        mc.mcdetector.Axis(0.0, 1.0e-3, 100, logscale=True),
        cosmin=cosmin
    )

Cartesian
---------

In the next example, we create a
:py:class:`~xopto.mcml.mcdetector.cartesian.Cartesian` detector that collects
photon packets through a grid of rectangular accumulators that span 
:math:`[-1, 1]` |nbsp| mm along the x axis and :math:`[-2, 2]` |nbsp| mm along
the y axis. The size of the rectangular accumulators is
:math:`(x,y)=(10,10)` |nbsp| μm. Note that parameters :code:`cosmin` and
:code:`direction` can be optionally used to limit the acceptance cone of the
detector.

.. code-block:: python

    from xopto.mcml import mc

    det = mc.mcdetector.Cartesian(
        xaxis = mc.mcdetector.Axis(-1.0, 1.0, 200),
        yaxis = mc.mcdetector.Axis(-2.0, 2.0, 400),
    )

If a Cartesian detector is created with a single axis, the axis configuration
is applied to the x and y axis.

.. code-block:: python

    from xopto.mcml import mc

    det = mc.mcdetector.Cartesian(mc.mcdetector.Axis(-1.0, 1.0, 200))

Total
-----

To collect the weights of photon packets that leave the sample into a single
accumulator, use the :py:class:`xopto.mcml.mcdetector.total.Total` detector.
Note that parameters :code:`cosmin` and :code:`direction` can be optionally
used to limit the acceptance cone of the detector.

.. code-block:: python

    from xopto.mcml import mc

    det = mc.mcdetector.Total()

Symmetric
---------

To collect photon packets symmetrically across the x axis, regardless of the
y coordinate use the :py:class:`xopto.mcml.mcdetector.symmetric.SymmetricX`
detector. This type of a detector is useful for sources that are tilted along
the x axis and produce a reflectance or transmittance response that is
symmetric along the x axis and the location along the y axis is not important.
The distribution of accumulators is defined through
:py:class:`xopto.mcbase.mcutil.axis.SymmetricAxis` that also supports logarithmic
spacing of the accumulators. Note that parameters :code:`cosmin` and
:code:`direction` can be optionally used to limit the acceptance cone of the
detector.
In the following example we create an instance of
:py:class:`~xopto.mcml.mcdetector.symmetric.SymmetricX` detector that spans the
range :math:`[-10,10]` |nbsp| mm along the axis with the accumulator size along
the set to 10 |nbsp| μm.

.. code-block:: python

    from xopto.mcml import mc

    det = mc.mcdetector.SymmetricX(mc.mcdetector.SymmetricAxis(0.0, 10.0, 1000))

Probe
-----

Probe detectors that utilize various layouts of optical fibers are intended for
use with complex surface layouts that break the radial symmetry of the
reflectance / transmittance at the sample surface. Consequently,
the :py:class:`~xopto.mcml.mcdetector.radial.Radial` detector cannot be used to
derive the reflectance collected through the individual fibers.

In the following example we create a six-around-one layout of optical fibers.
All the optical fibers have a fused silica core with a 400 |nbsp| μm diameter,
cladding diameter 420 |nbsp| μm, :math:`NA` 0.22 and are tightly
packed, i.e. the distance between the cores of the central and the surrounding
optical fibers is 420 |nbsp| μm. The refractive index of the fused silica is
taken from the :py:mod:`xopto.materials.ri` module.

.. code:: python

    from xopto.mcml import mc
    from xopto.mcml.mcutil import fiber
    from xopto.materials import ri

    fib = fiber.MultimodeFiber(
        400e-6,
        420e-6,
        ncore=ri.glass.fusedsilica.default(550e-9),
        na=0.22
    )
    detector = mc.mcdetector.SixAroundOne(fib)

Note that the geometry of the layout (:code:`spacing` of fibers) and the
properties of the optical fiber can be changed through accessing the related
class properties.

.. code:: python

    detector.spacing = 800e-6
    detector.fiber.na = 0.23
    detector.fiber.dcore = 200e-6
    detector.fiber.dcladding = 220e-6
    detector.fiber.ncore = ri.glass.fusedsilica.default(400e-9)

Use in Monte Carlo simulator
----------------------------

To use the detectors in Monte Carlo simulations, we need to populate an instance
of :py:class:`xopto.mcml.mcdetector.base.Detectors`, that can independently set
the detectors at the top and bottom sample surfaces and the detector for
photon packets that are specularly reflected at the medium-sample boundary
when launched. Note that it is not required to populate all the three detectors
(top, bottom and specular) as illustrated in this example.

.. code-block:: python

    from xopto.mcml import mc

    top = mc.mcdetector.Radial(mc.mcdetector.Axis(0.0, 1.0, 100))
    bottom = mc.mcdetector.Cartesian(mc.mcdetector.Axis(-1.0, 1.0, 200))
    specular = mc.mcdetector.Total()

    detectors = mc.mcdetector.Detectors(top=top, bottom=bottom, specular=specular)

After completing the Monte Carlo simulations, the returned instance of
:py:class:`~xopto.mcml.mcdetector.base.Detectors` can be used to access the
collected reflectance (also used for transmittance) or the accumulated raw
weight of the photon packets. Note that the reflectance / transmittance is
computed as the raw weight in the accumulator divided by the surface area of
the accumulator (weight/m :superscript:`2`) .

.. code-block:: python

    top_reflectance = detectors.top.reflectance
    top_raw weight = detectors.top.raw

    bottom_reflectance = detectors.bottom.reflectance
    bottom_raw weight = detectors.bottom.raw

    specular_reflectance = detectors.specular.reflectance
    specular_raw weight = detectors.specular.raw

