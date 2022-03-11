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

.. _mccyl-detector-label:

Detectors
=========

Detectors are used to collect the weights of photon packets that escape the
sample or to collect the weight of photon packets that are specularly
reflected when launching with the selected source.
The the two detector locations can be populated and passed to the Monte Carlo
simulator through the :py:class:`xopto.mccyl.mcdetector.base.Detectors`.
All the detector types that can be used in the Monte Carlo simulator
are implemented by subclassing :py:class:`xopto.mccyl.mcdetector.base.Detector`.
The :py:mod:`xopto.mccyl.mcdetector` module includes two different detectors:

* :py:class:`xopto.mccyl.mcdetector.total.Total` implements a 
  detector that collects the weight of photon packets into a single accumulator.

* :py:class:`xopto.mccyl.mcsdetector.fiz.FiZ` implements a 2D detector
  that collects the photon packets according to their polar angle
  :math:`\varphi` and :math:`z` location when leaving the sample and entering
  the surrounding medium`.

The individual detectors are conveniently imported into the
:py:mod:`xopto.mccyl.mc` and :py:mod:`xopto.mccyl.mcdetector` modules.

FiZ
------

The following example shows how to create a
:py:class:`~xopto.mccyl.mcdetector.fiz.FiZ` detector
that collects photon packets with polar angle :math:`\varphi` from -180° to
180° (1° step) and :math:`z` coordinate from -1 to 1 |nbsp| mm
(10 |nbsp| μm step). Distribution of accumulators along the two axis
is defined by :py:class:`~xopto.mcbase.mcutil.axis.Axis`.

.. code-block:: python

    from xopto.mccyl import mc
    import numpy as np

    det = mc.mcdetector.FiZ(
      fi=mc.mcdetector.Axis(np.pi, np.pi, 360),
      z=mc.mcdetector.Axis(-1.0e-3, 1.0e-3, 200)
    )

.. note::

    The valid range of polar angle :math:`\varphi` is from :math:`-\pi`
    (-180°) to :math:`\pi` (180°). The first and last accumulator along the
    polar :math:`\varphi` and :math:`z` axis also collect the weight of all the
    photon packets that exceed the lower and upper bounds of the range.

Total
-----

To collect the weights of photon packets that leave the sample into a single
accumulator, use the :py:class:`xopto.mccyl.mcdetector.total.Total` detector.
Note that parameters :code:`cosmin` and :code:`direction` can be optionally
used to limit the acceptance cone of the detector.

.. code-block:: python

    from xopto.mccyl import mc

    det = mc.mcdetector.Total()


Use in Monte Carlo simulator
----------------------------

To use the detectors in Monte Carlo simulations, we need to populate an instance
of :py:class:`xopto.mccyl.mcdetector.base.Detectors`, that can independently set
the detectors at the outer sample surface and the detector for
photon packets that are specularly reflected at the medium-sample boundary
when launched. Note that it is not required to populate the two detectors
(outer and specular) as illustrated in this example.

.. code-block:: python

    from xopto.mccyl import mc

    outer = mc.mcdetector.FiZ(
      fi=mc.mcdetector.Axis(-np.pi, np.pi, 360),
      z=mc.mcdetector.Axis(-1.0e-3, 1.0e-3, 200)
    )
    specular = mc.mcdetector.Total()

    detectors = mc.mcdetector.Detectors(outer=outer, specular=specular)

After completing the Monte Carlo simulations, the returned instance of
:py:class:`~xopto.mccyl.mcdetector.base.Detectors` can be used to access the
collected reflectance (also used for transmittance) or the accumulated raw
weight of the photon packets.

.. note::

    The reflectance / transmittance is computed as the raw weight in the
    accumulator divided by the number of launched photon packets and if
    applicable also divided by the surface area of the accumulator
    (weight/m :superscript:`2`).

.. code-block:: python

    outer_reflectance = detectors.outer.reflectance
    outer weight = detectors.outer.raw

    specular_reflectance = detectors.specular.reflectance
    specular_raw weight = detectors.specular.raw

