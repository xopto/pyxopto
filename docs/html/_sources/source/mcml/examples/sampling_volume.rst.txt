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

.. include:: ../../common.rst

.. _example-sampling_volume:

Sampling volume simulations
===========================

This example (available in `examples/mcml/sampling_volume`) covers all the necessary steps for simulating sampling volume utilizing multimode optical fibers as sources and detectors. Sampling volume gives us an understanding of which part of the turbid medium under investigation is primarily responsible for the given detected signal (Meglinsky, I. V., and S. J. Matcher, Med. Biol. Eng. Comput., 39(1), 44-50 (2001).). This example uses a similar approach as in :ref:`example-photon-packet-tracing`, which should be read first before proceeding with this example.

Importing the required modules and submodules
---------------------------------------------

Like in the example :ref:`example-photon-packet-tracing`, we import the necessary modules. In addition, we import the submodule :py:mod:`xopto.mcml.mcutil.axis`, which comprises helper classes for accumulator axis definitions.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as pp

    from xopto.mcml import mc
    from xopto.mcml.mcutil import fiber, axis
    from xopto.cl import clinfo

.. note::

    Sampling volume simulations can be done in the voxelated Monte Carlo model by using the :py:mod:`xopto.mcvox.mc` submodule.

Computational device
--------------------

Select the desired OpenCL computational device (see :ref:`opencl-devices-label`).

.. code-block:: python

    cl_device = clinfo.gpu(platform='nvidia')

.. note::

    In this example we have selected the first computational device listed under the Nvidia platform. The string should be changed according to the installed hardware devices.ed to an available computation device. It is sufficient to provide only a unique part that pertains to the desired computational device.

The layer stack
---------------

In this example we use a single-layered turbid medium of 1 |nbsp| cm thickness. The refractive index of the single layer was set to 1.33, while the surrounding medium refractive index is set to 1.45. The absorption and scattering coefficients of the single layer are set to 2 |nbsp| 1/cm and 250 |nbsp| 1/cm. Finally, the Henyey-Greenstein phase function is utilized with anisotropy factor of 0.9 in the single layer turbid medium.

.. code-block:: python

    d = 1e-2
    layers = mc.mclayer.Layers([
        mc.mclayer.Layer(d=0.0, n=1.45, mua=0.0, mus=0.0, pf=mc.mcpf.Hg(0.0)),
        mc.mclayer.Layer(d=d, n=1.33, mua=2e2, mus=250e2, pf=mc.mcpf.Hg(0.90)),
        mc.mclayer.Layer(d=0.0, n=1.45, mua=0.0, mus=0.0, pf=mc.mcpf.Hg(0.0)),
    ])

.. note::

    All quantities MUST be provided in appropriate units, i.e., distances in m, absorption and scattering coefficients in 1/m.

Source
------

The  multimode optical fiber source was defined similarly to the example :ref:`example-photon-packet-tracing` having a core diameter :code:`dcore` of 200 |nbsp| μm, combined diameter of core and cladding :code:`dcladding` of 220 |nbsp| μm, core refractive index :code:`ncore` of 1.45 and numerical aperture :code:`na` of 0.22. The multimode optical fiber source is also displaced along the :math:`x` coordinate by -220 |nbsp| μm and positioned normally to the turbid medium.

.. code-block:: python

    source = mc.mcsource.UniformFiber(
        fiber=fiber.MultimodeFiber(
            dcore=200e-6, dcladding=220e-6, ncore=1.45, na=0.22
        ),
        direction=(0.0, 0.0, 1.0),
        position=(-220e-6, 0.0, 0.0)
    )

Detector
--------

Unlike in the example :ref:`example-photon-packet-tracing`, we define the multimode optical fiber detector at the top of the turbid medium. For the multimode optical fiber detector, we use the same properties as for the multimode optical fiber source, while the displacement along the :math:`x` coordinate is in positive direction by 220 |nbsp| μm. 

.. code-block:: python

    detector_top = mc.mcdetector.FiberArray(
        [fiber.FiberLayout(
            fiber=fiber.MultimodeFiber(
            dcore=200e-6, dcladding=220e-6, ncore=1.45, na=0.22
            ),
            position=(220e-6, 0.0, 0.0),
            direction=(0.0, 0.0, 1.0)
        ),]
    )

    detectors = mc.mcdetector.Detectors(
        top=detector_top
    )

Trace
-----

For sampling volume simulations, the photon packet traces are required. The trace object is defined similarly to the example :ref:`example-photon-packet-tracing`. The only difference is in the geometrical properties, where the photon packet traces are selected only if the photon packet is detected in the multimode optical fiber detector. As such, the final :math:`z` coordinate should be within a small interval from the upper boundary :math:`z=0`. The final direction of the detected photon packet should be within the acceptance cone corresponding to the numerical aperture 0.22. Note that the detected photon packets are already propagated over the boundary and the surrounding medium is therefore fused silica or medium with refractive index approx. 1.45. Finally, the photon packets are detected only within a circle of radius 100 |nbsp| μm, which is displaced by the same amount as the multimode optical fiber detector.

.. code-block:: python

    nphotons = 100000 
    recorded_traces = 1000
    number_of_events = 500

    trace = mc.mctrace.Trace(
        maxlen=number_of_events, options=mc.mctrace.Trace.TRACE_ALL,
        filter=mc.mctrace.Filter(
            z=(-1e-9, 1e-9), 
            pz=(-1, -np.cos(np.arcsin(0.22/1.45))), 
            r=(0, 100e-6, (220e-6, 0))
        )
    )

Sampling volume
---------------

Here we construct the sampling volume object :py:class:`~xopto.mcbase.mcsv.SamplingVolume` which allows fast processing of the photon packet traces using the OpenCL parallelization. The constructor :py:class:`~xopto.mcbase.mcsv.SamplingVolume` accepts :math:`x`, :math:`y` and :math:`z` axis objects which define accumulator voxels along each axis into which traces are scored. The accumulators along :math:`x` and :math:`y` axes in our case span from -0.75 to 0.75 mm and are split into 200 intervals. The accumulator along the :math:`z` axis span from 0 to 1 mm and is split into 200 intervals.

.. code-block:: python

    nx = 200
    ny = 200
    nz = 200
    sv = mc.mcsv.SamplingVolume(
        xaxis=mc.mcsv..Axis(-0.75e-3, 0.75e-3, nx), 
        yaxis=mc.mcsv.Axis(-0.75e-3, 0.75e-3, ny),
        zaxis=mc.mcsv.Axis(0e-3, 1e-3, nz)
    )

The Monte Carlo simulator object
--------------------------------

The Monte Carlo simulator object :py:class:`~xopto.mcml.mc.Mc` is defined similarly to the example :ref:`example-photon-packet-tracing`. In this case the termination radius is set to a bigger value to allow the photon packets to propagate longer distances.

.. code-block:: python

    mc_obj = mc.Mc(
        layers=layers, 
        source=source, 
        detectors=detectors, 
        trace=trace, 
        cl_devices=cl_device
    )
    mc_obj.rmax = 1e-1

Running the Monte Carlo simulations
-----------------------------------

Similarly as provided in the example :ref:`example-photon-packet-tracing`, the simulations are run multiple times within the :code:`while` loop in order to detect a specified number of photon packet traces :code:`recorded_traces`. Subsequently, the sampling volume is computed using the method :py:meth:`~xopto.mcml.mc.Mc.sampling_volume`, which accepts the final trace object and the constructed sampling volume object :code:`sv` from before.

.. code-block:: python

    output = mc_obj.run(nphotons, verbose=True)
    while output is None or len(output[0]) < recorded_traces:
        output = mc_obj.run(nphotons, verbose=True, out=output)
        print('Photon packet traces collected:', len(output[0]))
    trace = output[0]
    sv = mc_obj.sampling_volume(trace, sv)

Sampling volume visualization
-----------------------------

The data within the sampling volume object :code:`sv` can be accessed via the :py:attr:`~xopto.mcml.mc.mcsv.data` attribute corresponding to a 3D numpy array with the 0 axis/indices corresponding to :math:`z` axis, 1 axis/indices corresponding to :math:`y` axis and 2 axis/indices corresponding to :math:`x` axis. In this example we visualize the sampling volume summed along the :math:`y` axis.

.. code-block:: python

    pp.figure()
    pp.title('Sampling volume')
    pp.imshow(sv.data.sum(axis=1), extent=(-0.75, 0.75, 0, 1))
    pp.xlabel('x (mm)')
    pp.ylabel('z (mm)')
    pp.show()

.. image:: ../../images/examples/mcml/sampling_volume/sampling_volume.svg
    :height: 480 px
    :width: 640 px
    :scale: 75 %
    :alt: Sampling volume
    :align: center

.. note::

    For a better signal-to-noise ratio a significantly larger amount of photon packet traces that satisfy the filter should be detected. In this example we have set the number of detected photon packet traces to 1000. However, increasing this number will also prolong the simulation times.

The complete example
--------------------
    
.. literalinclude:: ../../../../examples/mcml/sampling_volume.py

You can run this example from the root directory of the PyXOpto package as:

.. code-block:: bash

    python examples/mcml/sampling_volume.py
    