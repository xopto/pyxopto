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

.. _mcvox-sampling-volume-label:

.. include:: ../common.rst

Sampling volume
===============

Sampling volume is a tool for visualizing the path of photon packets through the
sample for a particular source and detector configuration. The analysis is
based on the traces / paths of photon packets recorded by the Monte Carlo
simulations. Briefly, each sampling volume accumulator that intersects with
the photon packet path is added a value that is the product of the terminal
photon packet weight and the length traveled in that particular accumulator.
For more details on how to run simulations that record / trace the paths of
photon packets see :ref:`Trace <mcvox-trace-label>`.

The recorded traces are returned by the Monte carlo simulator as a
:py:class:`xopto.mcbase.mctrace.Trace` object that can be passed to
the :py:meth:`xopto.mcvox.mc.Mc.sampling_volume` method of the Monte Carlo
simulator to process the traces. The traces are efficiently processed in
parallel using an OpenCL kernel.
The sampling volume utilizes axis-aligned voxelized 3D accumulators implemented
in :py:class:`xopto.mcbase.mcsv.SamplingVolume` and 
conveniently imported into the :py:mod:`xopto.mcvox.mcsv` and
:py:mod:`xopto.mcvox.mc` modules. 
The range and size of voxels along the main axes is set with three
:py:class:`~xopto.mcbase.mcuti.axis.Axis` objects. The sampling volume data are
stored in a 3D numpy array that can be accessed through the
:py:attr:`~xopto.mcbase.mcsv.SamplingVolume.data` property. The individual
voxels are addressed / sliced as :python:`data[z_slise, y_slice, x_slice]`.

In the following example we voxelize a volume that extends along the
x axis from -1 to 1 |nbsp| mm, along the y axis from -1 to 1 |nbsp| mm and along
the z axis from 0 to 2 |nbsp| mm. The size of voxels along all the three axes is
set to 0.01 |nbsp| mm.

.. code-block:: python

    from xopto.mcvox import mc

    sv = mc.mcsv.SamplingVolume(
        xaxis=mc.mcsv.Axis(-1.0e-3, 1.0e-3, 200),
        yaxis=mc.mcsv.Axis(-1.0e-3, 1.0e-3, 200),
        zaxis=mc.mcsv.Axis(0.0, 2.0e-3, 200)
    )

Given that we have an instance of the Monte Carlo simulator
:py:class:`xopto.mcvox.mc.Mc` (:code:`mc_obj`) and a trace
(:code:`trace_obj`) that was produced by this simulator instance, a call
to the :py:meth:`~xopto.mcvox.mc.Mc.sampling_volume` of the simulator
instance will process the trace and return the result as a 
:py:class:`~xopto.mcbase.mcsv.SamplingVolume` instance.

.. code:: python

    sv_res = mc.sampling_volume(trace, sv)

After processing the traces, the voxelized sampling volume accumulator can be
viewed using the :py:meth:`~xopto.mcbase.mcsv.SamplingVolume.plot` method. The
method allows slicing across the voxelized volume along one of the main
coordinate axis. Use the :code:`axis` parameter to set the coordinate
axis:

* :python:`axis='x'` slice along the x axis,
* :python:`axis=y` slice along the y axis,
* :python:`axis='z'` slice along the z axis.

To show integral projections along one of the main coordinate axis call with:

* :python:`axis='xproj'` for an integral projection along the x axis,
* :python:`axis=yproj` for an integral projection along the y axis,
* :python:`axis='zproj'` for an integral projection along the z axis.

The slices can be viewed in a linear or logarithmic scale. Call with
:python:`scale='lin'` to show the slices in a linear intensity scale.
By default the color coding of each slice is automatically scaled to the range
of weights (:python:`autoscale=True`). 

.. code-block:: python

    sv_res.plot(axis='z')    # slice along the z axis
    sv_res.plot(axis='y')    # slice along the y axis
    sv_res.plot(axis='x')    # slice along the x axis

    sv_res.plot(axis='zproj')    # integral projection / sum along the z axis
    sv_res.plot(axis='yproj')    # integral projection / sum along the y axis
    sv_res.plot(axis='xproj')    # integral projection / sum along the x axis