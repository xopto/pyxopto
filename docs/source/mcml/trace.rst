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

.. _mcml-trace-label:

Trace
=======

Trace allows collection of detailed data along the photon packet path. The
trace functionality is implemented in :py:class:`xopto.mcbase.mctrace.Trace`
and also conveniently imported into the :py:mod:`xopto.mcml.mctrace` and
:py:mod:`xopto.mcml.mc` modules.

The events that are traced by :py:class:`xopto.mcbase.mctrace.Trace`: can be
customized through the :code:`options` argument.

* :python:`options=Trace.TRACE_START` - only the start state of the photon
  packets is traced (this is the state of the photon packet after launching and
  subtraction the weight due to specular reflections at the source-medium or
  sample boundary).
* :python:`options=Trace.TRACE_END` - only the end state of the photon packet
  is traced (the photon packet is absorbed, leaves the sample through the top
  or botom surface or escapes the simulation radius).
* :python:`options=Trace.TRACE_START | Trace.TRACE_END` - the start and the end
  states of the photon packet are traced.
* :python:`options=Trace.TRACE_ALL` - the entire path of the photon packet is
  traced (default).

The following events are traced (in addition to the start and end state of
the photon packet):

* refraction or reflection at the layer/voxel boundary.
* crossing the layer or voxel boundary (even if the materials of the geometrical
  entities are the same)
* scattering and absorption events

Each trace entry includes eight floating-point numbers in the following order:

* x - x coordinate of the photon packet position (m).
* y - y coordinate of the photon packet position (m).
* z - z coordinate of the photon packet position (m).
* px - x component of the propagation direction vector (m).
* py - y component of the propagation direction vector (m).
* pz - z component of the propagation direction vector (m).
* w - weight of the photon packet.
* pl - total optical path length since the start event. Note that tracing of
       the optical path length can be tuned off by passing :python:`plon=False`
       to :py:meth:`~xopto.mcbase.mctrace.Trace`.

After running a Monte carlo simulation, the trace data are stored in
a structured numpy array of shape :python:`(num_packets, max_len)`,
where :python:`num_packets` is the number of
launched photon packets and :python:`maxlen` is the maximum length of the
trace that is set through the constructor parameter :python:`maxlen`.
The actual number of trace entries for each photon packet can be
obtained from the :py:attr:`~xopto.mcbase.mctrace.Trace.n` property, which is a
numpy vector of length :python:`num_packets`.
.

When tracing the entire path of the photon packets, it is important to select
an adequate maximum length of the trace. The value will depend on a number of
factors including the optical properties of the sample and the configuration
of the investigated source and detector. It is advised to iteratively determine
the maximum length of the trace by observing the number of trace events
:py:attr:`xopto.mcbase.mctrace.Trace.n`. If the maximum trace length is
exceeded during the Monte Carlo simulation, the last trace entry is overwritten
with the final state of the photon packet. 

Trace Filter
------------

Frequently, only the traces of photon packets that end with a particular state
are useful. A versatile filter :py:class:`~xopto.mcbase.mctrace.Filter` is
available for filtering the photon packet traces according to the terminal
position, propagation direction, and optical path length of the photon packet.
The trace filter :py:class:`~xopto.mcbase.mctrace.Filter` is also conveniently
imported into the :py:mod:`xopto.mcml.mctrace` and
:py:mod:`xopto.mcvox.mctrace` modules.

The traces can be filtered by the final state of the photon packet.
For a detailed list of available options see
:py:class:`~xopto.mcbase.mctrace.Filter`. If a filter object is
passed to the :py:class:`~xopto.mcbase.mctrace.Trace` constructor (parameter
:code:`filter`) the results of the Monte Carlo simulation will be automatically
filtered. Alternatively, filters can be easily applied to the existing
:py:class:`~xopto.mcbase.mctrace.Trace` instances by the
:py:meth:`xopto.mcbase.mctrace.Filter.__call__` method. By default, the filter
modifies the trace (does not create a new trace object for the filtered results)
and returns it. This behavior can be changed by passing :python:`update=False`
to the :py:meth:`~xopto.mcbase.mctrace.Filter.__call__` method. This will create
a new :py:class:`~xopto.mcbase.mctrace.Trace` instance for the filtered
trace. Note that the filter will drop all the traces that exceed the maximum
trace length, even if they satisfy the filter. The number of dropped traces is
returned as the second output of
:py:meth:`~xopto.mcbase.mctrace.Filter.__call__`.

In the following example we create a trace with a maximum number of
events / length set to 1000 and a filter that only selects traces of photon
packets that end at the top sample surface and are escaping within
20 :superscript:`o` of the surface normal. Note that the z component of
the propagation direction is negative for photon packets that escape the
sample through the top surface.

.. code-block:: python

    from xopto.mcml import mc
    import numpy as np

    filt = mc.mctrace.Filter(
        z=(-float('inf'), 0.0)
        pz=(-1.0, -np.cos(np.deg2rad(20.0)))
    )
    mc.mctrace.Trace(maxlen=1000, filter=filt)

    # apply the filter and modify the input trace
    some_trace, dropped = filt(some_trace)

    # apply the filter and create a new trace instance
    filtered_trace, dropped = filt(some_trace, False)

To visualize the paths of photon packets after completing a simulation,
use the :py:meth:`~xopto.mcbase.mctrace.Trace.plot` method. The traces can be\
plot in several different views:

* :python:`view='3d'` produces a 3D view of the traces (default).
* :python:`view='xy'` produces a 2D view of the traces projected onto the x-y plane.
* :python:`view='xz'` produces a 2D view of the traces projected onto the x-z plane.
* :python:`view='yz'` produces a 2D view of the traces projected onto the y-z plane.

The :code:`show` parameter controls if the plot window is shown before returning
from :py:meth:`~xopto.mcbase.mctrace.Trace.plot`. Note that the starting point
of the traced path is highlighted with a green marker and the final point of
the traced path with a red marker.

.. code:: python

    trace.plot(show=False) # produces a 3D view
    trace.plot(view='xy', show=False) # produces a 2D view in the x-y plane
    trace.plot(view='xz', show=False) # produces a 2D view in the x-z plane
    trace.plot(view='yz', show=True) # produces a 2D view in the y-z plane