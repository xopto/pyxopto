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

.. _mccyl-tracing-example-label:

Photon packet tracing
=====================

In this example we trace photon packets of a collimated Gaussian beam that
passes through an absorbing glass cylinder.

This example is available in the `examples/mccyl/tracing.py` file.

The functionality of cylindrical Monte Carlo (MC) is available through the
:py:mod:`xopto.mccyl.mc` module.

.. code-block:: python

    from xopto.mccyl import mc

The layer stack
---------------

The concentric layers of the medium can be defined through the
:py:mod:`xopto.mccyl.mclayer` submodule. The concentric layers are stacked
from the outermost to the innermost layer and centered at :math:`(x,y)=(0,0)`.

The outermost layer of the stack is used to describe the
medium that surrounds the sample. Therefore, at least two layers must be always
defined, namely the layer of the surrounding medium and one sample layer!

The diameter of the outermost layers will be automatically set to infinity
regardless of the specified layer diameter.

Note that all the layers in the stack must use the same scattering phase
function model. A variety of scattering phase function models is available
through the :py:mod:`xopto.mccyl.mcpf` submodule.

The outer diameter of the cylinder is set to 10 |nbsp| mm and the
refractive index to 1.45, the absorption coefficient
mua=1.0 |nbsp| cm :superscript:`-1` and the scattering coefficient to
mua=0.0 |nbsp| cm :superscript:`-1`

The glass cylinder is surrounded by a 20 |nbsp| mm thick layer of air with a
arefractive index of 1.0.

The sourrounding medium is also set to air with a arefractive index of 1.0.

.. code-block:: python

    layers = mc.mclayer.Layers(
        [
            mc.mclayer.Layer(n=1.00, d=np.inf,  mua=0.0e2, mus=0.0e2,  pf=mc.mcpf.Hg(0.0))
            mc.mclayer.Layer(n=1.00, d=50e-3,  mua=0.0e2, mus=0.0e2,  pf=mc.mcpf.Hg(0.0))
            mc.mclayer.Layer(n=1.45, d=10e-3,   mua=0.0e2, mus=0.0e2,  pf=mc.mcpf.Hg(0.0)),
        ]
    )

Note that the absorption coefficient :code:`mua`, scattering coefficient
:code:`mus` and the scattering phase function :code:`pf` of the
outermost layer are not used in the MC simulations, since the
photon packets are not propagated through the surrounding medium.
However, the refractive index :code:`n` of the two outermost
layer is used to properly refract/reflect the photon packet at the sample
surface when launched by the source or when escaping the sample.
The layer diameter :code:`d` should be given in m and the values of the
scattering coefficient :code:`mus` and absorption coefficient
:code:`mua` in m :superscript:`-1`.

The photon packet source
------------------------

Different sources of photon packets are available through the
:py:mod:`xopto.mccyl.mcsource` module. In this example we use a
collimated Gaussian source that enters the sample in the direction of the
:math:`x` axis. The standard deviation of the beam is set to
:math:`\sigma=1` |nbsp| mm.

.. code-block:: python

    source = mc.mcsource.GaussianBeam(1e-3)

The trace
------------

The traces of photon packets can be collected by a trace detector
:py:class:`~xopto.mcbase.mctrace.Trace` that is
available through the :py:mod:`xopto.mccyl.mctrace` module. The maximum
number of events that can be traced for each photon packet is limited to 128.

.. code-block:: python

    trace = mc.mctrace.Trace(128)

The OpenCL device
-----------------

The OpenCL device that will run the MC simulations can be selected through the
:py:mod:`xopto.cl.clinfo` module. In the following example we pick the first
available GPU device.

.. code-block:: python

    gpu = clinfo.gpu()


The Monte Carlo simulator
-------------------------

Next, we create a Monte Carlo simulator instance from the created layer stack,
photon packet source and trace detector.

.. code-block:: python

    mc_obj = mc.Mc(layers, source, None, trace, cl_devices=gpu)

Optionally, we can limit the maximum simulation radius that is measured from
the position of the photon packet source. In this example, we limit the
simulation radius to 100 |nbsp| mm, which exceeds the radius of the outer
sample boundary (50 |nbsp| mm)

.. code-block:: python

    mc_obj.rmax = 1000.0e-3

Finally, we can run the simulator instance with a given number of photon
packets (100 in this example) and collect the results. The simulator
returns three objects/results, namely the trace, fluence and detectors. Since
in this example we only use the trace, the remaining two
results (detectors and fluence) will be returned as :python:`None`. 

.. code-block:: python

    trace_res, fluence_res, detectors_res = mc_obj.run(100)

Visualization of results
------------------------

We can plot the simulation results using the
`matplotlib.pyplot <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.html>`_
module. 

The traces can be plot using the :py:meth:`xopto.mcbase.mctrace.Trace.plot` method.
In this example we plot the traces in the x-y plane, which can be selected
through the :python:`view` argument. The results are plot onto an existing
axis that is set through the :python:`ax` parameter. Finally, we also plot
the glass cylinder and the outer boundary of the sample, set the aspect ration
of the x and y axis to equal, add a legend and show the plot window.

.. code-block:: python

    # crate a plot window and a plot axis
    fig, ax = pp.subplots()

    # plot the packet traces projected onto the x-y plane
    trace_res.plot(view='xy', ax=ax)
    # plot the glass cylindr
    cylinder = pp.Circle((0.0, 0.0), 5.0e-3, color='blue', alpha=0.25,
                        fill=True, label='glass cylinder', zorder=2.0)
    ax.add_patch(cylinder)
    # plot the outer boundary of the sample
    detector = pp.Circle((0.0, 0.0), 25e-3, color='black',
                        fill=False, label='sample boundary')
    ax.add_patch(detector)

    # add axis labels
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')

    # make aspect ratio of the x and y axis equal
    ax.set_aspect('equal')

    # show plot legend
    pp.legend()

    # display the plot window
    pp.show()

Collected traces of the photon packets.

.. figure:: ../../images/examples/mccyl/tracing/trace.svg
    :align: center
    :width: 100 %

The complete example
--------------------

.. literalinclude:: ../../../../examples/mccyl/tracing.py

This example can be run from the root directory of the PyXOpto package as:

.. code-block:: bash

    python examples/mccyl/tracing.py
