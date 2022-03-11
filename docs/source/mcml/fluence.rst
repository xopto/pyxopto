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

.. _mcml-fluence-label:

.. include:: ../common.rst

Fluence
=======

Fluence detector is used to obtain the spatial distribution of the
energy (packet weights) deposited into the sample through the process of
absorption or optionally the fluence rate.
Three different implementations of the fluence detector are
available in the :py:mod:`xopto.mcbase.mcfluence` module:

* :py:class:`~xopto.mcbase.mcfluence.fluencerz.FluenceRz` - Radially symmetric implementation :math:`(r, z)`.
* :py:class:`~xopto.mcbase.mcfluence.fluencecyl.FluenceCyl` - Implementation in cylindrical coordinates :math:`(r, \varphi, z)`. 
* :py:class:`~xopto.mcbase.mcfluence.fluence.Fluence` - Implementation in Cartesian coordinates :math:`(x, y, z)`.

The following are time-resolved variants of the above fluence detectors:

* :py:class:`~xopto.mcbase.mcfluence.fluencerz.FluenceRzt` - Implementation in Cartesian coordinates wit additional temporal axis :math:`(x, y, z, t)`.
* :py:class:`~xopto.mcbase.mcfluence.fluencecyl.FluenceCylt` - Implementation in Cartesian coordinates wit additional temporal axis :math:`(x, y, z, t)`.
* :py:class:`~xopto.mcbase.mcfluence.fluencet.Fluencet` - Implementation in Cartesian coordinates wit additional temporal axis :math:`(x, y, z, t)`.

The observed volume voxelization is axis-aligned and allows independent
adjustment of the voxel size along the x, y and z axis for
:py:class:`~xopto.mcbase.mcfluence.fluence.Fluence`, along r, φ and z axis for .
:py:class:`~xopto.mcbase.mcfluence.fluencecyl.FluenceCyl` and along the r and z axis for
:py:class:`~xopto.mcbase.mcfluence.fluencerz.FluenceRz`.
The time-resolved variants of the fluence detectors have an additional
temporal (t) axis. The time along the temporal axis is measured in seconds.

The fluence accumulators utilize NumPy arrays and follow C-indexing order,
where the last index represents the most rapidly changing memory location:

* :py:class:`~xopto.mcbase.mcfluence.fluencerz.FluenceRz` - :code:`[z, r]`
* :py:class:`~xopto.mcbase.mcfluence.fluencecyl.FluenceCyl` - :code:`[z, φ, r]`
* :py:class:`~xopto.mcbase.mcfluence.fluence.Fluence` - :code:`[z, y, x]`
* :py:class:`~xopto.mcbase.mcfluence.fluencerz.FluenceRzt` - :code:`[z, r, t]`
* :py:class:`~xopto.mcbase.mcfluence.fluencecyl.FluenceCylt` - :code:`[z, φ, r, t]`
* :py:class:`~xopto.mcbase.mcfluence.fluencet.Fluencet` - :code:`[z, y, x, t]`


The :py:class:`~xopto.mcbase.mcfluence.fluence.Fluence`,
:py:class:`~xopto.mcbase.mcfluence.fluencet.Fluencet`,
:py:class:`~xopto.mcbase.mcfluence.fluencecyl.FluenceCyl`,
:py:class:`~xopto.mcbase.mcfluence.fluencecyl.FluenceCylt`,
:py:class:`~xopto.mcbase.mcfluence<.fluencerz.FluenceRz` and
:py:class:`~xopto.mcbase.mcfluence<.fluencerz.FluenceRzt`
classes are also conveniently
imported into the :py:mod:`xopto.mcml.mcfluence`, :py:mod:`xopto.mccyl.mcfluence`
:py:mod:`xopto.mcvox.mcfluence` modules.

The fluence detector allows collection of unscaled deposited / absorbed
packet weights or of the fluence rate where the deposited packet weights are
scaled by the inverse of the absorption coefficient. The collection method can
be selected through the :python:`mode` argument that should be set
to :python:`'deposition'` (default) for the raw weights and :python:`'fluence'`
for the fluence rate.

In the following example we voxelize a volume that extends along the
x axis from -5 to 5 |nbsp| mm, along the y axis from -5 to 5 |nbsp| mm and along
the z axis from 0 to 10 |nbsp| mm. The voxel size along all the three axis is
set to 0.1 |nbsp| mm.

.. code-block:: python

    from xopto.mcml import mc

    fluence = mc.mcfluence.Fluence(
        xaxis=mc.mcfluence.Axis(-5.0e-3, 5.0e-3, 100),
        yaxis=mc.mcfluence.Axis(-5.0e-3, 5.0e-3, 100),
        zaxis=mc.mcfluence.Axis(0.0, 10.0e-3, 100)
    )

The :py:class:`~xopto.mcbase.mcutil.axis.Axis` instances along the x, y and z
axis can be accessed through the 
:py:class:`~xopto.mcbase.mcfluence.fluence.Fluence.xaxis`,
:py:class:`~xopto.mcbase.mcfluence.fluence.Fluence.yaxis` and
:py:class:`~xopto.mcbase.mcfluence.fluence.Fluence.zaxis` properties. Likewise, the
centers of the voxels along the x, y and z axis can be accessed through the
:py:class:`~xopto.mcbase.mcfluence.fluence.Fluence.x`,
:py:class:`~xopto.mcbase.mcfluence.fluence.Fluence.y` and
:py:class:`~xopto.mcbase.mcfluence.fluence.Fluence.z` properties.

.. code-block:: python

    xaxis = fluence.xaxis
    yaxis = fluence.yaxis
    zaxis = fluence.zaxis

    x_centers = fluence.x    # or as xaxis.centers
    y_centers = fluence.y    # or as yaxis.centers
    z_centers = fluence.z    # or as zaxis.centers

The fluence data are stored in a 3D numpy array that can be accessed through the
:py:attr:`~xopto.mcbase.mcfluence.fluence.Fluence.raw` or
:py:attr:`~xopto.mcbase.mcfluence.fluence.Fluence.data` properties. The individual
voxels are addressed / sliced as :python:`data[z_slice, y_slice, x_slice]`.

.. note::

    The values returned by the :py:attr:`~xopto.mcbase.mcfluence.fluence.Fluence.raw`
    property are the unscaled sums of deposited / absorbed photon packet weights
    while the values returned by the :py:attr:`~xopto.mcbase.mcfluence.fluence.Fluence.data`
    property are scaled by the inverse of the product between the voxel volume and
    the number of launched packets. The time-resolved variants of fluence
    detectors scale the values returned by
    :py:attr:`~xopto.mcbase.mcfluence.fluencet.Fluencet.data`
    with the inverse of the product between the voxel volume, the number of
    launched packets and the corresponding time interval. 
    Also note that the voxel size is constant
    for the :py:class:`~xopto.mcbase.mcfluence.fluence.Fluence` but changes with location
    for the :py:class:`~xopto.mcbase.mcfluence.fluencerz.FluenceRz` and
    :py:class:`~xopto.mcbase.mcfluence.fluencecyl.FluenceCyl`.


After completing the Monte Carlo simulation, the fluence accumulator can be
viewed using the :py:meth:`~xopto.mcbase.mcfluence.fluence.Fluence.plot` method. The
method allows slicing across the voxelized volume along one of the main
coordinate axis. Use the :code:`axis` parameter to set the coordinate
axis:

* :python:`axis='x'` slice along the x axis,
* :python:`axis='y'` slice along the y axis,
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

    fluence.plot(axis='z')    # slice along the z axis
    fluence.plot(axis='y')    # slice along the y axis
    fluence.plot(axis='x')    # slice along the x axis

    fluence.plot(axis='zproj')    # integral projection / sum along the z axis
    fluence.plot(axis='yproj')    # integral projection / sum along the y axis
    fluence.plot(axis='xproj')    # integral projection / sum along the x axis