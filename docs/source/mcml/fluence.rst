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
absorption. The fluence detector is implemented in the
:py:class:`xopto.mcbase.mcfluence.Fluence` class. The observed volume
voxelization is axis-aligned and allows independent adjustment of the voxel
size along the x, y and z axis.

The :py:class:`~xopto.mcbase.mcfluence.Fluence` is also conveniently
imported into the :py:mod:`xopto.mcml.mcfluence` and
:py:mod:`xopto.mcvox.mcfluence` modules.

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
:py:class:`~xopto.mcbase.mcfluence.Fluence.xaxis`,
:py:class:`~xopto.mcbase.mcfluence.Fluence.yaxis` and
:py:class:`~xopto.mcbase.mcfluence.Fluence.zaxis` properties. Likewise, the
centers of the voxels along the x, y and z axis can be accessed through the
:py:class:`~xopto.mcbase.mcfluence.Fluence.x`,
:py:class:`~xopto.mcbase.mcfluence.Fluence.y` and
:py:class:`~xopto.mcbase.mcfluence.Fluence.z` properties.

.. code-block:: python

    xaxis = fluence.xaxis
    yaxis = fluence.yaxis
    zaxis = fluence.zaxis

    x_centers = fluence.x    # or as xaxis.centers
    y_centers = fluence.y    # or as yaxis.centers
    z_centers = fluence.z    # or as zaxis.centers

The fluence data are stored in a 3D numpy array that can be accessed through the
:py:attr:`~xopto.mcbase.mcfluence.Fluence.data` property. The individual voxels
are addressed / sliced as :python:`data[z_slice, y_slice, x_slice]`. Note that
the values in the fluence data array are the unscaled sums of
deposited / absorbed photon packet weights.

After completing the Monte Carlo simulation, the fluence accumulator can be
viewed using the :py:meth:`~xopto.mcbase.mcfluence.Fluence.plot` method. The
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

    fluence.plot(axis='z')    # slice along the z axis
    fluence.plot(axis='y')    # slice along the y axis
    fluence.plot(axis='x')    # slice along the x axis

    fluence.plot(axis='zproj')    # integral projection / sum along the z axis
    fluence.plot(axis='yproj')    # integral projection / sum along the y axis
    fluence.plot(axis='xproj')    # integral projection / sum along the x axis