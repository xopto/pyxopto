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

.. _mccyl-surface-label:

.. include:: ../common.rst

Surface layouts
===============

Surface layouts can be used to describe complex interfaces at the outer and
inner sample surfaces, such as thin absorbing or reflective coatings that break
the radial symmetry of the layer stack.
The surface layouts are subclasses of 
:py:class:`xopto.mccyl.mcsurface.base.SurfaceLayoutBase`. Note that the outer
and internal layer surfaces can be assigned a different surface layout. The
provided OpenCL code will be executed each time the photon packet hits the
outer sample surface or any of the internal sample layers.

The surface layuts are passed to a Simulator constructor
:py:meth:`xopto.mccyl.mc.Mc.__init__` using the :python:`surface` argument,
which can take an instance of the
:py:class:`xopto.mccyl.mcsurface.base.SurfaceLayouts` class.


