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

.. _mcvox-simulator-data-types-label:

.. include:: ../common.rst

Simulator data types
====================

Data types can be used to customize the OpenCL kernel of the Monte Carlo
simulator. The following 4 data types of the kernel can be customized:

#. Default integer types

    * :py:class:`~xopto.mcbase.mctypes.McInt32` - 32-bit integer
    * :py:class:`~xopto.mcbase.mctypes.McInt64` - 64-bit integer

#. Photon packet accumulator types

    * :py:class:`~xopto.mcbase.mctypes.McAccu32` - 32-bit unsigned integer
    * :py:class:`~xopto.mcbase.mctypes.McAccu64` - 64-bit unsigned integer

#. Floating point types and computational precision

    * :py:class:`~xopto.mcbase.mctypes.McFloat` - single precision floating-point
    * :py:class:`~xopto.mcbase.mctypes.McDouble` - double precision floating-point

#. Photon packet counter types

    * :py:class:`~xopto.mcbase.mctypes.McCnt32` - 32-bit unsigned integer counter
    * :py:class:`~xopto.mcbasem.mctypes.McCnt64` - 64-bit unsigned integer counter

The 4 data types are passed to the :py:meth:`xopto:mcvox.mc.Mc`
constructor as an instance of :py:class:`xopto.mcbase.mctypes.McDataTypes`.
Even though any combination of the 4 data types can be used, the preconfigured
sets of data types will be likely sufficient for the vast majority of
applications:

* :py:class:`xopto.mcbase.mctypes.McDataTypesSingle` - default set of types
  (single precision floating-point arithmetics, 32-bit photon packet counter)

    * :py:class:`~xopto.mcbase.mctypes.McInt32` - 32-bit integer
    * :py:class:`~xopto.mcbase.mctypes.McCnt32` - 32-bit unsigned photon packet counter
    * :py:class:`~xopto.mcbase.mctypes.McAccu64` - 64-bit unsigned weight accumulators
    * :py:class:`~xopto.mcbase.mctypes.McFloat` - single precision floating-point arithmetics

* :py:class:`xopto.mcbase.mctypes.McDataTypesSingleCnt64`
  (single precision floating-point arithmetics, 64-bit photon packet counter)

    * :py:class:`~xopto.mcbase.mctypes.McInt32` - 32-bit integer
    * :py:class:`~xopto.mcbase.mctypes.McCnt64` - 64-bit unsigned photon packet counter
    * :py:class:`~xopto.mcbase.mctypes.McAccu64` - 64-bit unsigned weight accumulators
    * :py:class:`~xopto.mcbase.mctypes.McFloat` - single precision floating-point arithmetics


* :py:class:`xopto.mcbase.mctypes.McDataTypesDouble`
  (double precision floating-point arithmetics, 32-bit photon packet counter)

    * :py:class:`~xopto.mcbase.mctypes.McInt32` - 32-bit integer
    * :py:class:`~xopto.mcbase.mctypes.McCnt32` - 32-bit unsigned photon packet counter
    * :py:class:`~xopto.mcbase.mctypes.McAccu64` - 64-bit unsigned weight accumulators
    * :py:class:`~xopto.mcbase.mctypes.McDouble` - double precision floating-point arithmetics

* :py:class:`xopto.mcbase.mctypes.McDataTypesDoubleCnt64`
  (double precision floating-point arithmetics, 64-bit photon packet counter)

    * :py:class:`~xopto.mcbase.mctypes.McInt32` - 32-bit integer
    * :py:class:`~xopto.mcbase.mctypes.McCnt64` - 32-bit unsigned photon packet counter
    * :py:class:`~xopto.mcbase.mctypes.McAccu64` - 64-bit unsigned weight accumulators
    * :py:class:`~xopto.mcbase.mctypes.McDouble` - double precision floating-point arithmetics

.. note::

    The 32-bit photon packet weight accumulators
    :py:class:`~xopto.mcbase.mctypes.McAccu32` should be avoided, since they can
    overflow even for simulations with a small number of photon packets.

Single precision floating-point arithmetics 
:py:class:`~xopto.mcbase.mctypes.McFloat` should be sufficient for most
application. Note that double precision floating-point arithmetics
:py:class:`~xopto.mcbase.mctypes.McDouble` is significantly slower
on GPU devices.

If more that 4,294,967,295 photon packets need to be simulated in a single run
of the Monte Carlo simulator, use a 64-bit photon packet counter
:py:class:`~xopto.mcbase.mctypes.McCnt64`. The performance penalty will depend
on the OpenCL device support for 64-bit atomic operations
(cl_khr_int64_base_atomics
OpenCL `extension <https://www.khronos.org/registry/OpenCL/sdk/1.2/docs/man/xhtml/EXTENSION.html>`_)