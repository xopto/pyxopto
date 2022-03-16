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

.. _mcvox-opencl-build-options-label:

.. include:: ../common.rst

OpenCL build options
====================

OpenCL build options can be used to further optimize the performance of
the simulator kernel. These options can be passed to the constructor
:py:func:`xopto.mcvox.mc.Mc.__init__` by populating
the :code:`cl_build_options` argument wit a list of OpenCL build options.
The build options can be passed as strings defined by the OpenCL
`standard <https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/clBuildProgram.html>`
or by using predefined options from the :py:mod:`xopto.cl.cloptions`.

The build options are conveniently imported into
the :py:mod:`xopto.mcvox.mc` as :py:mod:`~xopto.cl.cloptions`.

For a full list of options see :py:mod:`xopto.cl.cloptions`

.. note::
    Note that the most significant performance improvment is usually obtained
    by using the :code:`'-cl-fast-relaxed-math'`
    (:py:data:`~xopto.cl.cloptions.FastRelaxedMath`) build option.`