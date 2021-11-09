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

.. _mcvox-simulator-options-label:

.. include:: ../common.rst

Simulator options
=================

Simulator options are used to customize the OpenCL kernel. All kernel
options are derived from the :py:class:`xopto.mcbase.mcoptions.McOption` class.
The boolean options are implemented
by subclassing :py:class:`~xopto.mcbase.mcoptions.McBoolOption`,
the integer options
by subclassing :py:class:`~xopto.mcbase.mcoptions.McIntOption`,
the floating-point
options by subclassing :py:class:`~xopto.mcbase.mcoptions.McFloatOption` and
the data type options
by subclassing :py:class:`~xopto.mcbase.mcoptions.McTypeOption`.

The :py:mod:`xopto.mcbase.mcoptions` is also conveniently imported into the
:py:mod:`xopto.mcvox.mc` and :py:mod:`xopto.mcvox.mcoptions` modules.

The list of available options is as follows:

* :py:class:`~xopto.mcbase.mcoptions.McUseBallisticKernel`
  (default is :py:class:`~xopto.mcbase.mcoptions.McUseBallisticKernel.off`) -
  can be used to switch to ballistic implementation of the Monte Carlo kernel.
  The ballistic Monte Carlo kernel either scatters or fully absorbs the packet
  in each simulation step that does not intersect the boundaries of the current
  material (e.g. voxel or layer). The probability of the outcome depends on
  the values of the absorption coefficient :math:`\mu_a`, scattering
  coefficient :math:`\mu_s` and a random number :math:`\xi` sampled from a
  uniform distribution :math:`[0, 1]`:

    - if :math:`\frac{\mu_a}{\mu_a + \mu_s} \geq \xi` - fully absorb else scatter
      the packet according to the scattering phase function

  Note that in the ballistic Monte Carlo kernel the photon packets do not
  undergo termination through weight threshold and / or lottery and hence the
  values of options :py:class:`~xopto.mcbase.mcoptions.McMinimumPacketWeight`
  and :py:class:`~xopto.mcbase.mcoptions.McUseLottery` are not used by the
  ballistic kernel.
  In contrast, the regular Monte Carlo kernel terminates the photon packets
  through the mechanisms of weight threshold and lottery. In each simulation
  step that does not intersect with the boundaries of the current material, the
  weight :math:`w` of the photon packet is reduced by
  :math:`w\frac{\mu_a}{\mu_a + \mu_s}`.
  The ballistic and regular Monte Carlo kernel converge towards the same
  solution. While the ballistic kernel is faster it produces significantly
  more noisy results than the regular kernel. 

* :py:class:`~xopto.mcbase.mcoptions.McUseNativeMath`
  (default is :py:class:`~xopto.mcbase.mcoptions.McUseNativeMath.off`) -
  can be used to enable device-native math functions. Native math usually gives
  some performance benefit, but might not be fully compliant with the precision
  defined by the IEEE standards.

* :py:class:`~xopto.mcbase.mcoptions.McMaterialMemory`
  (default is :py:attr:`~xopto.mcbase.mcoptions.McMaterialMemory.global_mem`) -
  OpenCL memory type used to store the array of materials.
  Selecting :py:attr:`~xopto.mcbase.mcoptions.McIntLutMemory.constant_mem` memory
  will likely lead to a performance boost, in particular on
  older GPUs. However, note that the amount of available constant
  (:py:attr:`~xopto.mcbase.mcoptions.McIntLutMemory.constant_mem`)
  memory on GPUs is typically limited to about 64k. The constant memory is
  also used to hold some performance-critical simulator data and can be also
  used to hold lookup table data of integer or floating-point type.

* :py:class:`~xopto.mcbase.mcoptions.McIntLutMemory`
  (default is :py:attr:`~xopto.mcbase.mcoptions.McIntLutMemory.constant_mem`) -
  OpenCL memory type used to hold the integer type lookup table data.
  Selecting :py:attr:`~xopto.mcbase.mcoptions.McIntLutMemory.global_mem` memory
  will likely lead to a significant performance degradation, in particular on
  older GPUs. Note that the amount of available "constant_mem" memory on GPUs is
  typically limited to about 64k.

* :py:class:`~xopto.mcbase.mcoptions.McFloatLutMemory`
  (default is :py:attr:`~xopto.mcbase.mcoptions.McFloatLutMemory.constant_mem`) -
  OpenCL memory type used to hold the floating-point type lookup table data.
  Selecting :py:attr:`~xopto.mcbase.mcoptions.McFloatLutMemory.global_mem` memory
  will likely lead to a significant performance degradation, in particular on
  older GPUs. Note that the amount of available "constant_mem" memory on GPUs is
  typically limited to about 64 |nbsp| kB.

* :py:class:`~xopto.mcbase.mcoptions.McDebugMode`
  (default is :py:class:`~xopto.mcbase.mcoptions.McDebugMode.off`) -
  Can be used to enable kernel debug mode that will print information to
  the console output. Note that this mode requires printf functionality in
  the  OpenCL kernel that might not be supported by all the OpenCL devices.
  If this mode is turned on :py:class:`~xopto.mcbase.mcoptions.McDebugMode.on`,
  run the Monte Carlo simulator with one or a small number of photon packets.

* :py:class:`~xopto.mcbase.mcoptions.McUseSoft64Atomics`
  (default is :py:class:`~xopto.mcbase.mcoptions.McUseSoft64Atomics.off`) -
  Can be used to force a software-based implementation of 64-bit atomic
  operations.

* :py:class:`~xopto.mcbase.mcoptions.McUseLottery`
  (default is :py:class:`~xopto.mcbase.mcoptions.McUseSoft64Atomics.on`) -
  Can be used to disable termination of photon packets through lottery.

* :py:class:`~xopto.mcbase.mcoptions.McMinimumPacketWeight`
  (default is 10 :superscript:`-4`) -
  Sets the minimum photon packet weight allowed before starting packet
  termination or lottery.

* :py:class:`~xopto.mcbase.mcoptions.McPacketLotteryChance`
  (default is 0.1) -
  Terminate photon packet by lottery if the value of a uniform random
  number from [0, 1] exceeds this value.

* :py:class:`~xopto.mcbase.mcoptions.McUsePackedStructures`
  (default is :py:class:`~xopto.mcbase.mcoptions.McUsePackedStructures.off`) -
  Can be used to force the use of tightly packed structures in the OpenCL code.
  Note that packed structures can lead to significant performance
  degradation of the MonteCarlo kernel. This option is the last resort if
  the fields of the OpenCL and host structures cannot be properly aligned.
  When declaring a new OPenCL or host structure, always start the declaration
  with the largest data type and move towards smaller data types. Use data types
  that are of size no less than 4 bytes.
