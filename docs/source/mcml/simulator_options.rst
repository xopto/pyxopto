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

.. _mcml-simulator-options-label:

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
:py:mod:`xopto.mcml.mc` and :py:mod:`xopto.mcml.mcoptions` modules.

The list of available options is as follows:

* :py:class:`~xopto.mcbase.mcoptions.McMethod`
  (default is :py:class:`~xopto.mcbase.mcoptions.McMethod.albedo_weight`) -
  can be used to set the Monte Carlo method. A detailed description of the
  methods can be found in:

  #. A. Sassaroli and F. Martelli, *Equivalence of four Monte Carlo*
     *methods for photon migration in turbid media*,
     J Opt Soc Am A Opt Image Sci Vis, **29** (10), 2110-2117 (2012),
     https://doi.org/10.1364/JOSAA.29.002110.

  Three different photon packet stepping methods are available:

    - Albedo Rejection
      (:py:class:`~xopto.mcbase.mcoptions.McMethod.albedo_rejection`
      or :py:class:`~xopto.mcbase.mcoptions.McMethod.ar`)

      Propagation step :math:`s=-\frac{\ln(\xi_1)}{\mu_t}` is derived
      from a uniform random number :math:`\xi_1` from interval
      :math:`[0, 1]` and the total attenuation coefficient
      :math:`\mu_t=\mu_a + \mu_s`, which is the sum of the absorption
      :math:`\mu_a` and scattering coefficients :math:`\mu_s`.
      The packet is fully absorbed if no boundaries are
      hit along the step :math:`s` and a uniform random number
      :math:`\xi_2` from :math:`[0, 1]` fulfils
      :math:`\xi_2 \leq \frac{\mu_a}{\mu_t}`.  If the
      packet hits geometry boundaries, it is propagated to the boundary,
      where the boundary interactions are processed and a new step is
      started afterwards.

    - Albedo weight
      (:py:class:`~xopto.mcbase.mcoptions.McMethod.albedo_weight`
      or :py:class:`~xopto.mcbase.mcoptions.McMethod.aw`)

      Propagation step :math:`s=-\frac{\ln(\xi_1)}{\mu_t}` is derived
      from a uniform random number :math:`\xi_1` from interval
      :math:`[0, 1]` and the total attenuation coefficient
      :math:`\mu_t=\mu_a + \mu_s`, which is the sum of the absorption
      :math:`\mu_a` and scattering coefficients :math:`\mu_s`.
      The packet is scattered and partially absorbed if no boundaries are
      hit along the step :math:`s`. The fraction of the absorbed
      weight is computed as :math:`\frac{\mu_a}{\mu_t}`. If the
      packet hits geometry boundaries, it is propagated to the boundary,
      where the boundary interactions are processed and a new step is
      started afterwards.

    - Microscopic Beer-Lambert
      (:py:class:`~xopto.mcbase.mcoptions.McMethod.microscopic_beer_lambert`
      or :py:class:`~xopto.mcbase.mcoptions.McMethod.mbl`)

      Propagation step :math:`s=-\frac{\ln(\xi_1)}{\mu_s}` is derived
      from a uniform random number :math:`\xi_1` from interval
      :math:`[0, 1]` and the scattering coefficient :math:`\mu_s`.
      The packet is absorbed regardless if the boundaries are
      hit along the step.  If the packet hits a geometry boundary,
      it is propagated to the boundary, where the boundary
      interactions are processed and a new step is
      started afterwards. The fraction of the absorbed
      weight is computed as :math:`1 - \exp^{-\mu_a s_b}`, where
      :math:`s_b` is the distance to the geometry boundary or the
      full step if the boundary is not hit.

  The Albedo Rejection implementation of the Monte Carlo method is the fastest
  but produces noisy results. The Microscopic Beer-Lambert method produces
  higher quality Fluence or Energy Deposition simulations if the size of
  the voxels in the deposition / fluence grid is smaller than the
  mean free path of the packets, i.e. there is on average less than one
  absorption event per deposition/fluence voxel. The default Albedo Weight
  method is recommended for reflectance simulations and matches the
  performance of the Microscopic Beer-Lambert method if the mean free path
  of the packets is equal or smaller than the voxel size of the
  deposition / fluence grid.

* :py:class:`~xopto.mcbase.mcoptions.McUseNativeMath`
  (default is :py:class:`~xopto.mcbase.mcoptions.McUseNativeMath.off`) -
  can be used to enable device-native math functions. Native math usually gives
  some performance benefit, but might not be fully compliant with precisions
  defined by the IEEE standards.

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
  the  OpenCL kernel that might not be supported by all OpenCL devices.
  If this mode is turned on :py:class:`~xopto.mcbase.mcoptions.McDebugMode.on`,
  run the Monte Carlo simulator with one or a small number of photon packets.

* :py:class:`~xopto.mcbase.mcoptions.McUseSoft64Atomics`
  (default is :py:class:`~xopto.mcbase.mcoptions.McUseSoft64Atomics.off`) -
  Can be used to force software implementation of 64-bit atomic operations.

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
  When declaring OPenCL or host structures always start with the
  largest data type and move towards smaller data types. Use data types
  that are of size no less than 4 bytes.




  

  