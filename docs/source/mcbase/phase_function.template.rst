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

.. _{{ MC }}-pf-label:

.. include:: ../common.rst


Scattering phase functions
==========================

All the scattering phase functions that can be used in the Monte Carlo
simulator are implemented by subclassing
:py:class:`xopto.mcbase.mcpf.pfbase.PfBase`. The 
:py:mod:`xopto.mcbase.mcpf` module includes a number of commonly
used scattering phase functions:

* :py:class:`xopto.mcbase.mcpf.hg.Hg` implements the Henyey-Greenstein
  scattering phase function

* :py:class:`xopto.mcbase.mcpf.hg.Hga` implements the anisotropic
  Henyey-Greenstein scattering phase function

* :py:class:`xopto.mcbase.mcpf.hg.HgDir` implements the directional
  Henyey-Greenstein scattering phase function

* :py:class:`xopto.mcbase.mcpf.hg2.Hg2` implements a dual Henyey-Greenstein
  scattering phase function, i.e. a linear combination of two 
  :py:class:`~xopto.mcbase.mcpf.hg.Hg` scattering phase functions where
  one is scattering in the forward and the other in the backward direction

* :py:class:`xopto.mcbase.mcpf.mhg.MHg` implements the Modified 
  Henyey-Greenstein scattering phase function

* :py:class:`xopto.mcbase.mcpf.gk.Gk` implements the Gegenbauer kernel
  scattering phase function

* :py:class:`xopto.mcbase.mcpf.gk2.Gk2` implements a dual Gegenbauer kernel
  scattering phase function, i.e. a linear combination of two 
  :py:class:`~xopto.mcbase.mcpf.gk.Gk` scattering phase functions where
  one is scattering in the forward and the other in the backward direction

* :py:class:`xopto.mcbase.mcpf.mgk.MGk` implements the modified
  Gegenbauer kernel scattering phase function

* :py:class:`xopto.mcbase.mcpf.pc.Pc` implements the Power of
  cosine scattering phase function

* :py:class:`xopto.mcbase.mcpf.mpc.MPc` implements the Modified Power of
  cosine scattering phase function

* :py:class:`xopto.mcbase.mcpf.rayleigh.Rayleigh` implements the isotropic Rayleigh
  scattering phase function

* :py:class:`xopto.mcbase.mcpf.lut.Lut` and :py:class:`xopto.mcbase.mcpf.lut.LutEx`
  implement the lookup table-based scattering phase functions that need to
  be used when there is no analytical inverse of the cumulative probability
  density function. In this case, the scattering angle cannot be sampled
  analytically, instead an efficient numerical sampling scheme is required.

  .. note::

    The use of :py:class:`~xopto.mcbase.mcpf.lut.Lut`
    and :py:class:`~xopto.mcbase.mcpf.lut.LutEx` is tightly coupled to the
    :py:mod:`xopto.pf` module that implements numerous scattering phase functions
    with a number of numerical utilities for computing quantifiers such as
    :math:`\gamma`, :math:`\delta`, :math:`\sigma` and Legendre moments,
    computing the scattering cross sections of spherical particles, angular
    distribution of scattering probabilities, working with monodisperse,
    or various standard and custom size distributions, layered spherical
    particles, etc.  

The individual scattering phase functions are conveniently imported into
the :py:mod:`xopto.{{ MC }}.mc` and :py:mod:`xopto.{{ MC }}.mcpf` modules.

Analytical sampling
-------------------

The following example shows how to create a Henyey-Greenstein scattering
phase function with the anisotropy :math:`g` set to 0.8:

.. code-block:: python

    from xopto.{{ MC }} import mc

    pf = mc.mcpf.Hg(0.8)

The parameters of the scattering phase function can be updated at any time
through the class properties. The following examples changes the value
of parameter :code:`g` from 0.8 to 0.9.

.. code-block:: python

    pf.g = 0.9

The number and naming of parameters depends on the type of the scattering phase
function. The Gegenbauer kernel scattering phase function
:py:class:`xopto.mcbase.mcpf.gk.Gk` exposes two parameters,
namely :code:`gg` and :code:'a'.

.. code-block:: python

    pf = mc.mcpf.Gk(gg=0.9, a=0.5)

Any of the two parameters can be accessed through the class properties:

.. code-block:: python

    pf.gg = 0.7
    pf.a = 0.4

Lookup table-based numerical sampling
-------------------------------------

Complex scattering phase functions, such as the scattering phase functions of
spherical particles that can be computed with the Mie theory, can be used with
the Monte Carlo simulator through the lookup table scattering phase functions
:py:class:`xopto.mcbase.mcpf.lut.Lut` or
:py:class:`xopto.mcbase.mcpf.lut.LutEx`. In the first step we need to create
an instance of a scattering phase function that enables computation of various
scattering phase function quantifiers, such as the Legendre moments,
:math:`\gamma`, :math:`\delta`, :math:`\sigma`, etc. These
can be found in the :py:mod:`xopto.pf` that conveniently imports all
the implemented scattering phase functions.

The following example shows how to create a Monte Carlo simulator-compatible
scattering phase function :py:meth:`xopto.mcbase.mcpf.lut.Lut` for spherical
polystyrene particles of diameter 1.0 |nbsp| μm suspended in water and for
550 |nbsp| nm light. Note that we utilize the :py:mod:`xopto.materials` package
from which we import the refractive index module :py:mod:`xopto.materials.ri`.

.. code-block:: python

    from xopto import pf
    from xopto.{{ MC }} import mc
    from xopto.materials import ri

    mie = pf.Mie(ri.polystyrene.default(550e-9), ri.water.default(550e-9), 1.0e-6, 550e-9)

    mc_mie = mc.mcpf.Lut(*mie.mclut())

The same result can be accomplished by :py:meth:`xopto.mcbase.mcpf.lut.LutEx`
that takes the scattering phase function type and a list of arguments for the 
corresponding constructor (parameters of the scattering phase function).

.. code-block:: python

    from xopto import pf
    from xopto.{{ MC }} import mc
    from xopto.materials import ri

    params = [ri.polystyrene.default(550e-9), ri.water.default(550e-9), 1.0e-6, 550e-9]

    mc_mie = mc.mcpf.LutEx(pf.Mie, params)

The default lookup table size is set to 2000, which should yield an accurate
representation for the vast majority of the scattering phase functions.
However, for scattering phase functions with an extremely high anisotropy that
exceeds 0.95, a larger lookup table size might be required.
The size of the lookup table can be controlled with the :code:`lutsize`
parameter.

.. code-block:: python

    from xopto import pf
    from xopto.{{ MC }} import mc
    from xopto.materials import ri

    params = [ri.polystyrene.default(550e-9), ri.water.default(550e-9), 1.0e-6, 550e-9]
    mc_mie = mc.mcpf.LutEx(pf.Mie, params, lutsize=4000)

    mie = pf.Mie(ri.polystyrene.default(550e-9), ri.water.default(550e-9), 1.0e-6, 550e-9)
    mc_mie = mc.mcpf.Lut(*mie.mclut(lutsize=4000))

For scattering phase functions that come in a nonparametric form, such as when
measured with a goniometer, use the :py:mod:`xopto.pf.discrete.Discrete` that
can take values defined at discrete scattering angles. Then follow the above
examples to obtain a Monte Carlo simulator-compatible scattering phase
function with
:py:meth:`~xopto.mcbase.mcpf.lut.Lut` or
:py:meth:`~xopto.mcbase.mcpf.lut.LutEx`.

.. note::
  Any scattering phase function in :py:mod:`xopto.pf` can
  be converted into a lookup table-based Monte Carlo simulator-compatible
  scattering phase function.

The following example crates a lookup table-based
implementation of the Henyey-Greenstein (:py:class:`~xopto.mcbase.mcpf.mhg.Hg`)
scattering phase function with anisotropy :math:`g=0.8`:

.. code-block:: python

    from xopto import pf
    from xopto.{{ MC }} import mc

    params = [0.8]
    hg_lut = mc.mcpf.LutEx(pf.Hg, params)

The use of lookup-table based scattering phase functions in the Monte Carlo
simulations will not have a notable performance impact as long as all the
lookup table data can be kept in the constant memory of the OpenCL device,
the size of which is for a typical GPU around 64 |nbsp| kB.
