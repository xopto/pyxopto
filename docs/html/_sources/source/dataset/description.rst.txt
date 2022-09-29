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

.. _dataset-description-label:

Summary
=======

Dataset modules contain templates for rendering reference simulation scripts
for a number of different sample, source and detector configurations. The
reference scripts can be rendered by invoking the appropriate functions or
by running the modules. There are several parameters that can be used to
control the rendering process when running the modules from a command line:

* :sh:`-o`
    Option can be used to specify a target directory for the rendered
    scripts and the related datasets. The current working directory is used by
    default. The scripts will be rendered into :sh:`run` subdirectory and the
    dataset produced by the reference simulation scripts will be saved into
    the :sh:`data` directory. 

* :sh:`-d`
    Option can be used to specify a default OpenCL device that will
    be used to run the simulations. This value is matched against the vendor
    and device name (or part of it). The first available OpenCL-enabled GPU
    device is used by default. Note that the target OpenCL device can be also
    selected by setting the :sh:`CL_DEVICE` environment variable before running
    the rendered script.

* :sh:`-i`
    Option can be used to specify the index of the OpenCL device to use.
    This option is useful when multiple identical devices are installed.
    Note that the index is zero-based and defaults to 0. Note that the target
    OpenCL device index can be also selected by setting the :sh:`CL_INDEX`
    environment variable before running the rendered script.

* :sh:`-v`
    Flag enables verbose output.

* :sh:`-t`
    Flag performs a test run and implicitly enables verbose output.
    Scripts will be rendered but no files will be written to the target
    directory. The output will show the full paths of the generated scripts.

The following example will render all the available scripts into the
current working directory and produce a verbose output:

.. code-block:: sh

    python -m xopto.dataset.render.all -v

Details
=======

The following sections provide details on different reference datasets and
instructions for rendering / generating scripts that will compute these
datasets. At the end of each section there is a breakdown of the directory
tree and naming conventions that are used to organize the datasets. 

MCML comparison datasets
------------------------

Module :py:mod:`xopto.dataset.render.mcml_comparison` => :sh:`run/mcml_comparison/`

Renders scripts that will produce datasets for comparison with the MCML package.
These datasets are based on
:py:class:`~xopto.mcbase.mcpf.hg.Hg` scattering phase function, a
:py:class:`~xopto.mcml.mcsource.line.Line` source, a
:py:class:`~xopto.mcml.mcdetector.radial.Radial` reflectance and
transmittance detector and a :py:class:`~xopto.mcbase.mcfluence.FluenceRz`
fluence / deposition detector. The refractive index of the surrounding medium
is set to 1.0. The simulation termination radius is set to 1000 mm and each
simulation is run with 100 million photon packets.
The optical properties of the top sample layer are sapled from
:math:`\mu_a = [0.0, 2.5, 5.0]` |nbsp| cm :sup:`-1`,
:math:`\mu_s' = [5.0, 20.0, 35.0]` |nbsp| cm :sup:`-1` and
:math:`g = [0.1, 0.5, 0.9]`.


.. list-table:: Source, reflectance/transmittance detector and fluence/deposition detector configurations
    :header-rows: 1
    :widths: 50, 50, 100, 100

    * -
      - Parameter
      - Value
      - Description
    * - :py:class:`~xopto.mcml.mcsource.line.Line`
      -
      -
      - normal incidence
    * - :py:class:`~xopto.mcml.mcdetector.radial.Radial`
      - :math:`raxis`
      - :py:class:`~xopto.mcbase.mcutil.axis.Axis` (start=0, stop=0.005, n=500)
      - reflectance and transmittance accumulators
    * - :py:class:`~xopto.mcbase.mcfluence.FluenceRz`
      - :math:`raxis`
      - :py:class:`~xopto.mcbase.mcutil.axis.Axis` (start=0, stop=0.005, n=500)
      - radial :math:`r` axis
    * -
      - :math:`zaxis`
      - :py:class:`~xopto.mcbase.mcutil.axis.Axis` (start=0, stop=0.005, n=500)
      - :math:`z` axis
    * -
      - :math:`mode`
      - 'deposition'
      - weight deposition mode

* :sh:`run/mcml_comparison/1-layer-100mm.py`
    A single 100 |nbsp| mm thick sample layer.

    .. list-table:: A single layer sample of 100 |nbsp| mm thickness
        :header-rows: 1
        :widths: 50, 50, 100, 100

        * -
          - Parameter
          - Value
          - Description
        * - :py:class:`~xopto.mcml.mclayer.layer.Layer`
          - :math:`d`
          - 100 |nbsp| mm
          - thickness
        * -
          - :math:`\mu_a`
          - 0.0, 2.5 or 5.0 |nbsp| cm :sup:`-1`
          - absorption coefficient
        * -
          - :math:`\mu_s'`
          - 5.0, 20.0 or 35.0 |nbsp| cm :sup:`-1`
          - reduced scattering coefficient
        * -
          - :math:`n`
          - 1.337
          - refractive index
        * -
          - :math:`g`
          - 0.1, 0.5 or 0.9
          - :py:class:`~xopto.mcbase.mcpf.hg.Hg` (g)

* :sh:`run/mcml_comparison/1-layer-1mm.py`
    A single 1 |nbsp| mm thick sample layer.

    .. list-table:: A single layer sample of 1 |nbsp| mm thickness
        :header-rows: 1
        :widths: 50, 50, 100, 100

        * -
          - Parameter
          - Value
          - Description
        * - :py:class:`~xopto.mcml.mclayer.layer.Layer`
          - :math:`d`
          - 1 |nbsp| mm
          - thickness
        * -
          - :math:`\mu_a`
          - 0.0, 2.5 or 5.0 |nbsp| cm :sup:`-1`
          - absorption coefficient
        * -
          - :math:`\mu_s'`
          - 5.0, 20.0 or 35.0 |nbsp| cm :sup:`-1`
          - reduced scattering coefficient
        * -
          - :math:`n`
          - 1.337
          - refractive index
        * -
          - :math:`g`
          - 0.1, 0.5 or 0.9
          - :py:class:`~xopto.mcbase.mcpf.hg.Hg` (g)

* :sh:`run/mcml_comparison/2-layer-100um-1mm.py`
    Two sample layers, a 0.1 |nbsp| mm top and a 1.0 |nbsp| mm bottom layer.
    The optical properties are varied only in the top sample layer
    (the refractive index is fixed to 1.462). The optical properties of
    the bottom sample layer are fixed to
    :math:`\mu_a=0.5` |nbsp| cm :sup:`-1`,
    :math:`\mu_s'=20.0` |nbsp| cm :sup:`-1`, :math:`g=0.8` and
    :math:`n=1.337`.

    .. list-table:: A 2-layer sample
        :header-rows: 1
        :widths: 50, 50, 100, 100

        * -
          - Parameter
          - Value
          - Description
        * - :py:class:`~xopto.mcml.mclayer.layer.Layer`
          - :math:`d`
          - 0.1 |nbsp| mm
          - thickness
        * -
          - :math:`\mu_a`
          - 0.0, 2.5 or 5.0 |nbsp| cm :sup:`-1`
          - absorption coefficient
        * -
          - :math:`\mu_s'`
          - 5.0, 20.0 or 35.0 |nbsp| cm :sup:`-1`
          - reduced scattering coefficient
        * -
          - :math:`n`
          - 1.462
          - refractive index
        * -
          - :math:`g`
          - 0.1, 0.5 or 0.9
          - :py:class:`~xopto.mcbase.mcpf.hg.Hg` (g)
        * - :py:class:`~xopto.mcml.mclayer.layer.Layer`
          - :math:`d`
          - 1.0 |nbsp| mm
          - thickness
        * -
          - :math:`\mu_a`
          - 0.5 |nbsp| cm :sup:`-1`
          - absorption coefficient
        * -
          - :math:`\mu_s'`
          - 20 |nbsp| cm :sup:`-1`
          - reduced scattering coefficient
        * -
          - :math:`n`
          - 1.337
          - refractive index
        * -
          - :math:`g`
          - 0.8
          - :py:class:`~xopto.mcbase.mcpf.hg.Hg` (g)

**Dataset files**

Executing the generated scripts will save the simulation results / datasets
into compressed numpy data files that will be organized as follows:

.. code-block:: sh

    data/mcml_comparison/<sample>/line/radial/hg/g-<g>/mua-<mua>-musr-<musr>-invcm.npz

The values of placeholders <> are as follows:

* :sh:`<sample>` can take the following values:

    * :sh:`1-layer-1mm`
      A single layer 1 |nbsp| mm thick medium.

    * :sh:`1-layer-100mm`
      A single layer 100 |nbsp| mm thick medium.

    * :sh:`2-layer-100um-1mm`
      A two-layer medium with 0.1 |nbsp| top layer and 1 |nbsp| mm bottom layer

* :sh:`<g>` is the anisotropy formatted with two decimal digits and :sh:`_`
  as the decimal separator, e.g :sh:`0_10` for :math:`g=0.1`, 
  :sh:`0_50` for :math:`g=0.5` or :sh:`0_90` for :math:`g=0.9`.

* :sh:`<mua>` is the absorption coefficient in units of cm :sup:`-1` with
  two decimal digits and :sh:`_` as the decimal separator, e.g :sh:`2_50`
  for :math:`\mu_a=2.5` |nbsp| cm :sup:`-1` .

* :sh:`<musr>` is the reduced scattering coefficient in units of cm :sup:`-1`
  with two decimal digits and :sh:`_` as a decimal separator, e.g :sh:`20_00`
  for :math:`\mu_a=20.0` |nbsp| cm :sup:`-1` .

Layered media datasets
----------------------

Module :py:mod:`xopto.dataset.render.mcml` => :sh:`run/mcml/`

Renders scripts that produce datasets for the layered MC kernel. These will
include a number of different source, detector and sample scattering
phase function configurations. The simulation termination radius is set
to 25 mm, except for the spatial frequency domain(SFD) dataset that 
uses a 150 mm simulation termination radius. The datasets are produced with
100 million photon packets, except for the SFD dataset and all the datasets
that use an optical fiber probe with a linear layout of 6 fibers. These
datasets are run with 1000 million photon packets.
The optical properties of the sample are varied according to the values
in the following two tables.

.. list-table:: Absorption and reduced scattering coefficients
    :header-rows: 1
    :widths: 50, 50, 50, 50

    * - Parameter
      - From (cm :sup:`-1`)
      - To (cm :sup:`-1`)
      - Points
    * - :math:`\mu_a`
      - 0.0
      - 5.0
      - 21
    * - :math:`\mu_s'`
      - 5.0
      - 35.0
      - 21


.. list-table:: Scattering phase functions
    :header-rows: 1
    :widths: 50, 50, 100

    * -
      - Parameter
      - Values
    * - :py:class:`~xopto.mcbase.mcpf.hg.Hg` :math:`(g)`
      - :math:`g`
      - 0.1, 0.3, 0.5, 0.7, 0.9
    * - :py:class:`~xopto.mcbase.mcpf.mhg.MHg` :math:`(g, \beta)`
      - :math:`g`
      - 0.1, 0.3, 0.5, 0.7, 0.9
    * - 
      - :math:`\beta`
      - 0.0, 0.2, 0.4, 0.6, 0.8, 1.0
    * - :py:class:`~xopto.mcbase.mcpf.gk.Gk` :math:`(g, \alpha)`
      - :math:`g`
      - 0.1, 0.3, 0.5, 0.7, 0.9
    * -
      - :math:`\alpha`
      - -0.5, 0.0, 0.5, 1.0, 1.5
    * - :py:class:`~xopto.pf.mie.Mie` :math:`(n_1, n_2, d, \lambda)`
      - :math:`n_1`
      - 1.337
    * -
      - :math:`n_2`
      - 1.462 (*fused silica*)
    * -
      - :math:`diameter`
      - 0.25, 0.5, 1.0, 2.0, 4.0 |nbsp| μm
    * -
      - :math:`wavelength`
      - 500 |nbsp| nm
    * - :py:class:`~xopto.pf.mie.Mie` :math:`(n_1, n_2, d, \lambda)`
      - :math:`n_1`
      - 1.337
    * -
      - :math:`n_2`
      - 1.603 (*polystyrene*)
    * -
      - :math:`diameter`
      - 0.25, 0.5, 1.0, 2.0, 4.0 μm
    * -
      - :math:`wavelength`
      - 500 |nbsp| nm

The refractive index of the sample is set to 1.337.

Reference simulation scripts are rendered for the following basic sources
that use a laterally uniform boundary between the sample and the
surrounding medium.

.. list-table:: Basic sources with a uniform sample-source interface and related reflectance detectors
    :header-rows: 1
    :widths: 50, 50, 50, 100

    * - Source
      - Parameter
      - Value
      - Reflectance detector
    * - :py:class:`~xopto.mcml.mcsurce.line.Line`
      -
      -
      - :py:class:`~xopto.mcml.mcdetector.radial.Radial` (:py:class:`~xopto.mcbase.mcutil.axis.Axis` (start=0, stop=0.005, n=500))
    * - :py:class:`~xopto.mcml.mcsource.uniformbeam.UniformBeam`
      - :math:`diameter`
      - 200 |nbsp| μm
      - :py:class:`~xopto.mcml.mcdetector.radial.Radial` (:py:class:`~xopto.mcbase.mcutil.axis.Axis` (start=0, stop=0.005, n=500))
    * - :py:class:`~xopto.mcml.mcsource.gaussianbeam.GaussianBeam`
      - :math:`FWHM`
      - 100 |nbsp| μm
      - :py:class:`~xopto.mcml.mcdetector.radial.Radial` (:py:class:`~xopto.mcbase.mcutil.axis.Axis` (start=0, stop=0.005, n=500))
    * - :py:class:`~xopto.mcml.mcsource.fiber.UniformFiber`
      - :math:`dcore`
      - 200 |nbsp| μm
      - :py:class:`~xopto.mcml.mcdetector.radial.Radial` (:py:class:`~xopto.mcbase.mcutil.axis.Axis` (start=0, stop=0.005, n=500), cosmin=0.98637)
    * -
      - :math:`dcladding`
      - 220 |nbsp| μm
      -
    * -
      - :math:`ncore`
      - 1.462
      -
    * -
      - :math:`na`
      - 0.22
      -

The refractive index of the surrounding medium is set to 1.0 except when
using the :py:class:`~xopto.mcml.mcsource.fiber.UniformFiber` source,
when the refractive index of the surrounding medium follows the
refractive index of the fiber core 1.462.

The reflectance of basic sources is collected with a radial detector
:py:class:`~xopto.mcml.mcdetector.radial.Radial` |nbsp|
(:py:class:`~xopto.mcbase.mcutil.axis.Axis` (start=0, stop=0.005, n=500))
that extends from 0 to 5 |nbsp| mm with 5oo concentric ring accumulators,
each 5 |nbsp| μm wide.
The acceptance angle is unlimited, except for the 
:py:class:`~xopto.mcml.mcsource.fiber.UniformFiber` source, for which it is
limited by the NA of the fiber core. The acceptance angle within the
fiber core is computed as :math:`\cos \theta = \sqrt{1 - (NA/n_{core})^2}`.

Reference simulation scripts are also rendered for optical fiber probe
sources that use a surface layout to more accurately describe the interface
between the optical fiber probe tip and the sample.
All the probe sources launch the photon packets with the
:py:class:`~xopto.mcml.mcsource.fiber.UniformFiber` source.

.. list-table:: Optical fiber probe sources with a detailed sample-source interface and related reflectance detectors
    :header-rows: 1
    :widths: 50, 50, 50, 100, 100

    * - Probe
      - Parameter
      - Value
      - Description
      - Reflectance detector
    * - :py:class:`~xopto.mcml.mcsurface.probe.sixaroundone.SixAroundOne`
      - :math:`dcore`
      - 200 |nbsp| μm
      - six-around-one layout
      - :py:class:`~xopto.mcml.mcdetector.probe.sixaroundone.SixAroundOne`
    * -
      - :math:`dcladding`
      - 220 |nbsp| μm
      -
      -
    * -
      - :math:`ncore`
      - 1.462
      -
      -
    * -
      - :math:`na`
      - 0.22
      -
      -
    * -
      - :math:`spacing`
      - 220 |nbsp| μm
      -
      -
    * -
      - :math:`diameter`
      - 6.0 mm
      -
      -
    * -
      - :math:`reflectivity`
      - 0.6
      -
      -
    * - :py:class:`~xopto.mcml.mcsurface.probe.sixaroundone.SixAroundOne`
      - :math:`dcore`
      - 400 |nbsp| μm
      - six-around-one layout
      - :py:class:`~xopto.mcml.mcdetector.probe.sixaroundone.SixAroundOne`
    * -
      - :math:`dcladding`
      - 420 |nbsp| μm
      -
      -
    * -
      - :math:`ncore`
      - 1.462
      -
      -
    * -
      - :math:`na`
      - 0.22
      -
      -
    * -
      - :math:`spacing`
      - 420 |nbsp| μm
      -
      -
    * -
      - :math:`diameter`
      - 6.0 mm
      -
      -
    * -
      - :math:`reflectivity`
      - 0.6
      -
      -
    * - :py:class:`~xopto.mcml.mcsurface.probe.lineararray.LinearArray`
      - :math:`dcore`
      - 200 |nbsp| μm
      - linear layout of 6 fibers
      - :py:class:`~xopto.mcml.mcdetector.probe.lineararray.LinearArray`
    * -
      - :math:`dcladding`
      - 220 |nbsp| μm
      -
      -
    * -
      - :math:`ncore`
      - 1.462
      -
      -
    * -
      - :math:`na`
      - 0.22
      -
      -
    * -
      - :math:`n`
      - 6
      -
      -
    * -
      - :math:`spacing`
      - 220 |nbsp| μm
      -
      -
    * -
      - :math:`diameter`
      - 6.0 mm
      -
      -
    * -
      - :math:`reflectivity`
      - 0.6
      -
      -
    * - :py:class:`~xopto.mcml.mcsurface.probe.lineararray.LinearArray`
      - :math:`dcore`
      - 100 |nbsp| μm
      - single fiber layout
      - :py:class:`~xopto.mcml.mcdetector.probe.lineararray.LinearArray`
    * -
      - :math:`dcladding`
      - 120 |nbsp| μm
      -
      -
    * -
      - :math:`ncore`
      - 1.462
      -
      -
    * -
      - :math:`na`
      - 0.22
      -
      -
    * -
      - :math:`n`
      - 1
      -
      -
    * -
      - :math:`diameter`
      - 6.0 mm
      -
      -
    * -
      - :math:`reflectivity`
      - 0.6
      -
      -
    * - :py:class:`~xopto.mcml.mcsurface.probe.lineararray.LinearArray`
      - :math:`dcore`
      - 200 |nbsp| μm
      - single fiber layout
      - :py:class:`~xopto.mcml.mcdetector.probe.lineararray.LinearArray`
    * -
      - :math:`dcladding`
      - 220 |nbsp| μm
      -
      -
    * -
      - :math:`ncore`
      - 1.462
      -
      -
    * -
      - :math:`na`
      - 0.22
      -
      -
    * -
      - :math:`n`
      - 1
      -
      -
    * -
      - :math:`diameter`
      - 6.0 mm
      -
      -
    * -
      - :math:`reflectivity`
      - 0.6
      -
      -
    * - :py:class:`~xopto.mcml.mcsurface.probe.lineararray.LinearArray`
      - :math:`dcore`
      - 400 |nbsp| μm
      - single fiber layout
      - :py:class:`~xopto.mcml.mcdetector.probe.lineararray.LinearArray`
    * -
      - :math:`dcladding`
      - 420 |nbsp| μm
      -
      -
    * -
      - :math:`ncore`
      - 1.462
      -
      -
    * -
      - :math:`na`
      - 0.22
      -
      -
    * -
      - :math:`n`
      - 1
      -
      -
    * -
      - :math:`diameter`
      - 6.0 mm
      -
      -
    * -
      - :math:`reflectivity`
      - 0.6
      -
      -
    * - :py:class:`~xopto.mcml.mcsurface.probe.lineararray.LinearArray`
      - :math:`dcore`
      - 800 |nbsp| μm
      - single fiber layout
      - :py:class:`~xopto.mcml.mcdetector.probe.lineararray.LinearArray`
    * -
      - :math:`dcladding`
      - 820 |nbsp| μm
      -
      -
    * -
      - :math:`ncore`
      - 1.462
      -
      -
    * -
      - :math:`na`
      - 0.22
      -
      -
    * -
      - :math:`n`
      - 1
      -
      -
    * -
      - :math:`diameter`
      - 6.0 mm
      -
      -
    * -
      - :math:`reflectivity`
      - 0.6
      -
      -

The reflectance of optical fiber probe sources is collected only through
the individual optical fibers of the probe.

**SFD source-detector arrangements**

The SFD datasets are computed for two source-detector configurations and
include raw reflectance and the corresponding frequency-domain reflectance,
which is computed for spatial frequencies from 0.00 to 0.80 |nbsp| mm :sup:`-1`
in 0.01 |nbsp| mm :sup:`-1` steps:

* A normally incident :py:class:`~xopto.mcml.mcsource.line.Line` source and a
  radial :py:class:`~xopto.mcml.mcdetector.radial.Radial` |nbsp| 
  (:py:class:`~xopto.mcbase.mcutil.axis.Axis` |nbsp|
  (start=0, stop=0.15, n=4000, logscale=True), cosmin=0.98481)
  reflectance detector that uses 4000 logarithmically spaced concentric
  accumulators from 0 to 150 |nbsp| mm.
  The acceptance angle is limited to 10°. Hankel transform is used to
  compute the spatial frequency-domain reflectance. Note that this transform
  produces real values. 

* A normally incident :py:class:`~xopto.mcml.mcsource.line.Line` source and a
  tilted linear detector with 20° incidence (along the :math:`x` axis).
  The accumulators of the detector extend to infinity along the positive
  and negative :math:`y` axis and follow a logarithmic spacing along the
  positive and negative direction of the :math:`x` axis
  :py:class:`~xopto.mcml.mcdetector.symmetric.SymmetricX` |nbsp|
  (:py:class:`~xopto.mcbase.mcutil.axis.SymmetricAxis` |nbsp|
  (center=0, range=0.15, n_half=4000, logscale=True), cosmin=0.98480).
  The described detector uses 8000 (4000 in each direction along the
  :math:`x` axis) logarithmically spaced accumulators from :math:`x=-150`
  to :math:`x=150` |nbsp| mm. The acceptance angle of the detector
  is limited to 10° around the tilted detector axis. Fourier transform is used
  to compute the spatial frequency-domain reflectance. Note that this transform
  produces complex values with amplitude and phase information.


Note that the SFD datasets are run with 1000 million photon packets and that
the simulation termination radius is set to 150 |nbsp| mm.

**Dataset files**

Executing the generated scripts will save the simulation results / datasets
into compressed numpy data files that will be organized as follows:

.. code-block:: sh

    data/mcml/<sample>/<source>/<detector>/<pf>/<pf_param_1>/<pf_param_2>/mua-<mua>-musr-<musr>-invcm.npz

The values of placeholders <> are as follows:

* :sh:`<sample>` can take the following values:

    * :sh:`1-layer-semiinfinite`
        A single samplelayer of infinite thickness.

* :sh:`<source>` is the type of the photon packet source / probe used in the datasets:

    * :sh:`line` for a infinitely narrow line source
      :py:class:`~xopto.mcml.mcsource.line.Line`.

    * :sh:`collimated-200um` for a
      :py:class:`~xopto.mcml.mcsource.uniformbeam.UniformBeam` with
      a 200 |nbsp| µm beam diameter.

    * :sh:`gaussian-fwhm-100um` for a
      :py:class:`~xopto.mcml.mcsource.gaussianbeam.GaussianBeam` with
      100 |nbsp| µm beam FWHM.

    * :sh:`fiber-200um-0_22na` for a
      :py:class:`~xopto.mcml.mcsource.fiber.UniformFiber` with fiber core
      diameter of 200 |nbsp| µm and 0.22 NA.

    * :sh:`six-around-one-200um-0_22na` for a
      :py:class:`~xopto.mcml.mcsurface.probe.sixaroundone.SixAroundOne`
      probe surface layout with a central optical fiber source
      :py:class:`~xopto.mcml.mcsource.fiber.UniformFiber`. All optical fibers
      have a core diameter of 200 |nbsp| µm, cladding diameter of 220 |nbsp| µm,
      0.22 NA and are tightly packed.

    * :sh:`six-around-one-400um-0_22na` for a
      :py:class:`~xopto.mcml.mcsurface.probe.sixaroundone.SixAroundOne`
      probe surface layout with a central optical fiber source
      :py:class:`~xopto.mcml.mcsource.fiber.UniformFiber`. All optical fibers
      have a core diameter of 400 |nbsp| µm, cladding diameter of 420 |nbsp| µm,
      0.22 NA and are tightly packed.

    * :sh:`six-linear-array-200um-0_22na` for a
      :py:class:`~xopto.mcml.mcsurface.probe.lineararray.LinearArray`
      probe surface layout with 6 tightly packed optical fibers, a leftmost
      source fiber (negative :math:`x` coordinate)
      :py:class:`~xopto.mcml.mcsource.fiber.UniformFiber`. All optical fibers
      have a core diameter of 200 |nbsp| µm, cladding diameter of 220 |nbsp| µm,
      0.22 NA and are tightly packed.

    * :sh:`single-fiber-100um-0_22na` for a
      :py:class:`~xopto.mcml.mcsurface.probe.lineararray.LinearArray`
      probe surface layout with one central optical fiber, a fiber source
      :py:class:`~xopto.mcml.mcsource.fiber.UniformFiber`. The optical fiber
      has a core diameter of 100 |nbsp| µm, cladding diameter of 120 |nbsp| µm,
      0.22 NA.

    * :sh:`single-fiber-200um-0_22na` for a
      :py:class:`~xopto.mcml.mcsurface.probe.lineararray.LinearArray`
      probe surface layout with one central optical fiber, a fiber source
      :py:class:`~xopto.mcml.mcsource.fiber.UniformFiber`. The optical fiber
      has a core diameter of 200 |nbsp| µm, cladding diameter of 220 |nbsp| µm,
      0.22 NA.

    * :sh:`single-fiber-400um-0_22na` for a
      :py:class:`~xopto.mcml.mcsurface.probe.lineararray.LinearArray`
      probe surface layout with one central optical fiber, a fiber source
      :py:class:`~xopto.mcml.mcsource.fiber.UniformFiber`. The optical fiber
      has a core diameter of 400 |nbsp| µm, cladding diameter of 420 |nbsp| µm,
      0.22 NA.

    * :sh:`single-fiber-800um-0_22na` for a
      :py:class:`~xopto.mcml.mcsurface.probe.lineararray.LinearArray`
      probe surface layout with one central optical fiber, a fiber source
      :py:class:`~xopto.mcml.mcsource.fiber.UniformFiber`. The optical fiber
      has a core diameter of 800 |nbsp| µm, cladding diameter of 820 |nbsp| µm,
      0.22 NA.

* :sh:`<detector>` is the type of detector used by the datasets:
    
    * :sh:`radial` for simple sources with laterally uniform source-sample boundary,

    * :sh:`probe` for optical fiber probes with surface layout.

* :sh:`<pf>` is the type of scattering phase function used in the datasets:
    
    :sh:`hg` for :py:class:`~xopto.mcbase.mcpf.hg.Hg`.

    :sh:`mhg` for :py:class:`~xopto.mcbase.mcpf.mhg.MHg`.

    :sh:`gk` for :py:class:`~xopto.mcbase.mcpf.gk.Gk`.

    :sh:`mie-polystyrene` for :py:class:`~xopto.pf.mie.Mie` - a water
    suspension of polystyrene spheres.

    :sh:`mie-fusedsilica` for :py:class:`~xopto.pf.mie.Mie` - a water
    suspension of fused silica spheres.

* :sh:`pf_param_1`: is the first parameter of the scattering phase function
  formatted with two decimal digits and using :sh:`_` as the decimal separator:

    * :sh:`g-<g>` for :py:class:`~xopto.mcbase.mcpf.hg.Hg`, e.g. :sh:`g-0_10` for :math:`g=0.1`.

    * :sh:`g-<g>` for :py:class:`~xopto.mcbase.mcpf.mhg.MHg`, e.g. :sh:`g-0_50` for :math:`g=0.5`.

    * :sh:`g-<g>` for :py:class:`~xopto.mcbase.mcpf.gk.Gk`, e.g. :sh:`g-0_90` for :math:`g=0.9`.

    * :sh:`diameter-<g>um` for :py:class:`~xopto.pf.mie.Mie`, e.g. :sh:`diameter-0_25` for :math:`d=0.25` |nbsp| µm.

* :sh:`pf_param_2`: is the second parameter of the scattering phase function
  formatted with two decimal digits and using :sh:`_` as the decimal
  separator. An exception to this rule is the wavelength parameter of the
  :py:class:`~xopto.pf.mie.Mie` scattering phase function that is converted to
  nm and formatted as an integer.
  This placeholder is not used with the :py:class:`~xopto.mcbase.mcpf.hg.Hg`
  scattering phase function.

    * :sh:`b-<b>` for :py:class:`~xopto.mcbase.mcpf.mhg.MHg`, e.g. :sh:`b-0_60` for :math:`b=0.6`.

    * :sh:`a-<a>` for :py:class:`~xopto.mcbase.mcpf.gk.Gk`, e.g. :sh:`a-0_50` for :math:`a=0.5`.

    * :sh:`wavelength-<w>nm` for :py:class:`~xopto.pf.mie.Mie`, e.g. :sh:`wavelength-500nm` for :math:`wavelength=500` |nbsp| nm.

* :sh:`<mua>` is the absorption coefficient in units of cm :sup:`-1` with
  two decimal digits and :sh:`_` as the decimal separator, e.g :sh:`2_50`
  for :math:`\mu_a=2.5` |nbsp| cm :sup:`-1`.

* :sh:`<musr>` is the reduced scattering coefficient in units of cm :sup:`-1`
  with two decimal digits and :sh:`_` as a decimal separator, e.g :sh:`20_00`
  for :math:`\mu_a=20.0` |nbsp| cm :sup:`-1`.

Voxelized media datasets
------------------------

Module :py:mod:`xopto.dataset.render.mcvox` => :sh:`run/mcvox/`
  
Renders a dataset scripts for computing fluence / deposition
datasets with the MC kernel for voxelized media. A two-layer skin model
with an embedded blood vessel is used. The depth/position of the blood vessel
along the :math:`z` axis is varied from 0.2 to 0.8 |nbsp| mm in steps of
0.025 |nbsp| mm.
The refractive index of the surrounding medium is set to 1.337.
The simulations are run with 1000 million photon packets.

.. list-table:: A two-layer skin model with an embedded blood vessel
    :header-rows: 1
    :widths: 50, 50, 100, 100

    * -
      - Parameter
      - Value
      - Description
    * - :py:class:`~xopto.mcvox.mcsource.line.Line`
      -
      -
      - normally incident
    * - :py:class:`~xopto.mcvox.mcmaterial.Material`
      - :math:`\mu_a`
      - 16.5724 |nbsp| cm :sup:`-1`
      - epidermis
    * -
      - :math:`\mu_s`
      - 375.9398 |nbsp| cm :sup:`-1`
      -
    * -
      - :math:`n`
      - 1.337
      -
    * -
      - :math:`pf`
      - :py:class:`~xopto.mcbase.mcpf.hg.Hg` (0.9)
      -
    * - :py:class:`~xopto.mcvox.mcmaterial.Material`
      - :math:`\mu_a`
      - 45.85 |nbsp| cm :sup:`-1`
      - dermis
    * -
      - :math:`\mu_s`
      - 356.5406 |nbsp| cm :sup:`-1`
      -
    * -
      - :math:`n`
      - 1.337
      -
    * -
      - :math:`pf`
      - :py:class:`~xopto.mcbase.mcpf.hg.Hg` (0.9)
      -
    * - :py:class:`~xopto.mcvox.mcmaterial.Material`
      - :math:`\mu_a`
      - 230.5427 |nbsp| cm :sup:`-1`
      - blood vessel
    * -
      - :math:`\mu_s`
      - 93.985 |nbsp| cm :sup:`-1`
      -
    * -
      - :math:`n`
      - 1.337
      -
    * -
      - :math:`pf`
      - :py:class:`~xopto.mcbase.mcpf.hg.Hg` (0.9)
      -
    * - :py:class:`~xopto.mcbase.mcfluence.Fluence`
      - :math:`xaxis`
      - :py:class:`~xopto.mcbase.mcutil.axis.Axis` |nbsp| (start=-502.5e-6, stop=502.5e-6, n=201)
      - sample and deposition voxelization
    * -
      - :math:`yaxis`
      - :py:class:`~xopto.mcbase.mcutil.axis.Axis` |nbsp| (start=-502.5e-6, stop=502.5e-6, n=201)
      -
    * -
      - :math:`zaxis`
      - :py:class:`~xopto.mcbase.mcutil.axis.Axis` |nbsp| (start=0.0, stop=0.001, n=200)
      -
    * - :py:class:`~xopto.mcvox.mcdetector.cartesian.Cartesian`
      - :math:`xaxis`
      - :py:class:`~xopto.mcbase.mcutil.axis.Axis` |nbsp| (start=-502.5e-6, stop=502.5e-6, n=201)
      - reflectance and transmittance detectors
    * -
      - :math:`yaxis`
      - :py:class:`~xopto.mcbase.mcutil.axis.Axis` |nbsp| (start=-502.5e-6, stop=502.5e-6, n=201)
      -
    * - **Blood vessel**
      - :math:`diameter`
      - 200.0 |nbsp| μm
      - in dermis
    * -
      - :math:`position`
      - (0, 0, 0.0002-0.0008)
      - (x, y, z)
    * -
      - :math:`direction`
      - (0, 1, 0)
      - (x, y, z)
    * - **Epidermis**
      - :math:`thickness`
      - 100 |nbsp| μm
      -

**Dataset files**

Executing the generated scripts will save the simulation results / datasets
into compressed numpy data files that will be organized as follows:

.. code-block:: sh

    data/mcvox/fluence/2-layer-skin-<diameter>um-vessel-<depth>um-depth-deposition.npz

The values of placeholders <> are as follows:

* :sh:`<diameter>` is the diameter of the blood vessel in units of μm, formatted
  as an integer value, e.g :sh:`200` for a 200 |nbsp| μm blood vessel.

* :sh:`<depth>` is the :math:`z` coordinate (depth) of the blood vessel in units
  of μm, formatted as an integer value, e.g :sh:`500` for :math:`z=500`
  |nbsp| μm.

Sampling volume datasets
------------------------

Module :py:mod:`xopto.dataset.render.sv` => :sh:`run/sv/`
  
Renders a reference dataset script for computing sampling-volume
datasets. The sampling volume dataset is computed for
a semi-infinite homogeneous medium for an optical fiber probe with two
optical fibers placed at a distance of 0.5 |nbsp| mm. The refractive index
of the surrounding medium is set to 1.0. Simulations are run in batches
until 1,000,000 photon packet traces that reach the detector fiber are
collected and converted to sampling volume information. The trace capacity is
limited to 1000 events. The simulation termination radius is set to
25 |nbsp| mm. 

.. list-table:: Sampling volume for a probe with two optical fibers
    :header-rows: 1
    :widths: 50, 50, 100, 100

    * -
      - Parameter
      - Value
      - Description
    * - :py:class:`~xopto.mcml.mcsurface.probe.lineararray.LinearArray`
      - :math:`dcore`
      - 200 |nbsp| μm
      - linear layout of 2 fibers
    * -
      - :math:`dcladding`
      - 220 |nbsp| μm
      -
    * -
      - :math:`ncore`
      - 1.462
      -
    * -
      - :math:`na`
      - 0.22
      -
    * -
      - :math:`spacing`
      - 500 |nbsp| μm
      -
    * -
      - :math:`n`
      - 2
      -
    * -
      - :math:`diameter`
      - 6.0 |nbsp| mm
      -
    * -
      - :math:`reflectivity`
      - 0.6
      -
    * - :py:class:`~xopto.mcbase.mctrace.Trace`
      - :math:`maxlen`
      - 1000
      - packet trace configuration
    * - :py:class:`~xopto.mcbase.mcsv.SamplingVolume`
      - :math:`xaxis`
      - :py:class:`~xopto.mcbase.mcutil.axis.Axis` |nbsp| (start=-0.00075, stop=0.00075, n=300)
      - sampling volume voxelization
    * -
      - :math:`yaxis`
      - :py:class:`~xopto.mcbase.mcutil.axis.Axis` |nbsp| (start=-0.00075, stop=0.00075, n=300)
      -
    * -
      - :math:`zaxis`
      - :py:class:`~xopto.mcbase.mcutil.axis.Axis` |nbsp| (start=0.0, stop=0.001, n=200)
      -
    * - :py:class:`~xopto.mcml.mclayer.layer.Layer`
      - :math:`\mu_a`
      - 2.0 |nbsp| cm :sup:`-1`
      -  sample layer
    * -
      - :math:`\mu_s`
      - 500.0 |nbsp| cm :sup:`-1`
      - 
    * -
      - :math:`n`
      - 1.337
      - 
    * -
      - :math:`pf`
      - :py:class:`~xopto.mcbase.mcpf.hg.Hg` |nbsp| (0.95)
      - 

**Dataset files**

Executing the generated scripts will save the simulation results / datasets
into compressed numpy data files that will be organized as follows:

.. code-block:: sh

    data/mcml/1-layer-semiinfinite/sv/reflectance/fiber-200um-0_22na/sds-<sds>um/hg/g-<g>/um-mua-<mua>-musr-<musr>-invcm.npz

The values of placeholders <> are as follows:

* :sh:`<sds>` is the distance between the centers of the source and detector
  fibers in units of μm, formatted as an integer value, e.g :sh:`500` for a 
  500 |nbsp| μm distance.

* :sh:`<g>` is the anisotropy formatted with two decimal digits and :sh:`_`
  as the decimal separator, e.g :sh:`0_15` for :math:`g=0.15`.

* :sh:`<mua>` is the absorption coefficient in units of cm :sup:`-1` with
  two decimal digits and :sh:`_` as the decimal separator, e.g :sh:`2_50`
  for :math:`\mu_a=2.5` |nbsp| cm :sup:`-1` .

* :sh:`<musr>` is the reduced scattering coefficient in units of cm :sup:`-1`
  with two decimal digits and :sh:`_` as a decimal separator, e.g :sh:`20_00`
  for :math:`\mu_s'=20.0` |nbsp| cm :sup:`-1`.

All available datasets
----------------------

* :py:mod:`xopto.dataset.render.all`
  Renders all the available dataset scripts.
