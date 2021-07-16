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

.. _example-scattering-phase-function:

Reflectance spectrum dependence on scattering phase function
============================================================

This example covers the influence of the scattering phase function on the reflectance spectrum as acquired with optical fiber probes with small source-detector separations. The wavelength dependency of the scattering phase function and the scattering coefficient is calculated using Mie theory for a water suspension of 1 μm diameter polystyrene microspheres which are readily available from suppliers such as Polysciences, MicroParticles GmbH, Bangs Laboratories Inc., etc. The reflectance spectrum simulated using the scattering phase function obtained by Mie theory is compared to the most commonly used Henyey-Greenstein scattering phase function which takes into account only the anisotropy factor :math:`g`.

.. note::

    This example might run several minutes, depending on the available computational equipment.

Importing the required modules and submodules
---------------------------------------------
In this example, we import submodules from xopto, which enable an interface to Monte Carlo simulations (:py:mod:`xopto.mcml.mc`), definitions of optical fibers (:py:mod:`xopto.mcml.mcutil.fibers`), detection accumulators axes (:py:mod:`xopto.mcml.mcutil.axis`), utilities for mathematical integration of accumulators (:py:mod:`xopto.util.convolve`), functions and materials for scattering phase function calculations based on Mie theory (:py:mod:`xopto.pf.miepolystyrene`, :py:mod:`xopto.materials.ri.polystryene` and :py:mod:`xopto.materials.ri.water`) and finally submodule for computational device definitions (:py:mod:`xopto.cl.clinfo`). We also require the :py:mod:`numpy` and :py:mod:`matplotlib.pyplot` modules for mathematical functions and plotting.

.. code-block:: python

    from xopto.mcml import mc
    from xopto.mcml.mcutil import fiber, axis
    from xopto.util import convolve
    from xopto.pf import miepolystyrene
    from xopto.materials.ri import polystyrene, water
    from xopto.cl import clinfo

    import numpy as np
    import matplotlib.pyplot as pp

Calculating the spectrum of the scattering coefficient
------------------------------------------------------
The spectrum or wavelength dependency of the scattering coefficient from 1 μm diameter polystyrene microspheres suspended in water can be calculated using the Mie theory. For this, we first require the wavelength dependency of the refractive index of polystryene and water. The refractive index is available within the :py:mod:`xopto.materials.ri` package. We import the submodules :py:mod:`~xopto.materials.ri.polystyrene` and :py:mod:`~xopto.materials.ri.water` within which various refractive index models are available as measured by different authors. In this example, we utilize Nikolov's measurements of wavelength dependent refractive index of polystyrene and Daimon's measurements of wavelength dependent refractive index of water. It should be noted that variables :code:`ri_poly` and :code:`ri_water` are now instances of  :py:mod:`~xopto.materials.ri.polystyrene.Nikolov()` and :py:mod:`~xopto.materials.ri.water.Daimon()`. Initialization also accepts a temperature parameter :code:`t` in K. For this example we used the default value that corresponds to 20°C.

.. code-block:: python

    ri_poly = polystyrene.Nikolov()
    ri_water = water.Daimon()

Next, we choose a single wavelength (600 nm) and set the desired value of the reduced scattering coefficient at this wavelength (15 1/cm). This will help us determine a multiplicative factor that will be used to obtain the scattering coefficient at each wavelength from the scattering cross section. According to the input keyword parameters :code:`diameter`, :code:`wavelength˙` and refractive indices of polystyrene :code:`ripolystyrene` and water :code:`riliquid`, we construct a Mie theory model for polystryene microspheres suspended in water using :py:class:`~xopto.pf.miepolystyrene.MiePolystyrene`. The Mie theory model :py:class:`~xopto.pf.miepolystyrene.MiePolystyrene` allows us to calculate the anisotropy factor and scattering cross section at the chosen wavelength using methods :py:meth:`~xopto.pf.miepolystyrene.MiePolystyrene.g` and :py:meth:`~xopto.pf.miepolystyrene.MiePolystyrene.scs`. Note that input indices in method :py:meth:`~xopto.pf.miepolystyrene.MiePolystyrene.g` correspond to different scattering phase function moments. Input index 1 corresponds to the anisotropy factor.

In principle, the scattering coefficient can be calculated from the scattering cross section (:math:`\sigma_s`), if the number density :math:`[1/m^3]` of polystyrene microspheres is known:

.. math::

    \mu_s = n \; \sigma_s

However, in our case we are interested to obtain a reduced scattering coefficient of 15 1/cm at 600 nm, which is a common value for biological tissues. In this case we can calculate a multiplicative factor that takes the role of the number density. Firstly, we calculate the target scattering coefficient :code:`mus_target` from the target reduced scattering coefficient by dividing the latter with the :code:`1-mie.g(1)`. Subsequently, we can obtain the number density (:code:`n`) by dividing the target scattering coefficient :code:`mus_target` by the scattering cross section :code:`mie.scs()`.

.. code-block:: python

    wavelength = 600e-9
    musr_target = 15e2
    mie = miepolystyrene.MiePolystyrene(
        diameter=1.0e-6,
        wavelength=wavelength,
        ripolystyrene=ri_poly(wavelength),
        riliquid=ri_water(wavelength)
    )

    mus_target = musr_target/(1-mie.g(1))
    n = mus_target/mie.scs()

Finally, we can calculate the spectrum of the scattering coefficient. Firstly, we define the array of wavelength points spanning from 450 to 800 nm in 2 nm steps. We also define zero-initialized 1D arrays :code:`mus` and :code:`g` into which we will store the scattering coefficients and anisotropy factors as calculated by the Mie theory. Using the :code:`for` loop, we iterate through the wavelength points and at each iteration define a new object of Mie theory model with new input parameters that depend on the wavelength, i.e., the wavelength point :code:`wavelength` and refractive indices of polystyrene :code:`ripolystyrene` and water :code:`riliquid`.

Since the number density does not change with wavelength, the previously calculated multiplicative factor :code:`n` representing the number density can be used to multiply the scattering cross section at each wavelength, thus resulting in the scattering coefficient. The scattering coefficient and also the anisotropy factor are then stored in separate arrays :code:`mus` and :code:`g` at each iteration.

The obtained spectrum of the reduced scattering coefficient is shown in the figure below. Note that at 600 nm the reduced scattering is 15 1/cm as expected.

.. code-block:: python

    wavelengths = np.arange(450e-9, 801e-9, 2e-9)
    mus = np.zeros_like(wavelengths)
    g = np.zeros_like(wavelengths)
    for i, w in enumerate(wavelengths):

        mie = miepolystyrene.MiePolystyrene(
            diameter=1.0e-6,
            wavelength=w,
            ripolystyrene=ri_poly(w),
            riliquid=ri_water(w)
        )

        mus[i] =  factor * mie.scs()
        g[i] = mie.g(1)

    pp.figure()
    pp.plot(wavelengths, mus * (1-g))
    pp.show()

.. image:: ../../images/examples/mcml/scattering_phase_function/red_scat_spectrum_pf.svg
    :height: 480 px
    :width: 640 px
    :scale: 75 %
    :alt: Scattering phase function
    :align: center

Computational device
--------------------

In this section we select the desired OpenCL computational device (see also :ref:`opencl-devices-label`).

.. code-block:: python

    cl_device = clinfo.gpu(platform='nvidia')

.. note::

    In this example we have selected the first computational device listed under the Nvidia platform. The string should be changed according to the installed hardware devices.

The layer stack
---------------

We define two instances of the layer stack, the first instance :code:`layers_hg` corresponds to a simple case of using the Henyey-Greenstein scattering phase function, the second instance :code:`layers_mie` corresponds to the Mie scattering phase function. Other optical properties between the two instances of the layer stack are the same. Each instance is defined by the :py:class:`xopto.mcml.mclayer.layer.Layers` constructor, which accepts a list of :py:class:`xopto.mcml.mclayer.layer.Layer` objects in an ascending order along the z axis, which points into the medium. The top and bottom layers corresponding to the top and bottom surrounding medium should always be defined. While all the optical properties of the top and bottom layers are set, only the refractive index is used to properly refract/reflect the photon packet at the turbid medium boundary. In this example we define a single layer of turbid medium of thickness 1 cm to which we initially assign some arbitrary absorption and scattering coefficients. These are going to be dynamically changed at each simulation for a particular wavelength point.

The type of the scattering phase function must be the same for all of the layers, while the phase function specific parameter can be changed from layer to layer. For the first layer stack instance :code:`layers_hg` we specify the Henyey-Greenstein scattering phase function by using the constructor :py:class:`~xopto.mcbase.mcpf.hg.Hg`. The anisotropy factor will be changed at each wavelength interation according to the one stored in the array :code:`g`. For the second layer stack instance :code:`layers_mie`, we set a general lookup table based sampling scheme of a Mie scattering phase function cumulative distribution function using the constructor :py:class:`~xopto.mcbase.mcpf.Lut`. The lookup table based sampling is required since the Mie scattering phase function does not exhibit an analytical inverse of the cumulative distribution function. The constructor :py:class:`~xopto.mcbase.mcpf.Lut` accepts parameters that are returned by the :py:meth:`~xopto.pf.miepolystyrene.MiePolystyrene.mclut` method.

.. code-block:: python

    d = 1e-2
    layers_hg = mc.mclayer.Layers([
        mc.mclayer.Layer(d=0.0, n=1.45, mua=0.0, mus=0.0, pf=mc.mcpf.Hg(0.0)),
        mc.mclayer.Layer(d=d, n=1.33, mua=0.0, mus=0.0, pf=mc.mcpf.Hg(0.8)),
        mc.mclayer.Layer(d=0.0, n=1.45, mua=0.0, mus=0.0, pf=mc.mcpf.Hg(0.0)),
    ])

    mie = miepolystyrene.MiePolystyrene(
        diameter=1.0e-6,
        wavelength=w,
        ripolystyrene=ri_poly(w),
        riliquid=ri_water(w)
    )

    mie_pf = mc.mcpf.Lut(*mie.mclut())
    layers_mie = mc.mclayer.Layers([
        mc.mclayer.Layer(d=0.0, n=1.45, mua=0.0, mus=0.0, pf=mie_pf),
        mc.mclayer.Layer(d=d, n=1.33, mua=0.0, mus=0.0, pf=mie_pf),
        mc.mclayer.Layer(d=0.0, n=1.45, mua=0.0, mus=0.0, pf=mie_pf),
    ])

Source
------

We define a multimode optical fiber source having a core diameter :code:`dcore` of 200 |nbsp| μm, combined diameter of core and cladding :code:`dcladding` of 220 |nbsp| μm, core refractive index :code:`ncore` of 1.45 and numerical aperture :code:`na` of 0.22. The multimode optical fiber is constructed by :py:class:`~xopto.mcbase.mcutil.fiber.MultimodeFiber`. It is then passed to the optical fiber source, which launches photon packets uniformly into the solid angle :py:class:`~xopto.mcml.mcsource.fiber.UniformFiber`. The optical fiber source is positioned normally to the turbid medium as assigned by the :code:`direction` and :code:`position` keyword arguments. Note that these are the default values of :code:`direction` and :code:`position`.

.. code-block:: python

    source = mc.mcsource.UniformFiber(
        fiber=fiber.MultimodeFiber(
            dcore=200e-6, dcladding=220e-6, ncore=1.45, na=0.22
        ),
        direction=(0.0, 0.0, 1.0),
        position=(0.0, 0.0, 0.0)
    )

Detector
--------

As a detector, we define a radial accumulator spanning from 0 to 1 mm with 250 bins. We will utilize this detector to integrate over the annular accumulator rings that overlap with a virtual multimode fiber detector. Such an approach effectively collects all of the photon packets that remit within a certain radius interval and thus significantly improves the signal-to-noise ratio in comparison to a single multimode optical fiber detector placed tightly to the source fiber. The radial accumulator is defined using an instance of :py:class:`~xopto.mcml.mcdetector.radial.Radial` which accepts the radial axis instance :py:class:`~xopto.mcbase.mcutil.axis.RadialAxis` and a :code:`cosmin` parameter, which defines the minimal cosine of the polar angle within which the photon packets are accepted. Note that this parameter is related to the numerical aperture of the fiber.

.. code-block:: python

    detector_top_radial = mc.mcdetector.Radial(
        raxis=axis.RadialAxis(
            start=0.0, stop=1e-3, n=250
        ),
        cosmin=np.cos(np.arcsin(0.22/1.45))
    )

The Monte Carlo simulator objects
---------------------------------

We define two :py:class:`~xopto.mcml.mc.MC` simulator objects, which accept the defined layer stacks, sources, detectors and OpenCL ready device context. The :code:`detector` keyword argument accepts an instance of :py:class:`~xopto.mcml.mcdetector.base.Detectors`, which in turn accepts a top, bottom and specular detectors. In this example we have set only the top detector :code:`detector_top_radial`, which we have defined above.

In both cases we also define a maximum simulation radius :py:attr:`xopto.mcml.mc.Mc.rmax` beyond which the photon packets are not propagated since they likely do not contribute to the reflectance signal.

.. code-block:: python

    mc_obj_hg = mc.Mc(
        layers=layers_hg, 
        source=source, 
        detectors=mc.mcdetector.Detectors(top=detector_top_radial), 
        cl_devices=cl_device
    )

    mc_obj_mie = mc.Mc(
        layers=layers_mie, 
        source=source, 
        detectors=mc.mcdetector.Detectors(top=detector_top_radial), 
        cl_devices=cl_device
    )

    mc_obj_hg.rmax = 1e-2
    mc_obj_mie.rmax = 1e-2

Running the Monte Carlo simulations and reflectance processing
--------------------------------------------------------------

Before the simulations are run, we define two arrays that will store reflectance :code:`reflectance_hg` and :code:`reflectance_mie` and correspond to the case with Henyey-Greenstein and Mie scattering phase functions, respectively. Next, using the :code:`for` loop, we iterate over the wavelengths and at each iteration change the scattering coefficient and the scattering phase function of the single layer stored within the  :py:class:`~xopto.mcml.mclayer.layer.Layers` object returned by the attribute :py:attr:`~xopto.mcml.mc.Mc.layers`. Note that the turbid medium single layer :py:class:`~xopto.mcml.mclayer.layer.Layer` is accessed with an index 1. In the case of Henyey-Greenstein phase function only the anisotropy factor is changed according to the one calculated by the Mie theory. In the case of Mie scattering phase function the lookup table for sampling the cumulative distribution function has to be recalculated at each wavelength point, since the shape of the scattering phase function changes. Each simulation is run with :math:`10^8` launched photon packets. Note that we only access the last returned parameter which corresponds to the :py:class:`~xopto.mcml.mcdetector.base.Detectors` object and separately stores top and bottom detectors, if defined. In this example we have defined only the top detector :py:class:`~xopto.mcml.mcdetector.radial.Radial` that can be accessed via the attribute :py:attr:`~xopto.mcml.mcdetector.base.Detectors.top`. The acquired reflectance and radial positions of the :py:class:`~xopto.mcml.mcdetector.radial.Radial` accumulator bins can be acessed via, e.g., :code:`detector_hg.top.r` and :code:`detector_hg.top.r`, respectively. Since we are interested in the reflectance as detected by a multimode optical fiber detector located tightly next to the source fiber, we utilize the function :py:func:`~xopto.util.convolve.fiber_reflectance`, which accepts the bin positions, acquired signal at each corresponding radial accumulator, source-detector separation (:code:`sds`) and diameter of the multimode fiber detector core (:code:`dcore`). The returned reflectance is then stored to the :code:`reflectance_hg` and :code:`reflectance_mie` arrays. It should be noted that :py:func:`~xopto.util.convolve.fiber_reflectance` function returns a 2D array of shape (number of radial reflectance accumulator arrays, number of detector fibers). Since a single radial reflectance accumulator and single detector fibers is used, we access the desired calculated reflectance via :code:`[0,0]`.

.. code-block:: python

    reflectance_hg = np.zeros_like(wavelengths)
    reflectance_mie = np.zeros_like(wavelengths)

    nphotons = 1e8
    for i, w in enumerate(wavelengths):

        mc_obj_hg.layers[1].mus = mus[i]
        mc_obj_hg.layers[1].pf = mc.mcpf.Hg(g[i])

        mc_obj_mie.layers[1].mus = mus[i]

        mie = miepolystyrene.MiePolystyrene(
            diameter=1.0e-6,
            wavelength=w,
            ripolystyrene=ri_poly(w),
            riliquid=ri_water(w)
        )
        mie_pf = mc.mcpf.Lut(*mie.mclut())
        mc_obj_mie.layers[1].pf = mie_pf

        detector_hg = mc_obj_hg.run(nphotons, verbose=True, wgsize=256)[-1]
        detector_mie = mc_obj_mie.run(nphotons, verbose=True, wgsize=256)[-1]

        reflectance_hg[i] = convolve.fiber_reflectance( 
            detector_hg.top.r,
            detector_hg.top.reflectance,
            sds=220e-6,
            dcore=200e-6
        )[0,0]

        reflectance_mie[i] = convolve.fiber_reflectance( 
            detector_mie.top.r,
            detector_mie.top.reflectance,
            sds=220e-6,
            dcore=200e-6
        )[0,0]


Comparing the obtained reflectance
----------------------------------

Finally, we compare the reflectance acquired with a multimode fiber when the two underlying turbid media exhibit the same scattering coefficients and anisotropy factors and thus the same reduced scattering coefficient, while the only difference is the utilized scattering phase function. As can be observed in the plot below, the two reflectance spectra differ significantly. This is due to the small source-detector separations that render reflectance highly sensitive to the scattering phase function. By increasing the source-detector separations, the different should become less obvious, especially when the separations would be in the so called diffusion regime, where reflectance only depends on the absorption and reduced scattering coefficients, which are the same for both media.

.. code-block:: python

    pp.figure()
    pp.plot(1e9*wavelengths, 100*reflectance_hg, label='Henyey-Greenstein pf')
    pp.plot(1e9*wavelengths, 100*reflectance_mie, label='Mie pf')
    pp.legend()
    pp.xlabel('Wavelength (nm)')
    pp.ylabel('Reflectance (%)')
    pp.show()

.. image:: ../../images/examples/mcml/scattering_phase_function/reflectance_phase_function.svg
    :height: 480 px
    :width: 640 px
    :scale: 75 %
    :alt: Traces of filtered photon packets
    :align: center


The complete example
--------------------

.. literalinclude:: ../../../../examples/mcml/reflectance_phase_function.py

You can run this example from the root directory of the PyXOpto package as:

.. code-block:: bash

    python examples/mcml/reflectance_phase_function.py