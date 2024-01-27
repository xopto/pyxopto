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

.. _mcml-layer-stack-label:

The layer stack
===============

The sample and the surrounding medium are described by a stack of layers.
The layers are stacked along the positive z-axis. The first layer
in the stack is used to describe the medium above the top sample
surface and the last layer in the stack is used to describe the
medium bellow the bottom sample surfaces.

The classes that allow definition of individual layers and forming of
layer stacks can be found in the :py:mod:`xopto.mcml.mclayer.layer` module. The
module is for convenience also imported into the main simulator module
:py:mod:`xopto.mcml.mc`.

An isotropic layer is an instance of :py:class:`xopto.mcml.mclayer.layer.Layer`
while an anisotropic layer is an instance of
:py:class:`xopto.mcml.mclayer.layer.AnisotropicLayer`.
The constructor takes several parameters that define the geometry and
optical properties of the layer material.

1. The layer thickness is defined by parameter :code:`d` (m).
2. The refractive index is defined by parameter :code:`n`.
3. The scattering coefficient is defined by parameter :code:`mus` (m :superscript:`-1`).
   For isotropic layer instances (:py:class:`xopto.mcml.mclayer.layer.Layer`) this is
   a scalar value. For anisotropic layer instances
   (:py:class:`xopto.mcml.mclayer.layer.AnisotropicLayer`) this is a 3x3 tensor.
   The scattering coefficient tensor can be initialized from a scalar value
   (all the diagonal elements are set to this value and the scattering becomes 
   isotropic) or a vector of three values that are assigned to the diagonal
   elements of the tensor.
4. The absorption coefficient is defined by parameter :code:`mua` (m :superscript:`-1`).
   For isotropic layer instances (:py:class:`xopto.mcml.mclayer.layer.Layer`) this is
   a scalar value. For anisotropic layer instances
   (:py:class:`xopto.mcml.mclayer.layer.AnisotropicLayer`) this is a 3x3 tensor.
   The absorption coefficient tensor can be initialized from a scalar value
   (all the diagonal elements are set to this value and the absorption becomes 
   isotropic) or a vector of three values that are assigned to the diagonal
   elements of the tensor.
5. The scattering phase function is defined by parameter :code:`pf` that can
   be any instance of :py:class:`xopto.mcbase.mcpf.pfbase.PfBase`.

The following example creates a layer with a refractive index 1.33,
absorption coefficient 1 |nbsp| cm :superscript:`-1` , scattering coefficient
50 |nbsp| cm :superscript:`-1`, thickness 10 |nbsp| mm and Henyey-Greenstein
scattering phase function with anisotropy 0.9.

.. code-block:: python

    from xopto.mcml import mc

    layer = mc.mclayer.Layer(d=0.01, n=1.33, mua=1.0e2, mus=50.0e2, pf=mc.mcpf.Hg(0.9))

All the optical and geometrical properties of a layer that were set
through the constructor :py:meth:`xopto.mcml.mclayer.layer.Layer` can be
later changed through accessing the instance properties. In case of
the scattering phase function, first access the 
:py:attr:`~xopto.mcml.mclayer.layer.Layer.pf` layer property
and from there any of the properties implemented by the scattering
phase function model. Note that the Henyey-Greenstein scattering phase
function from this example exposes only the anisotropy
:py:attr:`~xopto.mcbase.mcpf.hg.Hg.g`.

.. code-block:: python

    layer.mua = 0.5e2
    layer.mus = 40.0e2
    layer.d = 0.015
    layer.n = 1.452
    layer.pf.g = 0.9

The individual layers are then combined into a stack through the
:py:meth:`~xopto.mcml.mclayer.layer.Layers`. The constructor takes a list
of :py:class:`~xopto.mcml.mclayer.layer.Layer`. The created instance
manages the transfer of data between the host and the OpenCL device.

.. note::
    All layer instances that are combined into a layer stack must be of the
    same type (:py:class:`xopto.mcml.mclayer.layer.Layer` or
    :py:class:`xopto.mcml.mclayer.layer.AnisotropicLayer`).

.. code-block:: python

    layers = mc.mclayer.Layers(
        [
            mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0e2, mus=0.0e2,  pf=mc.mcpf.Hg(0.0)),
            mc.mclayer.Layer(d=0.1, n=1.3, mua=1.0e2, mus=50.0e2, pf=mc.mcpf.Hg(0.9)),
            mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0e2, mus=0.0e2,  pf=mc.mcpf.Hg(0.0))
        ]
    )

.. note::
    The first and last layer in the stack represent the surrounding
    medium. The absorption coefficient :code:`mua`, scattering coefficient
    :code:`mus` and the scattering phase function :code:`pf` of the
    topmost and bottommost layers are not used in the MC simulations, since the
    photon packets are not propagated through the surrounding medium.
    However, the refractive index :code:`n` of the two outermost
    layers is used to properly refract/reflect the photon packet at the layer
    boundaries when launched by the source or when escaping the sample.

The properties of layers in the layer stack can be modified at any time
through accessing the individual layers and from there the layer
properties. The number of layers (including the two layers of the surrounding
medium) can be determined through the builtin :py:func:`len`. The individual
layers can be accessed using the :code:`[]` operator or
:py:meth:`~xopto.mcml.mclayer.layer.Layers.layer` method.

.. code-block:: python

    layers[1].mua = 0.5e2
    layers.layer(1).mua = 0.5e2

    num_layers = len(layers)

The :py:class:`~xopto.mcml.mclayer.layer.Layers` and
:py:class:`~xopto.mcml.mclayer.layer.Layer` produce an informative
human-readable output when used with the :py:func:`print` builtin.

Applying :py:func:`print` to an individual layer (instances of
:py:class:`~xopto.mcml.mclayer.layer.Layer`) will produce the following output:

.. code-block:: python

    print(layers[1])

.. code-block:: text

    Layer(d=0.1, n=1.3, mua=100.0, mus=5000.0, pf=Hg(g=0.9)) # id 0x7F15B4CEF820.

Applying :py:func:`print` to a layer stack  (instance of
:py:class:`~xopto.mcml.mclayer.layer.Layers`) will produce the following output:

.. code-block:: python

    print(layers)

.. code-block:: text

    Layers([
        Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=Hg(g=0.0)),  # medium above the sample,
        Layer(d=0.1, n=1.3, mua=100.0, mus=5000.0, pf=Hg(g=0.9)),
        Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=Hg(g=0.0))   # medium bellow the sample
    ]) # id 0x7F15B4CEFD30.
