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

.. _mcvox-material-label:

Material
========

Each voxel of the volume can be assigned a different material. Materials are
defined as instances of :py:class:`xopto.mcvox.mcmaterial.material.Material`.
The constructor takes several parameters that define the optical properties of
the material.

1. The refractive index is defined by parameter :code:`n`.
2. The scattering coefficient is defined by parameter :code:`mus` (m :superscript:`-1`).
3. The absorption coefficient is defined by parameter :code:`mua` (m :superscript:`-1`).
4. The scattering phase function is defined by parameter :code:`pf` that can
   be any instance of :py:class:`xopto.mcbase.mcpf.pfbase.PfBase`.

The following example creates a material with a refractive index 1.33,
absorption coefficient 1 |nbsp| cm :superscript:`-1` , scattering coefficient
50 |nbsp| cm :superscript:`-1` and Henyey-Greenstein scattering phase function
with anisotropy 0.9.

.. code-block:: python

    from xopto.mcvox import mc

    material = mc.mcmaterial.Material(n=1.33, mua=1.0e2, mus=50.0e2, pf=mc.mcpf.Hg(0.9))

All the optical properties of a material that were set
through the constructor :py:meth:`~xopto.mcvox.mcmaterial.material.Material`
can be later changed through accessing the instance properties. In case of
the scattering phase function, first access the 
:py:attr:`~xopto.mcvox.mcmaterial.material.Material.pf` material property
and from there any of the properties implemented by the scattering
phase function model. Note that the Henyey-Greenstein scattering phase
function from this example exposes only the anisotropy
:py:attr:`~xopto.mcbase.mcpf.hg.Hg.g`.

.. code-block:: python

    material.mua = 0.5e2
    material.mus = 40.0e2
    material.n = 1.452
    material.pf.g = 0.9

The individual materials are then combined into a list through
:py:meth:`~xopto.mcvox.mcmaterial.material.Materials`. The constructor takes
a list of :py:class:`~xopto.mcvox.mcmaterial.material.Material`. The created
instance manages the transfer of data between the host and the OpenCL device.

.. code-block:: python

    materials = mc.mcmaterial.Materials(
        [
            mc.mcmaterial.Material(n=1.0, mua=0.0e2, mus=0.0e2,  pf=mc.mcpf.Hg(0.0)),
            mc.mcmaterial.Material(n=1.3, mua=1.0e2, mus=50.0e2, pf=mc.mcpf.Hg(0.9)),
            mc.mcmaterial.Material(n=1.2, mua=2.0e2, mus=10.0e2, pf=mc.mcpf.Hg(0.5))
        ]
    )

Note that the first material in the list represent the surrounding
medium. The absorption coefficient :code:`mua`, scattering coefficient
:code:`mus` and the scattering phase function :code:`pf` of the
first material are not used in the MC simulations, since the
photon packets are not propagated through the surrounding medium.
However, the refractive index :code:`n` of the surrounding medium
is used to properly refract/reflect the photon packet at the sample surface
when launched by the source or when escaping the sample.

Voxels of the samples can be labeled by a particular material through the
corresponding index of the material in the list. A label 0 will set the
voxel material to the same material as is used for the surrounding medium.

The properties of materials in the list can be modified at any time
through accessing the individual materials and from there the material
properties. The number of materials (including the material of the surrounding
medium) can be determined through the builtin :py:func:`len`. The individual
materials can be accessed using the :code:`[]` operator or
:py:meth:`~xopto.mcvox.mcmaterial.material.Materials.material` method.

.. code-block:: python

    materials[1].mua = 0.5e2
    materials.material(1).mua = 0.5e2

    num_materials = len(materials)

The :py:class:`~xopto.mcvox.mcmaterial.material.Material` and
:py:class:`~xopto.mcvox.mcmaterial.material.Materials` produce an informative
human-readable output when used with the :py:func:`print` builtin.

Applying :py:func:`print` to a material (instance of
:py:class:`~xopto.mcvox.mcmaterial.material.Material`) will produce the
following output:

.. code-block:: python

    print(materials[1])

.. code-block:: text

    Material(n=1.3, mua=100.0, mus=5000.0, pf=Hg(g=0.9)) # id 0x7F15B4CEF820.

Applying :py:func:`print` to a list of materials (instance of
:py:class:`~xopto.mcvox.mcmaterial.material.Materials`) will produce the
following output:

.. code-block:: python

    print(materials)

.. code-block:: text

    Materials([
        Material(n=1.0, mua=0.0, mus=0.0, pf=Hg(g=0.0)), # surrounding medium,
        Material(n=1.3, mua=100.0, mus=5000.0, pf=Hg(g=0.9)),
        Material(n=1.2, mua=200.0, mus=1000.0, pf=Hg(g=0.5))
    ]) # id 0x7F15B4CEFD30.
