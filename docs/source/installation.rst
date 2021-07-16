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

.. _installation-label:

.. include:: common.rst

Installation
============

.. _installation-opencl-label:

OpenCL
------

OpenCL programming language and development tools are available for a number
of platforms.

Nvidia
^^^^^^

Windows OS


On Windows OS Nvidia distributes OpenCL with the CUDA Toolkit. Note that that
`Cuda Tolkit 10.2 <https://developer.nvidia.com/cuda-10.2-download-archive?target_os=Windows&target_arch=x86_64>`_
is the last version that supports older Windows OS,
such as Windows 7, Windows 8 and Windows 8.1.
`The latest <https://developer.nvidia.com/cuda-downloads?target_os=Windows&target_arch=x86_64>`_
version of CUDA toolkit only supports Windows 10. Follow the above links to
download the CUDA Toolkit installer executable for your Windows. Run the
installer and follow the dialogs to install the CUDA Toolkit.

Linux OS


You can install the CUDA toolkit by downloading the distribution from
`NVIDIA Developer site <https://developer.nvidia.com/cuda-downloads?target_os=Linux>`_
and  follow the `instructions <https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#package-manager-installation>`_.

Some distributions, such as Ubuntu, also have good native support for installing
the NVIDIA drivers using :code:`ubuntu-drivers` and the autoinstall option that
will install the appropriate driver version for the installed hardware.

.. code:: bash

    sudo ubuntu-drivers autoinstall

Cuda Toolkit can be installed using the package manager:

.. code:: bash

    sudo apt update
    sudo apt install nvidia-cuda-toolkit

To enable OpenCL, you will also have to install :code:`ocl-icd-libopencl1` and
optionally a command line tool :code:`clinfo` that can be used to display all
the installed OpenCL devices.

.. code:: bash

    sudo apt update
    sudo apt install ocl-icd-libopencl1
    sudo apt install clinfo

.. _installation-python-label:

AMD
^^^

Drivers and the OpenCL runtime for AMD devices along with installation instructions
can be downloaded  from https://www.amd.com/en/support.

Intel
^^^^^

Instructions for installing Intel OpenCL runtime can be found
`here <https://software.intel.com/content/www/us/en/develop/articles/opencl-drivers.html>`__.
Note that some distributions of Linux have native support for installing the
required packages through the system package manager. On Ubuntu distribution
the drivers and OpenCL runtime are available through the
:code:`intel-opencl-icd` and :code:`beignet-opencl-icd` packages.

Python
------

PyXOpto requires `Python 3 <https://www.python.org/>`_ (version >= 3.7)
programming language and at least the following extension packages:

* `Scipy <https://scipy.org/>`_ version >= 1.4.1
* `Matplotlib <https://matplotlib.org>`_ version >= 1.18.1
* `Numpy <www.nupy.org>`_ version >= 1.18.4
* `PyOpencl <https://documen.tician.de/pyopencl>`_ version >= 2019.1.2
* `Shapely <https://shapely.readthedocs.io/en/latest/>`_ version >= 1.6.4
* `Jinja2 <https://jinja2docs.readthedocs.io/en/stable/>`_ version >= 2.9.8
* `Numba <https://numba.pydata.org/>`_ version >= 0.34.0

Additional packages are required to use modules that enable scattering phase
function calculations for layered spherical particles:

* `Scattnlay <https://github.com/ovidiopr/scattnlay>`_ version >= 2.2

Note that the distribution is tested with the listed versions.

Windows OS
^^^^^^^^^^

On Windows OS you can conveniently use
`WinPython <https://sourceforge.net/projects/winpython/files/>`_,
a portable distribution of the Python programming language for Windows 7/8/10.
The `WinPython <https://sourceforge.net/projects/winpython/files/>`_
distribution comes with many commonly used Python extension packages including
the ones required by the PyXOpto package. The installation also comes with
`Spyder <https://www.spyder-ide.org/>`_ and
`Visual Studio Code <https://code.visualstudio.com/>`_ integrated development
environments.

Alternatively, use the `Anaconda <https://anaconda.org/anaconda/python>`_ 
Python distribution that is available for a wide range of operating systems.
Follow the `instructions <https://docs.anaconda.com/anaconda/navigator/tutorials/manage-packages/>`__
to install the required extension packages.

`Here <https://www.lfd.uci.edu/~gohlke/pythonlibs/>`_ you can find many 32-
and 64-bit Windows binaries of open-source Python extension packages.

The following two examples should install the Python 3 programming language and the required packages on Windows OS.

Anaconda Python Distribution:


1. Install the `Anaconda Individual Edition <​https://www.anaconda.com/products/individual-b>`_.
2. Run Anaconda Prompt.
3. Add *conda-forge* channel:

.. code-block:: bash

    conda config -add channels conda-forge

4. Create a new conda environment named *custom_env*:

.. code-block:: bash

    conda create -n custom_env

1. Install the required Python packages:

.. code-block:: bash

    conda install -n custom_env scipy matplotlib numpy pyopencl shapely jinja2 numba

6. The custom environment *custom_env* can be activated:

.. code-block:: bash

    activate custom_env

WinPython distribution:


1. Download and extract the WinPython distribution (available at https://winpython.github.io/ and make sure that it contains preinstalled packages scipy, matplotlib, numpy, shapely, jinja2 and numba)
2. The package *pyopencl* is not included and has to be downloaded separately from Gohlke’s website https://www.lfd.uci.edu/~gohlke/pythonlibs/#pyopencl. Download the latest version that corresponds to the extracted Python version and architecture.
3. Run WinPython Control Panel from the extracted WinPython folder.
4. Add and install the package *pyopencl* by locating the downloaded .whl file.

Linux OS
^^^^^^^^

On Linux operating systems, Python is very likely already installed.
If not, use the OS package manager to install Python 3.
Individual Python packages can be installed through the Python package manager
`pip <https://pypi.org/project/pip/>`_. 
Note that on some Linux distribution the Python 3 and related Pip executables
are named python3 and pip3. In this case, the python and pip executables are
reserved for the deprecated Python 2 programming language.

The following example should install the latest Python 3 programming language
and related packages on Ubuntu or Debian distributions:

.. code-block:: bash

    sudo apt update
    sudo apt install python3
    sudo apt install python3-pip

    sudo pip3 install scipy
    sudo pip3 install matplotlib
    sudo pip3 install numpy
    sudo pip3 install Shapely
    sudo pip3 install numba
    sudo pip3 install pyopencl
    sudo pip3 install Jinja2
    sudo pip3 install python-scattnlay

Alternatively, install one of the supported Python distributions such
as `Anaconda <https://anaconda.org/anaconda/python>`_ and follow
the `instructions <https://docs.anaconda.com/anaconda/navigator/tutorials/manage-packages/>`__
to install the required extension packages.

Integrated development environments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
There are numerous Integrated development environments (IDE) that support the
Python programming language. As already pointed out,
`Winpython <https://sourceforge.net/projects/winpython/files/>`_
distribution comes with `Spyder <https://www.spyder-ide.org/>`_ IDE
and `Visual Studio Code <https://code.visualstudio.com/>`_ that are also
supported by the Anaconda distribution.
Other popular full-featured IDEs supported across different operating systems:

* Eclipse with `PyDev <http://www.pydev.org/>`_
* `PyCharm <https://www.jetbrains.com/pycharm/>`_

PyXOpto
-------
First, download or clone the PyXOpto source repository to a local directory.
The source code can be installed as a Python package or used directly from
the downloaded source.

As a Python package
^^^^^^^^^^^^^^^^^^^

PyXOpto can be installed as a package using the setup.py file. Run the
following command from the root directory of PyXOpto (the one with the
setup.py file).

.. code-block:: bash

    python setup.py install

This will also install the dependencies that include several popular
Python packages (scipy, matplotlib, numpy, pyopencl, shapely, numba,
and jinja2).

.. note::

    When using Windows OS and Anaconda Pytohon distribution, PyXOpto should be installed by running the Anaconda Prompt and activating the custom environment created in the Python installation section. When using the WinPython distribution, the WinPython Command Prompt found in the extracted WinPython folder should be run instead.

The PyXOpto installation can be verified by running an example found in *examples* folder, e.g., examples/mcml/basic.py:

.. code-block:: bash

    python basic.py


Using from source
^^^^^^^^^^^^^^^^^

To use the PyXOpto package from source, you will have to manually install all
the Python packages listed in the setup.py file. Follow the
the installation instructions for :ref:`Python<installation-python-label>`.

You will also have to manually include the root directory of the PyXOpto
package into the Python search path. This can be conveniently 
accomplished through setting the :code:`PYTHONPATH` environment variable.

On Linux operating system use:

.. code-block:: bash

    export PTYTHONPATH=path/to/pyxopto:$PYTHONPATH

On Windows operating systems use:

.. code-block:: bat

    set PTYTHONPATH=path\to\pyxopto;%PYTHONPATH%

After installing the dependencies and setting the environment variable
:code:`PYTHONPATH`, you should be able to import PyXOpto.

Binaries
--------

PyXOpto does not require compilation of C/C++ sources. However, performance
critical parts of the random number generator module :py:mod:`xopto.cl.clrng`
can optionally benefit from a binary library. The build process is fully
automated and can be executed by running the :py:func:`xopto.rebuild` function.

.. code-block:: python

    import xopto
    xopto.rebuild()

The build process requires `GCC <https://gcc.gnu.org/>`_ or `The Microsoft C++ 
(MSVC) compiler toolset <https://docs.microsoft.com/en-us/cpp/build/building-on-the-command-line?view=msvc-160>`_.

Local PyXOpto data
------------------

PyXOpto will save a small amount of data (built libraries, temporary files and
performance-critical precalculated data) in the home directory
of the current user :code:`~/.xopto`. This location can be customized through
the :code:`PYXOPTO_USER_PATH` environment variable.

Environment variables
---------------------

Some properties of PyXOpto can be customized through the use of environment
variables:

  - :code:`PYXOPTO_VERBOSE=0` -
    Turns on (:code:`PYXOPTO_VERBOSE=1`) or off (:code:`PYXOPTO_VERBOSE=1`)
    the verbose initialization mode.
  - :code:`PYXOPTO_USER_PATH="~/.xopto"` -
    Location for storing temporary files, built libraries and
    performance-critical precalculated data.

Citing PyXOpto
--------------

We, the authors of PyXOpto, expect that the package is used in accordance with
the `GPL3+ <https://www.gnu.org/licenses/gpl-3.0-standalone.html>`_
license and that any work using the PyXOpto package also cites the project
and at least one of the following references:

#. P. Naglič, F. Pernuš, B. Likar, and M. Bürmen, *Limitations of the commonly*
   *used simplified laterally uniform optical fiber probe-tissue interface in*
   *Monte Carlo simulations of diffuse reflectance*, Biomed. Opt. Expres,
   **6** (10), 3973-3988 (2015), https://doi.org/10.1364/BOE.6.003973.

#. P. Naglič, F. Pernuš, B. Likar, and M. Bürmen, *Lookup table-based sampling*
   *of the phase function for Monte Carlo simulations of light propagation in*
   *turbid media*, Biomed. Opt. Expres, **8** (3), 1895-1910 (2017),
   https://doi.org/10.1364/BOE.8.001895.

For commercial use or alternative licensing of PyXOpto please contact us at
info@xopto.eu.
