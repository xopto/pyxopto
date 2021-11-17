# PyXOpto

PyXOpto is a collection of python tools for performing Monte Carlo simulations
of light propagation in turbid media using massively parallel processing on a wide range of OpenCL-enabled devices. The tools allow steady-state and time-resolved simulations of light propagation, deposition and fluence simulations, tracing and filtering of photon packet paths, computation of sampling volumes for a number of source-detector configurations, support arbitrary scattering phase functions and are easy to customize or extend.

![Deposit simulations](/docs/source/static/animation/xopto_fluence.gif)
<br>*Energy deposit simulations for a voxelized volume with a laterally moving Gaussian beam.*

![Sampling volume trajectories](/docs/source/static/animation/paths_transmittance_128.gif)
<br>*Transmittance trajectories of photon packets on the way from the source to the laterally displaced detector optical fiber.*

![Sampling volume in transmittance configuration - source-detector separation](/docs/source/static/animation/sds_transmittance_10000.gif)
<br>*Sampling volume in transmittance configuration as a function of the lateral displacement between the source and detector optical fiber.*

![Sampling volume in reflectance configuration - source-detector separation](/docs/source/static/animation/sds_reflectance_10000.gif)
<br>*Sampling volume in reflectance configuration as a function of the distance between the source and detector optical fiber.*

![Sampling volume in reflectance configuration - incidence angle](/docs/source/static/animation/incidence_angle_100000.gif)
<br>*Sampling volume in reflectance configuration as a function of the incidence angle of the source optical fiber.*

# Documentation

Full documentation of PyXOpto is available [here](https://xopto.github.io/pyxopto/).

## Installation

Detailed instructions are available [here](https://xopto.github.io/pyxopto/docs/html/source/installation.html).

### Python
PyXOpto requires a Python 3 installation. Most of the Linux OS distributions will come with a preinstalled Python 3. On Windows OS, the easiest way to install the Python 3 programming language is to use the [WinPython](https://github.com/winpython) or [Anaconda](https://docs.anaconda.com/anaconda/install/windows/) distributions. There are numerous integrated development environments that work with Python, among these [Visual Studio Code](https://code.visualstudio.com) and [PyCharm](https://www.jetbrains.com/pycharm/) are two popular cross-platform options. The [WinPython](https://sourceforge.net/projects/winpython/files/WinPython_3.9/3.9.4.0/) distributions can be downloaded with an embedded and preconfigured Visual Studio Code (e.g. [Winpython64-3.9.4.0cod.exe](https://sourceforge.net/projects/winpython/files/WinPython_3.9/3.9.4.0/Winpython64-3.9.4.0cod.exe/download)).

### PyXOpto
First, download or clone the PyXOpto source repository to a local directory. The source code can be installed as python package or used independently from the dowloaded source.

#### As a python package
PyXOpto can be installed as a package using the setup.py file. Run the
following command from the root directory of PyXOpto (the one with the setup.py file).
```bash
python setup.py install
```
This will also install the dependencies that include several popular Python packages (scipy, matplotlib, numpy, pyopencl, shapely, numba,
and jinja2) dependencies.

### Using from source
To use the PyXOpto package from source, you will have to manually install all the python dependencies listed in the setup.py file (scipy, matplotlib, numpy, pyopencl, shapely, numba, jinja2, ...). The easiest way to install the dependencies is to use the Python package installer [pip](https://pypi.org/project/pip/). Note that the WinPython distribution will likely come with many if not all the dependencies already installed. Also note that on some Linux distributions, the Python 3 executable is named `python3` and `python` is used for the deprecated Python 2.
You will also have to manually include the root directory of the PyXOpto package into the Python search path. This can be conveniently 
accomplished through setting the `PYTHONPATH` environment variable.
On Linux operating system use:
```bash
export PTYTHONPATH=path/to/pyxopto:$PYTHONPATH
```
On Windows operating systems use:
```cmd
set PTYTHONPATH=path\to\pyxopto;%PYTHONPATH%
```
After installing the dependencies and setting the environment variable
`PYTHONPATH`, you should be able to import PyXOpto.

## Basic Monte Carlo simulations of a layered medium
This basic example is available in the `examples/mcml/basic_example.py` file.

The functionality of layered Monte Carlo (MC) is accessible through the
`xopto.mcml.mc` module. First, create an empty file, 
e.g. `basic_example.py` and import the `xopto.mcml.mc` module.
```py
from xopto.mcml import mc
```

### The layer stack
The layers of the medium can be defined through the `mclayer` submodule. The
layers are stacked along the positive direction of the z coordinate axis. 

The topmost and bottommost layers of the stack are used to describe the
medium that surrounds the sample at the top and at the bottom surfaces, respectively.
Therefore, at least three layers must be always defined,
namely the two layers of the surrounding medium and one sample layer!

The bottom surface of the topmost layer (the surrounding medium) is
located at coordinate z=0. The positive direction of the z axis points in the
direction of the sample layer stack.

The thicknesses of the topmost and bottommost layers will be
automatically set to infinity regardless of the specified layer thickness.

Note that all the layers in the stack must use the same scattering phase
function model. A variety of scattering phase function models is available through the
`mcpf` submodule.

An example of a basic turbid sample of thickness d=10.0&nbsp;mm, with an absorption
coefficient mua=1.0&nbsp;1/cm, scattering coefficient mus=50.0&nbsp;1/cm,
a Henyey-Greenstein scattering phase function (`mc.mcph.Hg`) with an anisotropy g=0.8 and
refractive index 1.33 is as follows.
```py
layers = mc.mclayer.Layers(
    [
        mc.mclayer.Layer(n=1.00, d=np.inf,  mua=1.0e2, mus=50.0e2, pf=mc.mcpf.Hg(0.0)),
        mc.mclayer.Layer(n=1.33, d=10.0e-3, mua=1.0e2, mus=50.0e2, pf=mc.mcpf.Hg(0.8)),
        mc.mclayer.Layer(n=1.00, d=np.inf,  mua=1.0e2, mus=50.0e2, pf=mc.mcpf.Hg(0.0))
    ]
)
```
Note that the absorption coefficient `mua`, scattering coefficient `mus` and the
scattering phase function `pf` of the topmost and bottommost layers are not
used in the MC simulations, since the photon packets are not propagated through
the surrounding medium. However, the refractive index `n` of the two outermost
layers is used to properly refract/reflect the photon packet at the layer
boundaries when launched by the source or when escaping the sample.
The value of the layer thickness `d` should be given in m and
the values of the scattering `mus` and absorption `mua` coefficient in 1/m.

### The photon packet source
Different sources of photon packets are available through the `mcsource`
submodule. The following example creates a basic line source (infinitely thin) at the top sample surface
(x, y, z)=(0, 0, 0) with a perpendicular incidence (0, 0, 1).
```
source = mc.mcsource.Line()
```

### The detectors
The photon packets can be collected by a surface detector after exiting the
top or bottom sample surface. Different types of surface detectors are
available through the `mcdetector` submodule. Note that the top and bottom
sample surface can use different configurations and/or types of detectors.
```py
detectors = mc.mcdetector.Detectors(
    top = mc.mcdetector.Radial(
        mc.mcdetector.Axis(0.0, 10.0e-3, 1000, cosmin=np.deg2rad(20))
    ),
    bottom = mc.mcdetector.Radial(
        mc.mcdetector.Axis(0.0, 10.0e-3, 100)
    )
)
```
In the above example, we create two radial detectors one at the top and one at the bottom sample surface. The spacing between the
concentric accumulators of the radial detector at the top sample surface is set to 10&nbsp;μm, while the spacing of the concentric accumulators at the bottom sample surface is set to 100&nbsp;μm. Both detectors are accumulating photon packets from 0.0&nbsp;mm to 10.0&nbsp;mm. The detector at the top sample
surface only collects photon packets that exit the sample within 20&deg; of
the surface normal, while the detector at the bottom sample surface collects
all the photon packets that exit the sample.

### The OpenCL device
The OpenCL device that will run the MC simulations can be selected through the
`xopto.cl.clinfo` module. In the following example we pick the first available
GPU device.
```py
gpu = clinfo.gpu()
```

### The Monte Carlo simulator
Next, we create a Monte Carlo simulator instance from the defined layers,
photon packet source and detectors.
```py
mc_obj = mc.Mc(layers, source, detectors, cl_devices=gpu)
``` 
Optionally, we can limit the maximum simulation radius that is measured from
the position of the photon packet source. In this example, we limit the
simulation radius to 25&nbsp;mm.
```
mc_obj.rmax = 25.0e-3
```
Finally, we can run the simulator instance with a given number of photon
packets (10,000,000 in this example) and collect the results. The simulator
returns three objects/results, namely the trace, fluence and detectors. Since
in this basic example we only use the surface detectors, the remaining two
results (fluence and trace) will be set to `None`. 
```
trace_res, fluence_res, detectors_res = mc_obj.run(10e6)
```
Note that the photon packets that exit the sample within the acceptance cone
but at a distance/radius that exceeds the maximum radius of the
detector will be accumulated in the last concentric ring.

### Visualizing the results
We can plot the simulation results using the `mtplotlib.pyplot` module. For a better visualization of the reflectance/transmittance a logarithmic scale is used in the y axis of the plots.
```
fig, (ax1, ax2) = pp.subplots(2, 1)

ax1.semilogy(detectors_res.top.r*1e3, detectors_res.top.reflectance)
ax1.set_xlabel('Distance from source (mm)')
ax1.set_ylabel('Reflectance')

ax2.semilogy(detectors_res.bottom.r*1e3, detectors_res.bottom.reflectance)
ax2.set_xlabel('Distance from source (mm)')
ax2.set_ylabel('Reflectance')

pp.show()
```
Reflectance and transmittance collected by the surface detectors.
![](/docs/source/images/examples/mcml/basic/detector.svg)

### The complete example

[basic.py](/examples/mcml/basic.py)

You can run this example from the root directory of the PyXopto package as:
```bash
python examples/mcml/basic.py
```

# MC Dataset

Precomputed datasets of reflectance, transmittance, energy deposition and sampling volume for a number of different source, detector and sample
configurations are available through a separate repository
[MC Dataset](https://github.com/xopto/mcdataset).

# Citing PyXOpto

We, the authors of PyXOpto, expect that the package is used in accordance with
the [GPL3+](https://www.gnu.org/licenses/gpl-3.0-standalone.html)
license and that any work using the PyXOpto package also cites the project
and at least one of the following references:

 - P. Naglič, F. Pernuš, B. Likar, and M. Bürmen, *Limitations of the commonly*
   *used simplified laterally uniform optical fiber probe-tissue interface in*
   *Monte Carlo simulations of diffuse reflectance*, Biomed. Opt. Expres,
   **6** (10), 3973-3988 (2015), https://doi.org/10.1364/BOE.6.003973.

 - P. Naglič, F. Pernuš, B. Likar, and M. Bürmen, *Lookup table-based sampling*
   *of the phase function for Monte Carlo simulations of light propagation in*
   *turbid media*, Biomed. Opt. Expres, **8** (3), 1895-1910 (2017),
   https://doi.org/10.1364/BOE.8.001895.

For alternative licensing options of PyXOpto please contact us at info@xopto.eu.
