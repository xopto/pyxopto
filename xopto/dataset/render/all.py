

# -*- coding: utf-8 -*-
################################ Begin license #################################
# Copyright (C) Laboratory of Imaging technologies,
#               Faculty of Electrical Engineering,
#               University of Ljubljana.
#
# This file is part of PyXOpto.
#
# PyXOpto is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PyXOpto is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PyXOpto. If not, see <https://www.gnu.org/licenses/>.
################################# End license ##################################

from .mcml_comparison import render_mcml_comparison
from .mcml import render_mcml
from .mcvox import render_mcvox_fluence
from .sv import render_sv_reflectance
from .common import prepare_cli, process_cli


def render_all(target_dir: str = None,
               cl_device: str = None, cl_index: int = 0,
               test: bool = False, verbose: bool = False):
    '''
    Render all templates using default configuration.

    Parameters
    ----------
    target_dir: str
        Root directory for the dataset. The scripts will be rendered into
        "run" subdirectory and the dataset data will be saved into the data
        subdirectory. If None, the parent directory of this file will serve
        as the root directory. 
    cl_device: str
        Default OpenCL device name or None. The value can be also set through
        the CL_DEVICE environment variable.
    cl_index: int
        OpenCL device index (if multiple OpenCL devices of the same kind
        are installed). The value can be also set through the CL_INDEX
        environment variable.
    test: bool
        Do a test run. The run scripts will be rendered but not saved. This
        option will automatically enable the verbose mode.
    verbose: bool
        Enables verbose reporting.
    '''
    kwargs = {'target_dir': target_dir, 'test': test, 'verbose': verbose,
              'cl_device': cl_device, 'cl_index': cl_index}
    render_mcml_comparison(**kwargs)
    render_mcml(**kwargs)
    render_sv_reflectance(**kwargs)
    render_mcvox_fluence(**kwargs)

if __name__ == '__main__':
    parser = prepare_cli('Render run scripts for all the datasets')
    # no additional command line arguments are required
    kwargs = process_cli(parser)
    render_all(**kwargs)
