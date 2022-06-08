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

import os
import sys

__version__ = '0.2.0'


VERBOSE = False
''' Verbose reporting. '''

if 'PYXOPTO_VERBOSE' in os.environ and os.environ['PYXOPTO_VERBOSE']:
    VERBOSE = bool(os.environ['PYXOPTO_VERBOSE'])

if VERBOSE:
    print('Initializing PyXOpto {}'.format(__version__))
    os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'

ROOT_PATH = os.path.dirname(os.path.abspath(__file__))
''' The package directory. '''

PICKLE_PROTOCOL = 4
''' Default pickle protocol of the package. '''

DATA_PATH = os.path.join(ROOT_PATH, 'data')
''' Main data storage directory. '''

BIN_PATH = os.path.join(ROOT_PATH, 'bin')
''' Directory of binary/executable files. '''

PRIMES_PATH = os.path.join(DATA_PATH, 'primes')
''' Location of OpenCL random number generator seeds. '''


SRC_PATH = os.path.join(ROOT_PATH, 'src')
''' Location of compilable source code. '''

USER_TMP_PATH = os.path.join(ROOT_PATH, 'tmp')
''' Temporary data path. '''

MCBASE_PATH = os.path.join(ROOT_PATH, 'mcbase')
''' Root directory of MC base. '''

MCCYL_PATH = os.path.join(ROOT_PATH, 'mccyl')
''' Root directory of MC Cyl. '''

MCML_PATH = os.path.join(ROOT_PATH, 'mcml')
''' Root directory of MC ML. '''

MCVOX_PATH = os.path.join(ROOT_PATH, 'mcvox')
''' Root directory of MC VOX. '''

UTIL_PATH = os.path.join(ROOT_PATH, 'util')
''' Root directory of MC utilities. '''

USER_PATH = os.path.join(os.path.expanduser('~'), '.xopto', 'pyxopto')
''' Root path of user data storage. '''
if 'PYXOPTO_USER_PATH' in os.environ and os.environ['PYXOPTO_USER_PATH']:
    USER_PATH = str(os.environ['PYXOPTO_USER_PATH'])

USER_DATA_PATH = os.path.join(USER_PATH, 'data')
''' User directory for precalculated data. '''

USER_BIN_PATH = os.path.join(USER_PATH, 'bin')
''' User directory for built libraries or executables. '''

USER_TMP_PATH = os.path.join(USER_PATH, 'tmp')
''' User directory for temporary data. '''

USER_BUILD_PATH = os.path.join(USER_PATH, 'build')
''' User directory for build. '''

USER_PRIMES_PATH = os.path.join(USER_DATA_PATH, 'primes')
''' User directory for RNG primes. '''

def fulfills_min(requirement: str) -> bool:
    '''
    Checks if this version of PyXOpto fulfills the minimum requiremed version.

    Parameters
    ----------
    required: str
        Required version as a string, e.g. "0.2.1".

    Returns
    -------
    fulfills: bool
        Returns True if PyXOpto fullfills the requirement, else False.

    Note
    ----
    PyXOpto version is composed of three numbers (major, minor, patch), e.g.
    "0.2.0".
    '''
    this_version = [int(item) for item in __version__.split('.')]
    for index, item in enumerate(requirement.split('.')):
        item = item.strip()
        if index >= len(this_version):
            break
        if item == '*':
            continue
        if this_version[index] < int(item):
            return False

    return True


def make_user_dirs():
    '''
    Create all the user directories for data, binary and temporary files.
    '''
    for dir in (USER_DATA_PATH, USER_BIN_PATH, USER_TMP_PATH,
                USER_BUILD_PATH, USER_PRIMES_PATH):
        if not os.path.isdir(dir):
            try:
                os.makedirs(dir)
            except OSError:
                pass

make_user_dirs()

def rebuild(envvars: dict = None, verbose: bool = False):
    '''
    Rebuilds all the external dependencies of xopto (rng).

    Parameters
    ----------
    envvars: dict
        Additional environmental variables for the build process.

    verbose: bool
        Enables verbose output.

    Examples
    --------
    Build with default settings and external library installation paths.

    >>> import xopto
    >>> xopto.rebuild()
            Compiling rng.cpp ...
            Linking rng.o ...
            Moving rng64.so to xopto/bin ...
            Cleanup ...
            All done.
    >>>
    '''
    from xopto.util import build
    import platform

    if not os.path.isdir(USER_BIN_PATH):
        try:
            os.makedirs(USER_BIN_PATH)
        except OSError:
            pass

    if not os.path.isdir(USER_BUILD_PATH):
        try:
            os.makedirs(USER_BUILD_PATH)
        except OSError:
            pass

    env = None
    if envvars is not None:
        for key in envvars:
            os.environ[key] = envvars[key]
        env = dict(os.environ)

    srclist=[os.path.join(ROOT_PATH, 'src', 'rng', 'rng.cpp')]
    inclist = [os.path.join(ROOT_PATH, 'src', 'rng')]
    moveto = USER_BIN_PATH
    os_num_bits = build.os_bits()

    if platform.system() == 'Linux':
        build_target = 'rng{}.so'.format(os_num_bits)

        build.build_so_dll(
            srclist=srclist,
            inclist=inclist,
            target=build_target, moveto=moveto,
            envvars=envvars, cwd=USER_BUILD_PATH,
            verbose=verbose
        )

    elif platform.system() == 'Windows':
        build_target = 'rng{}.dll'.format(os_num_bits)

        build.build_so_dll(
            srclist=srclist,
            inclist=inclist,
            target=build_target, moveto=moveto,
            envvars=envvars, cwd=USER_BUILD_PATH,
            verbose=verbose
        )
    else:
        raise RuntimeError('Rebuild supported only on Linux and Windows OS!')

def rebuild_(envvars: dict = None):
    '''
    Rebuilds all the external dependencies of xopto (rng).

    Parameters
    ----------
    envvars: dict
        Additional environmental variables for the build process.

    Examples
    --------
    Build with default settings and external library installation paths.

    >>> import xopto
    >>> xopto.rebuild()
            Compiling rng.cpp ...
            Linking rng.o ...
            Moving rng64.so to xopto/bin ...
            Cleanup ...
            All done.
    >>>
    '''
    import subprocess
    import platform
    import math
    import sys

    if os.path.isdir(BIN_PATH):
        try:
            os.makedirs(BIN_PATH)
        except OSError:
            pass

    env = None
    if envvars is not None:
        for key in envvars:
            os.environ[key] = envvars[key]
        env = dict(os.environ)

    if platform.system() == 'Linux':
        cwd = os.path.join(ROOT_PATH, SRC_PATH, 'rng')
        if math.log(sys.maxsize)/math.log(2) > 32:
            proc = subprocess.Popen(
                ['./build64.sh'], cwd=cwd, env=env,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        else:
            proc = subprocess.Popen(
                ['./build32.sh'], cwd=cwd, env=env,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in proc.stdout:
            if isinstance(line, bytes):
                line = line.decode('utf8')
            sys.stdout.write(line)
    elif platform.system() == 'Windows':
        cwd = os.path.join(ROOT_PATH, SRC_PATH, 'rng')
        cwd = cwd.replace('/', '\\')
        if math.log(sys.maxsize)/math.log(2) > 32:
            proc = subprocess.Popen(
                ['build64.bat', 'exit'], cwd=cwd, shell=True,
                env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        else:
            proc = subprocess.Popen(
                ['build32.bat', 'exit'], cwd=cwd, shell=True,
                env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in proc.stdout:
            if isinstance(line, bytes):
                line = line.decode('utf8')
            sys.stdout.write(line)
    else:
        raise RuntimeError('Rebuild supported only on Linux and Windows OS!')
