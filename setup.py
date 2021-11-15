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

from setuptools import setup, find_packages

setup(
    name='PyXOpto',
    version='0.1.0',
    description='Python Monte Carlo light propagation tools',
    url='https://github.com/xopto/pyxopto',
    author='XOpto team',
    author_email='info@xopto.eu',
    license='GPLv3+',
    packages=find_packages(
        exclude=['PyXOpto.egg-info*', 'build*', 'dist*', 'maintenance*', 'docs*']),

    package_data={
        '': ['*.c', '*.cpp', '*.h', '*.npz', '*.pkl'],
    },
    #include_package_data=True,

    install_requires=[
        'scipy >= 1.4.1',
        'numpy >= 1.18.4',
        'matplotlib >= 1.18.1',
        'pyopencl >= 2019.1.2',
        'Shapely >= 1.6.4',
        'Jinja2 >= 2.9.6',
        'numba'
    ],

    classifiers=[
        'Development Status :: 1 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
