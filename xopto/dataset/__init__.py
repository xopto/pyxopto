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

from xopto import ROOT_PATH

DATASET_PATH = os.path.join(ROOT_PATH, 'dataset')

TEMPLATE_PATH = os.path.join(DATASET_PATH, 'template')
''' Directory with run templates. '''

MCML_TEMPLATE_PATH = os.path.join(TEMPLATE_PATH, 'mcml')
''' Directory with run templates for the layered MC. '''

MCVOX_TEMPLATE_PATH = os.path.join(TEMPLATE_PATH, 'mcvox')
''' Directory with run templates for the voxelized MC. '''

SV_TEMPLATE_PATH = os.path.join(TEMPLATE_PATH, 'sv')
''' Directory with run templates for sampling volume computations. '''
