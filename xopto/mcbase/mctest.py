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

import numpy as np

from xopto.mcbase import mcobject

class McTest:
    def __init__(self, mc_obj: mcobject.McObject,
                 reference: np.ndarray, description: str):
        '''
        Base class of Simulator tests.
    
        Parameters
        ----------
        mc: McObject
            Simulator instance.
        reference: np.ndarray
            Reference data expected to be produced by the test.
        description: str
            Description of the test.
        '''
        self._mc = mc_obj
        self._reference = np.asarray(reference)
        self._error = None
        self._description = str(description)

    def _get_mc(self) -> mcobject.McObject:
        return self._mc
    mc = property(_get_mc, None, None, 'Simulator instance.')

    def _get_description(self) -> str:
        return self._description
    def _set_description(self, description: str):
        self._description = str(description)
    description = property(_get_description, _set_description, None,
                           'Test description.')

    def _get_error(self) -> np.ndarray:
        return self._error
    def _set_error(self, error: np.ndarray):
        self._error = np.asarray(error)
    error = property(_get_error, _set_error, None, 'Test error.')

    def _get_reference(self) -> np.ndarray:
        return self._reference
    def _set_reference(self, data: np.ndarray):
        self._reference = data
    reference = property(_get_reference, _set_reference, None,
                         'Reference/expected results.')

    def _get_simulated(self) -> np.ndarray:
        return self._simulated
    def _set_simulated(self, result:np.ndarray):
        self._simulated = np.asarray(result)
    simulated = property(_get_simulated, _set_simulated, None,
                         'Simulated results.')

    def passed(self) -> bool:
        '''
        Implement this method in the derived calass. Return True if the
        test is successful or False if it fails.
        '''
        return False

    def progress_str(self, ndone: int, ntotal: int, tjob: float) -> str:
        '''
        Create a progress message.

        Parameters
        ----------
        ndone: int
            The number of completed tasks.
        ntotal: int
            The total number of tasks.
        job: float
            Time since start of the job.

        Returns
        report: str
            Progress report string.
        '''
        percentage = ndone/ntotal*100.0
        return '{:s}: {}/{} ({:.1f}%) in {:.3f} s'.format(
            self.description, ndone, ntotal, percentage, tjob)

