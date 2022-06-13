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

from typing import Type
import pickle
from xopto.pf import PfBase
from xopto import PICKLE_PROTOCOL
from xopto.mcbase.mcpf import LutEx


class ObjCache:
    @staticmethod
    def load(cachefile, **kwargs) -> 'ObjCache':
        '''
        Load ObjCache object from a file.

        Parameters
        ----------
        cachefile: str, file
            File location or a file like object.
        kwargs: dict
            Keyword arguments are passed to the :py:meth:`ObjCache`
            constructor.

        Returns
        -------
        cache: ObjCache
            ObjCache instance initialized with data from the specified file.
        '''
        return ObjCache(cachefile=cachefile, **kwargs)

    def __init__(self, cachefile: str = None, verbose: bool = False):
        '''
        ObjCache constructor.

        Parameters
        ----------
        cachefile: str, file
            File location or a binary file like object.
        verbose: bool
            Print information to stdout.
        '''
        self._cache = {}
        self._verbose = bool(verbose)

        if cachefile is not None:
            if isinstance(cachefile, str):
                with open(cachefile, 'rb') as fid:
                    cache = pickle.load(fid)
            else:
                cache = pickle.load(cachefile)

            if not isinstance(cache, dict):
                raise ValueError('The loaded object cache file is not valid!')
            self._cache = cache

    def save(self, cachefile: str):
        '''
        Save cache data to a file.

        Parameters
        ----------
        cachefile: str, file
            File location or a binary file like object.
        '''
        if isinstance(cachefile, str):
            with open(cachefile, 'wb') as fid:
                pickle.dump(self._cache, fid, protocol=PICKLE_PROTOCOL)
        else:
            pickle.dump(self._cache, cachefile, protocol=PICKLE_PROTOCOL)

    def verbose(self) -> bool:
        '''
        Returns the value of the verbose flag that was passed
        to the constructor.
        '''
        return self._verbose

    def get(self, obj_type: type, *args) -> object:
        '''
        Request an object of type obj_type and constructor arguments args
        from the cache. Object is created if one is not found in the cache.

        Parameters
        ----------
        obj_type: object type or callable
            Any python class or callable.
        args: list
            Parameters passed to the object constructor or callable.

        Returns
        -------
        obj: object
            An object created as :code:`obj_type(*args)` or an existing object
            from cache.
        '''
        obj_key = (obj_type, *args)

        obj = self._cache.get(obj_key)
        if obj is None:
            obj = obj_type(*args)
            self._cache[obj_key] = obj

        return obj

    def contains(self, obj_type: type, *args) -> bool:
        '''
        Checks if the specified object is in the cache.

        Parameters
        ----------
        obj_type: object type or callable
            Any python class or callable.
        args: list
            Parameters passed to the object constructor or callable.

        Returns
        -------
        found: bool
            Returns True if the specified object is found in the cache.
        '''
        obj_key = (obj_type, *args)
        return obj_key in self._cache

    def insert(self, value: object, obj_type: type, *args):
        '''
        Inserts the specified object into the cache if an existing
        entry is not found.

        Parameters
        ----------
        value: object
            Value that is attained by calling obj_type(*args)
        obj_type: object type or callable
            Any python class or callable.
        args: list
            Parameters passed to the object constructor or callable.

        Returns
        -------
        found: bool
            Returns True if the specified object was added to the cache or
            False if an existing object has been found in the cache.
        '''
        obj_key = (obj_type, *args)
        if obj_key not in self._cache:
            obj = obj_type(*args)
            self._cache[obj_key] = obj
            return True

        return False


class LutCache(ObjCache):
    @staticmethod
    def load(cachefile: str, **kwargs) -> 'LutCache':
        '''
        Load LutCache object from a file.

        Parameters
        ----------
        cachefile: str, file
            File location or a file like object.
        kwargs: dict
            Keyword arguments are passed to the LutCache constructor.

        Returns
        -------
        cache: ObjCache
            ObjCache instance initialized with data from the specified file.
        '''
        return LutCache(cachefile=cachefile, **kwargs)

    def __init__(self, lutsize: int = 2000, verbose: bool = False,
                 cachefile: str = None):
        '''
        Initializes a lookup table-based Monte Carlo scattering phase function
        cache (mcpf.LutPhaseFunction objects). This object internally saves
        all created instances of lookup table-based scattering phase
        functions and runs computations only if the requested scattering
        phase function cannot be found in the cache.

        Parameters
        ----------
        lutsize:int
            The lookup table size.
        verbose: bool
            Print information to stdout isf set to True.
        cachefile: str or file
            File location or a binary file like object from which to
            load LutCache.
        '''
        ObjCache.__init__(self, verbose=verbose, cachefile=cachefile)
        self._lutsize = int(lutsize)

        if cachefile is not None:
            if isinstance(cachefile, str):
                with open(cachefile, 'rb') as fid:
                    pickle.load(fid) # skip first
                    lut_options = pickle.load(fid)
            else:
                lut_options = pickle.load(cachefile)

            self._lutsize = int(lut_options['lutsize'])

    def save(self, cachefile: str):
        '''
        Save LutCache object to a file.

        Parameters
        ----------
        cachefile: str, file
            File location or a binary file like object.
        '''
        lut_options = {'lutsize': self._lutsize, 'type_name':'LutCache'}
        if isinstance(cachefile, str):
            with open(cachefile, 'wb') as fid:
                ObjCache.save(self, fid)
                pickle.dump(lut_options, fid, protocol=PICKLE_PROTOCOL)
        else:
            ObjCache.save(self, cachefile)
            pickle.dump(lut_options, cachefile, protocol=PICKLE_PROTOCOL)

    def get(self, pftype: Type[PfBase], pfargs: tuple, lutsize: int = None) \
            -> object:
        '''
        Prepare a new lookup table-based scattering phase function. Computations
        are made only if the requested scattering phase function is not found
        in the cache. This method can be used to pre-populate the cache with
        scattering phase functions.

        Parameters
        ----------
        pftype: Type[xopto.pf.PfBase]
            Any type of a scattering phase function that inherits from
            the :py:class:`xopto.pf.PfBase` class.
        pfargs: tuple
            Parameters passed to the scattering phase function constructor.
        lutsize: int
            Size of the lookup table passed to the
            py:meth:`xopto.mcbase.pf.lut.LutEx`.
            If None, the value passed to the constructor is used.

        Returns
        -------
        lut_pf: LutPhaseFunction
            Created or from cache loaded lookup table scattering phase function.
        '''
        if lutsize is not None:
            lutsize = int(lutsize)
        else:
            lutsize = self._lutsize

        pfargs = tuple(pfargs)

        if self.verbose():
            if ObjCache.contains(
                    self, LutEx, pftype, pfargs, lutsize):
                print(
                    'LUT scattering phase function:\n'\
                    '\t{}{},\n'\
                    '\twith lookup table size {} FOUND in cache!'.format(
                        pftype.__name__, str(pfargs), lutsize)
                )
            else:
                print(
                    'LUT scattering phase function:\n'\
                    '\t{}({}),\n'\
                    '\twith lookup table size {} NOT found in cache!'.format(
                        pftype.__name__, str(pfargs), lutsize)
                )

        obj = ObjCache.get(self, LutEx, pftype, pfargs, lutsize)

        return obj

    def insert(self, lut: LutEx, pftype: Type[PfBase], pfargs: tuple,
               lutsize: int = None) -> object:
        '''
        Inserts an existing lookup table-based scattering phase function into
        the cache.

        Parameters
        ----------
        lut: LutEx
            Instance of :py:class:`~xopto.mcbase.mcpf.LutEx` created with
            the specified scattering phase function type, arguments and
            lookup table size.
        pftype: Type[xopto.pf.PfBase]
            Any type of a scattering phase function that inherits from
            the :py:class:`xopto.pf.PfBase` class.
        pfargs: tuple
            Parameters passed to the scattering phase function constructor.
        lutsize: int
            Size of the lookup table passed to the
            py:meth:`xopto.mcbase.pf.lut.LutEx`.
            If None, the value passed to the constructor is used.

        Returns
        -------
        inserted: bool
            Return True if the scattering phase function was inserted into
            the cache and False if existing scattering phase function
            was found in the lookup table.
        '''
        if lutsize is not None:
            lutsize = int(lutsize)
        else:
            lutsize = self._lutsize

        pfargs = tuple(pfargs)

        return ObjCache.insert(self, lut, LutEx, pftype, pfargs, lutsize)
