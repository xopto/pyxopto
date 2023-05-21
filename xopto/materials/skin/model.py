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

from typing import Callable, Tuple, List

import numpy as np
from scipy.interpolate import interp1d

from xopto.materials import absorption, ri
from xopto.mcml import mc


class TissueDatabase:
    '''
    Default absorption coefficients (1/m) of tissue chromophores as a function
    of wavelength of light (m).
    '''

    oxyhem_mua = absorption.oxyhem.default
    ''' Absorption coefficient (1/m) of oxyhemoglobin as a function of wavelength.'''
    deoxyhem_mua = absorption.deoxyhem.default
    ''' Absorption coefficient (1/m) of deoxyhemoglobin as a function of wavelength.'''
    melanin_mua = absorption.melanin.default
    ''' Absorption coefficient (1/m) of melanin as a function of wavelength.'''
    carotenoid_mua = absorption.carotenoid.default
    ''' Absorption coefficient (1/m) of carotenoid as a function of wavelength.'''
    water_mua = absorption.water.default
    ''' Absorption coefficient (1/m) of water as a function of wavelength.'''
    fat_mua = absorption.fat.default
    ''' Absorption coefficient (1/m) of fat as a function of wavelength.'''


class ConstantValue:
    def __init__(self, value: float or int):
        '''
        A constant value that is internally represented as a function
        of wavelength. 
        '''
        self._value = value

    def __call__(self, wavelength):
        if hasattr(wavelength, '__len__'):
            value = np.zeros_like(wavelength)
            value.fill(self._value)
        else:
            value = self._value

        return value

    def _get_value(self) -> float or int:
        return self._value
    def _set_value(self, value: int or float):
        self._value = value
    value = property(_get_value, _set_value, None, 'Constant value.')

    def todict(self):
        return self._value

    @classmethod
    def fromdict(cls, data: dict or float):
        if isinstance(data, (float, int)):
            return cls(data)

        data_ = dict(data)
        T = data_.pop('type')
        if T != cls.__name__:
            raise TypeError('Expected data for type {} but got {}!'.format(
                cls.__name__, T))

        return cls(**data_)

    def __float__(self):
        return float(self._value)

    def __int__(self):
        return int(self._value)

    def __str__(self):
        return '{} # {}'.format(repr(self), id(self))

    def __repr__(self):
        return 'ConstantValue({})'.format(self._value)


class Polynomial:
    def __init__(self, *args, deg: int = None):
        '''
        Polynomial model of one parameter.
        Passing a single float value to this constructor will result
        in a constant scalar return value.

        Parameters
        ----------
        *args: float or Tuple[float, float]
            Polynomial coefficients as args[-1] + args[-2]*x + ... or
            coordinates (x, y) of points used to construct a polynomial.
        deg: int
            Degree of the polynomial to use when the arguments are data
            coordinates. If None, the degree is determined by the number
            of points.
        '''
        if isinstance(args[0], Polynomial):
            obj = args[0]
            args = obj.pk
        elif isinstance(args[0], (tuple, list)):
            tmp = np.stack(args)
            x, y = tmp[:, 0], tmp[:, 1]
            if deg is None:
                deg = len(x) - 1
            args = np.polyfit(x, y, deg).tolist()

        self._pk = None
        self._set_pk(args)

    def _set_pk(self, pk):
        if not hasattr(pk, '__len__'):
            pk = [pk,]
        self._pk = np.array(pk, dtype=np.float64).flatten()
    def _get_pk(self):
        return self._pk.tolist()
    pk = property(_get_pk, _set_pk, None,
                  'Polynomial coefficients as pk[-1] + pk[-2]*x + ...')

    def __call__(self, x: float or np.ndarray) -> float or np.ndarray:
        '''
        Evaluate the polynomial function at the given argument x.

        Parameters
        ----------
        x: float or np.ndarray
            Argument of the polynomial function.
        '''
        return np.polyval(self._pk, x)

    def todict(self) -> dict:
        '''
        Export model parameters to a dict.

        Returns
        -------
        model: dict
            Model as a dict with parameters named as in the constructor.
        '''
        return {'type': self.__class__.__name__, 'pk':self.pk}

    @classmethod
    def fromdict(cls, data: dict) -> 'Polynomial':
        '''
        Create a new instance from data in a dict.

        Parameters
        ----------
        data: dict
            Data for the new instance in a dict (as returned by the
            :py:meth:`Polynomial.todict` method).

        Returns
        -------
        obj: Polynomial
            A new instance.
        '''
        data_ = dict(data)
        T = data_.pop('type')
        if T != cls.__name__:
            raise ValueError('Expected data for type {} but got {}!'.format(
                cls.__name__, T))

        return cls(**data_)

    def update(self, data: dict):
        '''
        Update model parameters from dict.

        Parameters
        ----------
        data: dict
            Model parameters in a dict named as in the constructor. Parameters
            of the model that are not defined in the dict are left unchanged.
        '''
        if data.get('pk') is not None:
            self._set_pk(data.get('pk'))

    def __str__(self):
        n = self._pk.size
        parts = []
        for index, item in enumerate(self._pk):
            if n - index - 1 > 1:
                parts.append('{:g} x^{:d}'.format(item, n - index - 1))
            elif n - index - 1 == 1:
                parts.append('{:g} x'.format(item))
            else:
                parts.append('{:g}'.format(item))

        return 'p(x) = ' + ' + '.join(parts)

    def __repr__(self):
        args = ', '.join(['{:f}'.format(item) for item in self._pk])
        return 'Polynomial({:s}) # id {}'.format(args, id(self))


class ReducedScatteringPower:
    def __init__(self, a: float = 0.0, b: float = 0.0,
                 wavelength: float = 550e-9):
        '''
        Model of the  scattering coefficient that follows equation
        :math:`a(\lambda/\lambda_0)^{-b}`,
        where :math:`a` is the reduced scattering coefficient at
        wavelength :math:`\lambda_0` and :math:`b` is the scattering
        power.
        '''
        if isinstance(a, ReducedScatteringPower):
            obj = a
            a = obj.a
            b = obj.b
            wavelength = obj.wavelength

        self._a = self._b = self._wavelength = None
        self._set_a(a)
        self._set_b(b)
        self._set_wavelength(wavelength)

    def _set_a(self, a: float):
        self._a = float(a)
    def _get_a(self) -> float:
        return self._a
    a = property(
        _get_a, _set_a, None,
        'Reduced scattering coefficient (1/m) at the reference wavelength.')

    def _set_b(self, b: float):
        self._b = float(b)
    def _get_b(self) -> float:
        return self._b
    b = property(_get_b, _set_b, None,
                 'Power of the reduced scattering model.')

    def _set_wavelength(self, wavelength: float):
        self._wavelength = float(wavelength)
    def _get_wavelength(self) -> float:
        return self._wavelength
    wavelength = property(
        _get_wavelength, _set_wavelength, None,
        'Reference wavelength (m) in the reduced scattering model.')

    def todict(self) -> dict:
        '''
        Export model parameters to dict.

        Returns
        -------
        model: dict
            Model as a dict with parameters named as in the constructor.
        '''
        return {'type': self.__class__.__name__,
                'a':self.a, 'b':self.b, 'wavelength':self.wavelength}

    @classmethod
    def fromdict(cls, data: dict) -> 'ReducedScatteringPower':
        '''
        Create new model instance from data in the dict.

        Parameters
        ----------
        data: dict
            Model parameters as returned by the
            :py:meth:`ReducedScatteringPower.todict` method.

        Returns
        -------
        obj: ReducedScatteringPower
            A new instance of :py:class:`ReducedScatteringPower`.
        '''
        data_ = dict(data)
        T = data.pop('type')
        if T != cls.__name__:
            raise TypeError('Unexpected data for type {} but got {}!'.format(
                cls.__name__, T))

        return cls(**data_)

    def update(self, data: dict):
        '''
        Update model parameters from a dict. Only parameters
        that can be found in the dict are updated.

        Parameters
        ----------
        data: dict
            Values of model parameters that will be update. Note that only
            the parameters present in the dict will be updated.
            Model parameters in a dict. Parameters should be named as in the
            constructor.
        '''
        if data.get('a') is not None:
            self._set_a(data.get('a'))
        if data.get('b') is not None:
            self._set_b(data.get('b'))
        if data.get('wavelength') is not None:
            self._set_wavelength(data.get('wavelength'))

    def __call__(self, wavelength: float or np.ndarray) -> float or np.ndarray:
        '''
        Alias for the :py:meth:`ReducedScatteringPower.musr` method.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength of light (m) at which to compute the reduced scattering
            coefficient.

        Returns
        -------
        musr: float or np.ndarray
            The value of reduced scattering coefficient (1/m) at the
            specified wavelength.
        '''
        return self.musr(wavelength)

    def musr(self, wavelength: float or np.ndarray) -> float or np.ndarray:
        '''
        Reduced scattering coefficient as a function of wavelength.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength of light (m) at which to compute the reduced scattering
            coefficient.

        Returns
        -------
        musr: float or np.ndarray
            The value of reduced scattering coefficient (1/m) at the
            specified wavelength.
        '''
        return self._a*(wavelength/self._wavelength)**(-self._b)

    def __str__(self):
        return 'ReducedScatteringPower(a={:f}, b={:f}, wavelength={:g})'.format(
            self.a, self.b, self.wavelength)

    def __repr__(self):
        return '{:s} # id {}'.format(self.__str__(), id(self))


class DermisBaselineAbsorption:
    def __init__(self, kwater: float = 1.0, gamma: float = 0.5):
        '''
        Dermis baseline absorption excluding all standard chromophores, such
        as melanin, carotenoid, water, fat and blood.

        Parameters
        ----------
        kwater: float
            Correction factor related to the water content.
            Default value is set to 1.0.
        gamma: float
            Default value is set to 0.5.
        '''
        self._kwater = self._gamma = None
        self._set_gamma(gamma)
        self._set_kwater(kwater)

    def _set_gamma(self, gamma: float):
        self._gamma = float(gamma)
    def _get_gamma(self) -> float:
        return self._gamma
    gamma = property(_get_gamma, _set_gamma, None,
                     'Baseline absorption parameter gamma.')

    def _set_kwater(self, k: float):
        self._kwater = float(k)
    def _get_kwater(self) -> float:
        return self._kwater
    kwater = property(_get_kwater, _set_kwater, None,
                      'Baseline absorption absorption parameter water/water0.')

    def __call__(self, wavelength: float or np.ndarray) -> float or np.ndarray:
        '''
        Alias for the :py:meth:`DermisBaselineAbsorption.mua` method.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength of light (m).

        Returns
        -------
        mua: float or np.ndarray
            Absorption coefficient (1/m) at the specified wavelength.
        '''
        return self.mua(wavelength)

    def mua(self, wavelength: float or np.ndarray) -> float or np.ndarray:
        '''
        Absorption coefficient (1/m) as a function of wavelength of light (m).

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength of light (m)

        Returns
        -------
        mua: float or np.ndarray
            Absorption coefficient (1/m) at the specified wavelength.
        '''
        return 1e2*self._gamma*self._kwater* \
            (0.244 + 16.82*np.exp(-(1e9*wavelength - 400) / 80.5))

    def todict(self) -> dict:
        '''
        Export model parameters to dict.

        Returns
        -------
        data: dict
            Model parameters as passed to the constructor in a dict.
        '''
        return {'type': self.__class__.__name__,
                'gamma':self._gamma, 'kwater':self._kwater}

    @classmethod
    def fromdict(cls, data: dict) -> 'DermisBaselineAbsorption':
        '''
        Creates a new model instance from data in the dict.

        Parameters
        ----------
        data: dict
            Model data as returned by the
            :py:meth:`DermisBaselineAbsorption.todict` method.

        Returns
        -------
        obj:  DermisBaselineAbsorption
            A new instance of :py:class:`DermisBaselineAbsorption`.
        '''
        data_ = dict(data)
        T = data_.pop('type')
        if T != cls.__name__:
            raise TypeError('Unexpected data for type {} but got {}!'.format(
                cls.__name__, T))

        return cls(**data_)

    def update(self, data:dict):
        '''
        Update model parameters from a dict. Only parameters
        that can be found in the dict are updated.

        Parameters
        ----------
        data: dict
            Model parameters in a dict. Parameters should be named as in the
            constructor.
        '''
        if data.get('gamma') is not None:
            self._set_gamma(data.get('gamma'))
        if data.get('kwater') is not None:
            self._set_kwater(data.get('kwater'))

    def __str__(self):
        return 'DermisBaselineAbsorption(kwater={:f}, gamma={:f})'.format(
            self._kwater, self._gamma)

    def __repr__(self):
        return '{:s} # id {}'.format(self.__str__(), id(self))


class EpidermisBaselineAbsorption:
    def __init__(self, gamma: float = 0.5):
        '''
        Epidermis baseline absorption excluding all standard chromophores, such
        as melanin, carotenoid, water, fat and blood.

        Parameters
        ----------
        gamma: float
            Default value is set to 0.5.
        '''
        self._gamma = None
        self._set_gamma(gamma)

    def _set_gamma(self, gamma: float):
        self._gamma = float(gamma)
    def _get_gamma(self) -> float:
        return self._gamma
    gamma = property(_get_gamma, _set_gamma, None,
                     'Baseline absorption parameter gamma.')

    def __call__(self, wavelength: float) -> float:
        '''
        Alias for the mua method.
        '''
        return self.mua(wavelength)

    def mua(self, wavelength: float) -> float:
        '''
        Absorption coefficient (1/m) as a function of wavelength of light (m).

        Parameters
        ----------
        wavelength: float, np.ndarray
            Wavelength of light (m).

        Returns
        -------
        mua: float, np.ndarray
            Absorption coefficient (1/m) at the specified wavelength.
        '''
        return 1e2*self._gamma*\
            (0.244 + 85.3*np.exp(-(1e9*wavelength - 154)/66.2))

    def todict(self) -> dict:
        '''
        Export model parameters to dict.

        Returns
        -------
        data: dict
            Model parameters as passed to the constructor in a dict.
        '''
        return {'type': self.__class__.__name__, 'gamma':self._gamma}


    @classmethod
    def fromdict(cls, data: dict) -> 'EpidermisBaselineAbsorption':
        '''
        Creates a new model instance from data in the dict.

        Parameters
        ----------
        data: dict
            Model data as returned by the
            :py:meth:`EpidermisBaselineAbsorption.todict` method.

        Returns
        -------
        obj:  EpidermisBaselineAbsorption
            A new instance of :py:class:`EpidermisBaselineAbsorption`.
        '''
        data_ = dict(data)
        T = data_.pop('type')
        if T != cls.__name__:
            raise TypeError('Unexpected data for type {} but got {}!'.format(
                cls.__name__, T))

        return cls(**data_)

    def update(self, data: dict):
        '''
        Update model parameters from a dict. Only parameters
        that can be found in the dict are updated

        Parameters
        ----------
        data: dict
            Model parameters a dict. Parameters should be named as in the
            constructor.
        '''
        if data.get('gamma') is not None:
            self._set_gamma(data.get('gamma'))


    def __str__(self):
        return 'EpidermisBaselineAbsorption(gamma={:f})'.format(self._gamma)

    def __repr__(self):
        return '{:s} # id {}'.format(self.__str__(), id(self))


class SkinLayer:
    def __init__(self, d: float = 1e-3,
                 n: None or float or Callable[[float], float] = None,
                 musr: float or Callable[[float], float] = 15.0e2,
                 g: float or Callable[[float], float] = 0.8,
                 baseline_absorption: float or EpidermisBaselineAbsorption or
                                      DermisBaselineAbsorption = 0.0,
                 melanin: float = 0.0,
                 carotenoid: float = 0.0,
                 water: float = 0.0,
                 fat: float = 0.0,
                 blood: float = 0.0,
                 spo2: float = 0.97,
                 database: TissueDatabase = None):
        '''
        Model of a single skin layer.

        Parameters
        ----------
        d: float
            Layer thickness (m).
            Set to 1 mm by default.
        n: None or float or Callable[[float], float]
            Layer refractive index as a callable that takes one parameter,
            the wavelength of light (m), and returns the layer refractive index.
            Set to :py:func:`xopto.materials.ri.skin.default` by default.
        musr: float or Callable[[float], float]
            Reduced scattering coefficient as constant wavelength independent
            value or callable that takes the wavelength of light (m) and
            returns the reduced scattering coefficient (1/m).
            Instance of :py:class:`ReducedScatteringPower` can be used
            to define the reduced scattering coefficient wavelength dependence
            by a common power model.
        g: float or Callable[[float], float]
            A constant wavelength independent value of the scattering phase
            function anisotropy or a callable that takes the wavelength of
            light (m) and returns the scattering anisotropy g.
        baseline_absorption: float or Callable[[float], float]
            A constant wavelength independent value of the baseline absorption
            coefficient (1/m) or a callable that takes the wavelength of
            light (m) and returns the value of absorption coefficient (1/m).
            The volume fraction of baseline absorption is the remaining volume
            fraction that is not filled by the standard chromophores.
            Standard baseline absorption model for epidermis can be found in
            :py:attr:`xopto.materials.absorption.epidermis.default` and for
            dermis in :py:attr:`xopto.materials.absorption.dermis.default`.
        melanin: float
            Melanin volume fraction in the layer. Valid values are
            from 0.0 to 1.0.
        carotenoid: float
            Carotenoid volume fraction in the layer. Valid values are
            from 0.0 to 1.0.
        water: float
            Water volume fraction in the layer. Valid values are
            from 0.0 to 1.0.
        fat: float
            Fat volume fraction in the layer. Valid values are from 0.0 to 1.0.
        blood: float
            Blood volume fraction in the layer. Valid values are from
            0.0 to 1.0.
        spo2: float
            Oxygen saturation of the blood. Valid values are from 0.0 to 1.0.
        database: object
            Object with callables that take the wavelength of light (m)
            and return the absorption coefficient (1/m) of the chromophore.
            The object must have the following callable methods:

            - oxyhem_mua,
            - deoxyhem_mua,
            - melanin_mua,
            - carotenoid_mua,
            - water_mua,
            - fat_mua.

            Set to :py:class:`TissueDatabase` by default.
        '''
        if isinstance(d, SkinLayer):
            obj = d

            n = obj.n
            d = obj.d

            musr = obj.musr
            g = obj.g

            baseline_absorption = obj.baseline_absorption
            melanin = obj.melanin
            carotenoid = obj.carotenoid
            water = obj.water
            fat = obj.fat

            blood = obj.blood
            spo2 = obj.spo2

        if n is None:
            n = ri.skin

        if isinstance(n, (float, int)):
            n = ConstantValue(n)
        if isinstance(musr, (float, int)):
            musr = ConstantValue(musr)
        if isinstance(g, (float, int)):
            g = ConstantValue(g)
        if isinstance(baseline_absorption, (float, int)):
            baseline_absorption = ConstantValue(baseline_absorption)

        if database is None:
            database = TissueDatabase

        self._d = self._n = self._musr = self._g = None
        self._baseline_absorption = None
        self._melanin = self._carotenoid = self._water = self._fat = None
        self._blood = self._spo2 = None

        self._set_d(d)
        self._set_n(n)
        self._set_musr(musr)
        self._set_g(g)

        self._set_baseline_absorption(baseline_absorption)
        self._set_melanin(melanin)
        self._set_carotenoid(carotenoid)
        self._set_water(water)
        self._set_fat(fat)

        self._set_blood(blood)
        self._set_spo2(spo2)

        self._database = database

    def _set_d(self, d: float):
        self._d = max(float(d), 0.0)
    def _get_d(self) -> float:
        return self._d
    d = property(
        _get_d, _set_d, None,
        'Layer thickness (m).'
    )

    def _set_n(self, n: float or Callable[[float], float]):
        if isinstance(n, (float, int)):
            n = ConstantValue(n)
        self._n = n
    def _get_n(self) -> Callable[[float], float]:
        return self._n
    n = property(
        _get_n, _set_n, None,
        'Refractive index as a function of wavelength of light (m).'
    )

    def _get_musr(self) -> Callable[[float], float]:
        return self._musr
    def _set_musr(self, musr: float or Callable[[float], float]):
        if isinstance(musr, (float, int)):
            musr = ConstantValue(float(musr))
        self._musr = musr
    musr = property(_get_musr, _set_musr, None,
                     'Reduced scattering coefficient (1/m) as a function '
                     'of wavelength of light (m).')

    def _set_g(self, g: float or Callable[[float], float]):
        if isinstance(g, (float, int)):
            if g < -1.0 or g > 1.0:
                raise ValueError('Scattering anisotropy g out of valid '
                                 'range :math:`[-1, 1]`!')
            g = ConstantValue(float(g))

        self._g = g
    def _get_g(self) -> Callable[[float], float]:
        return self._g
    g = property(
        _get_g, _set_g, None,
        'Scattering anisotropy g as a function of wavelength of light (m).'
    )

    def _set_baseline_absorption(self, mua: float or Callable[[float], float]):
        if isinstance(mua, (int, float)):
            mua = ConstantValue(float(mua))
        self._baseline_absorption = mua
    def _get_baseline_absorption(self) -> Callable[[float], float]:
        return self._baseline_absorption
    baseline_absorption = property(
        _get_baseline_absorption, _set_baseline_absorption, None,
        'Baseline absorption coefficient (1/m) of the layer '
        'as a function of wavelength (m).')

    def _set_water(self, f: float):
        self._water = min(max(float(f), 0.0), 1.0)
    def _get_water(self) -> float:
        return self._water
    water = property(
        _get_water, _set_water, None,
        'Water volume fraction from 0.0 to 1.0.'
    )

    def _set_melanin(self, f: float):
        self._melanin = min(max(float(f), 0.0), 1.0)
    def _get_melanin(self) -> float:
        return self._melanin
    melanin = property(
        _get_melanin, _set_melanin, None,
        'Melanin volume fraction from 0.0 to 1.0.'
    )

    def _set_carotenoid(self, f: float):
        self._carotenoid = min(max(float(f), 0.0), 1.0)
    def _get_carotenoid(self) -> float:
        return self._carotenoid
    carotenoid = property(
        _get_carotenoid, _set_carotenoid, None,
        'Carotenoid volume fraction from 0.0 to 1.0.'
    )

    def _set_fat(self, f: float):
        self._fat = min(max(float(f), 0.0), 1.0)
    def _get_fat(self) -> float:
        return self._fat
    fat = property(
        _get_fat, _set_fat, None,
        'Fat volume fraction from 0.0 to 1.0.'
    )

    def _set_blood(self, f: float):
        self._blood = min(max(float(f), 0.0), 1.0)
    def _get_blood(self) -> float:
        return self._blood
    blood = property(
        _get_blood, _set_blood, None,
        'Blood volume fraction from 0.0 to 1.0.'
    )

    def _set_spo2(self, value: float):
        self._spo2 = min(max(float(value), 0.0), 1.0)
    def _get_spo2(self) -> float:
        return self._spo2
    spo2 = property(
        _get_spo2, _set_spo2, None,
        'Arterial blood oxygenation fraction from 0.0 to 1.0.')


    def mus(self, wavelength: float or np.ndarray) -> float or np.ndarray:
        '''
        Scattering coefficient (1/m) of the layer as a function of the
        wavelength of light (m). The scattering coefficient is computed
        by dividing the reduced scattering coefficient with 1 - g, where g is
        the wavelength dependent scattering anisotropy of the layer.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength of light (m)

        Returns
        -------
        musr: float or np.ndarray
            Scattering coefficient (1/m) of the layer at the specified
            wavelength of light (m).
        '''
        return self.musr(wavelength)/(1.0 - self.g(wavelength))

    def mua(self, wavelength: float or np.ndarray) -> float or np.ndarray:
        '''
        Absorption coefficient (1/m) of the layer as a function of the
        wavelength of light (m). The absorption coefficient can include
        absorption of the layer baseline, melanin absorption, water absorption,
        carotenoid absorption, fat absorption and blood absorption. The
        absorption of blood is oxygen saturation.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength of light (m) at which to compute the absorption
            coefficient (1/m).

        Returns
        -------
        mua: float or np.ndarray
            Absorption coefficient (1/m) at the specified
            wavelength (m) of light.
        '''
        mua_baseline = self.baseline_absorption(wavelength)

        mua_melanin = mua_carotenoid = mua_water = mua_fat = \
            mua_oxy = mua_deoxy = mua_blood = np.zeros_like(wavelength)
        if self._melanin > 0.0:
            mua_melanin = self._database.melanin_mua(wavelength)
        if self._carotenoid > 0.0:
            mua_carotenoid = self._database.carotenoid_mua(wavelength)
        if self._water > 0.0:
            mua_water = self._database.water_mua(wavelength)
        if self._fat > 0.0:
            mua_fat = self._database.fat_mua(wavelength)
        if self._blood > 0.0:
            mua_oxy = self._database.oxyhem_mua(wavelength)
            mua_deoxy = self._database.deoxyhem_mua(wavelength)
            mua_blood = self._spo2*mua_oxy + (1.0 - self._spo2)*mua_deoxy

        k_base = 1.0 - self._melanin - self._carotenoid - self._water - \
                 self._fat - self._blood

        mua = self._water*mua_water + self._melanin*mua_melanin + \
            self._carotenoid*mua_carotenoid + self._fat*mua_fat + \
            + self._blood*mua_blood + k_base*mua_baseline

        return mua

    def todict(self) -> dict:
        '''
        Export the skin layer model to a dict.

        Returns
        -------
        data: dict
            Model data exported to a dict.
        '''
        return {
            'type': self.__class__.__name__,
            'd': self.d,
            'n': self.n.value if isinstance(self.n, ConstantValue)
                              else self.n,
            'musr': self.musr.value if isinstance(self.musr, ConstantValue)
                                    else self.musr,
            'g': self.g.value if isinstance(self.g, ConstantValue)
                              else self.g,
            'baseline': self.baseline_absorption.value
                        if isinstance(self.baseline_absorption, ConstantValue)
                        else self.baseline_absorption,
            'melanin':self.melanin,
            'carotenoid':self.carotenoid,
            'water':self.water,
            'fat':self.fat,
            'blood':self.blood,
            'spo2':self.spo2,
        }

    def update(self, data: dict):
        '''
        Update the model parameters from a dict. Only the parameters that can
        be found in the supplied data are updated.

        data: dict
            Data for updating the model parameters.
        '''
        if data.get('d') is not None:
            self._set_d(data.get('d'))
        if data.get('n') is not None:
            self._set_n(data.get('n'))

        if data.get('musr') is not None:
            self._set_musr(data.get('musr'))
        if data.get('g') is not None:
            self._set_g(data.get('g'))

        if data.get('baseline_absorption') is not None:
            self._set_baseline_absorption(data.get('baseline_absorption'))
        if data.get('melanin') is not None:
            self._set_melanin(data.get('melanin'))
        if data.get('carotenoid') is not None:
            self._set_carotenoid(data.get('carotenoid'))
        if data.get('water') is not None:
            self._set_water(data.get('water'))
        if data.get('fat') is not None:
            self._set_fat(data.get('fat'))

        if data.get('blood') is not None:
            self._set_blood(data.get('blood'))
        if data.get('spo2') is not None:
            self._set_spo2(data.get('spo2'))

        return self

    def plot(self, wavelengths=None, *args, **kwargs):
        '''
        Plot the absorption and reduced scattering coefficients of the layer
        for the given wavelength range. The optional positional and keyword
        arguments are passed to the pyplot plot function.

        Parameters
        ----------
        wavelengths: np.ndarray
            Wavelengths at which to plot the absorption and reduced scattering
            coefficients of the layer.
        *args: list
            Optional positional arguments passed to the pyplot plot function.
        *kwargs: list
            Optional keyword arguments passed to the pyplot plot function.
        '''
        from matplotlib import pyplot as pp

        if wavelengths is None:
            wavelengths = np.linspace(400e-9, 1000e-9, 601)

        fig, ax = pp.subplots(1, 2)
        ax[0].plot(wavelengths*1e9, self.mua(wavelengths)*1e-2, *args, **kwargs)
        ax[0].grid()
        ax[0].set_ylabel('Absorption coefficient (1/cm)')
        ax[0].set_xlabel('Wavelength (nm)')

        ax[1].plot(wavelengths*1e9, self.musr(wavelengths)*1e-2, *args, **kwargs)
        ax[1].grid()
        ax[1].set_ylabel('Red. scatt. coefficient (1/cm)')
        ax[1].set_xlabel('Wavelength (nm)')

    def apply_to_mc_layer(self, mc_layer: mc.mclayer.Layer, wavelength: float):
        '''
        Apply the optical properties of skin computed at the specified
        wavelength to a Monte Carlo simulator layer instance.

        Parameters
        ----------
        layer: xopto.mcml.mclayer.layer.Layer
            Layer of a Monte Carlo simulator instance.
        wavelength: float
            Wavelength of light at which to compute the optical properties
            of the skin layer.
        '''
        if not isinstance(mc_layer.pf, mc.mcpf.Hg):
            TypeError('Monte Carlo simulator layers must use the HG '
                      'scattering phase function!')

        mc_layer.d = self.d
        mc_layer.n = self.n(wavelength)
        mc_layer.mua = self.mua(wavelength)
        mc_layer.mus = self.mus(wavelength)
        mc_layer.pf.g = self.g(wavelength)

    def __str__(self):
        return 'SkinLayer(d={:f}, n={}, '\
               'musr={}, g={}, '\
               'baseline_absorption={}, '\
               'melanin={:f}, carotenoid={:f}, water={:f}, fat={:f}, '\
               'blood={:f}, spo2={:f}'.format(
                   self.d, self.n, self.musr, self.g,
                   self.baseline_absorption,
                   self.melanin, self.carotenoid, self.water, self.fat,
                   self.blood, self.spo2)

    def __repr__(self):
        return '{:s} # id {}'.format(str(self), id(self))


class Skin3:
    def __init__(self, layers: None or Tuple[SkinLayer] = None):
        '''
        Initializes a 3-layer skin model that includes epidermis, dermis and 
        subcutaneous tissue layers.

        Parameters
        ----------
        layers: None or Tuple[SkinLayers]
            Stack of skin 3 skin layers (epidermis, dermis, subcutis).
            If None, a new default stack of 3 skin layers is created.
        '''
        if layers is None:
            layers = (
                SkinLayer(
                    d=0.1e-3,
                    musr=ReducedScatteringPower(15e2, 0.1, 580e-9),
                    n=ri.water.default,
                    g=0.9,
                    baseline_absorption=EpidermisBaselineAbsorption(),
                    melanin=0.05, water=0.02),
                SkinLayer(
                    d=2.0e-3,
                    musr=ReducedScatteringPower(20e2, 0.15, 580e-9),
                    n=ri.skin.default,
                    g=0.8,
                    baseline_absorption=DermisBaselineAbsorption(),
                    water=0.65, blood=0.04, spo2=0.97),
                SkinLayer(
                    d=np.inf,
                    musr=ReducedScatteringPower(10e2, 0.15, 580e-9),
                    n=ri.skin.default,
                    g=0.8,
                    baseline_absorption=DermisBaselineAbsorption(),
                    water=0.05, fat=0.4, blood=0.03, spo2=0.97,),
            )
        else:
            if len(layers) != 3:
                raise ValueError('The layer stack of skin3, must have '
                                 'exactly 3 SkinLayer instances!')
            if not all(isinstance(item, SkinLayer) for item in layer):
                raise TypeError('The layer stack of Skin3 must be formed by '
                                'SkinLayer instances!')

        self._layers = layers
        self._layer_names = ('epidermis', 'dermis', 'subcutis')

    def _get_epidermis(self) -> SkinLayer:
        return self._layer[0]
    epidermis = property(_get_epidermis, None, None, 'Epidermis skin layer')

    def _get_dermis(self) -> SkinLayer:
        return self._layer[1]
    dermis = property(_get_dermis, None, None, 'Dermis skin layer.')

    def _get_subcutis(self) -> SkinLayer:
        return self._layer[2]
    subcutis = property(_get_subcutis, None, None, 'Subcutis skin layer.')

    def _get_layers(self) -> Tuple[SkinLayer]:
        return self._layers
    layers = property(_get_layers, None, None,
                      'Skin layers as a tuple (epidermis, dermis, subcutis).')

    def update_mc_layers(self, layers: mc.mclayer.Layers, wavelength: float):
        '''
        Update the layer stack of the MC simulator instance. The simulator
        layer stack must contain exactly 5 layers, i.e. 3 inner layers for
        the skin and 2 outer layers for the surrounding medium.

        mc: xopto.mcml.mclayer.layer.Layers
            A Monte Carlo simulator layer stack.

        wavelength: float
            Wavelength of light that will be used to compute the optical
            properties of the skin layers.

        Note
        ----
        The surrounding medium layers (the first and the last layer in the
        stack) are not updated/changed by this call.
        '''
        if len(layers ) != 5:
            raise ValueError('The simulator instance must contain '
                             'exactly 5 layers!')
        
        for mc_layer, skin_layer in zip(layers[1:4], self._layers):
            skin_layer.apply_to_mc_layer(mc_layer, wavelength)

    def create_mc_layers(self, wavelength: float) -> mc.mclayer.Layers:
        '''
        Create and initialize a layer stack for Monte Carlo simulations.
        The refractive index of the surrounding medium is set to 1.0
        (air).

        Parameters
        ----------
        wavelength: float
            Wavelength of light that will be used to compute the optical
            properties of the skin layers.

        Returns
        -------
        layers: xopto.mcml.mclayer.layer.Layers
            A stack of monte carlo layers representing the skin model and
            two layers of the surrounding medium (the first and the last
            layer in the layer stack).
        '''
        layers = []
        for i in range(len(self._layers) + 2):
            layers.append(
                mc.mclayer.Layer(d=np.inf, mua=0.0, mus=0.0,
                                 n=1.0, pf=mc.mcpf.Hg(0.0))
            )
        self.update_mc_layers(layers, wavelength)

        return layers

    def todict(self) -> dict:
        '''
        Export skin model to a dict.

        Returns
        -------
        data: dict
            Skin model exported to dict.
        '''
        return {
            'type': self.__class__.__name__,
            'layers': [layer.todict() for layer in self._layers]
        }

    def _get_layers(self) -> List[SkinLayer]:
        return self._layers
    layers = property(_get_layers, None, None, 'The skin layer stack')

    def __getitem__(self, index) -> SkinLayer:
        return self._layers[index]

    def __len__(self) -> int:
        return len(self._layers)

    def __str__(self):
        return 'Skin3\n' + \
               '\n'.join([str(layer) for layer in self._layers])

    def __repr__(self):
        return str(self)


if __name__ == '__main__':
    from matplotlib import pyplot as pp

    layer = SkinLayer(
        d=0.8, n=ri.skin.default,
        musr=ReducedScatteringPower(20.0e2, 0.1, 550e-9), g=0.9,
        melanin=0.05, water=0.02,
        baseline_absorption=EpidermisBaselineAbsorption(),
    )
    layer.plot()
    pp.legend()
    pp.show()

    from xopto.util.sampler import UniformSampler, ConstantSampler
