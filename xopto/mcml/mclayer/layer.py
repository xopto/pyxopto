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

from typing import List
from xopto.mcbase.mcpf.pfbase import PfBase

import numpy as np

from xopto.mcml import mcobject
from xopto.mcml import mctypes
from xopto.mcml.mcutil import boundary
from xopto.mcml import cltypes
from xopto.mcml import mcpf


class Layer(mcobject.McObject):
    '''
    Class that represents a single sample layer.
    '''

    def cl_type(self, mc: mcobject.McObject) -> cltypes.Structure:
        '''
        Returns a structure data type that is used to represent one Layer
        instance in the OpenCL kernel of the Monte Carlo simulator.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.

        Returns
        -------
        opencl_t: ClLayer
            OpenCL Structure that represents a layer. 
        '''
        T = mc.types
        class ClLayer(cltypes.Structure):
            _fields_ = [
                ('thickness', T.mc_fp_t),
                ('top', T.mc_fp_t),
                ('bottom', T.mc_fp_t),
                ('n', T.mc_fp_t),
                ('cos_critical_top', T.mc_fp_t),
                ('cos_critical_bottom', T.mc_fp_t),
                ('mus', T.mc_fp_t),
                ('mua', T.mc_fp_t),
                ('inv_mut', T.mc_fp_t),
                ('mua_inv_mut', T.mc_fp_t),
                ('pf', self.pf.fetch_cl_type(mc))
            ]

        return ClLayer

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        '''
        Structure and related API that defines a layer in the Monte Carlo simulator.
        '''
        return '\n'.join((
            '/**',
            ' * @brief Data type describing a single sample layer.',
            ' * @note The members of this object are constant and do not change during the simulation.',
            ' * @{',
            ' */',
            'struct MC_STRUCT_ATTRIBUTES McLayer {',
            '	mc_fp_t thickness;                  /**< Layer thickness. */',
            '	mc_fp_t top;                        /**< Z coordinate of the layer top surface (z coordinate increases with the layer index). */',
            '	mc_fp_t bottom;                     /**< Z coordinate of the layer bottom surface (z coordinate increases with the layer index). */',
            '	mc_fp_t n;                          /**< Layer index of refraction. */',
            '	mc_fp_t cos_critical_top;           /**< Total internal reflection angle cosine for the above layer. */',
            '	mc_fp_t cos_critical_bottom;        /**< Total internal reflection angle cosine for the bellow layer. */',
            '	mc_fp_t mus;                        /**< Scattering coefficient. */',
            '	mc_fp_t mua;                        /**< Absorption coefficient. */',
            '	mc_fp_t inv_mut;                    /**< Reciprocal value of the total attenuation coefficient. */',
            '	mc_fp_t mua_inv_mut;                /**< Reciprocal value of the total attenuation coefficient multiplied by the absorption coefficient. */',
            '	McPf pf;                            /**< Scattering phase function parameters. */',
            '};',
            '/**',
            ' * @}',
            ' */',
            '',
            '/**',
            ' * @brief Evaluates to the layer thickness.',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define mc_layer_thickness(player) ((player)->thickness)',
            '',
            '/**',
            ' * @brief Evaluates to the z coordinate of the layer top surface.',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define mc_layer_top(player) ((player)->top)',
            '',
            '/**',
            ' * @brief Evaluates to the z coordinate of the layer bottom surface.',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define mc_layer_bottom(player) ((player)->bottom)',
            '',
            '/**',
            ' * @brief Evaluates to the refractive index of the layer.',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define mc_layer_n(player) ((player)->n)',
            '',
            '/**',
            ' * @brief Evaluates to the critical cosine (total internal reflection)',
            ' *        at the top layer boundary.',
            ' * @details If the absolute cosine of the angle of incidence',
            ' *          (with respect to z axis) is less than the critical cosine,',
            ' *          the incident packet is reflected at the boundary.',
            ' * @param[in] player Pointer to a layer object.',
            ' */',
            '#define mc_layer_cc_top(player) ((player)->cos_critical_top)',
            '',
            '/**',
            ' * @brief Evaluates to the critical cosine (total internal reflection)',
            ' *        at the bottom layer boundary.',
            ' * @details If the absolute cosine of the angle of incidence',
            ' *          (with respect to z axis) is less than the critical cosine, the',
            ' *          incident packet is reflected from the boundary.',
            ' * @param[in] player Pointer to a layer object.',
            ' */',
            '#define mc_layer_cc_bottom(player) ((player)->cos_critical_bottom)',
            '',
            '/**',
            ' * @brief Evaluates to the scattering coefficient along',
            ' *        the given propagation direction.',
            ' * @param[in] player Pointer to a layer instance.',
            ' * @param[in] pdir   Propagation direction vector.',
            ' */',
            '#define mc_layer_mus(player, pdir) ((player)->mus)',
            '',
            '/**',
            '* @brief Evaluates to the absorption coefficient along',
            ' *       the given propagation direction.',
            '* @param[in] player Pointer to a layer instance.',
            '* @param[in] pdir   Propagation direction vector.',
            '*/',
            '#define mc_layer_mua(player, pdir) ((player)->mua)',
            '',
            '/**',
            ' * @brief Evaluates to the reciprocal value of the total ',
            ' *        attenuation coefficient.',
            ' * @param[in] player Pointer to a layer instance.',
            ' * @param[in] pdir   Propagation direction vector.',
            ' */',
            '#define mc_layer_inv_mut(player, pdir) ((player)->inv_mut)',
            '',
            '/**',
            ' * @brief Evaluates to the reciprocal value of the total ',
            ' *        attenuation coefficient multiplied by the absorption coefficient.',
            ' * @param[in] player Pointer to a layer instance.',
            ' * @param[in] pdir   Propagation direction vector.',
            ' */',
            '#define mc_layer_mua_inv_mut(player, pdir) ((player)->mua_inv_mut)',
            '',
            '/**',
            ' * @brief Evaluates to the absorption coefficient of the layer multiplied',
            ' *        by the reciprocal of the total attenuation coefficient along',
            ' *        the given propagation direction.',
            ' * @param[in] player Pointer to a layer object.',
            ' * @param[in] pdir   Propagation direction vector.',
            ' *',
            ' * @returns   Absorption coefficient multiplied by the reciprocal',
            ' *            value of the total attenuation coefficient along the',
            ' *            given propagation direction.',
            ' */',
            '',
            '/**',
            ' * @brief Evaluates to a pointer to scattering phase function',
            ' *       of the layer.',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define mc_layer_pf(player) ((player)->pf)',
            '',
            '#if MC_ENABLE_DEBUG || defined(__DOXYGEN__)',
            '/**',
            ' * @brief Print one sample layer.',
            ' * param[in] prefix Can be used to pass indent string for the layer parameters."',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define dbg_print_layer(player, prefix) \\',
            '	printf(prefix "d: %.9f\\n" \\',
            '			prefix "top: %.9f\\n" \\',
            '			prefix "bottom: %.9f\\n" \\',
            '			prefix "n: %.9f\\n" \\',
            '			prefix "cctop: %.9f\\n" \\',
            '			prefix "ccbottom: %.9f\\n" \\',
            '			prefix "mua: %.9f\\n" \\',
            '			prefix "mus: %.9f\\n" \\',
            '			prefix "inv_mut: %.9f\\n" \\',
            '			prefix "mua_inv_mut: %.9f\\n", \\',
            '					(player)->thickness, (player)->top, (player)->bottom, (player)->n, \\',
            '					(player)->cos_critical_top, (player)->cos_critical_bottom, \\',
            '					(player)->mua, (player)->mus, \\',
            '					(player)->inv_mut, (player)->mua_inv_mut); \\',
            '			{ McPf const _dbg_pf=(player)->pf; dbg_print_pf(&_dbg_pf); };',
            '#else',
            '#define dbg_print_layer(player, label) ;',
            '#endif',
            ))

    def __init__(self, d: float, n: float, mua: float, mus: float,
                 pf: mcpf.PfBase):
        '''
        Layer object constructor.

        Parameters
        ----------
        d: float
            Layer thickness (m).
        n: float
            Index of refraction.
        mua: float
            Absorption coefficient (1/m).
        mus: float
            Scattering (NOT reduced) coefficient (1/m).
        pf: mcpf.PfBase
            Scattering phase function object that is derived from the
            :py:class:`xopto.mcbase.mcpf.pfbase.PfBase` class.

        The physical properties of the layer can be read or changed through
        member properties:

            - d: float - 
              Layer thickness (m).
            - n: float - 
              Index of refraction.
            - mua: float - 
              Absorption coefficient (1/m).
            - mus: float - 
              Scattering (NOT reduced) coefficient (1/m).
            - pf: mcpf.PfBase - 
              Scattering phase function object that is derived from the
              :py:class:`xopto.mcbase.mcpf.pfbase.PfBase` class.

        Note
        ----
        The layer boundaries in the z direction are increasing in the
        direction of the layer stack.
        Z coordinate 0.0 belongs to the top surface of the first sample layer!
        '''
        self._d = float(d)
        self._n = float(n)
        self._mua = float(mua)
        self._mus = float(mus)
        self._pf = pf

    def _set_d(self, d: float):
        self._d = float(d)
    def _get_d(self) -> float:
        return self._d
    d = property(_get_d, _set_d, None, 'Layer thickens (m).')

    def _set_n(self, n: float):
        self._n = float(n)
    def _get_n(self) -> float:
        return self._n
    n = property(_get_n, _set_n, None, 'Refractive index of the layer.')

    def _set_mua(self, mua: float):
        self._mua = float(mua)
    def _get_mua(self) -> float:
        return self._mua
    mua = property(_get_mua, _set_mua, None,
                   'Absorption coefficient of the layer (1/m).')

    def _set_mus(self, mus: float):
        self._mus = float(mus)
    def _get_mus(self) -> float:
        return self._mus
    mus = property(_get_mus, _set_mus, None,
                   'Scattering coefficient of the layer (1/m).')

    def _get_pf(self) -> mcpf.PfBase:
        return self._pf
    def _set_pf(self, pf: mcpf.PfBase):
        if type(self._pf) != type(pf):
            raise ValueError('The scattering phase function type '\
                             'of the layer must not change!')
        self._pf = pf
    pf = property(_get_pf, _set_pf, None, 'Phase function object.')

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Pack the layers into an OpenCL data type. The OpenCL data
        type is returned by the :py:meth:`Layers.cl_type` method.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: cltypes.Structure
            A structure representing a layer in the MC simulator.

        Returns
        -------
        target: cltypes.Structure
            Structure received as an input argument or a new
            instance if the input argument target is None.

        Note
        ----
        This method only packs the fields that do not depend on the
        properties of the other/neighboring layers in the stack
        (PACKS thickness, n, mua, mus, inv_mut, and mua_inv_mut, but
        DOES NOT PACK top, bottom, cos_critical_top, and cos_critical_bottom).
        '''
        if target is None:
            target_type = self.fetch_cl_type(mc)
            target = target_type()

        mut = self.mua + self.mus
        if mut > 0.0:
            inv_mut = 1.0/mut
        else:
            inv_mut = float('inf')

        if self.mus == 0.0:
            mua_inv_mut = 1.0
        else:
            mua_inv_mut = self.mua*inv_mut

        target.thickness = self.d
        target.n = self.n
        target.mua = self.mua
        target.mus = self.mus
        target.inv_mut = inv_mut
        target.mua_inv_mut = mua_inv_mut

        self.pf.cl_pack(mc, target.pf)

        return target

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'d':self._d, 'n':self._n, 'mua':self._mua, 'mus':self._mus,
                'pf':self._pf.todict(), 'type':'Layer'}

    @classmethod
    def fromdict(cls, data: dict) -> 'Layer':
        '''
        Create a new object from dict. The dict keys must match
        the parameter names defined by the constructor.
        '''
        data_ = dict(data)
        t = data_.pop('type')
        if t != 'Layer':
            raise ValueError('Cannot create a Layer instance from the data!')
        pf_data = data_.pop('pf')
        if not hasattr(mcpf, pf_data['type']):
            raise TypeError('Scattering phase function "{}" '
                            'not implemented'.format(pf_data['type']))
        pf_type = getattr(mcpf, pf_data['type'])
        if pf_type is None:
            raise TypeError('Scattering phase function type "{}" not '
                            'found!'.format(t))
        return cls(pf=pf_type.fromdict(pf_data), **data_)

    def __str__(self):
        return 'Layer(d={}, n={}, mua={}, mus={}, pf={})'.format(
            self._d, self._n, self._mua, self._mus, self._pf)

    def __repr__(self):
        return  '{:s} # id 0x{:>08X}.'.format(self.__str__(), id(self))


class AnisotropicLayer(mcobject.McObject):
    '''
    Class that represents a single anisotropic sample layer.
    '''

    def cl_type(self, mc: mcobject.McObject) -> cltypes.Structure:
        '''
        Returns a structure data type that is used to represent one Layer
        instance in the OpenCL kernel of the Monte Carlo simulator.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.

        Returns
        -------
        opencl_t: ClAnisotropicLayer
            OpenCL Structure that represents a layer. 
        '''
        T = mc.types
        class ClAnisotropicLayer(cltypes.Structure):
            _fields_ = [
                ('thickness', T.mc_fp_t),
                ('top', T.mc_fp_t),
                ('bottom', T.mc_fp_t),
                ('n', T.mc_fp_t),
                ('cos_critical_top', T.mc_fp_t),
                ('cos_critical_bottom', T.mc_fp_t),
                ('mus', T.mc_matrix3f_t),
                ('mua', T.mc_matrix3f_t),
                ('mut', T.mc_matrix3f_t),
                ('pf', self.pf.fetch_cl_type(mc))
            ]

        return ClAnisotropicLayer

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        '''
        Structure that defines the layer in the Monte Carlo simulator.
        '''
        return '\n'.join((
            '/**',
            ' * @brief Data type describing a single sample layer.',
            ' * @note The members of this object are constant and do not change during the simulation.',
            ' * @{',
            ' */',
            'struct MC_STRUCT_ATTRIBUTES McLayer {',
            '	mc_fp_t thickness;                  /**< Layer thickness. */',
            '	mc_fp_t top;                        /**< Z coordinate of the layer top surface (z coordinate increases with the layer index). */',
            '	mc_fp_t bottom;                     /**< Z coordinate of the layer bottom surface (z coordinate increases with the layer index). */',
            '	mc_fp_t n;                          /**< Layer index of refraction. */',
            '	mc_fp_t cos_critical_top;           /**< Total internal reflection angle cosine for the above layer. */',
            '	mc_fp_t cos_critical_bottom;        /**< Total internal reflection angle cosine for the bellow layer. */',
            '	mc_matrix3f_t mus_tensor;           /**< Scattering coefficient tensor. */',
            '	mc_matrix3f_t mua_tensor;           /**< Absorption coefficient tensor. */',
            '	mc_matrix3f_t mut_tensor;           /**< Total attenuation coefficient tensor. */',
            '	McPf pf;                            /**< Scattering phase function parameters. */',
            '};',
            '/**',
            ' * @}',
            ' */',
            '',
            '/**',
            ' * @brief Evaluates to the layer thickness.',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define mc_layer_thickness(player) ((player)->thickness)',
            '',
            '/**',
            ' * @brief Evaluates to the z coordinate of the layer top surface.',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define mc_layer_top(player) ((player)->top)',
            '',
            '/**',
            ' * @brief Evaluates to the z coordinate of the layer bottom surface.',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define mc_layer_bottom(player) ((player)->bottom)',
            '',
            '/**',
            ' * @brief Evaluates to the refractive index of the layer.',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define mc_layer_n(player) ((player)->n)',
            '',
            '/**',
            ' * @brief Evaluates to the critical cosine (total internal reflection)',
            ' *        at the top layer boundary.',
            ' * @details If the absolute cosine of the angle of incidence',
            ' *          (with respect to z axis) is less than the critical cosine,',
            ' *          the incident packet is reflected at the boundary.',
            ' * @param[in] player Pointer to a layer object.',
            ' */',
            '#define mc_layer_cc_top(player) ((player)->cos_critical_top)',
            '',
            '/**',
            ' * @brief Evaluates to the critical cosine (total internal reflection)',
            ' *        at the bottom layer boundary.',
            ' * @details If the absolute cosine of the angle of incidence',
            ' *          (with respect to z axis) is less than the critical cosine, the',
            ' *          incident packet is reflected from the boundary.',
            ' * @param[in] player Pointer to a layer object.',
            ' */',
            '#define mc_layer_cc_bottom(player) ((player)->cos_critical_bottom)',
            '',
            '/**',
            ' * @brief Evaluates to a pointer to the scattering coefficient tensor.',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define mc_layer_mus_tensor(player) (&(player)->mus_tensor)',
            '',
            '/**',
            ' * @brief Evaluates to a pointer to the absorption coefficient tensor.',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define mc_layer_mua_tensor(player) (&(player)->mua_tensor)',
            '',
            '/**',
            ' * @brief Evaluates to a pointer to the total attenuation coefficient tensor.',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define mc_layer_mut_tensor(player) (&(player)->mut_tensor)',
            '',
            '/**',
            ' * @brief Evaluates to the scattering coefficient along',
            ' *        the given propagation direction.',
            ' * @param[in] player Pointer to a layer instance.',
            ' * @param[in] pdir   Propagation direction vector.',
            ' */',
            '#define mc_layer_mus(player, pdir) tensor3f_project(mc_layer_mus_tensor(player), pdir)',
            '',
            '/**',
            '* @brief Evaluates to the absorption coefficient along',
            ' *       the given propagation direction.',
            '* @param[in] player Pointer to a layer instance.',
            '* @param[in] pdir   Propagation direction vector.',
            '*/',
            '#define mc_layer_mua(player, pdir) tensor3f_project(mc_layer_mua_tensor(player), pdir)',
            '',
            '/**',
            ' * @brief Evaluates to the total attenuation coefficient along',
            ' *        the given propagation direction.',
            ' * @param[in] player Pointer to a layer instance.',
            ' * @param[in] pdir   Propagation direction vector.',
            ' */',
            '#define mc_layer_mut(player, pdir) tensor3f_project(mc_layer_mut_tensor(player), pdir)',
            '',
            '/**',
            ' * @brief Evaluates to the absorption coefficient of the layer multiplied',
            ' *        by the reciprocal of the total attenuation coefficient along',
            ' *        the given propagation direction.',
            ' * @param[in] player Pointer to a layer object.',
            ' * @param[in] pdir   Propagation direction vector.',
            ' *',
            ' * @returns   Absorption coefficient multiplied by the reciprocal',
            ' *            value of the total attenuation coefficient along the',
            ' *            given propagation direction.',
            ' */',
            'static inline mc_fp_t mc_layer_mua_inv_mut(__constant McLayer const *player, mc_point3f_t const *pdir) {',
            '	mc_fp_t mua = mc_layer_mua(player, pdir);',
            '	mc_fp_t mut = mc_layer_mut(player, pdir);',
            '',
            '	return (mua != FP_0) ? ((mut != FP_0) ? mc_fdiv(mua, mut) : FP_INF) : FP_0;',
            '};',
            '',
            '/**',
            ' * @brief Evaluates to the reciprocal value of the total attenuation',
            ' *        coefficient along the given propagation direction.',
            ' * @param[in] player Pointer to a layer object.',
            ' * @param[in] pdir   Propagation direction vector.',
            ' *',
            ' * @returns   Reciprocal value of the total attenuation coefficient',
            ' *            along the given propagation direction.',
            ' */',
            'static inline mc_fp_t mc_layer_inv_mut(__constant McLayer const *player, mc_point3f_t const *pdir) {',
            '	mc_fp_t mut = mc_layer_mut(player, pdir);',
            '',
            '	return (mut != FP_0) ? mc_reciprocal(mut) : FP_INF;',
            '};',
            '',
            '/**',
            ' * @brief Evaluates to a pointer to scattering phase function',
            ' *       of the layer.',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define mc_layer_pf(player) ((player)->pf)',
            '',
            '#if MC_ENABLE_DEBUG || defined(__DOXYGEN__)',
            '/**',
            ' * @brief Print one sample layer.',
            ' * param[in] prefix Can be used to pass indent string for the layer parameters."',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define dbg_print_layer(player, prefix) \\',
            '	printf(prefix "d: %.9f\\n" \\',
            '			prefix "top: %.9f\\n" \\',
            '			prefix "bottom: %.9f\\n" \\',
            '			prefix "n: %.9f\\n" \\',
            '			prefix "cctop: %.9f\\n" \\',
            '			prefix "ccbottom: %.9f\\n", \\',
            '					(player)->thickness, (player)->top, (player)->bottom, (player)->n, \\',
            '					(player)->cos_critical_top, (player)->cos_critical_bottom); \\',
            '	dbg_print_matrix3f(prefix "mua:", &(player)->mua_tensor); \\',
            '	dbg_print_matrix3f(prefix "mus:", &(player)->mus_tensor); \\',
            '	dbg_print_matrix3f(prefix "mut:", &(player)->mut_tensor); \\',
            '	{ McPf const _dbg_pf=(player)->pf; dbg_print_pf(&_dbg_pf); };',
            '#else',
            '#define dbg_print_layer(player, label) ;',
            '#endif',
        ))

    def __init__(self, d: float, n: float,
                 mua: float or np.ndarray, mus: float or np.ndarray,
                 pf: mcpf.PfBase):
        '''
        Anisotropic layer object constructor.

        Parameters
        ----------
        d: float
            Layer thickness (m).
        n: float
            Index of refraction.
        mua: float or np.ndarray
            Absorption coefficient tensor (1/m). A scalar float value for an
            isotropic material. A vector of length 3 for the diagonal elements
            of the tensor (non-diagonal elements are set to 0).
            A numpy array of shape (3, 3) for the complete tensor.
        mus: float
            Scattering coefficient tensor (1/m). A scalar float value for
            an isotropic material. A vector of length 3 if only the diagonal
            elements of the tensor are nonzero (non-diagonal elements are set to 0).
            A numpy array of shape (3, 3) for the complete tensor.
        pf: mcpf.PfBase
            Scattering phase function object that is derived from
            :py:class:`xopto.mcbase.mcpf.pfbase.PfBase` class.

        The physical properties of the layer can be read or changed through
        properties:

            - d: float - 
              Layer thickness (m).
            - n: float - 
              Index of refraction.
            - mua: float or np.ndarray - 
              Absorption coefficient tensor (1/m).
            - mus: float np.ndarray - 
              Scattering (NOT reduced) coefficient tensor (1/m).
            - pf: mcpf.PfBase - 
              Scattering phase function object that is derived from
              :py:class:`xopto.mcbase.mcpf.pfbase.PfBase` class.

        Note
        ----
        The layer boundaries in the z direction are increasing in the
        direction of the layer stack.
        Z coordinate 0.0 belongs to the top surface of the first sample layer!
        '''
        self._d = float(d)
        self._n = float(n)
        self._mua = np.zeros((3, 3))
        self._mus = np.zeros((3, 3))
        self._pf = pf

        self._set_mua(mua)
        self._set_mus(mus)

    def _set_d(self, d: float):
        self._d = float(d)
    def _get_d(self) -> float:
        return self._d
    d = property(_get_d, _set_d, None, 'Layer thickness (m).')

    def _set_n(self, n: float):
        self._n = float(n)
    def _get_n(self) -> float:
        return self._n
    n = property(_get_n, _set_n, None, 'Refractive index of the layer.')

    def _set_mua(self, mua: float or np.ndarray):
        if isinstance(mua, (float, int)):
            self._mua[0, 0] = mua
            self._mua[1, 1] = mua
            self._mua[2, 2] = mua
        else:
            mua = np.asarray(mua, dtype=float)
            if mua.size == 3:
                self._mua[0, 0] = mua[0]
                self._mua[1, 1] = mua[1]
                self._mua[2, 2] = mua[2]
            else:
                self._mua[:] = mua
 
    def _get_mua(self) -> np.ndarray:
        return self._mua
    mua = property(_get_mua, _set_mua, None,
                   'Absorption coefficient tensor (3x3) of the layer (1/m).')

    def _set_mus(self, mus: float or np.ndarray):
        if isinstance(mus, (float, int)):
            self._mus[0, 0] = mus
            self._mus[1, 1] = mus
            self._mus[2, 2] = mus
        else:
            mus = np.asarray(mus, dtype=float)
            if mus.size == 3:
                self._mus[0, 0] = mus[0]
                self._mus[1, 1] = mus[1]
                self._mus[2, 2] = mus[2]
            else:
                self._mus[:] = mus
    def _get_mus(self) -> np.ndarray:
        return self._mus
    mus = property(_get_mus, _set_mus, None,
                   'Scattering coefficient tensor (3x3) of the layer (1/m).')

    def _get_pf(self) -> mcpf.PfBase:
        return self._pf
    def _set_pf(self, pf: mcpf.PfBase):
        if type(self._pf) != type(pf):
            raise ValueError('The scattering phase function type '\
                             'of the layer must not change!')
        self._pf = pf
    pf = property(_get_pf, _set_pf, None,
                  'Scattering phase function object.')

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Pack the layers into an OpenCL data type. The OpenCL data
        type is returned by the :py:meth:`Layers.cl_type` method.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: cltypes.Structure
            A structure representing a layer in the MC simulator.

        Returns
        -------
        target: cltypes.Structure
            Structure received as an input argument or a new
            instance if the input argument target is None.

        Note
        ----
        This method only packs the fields that do not depend on the
        properties of the other/neighboring layers in the stack
        (PACKS thickness, n, mua, mus, and mut, but
        DOES NOT PACK top, bottom, cos_critical_top, and cos_critical_bottom).
        '''
        if target is None:
            target_type = self.fetch_cl_type(mc)
            target = target_type()

        target.thickness = self.d
        target.n = self.n
        target.mua.fromarray(self.mua)
        target.mus.fromarray(self.mus)
        target.mut.fromarray(self.mua + self.mus)

        self.pf.cl_pack(mc, target.pf)

        return target

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'d':self._d, 'n':self._n, 'mua':self._mua, 'mus':self._mus,
                'pf':self._pf.todict(), 'type': 'AnisotropicLayer'}

    @classmethod
    def fromdict(cls, data: dict) -> 'Layer':
        '''
        Create a new object from dict. The dict keys must match
        the parameter names defined by the constructor.
        '''
        data_ = dict(data)
        t = data_.pop('type')
        if t != 'AnisotropicLayer':
            raise ValueError(
                'Cannot create an AnisotropicLayer instance from the data!')
        pf_data = data_.pop('pf')
        if not hasattr(mcpf, pf_data['type']):
            raise TypeError('Scattering phase function "{}" '
                            'not implemented'.format(pf_data['type']))
        pf_type = getattr(mcpf, pf_data['type'])
        return cls(pf=pf_type.fromdict(pf_data), **data_)

    def __str__(self):
        return 'AnisotropicLayer(d={}, n={}, mua={}, mus={}, pf={})'.format(
            self._d, self._n, self._mua, self._mus, self._pf)

    def __repr__(self):
        return  '{:s} # id 0x{:>08X}.'.format(self.__str__(), id(self))


class Layers(mcobject.McObject):
    '''
    Class that represents a stack of layers forming the sample.

    Note
    ----
    The topmost and bottommost layers of the stack are used to describe the
    medium that surrounds the sample top and bottom surfaces, respectively.
    Therefore, at least three layers must be always specified,
    namely two layers of the surrounding medium and one sample layer!
    The thicknesses of the topmost and bottommost layers will be automatically
    set to infinity regardless of the layer thickness set by the user.
    '''
    def __init__(self, layers: List[Layer or AnisotropicLayer] or 'Layers'):
        '''
        Constructs a managed sample layer stack from a list of sample layers.

        Note
        ----
        The topmost and bottommost layers of the stack are used to describe the
        medium that surrounds the sample top and bottom surfaces, respectively.
        Therefore, at least three layers must be always specified,
        namely two layers of the surrounding medium and one sample layer!

        The bottom surface of the topmost layer (the surrounding medium) is
        located at z=0. The positive direction of the z axis points in the
        direction of the layer stack.

        The thicknesses of the topmost and bottommost layers will be
        automatically set to infinity when passed to the OpenCL kernel
        (regardless of the layer thickness set by the user).

        Note that all layers must use the same scattering phase function
        model.

        Parameters
        ----------
        layers: List[Layer or AnisotropicLayer] or Layers
            A list of sample layers. Requires at least 3 items!
        '''
        self._layer_type = self._pf_type = None

        if isinstance(layers, Layers):
            self._layers = layers.tolist()
            self._pf_type = type(self._layers[0].pf)
        else:
            self._layers = list(layers)
            self.check()

    def check(self):
        '''
        Check if the layers are consistent and using a single
        scattering phase function type.
        Raises exception on error.
        '''
        if len(self._layers) < 3:
                ValueError('At least three layers are required, '
                           'but got only {:d}!'.format(len(self._layers)))

        if self._pf_type is None:
            self._pf_type = type(self._layers[1].pf)
        if self._layer_type is None:
            self._layer_type = type(self._layers[0])

        for layer in self._layers:
            if not isinstance(layer, (Layer, AnisotropicLayer)):
                raise TypeError(
                    'All the sample layers must be instances of Layer or '
                    'AnisotropicLayer but found {:s}!'.format(
                        type(layer).__name__))
            if self._layer_type != type(layer):
                raise TypeError(
                    'All the sample layers must use the same type!'
                    'Found {} and {}!'.format(
                        self._layer_type.__name__, layer.__class__.__name__))

            if type(layer.pf) != self._pf_type:
                raise TypeError(
                    'All the sample layers must use the same scattering phase '
                    'function model! Found {} and {}!'.format(
                        self._pf_type.__name__, type(layer.pf).__name__))


    def layer(self, index: int) -> Layer or AnisotropicLayer:
        '''
        Returns layer at the specified index. Note that the first layer
        (index 0) and the last layer (index -1) represent the medium
        surrounding the top and bottom surfaces of the sample, respectively
        '''
        return self._layers[index]

    def layer_index(self, z: float) -> int:
        '''
        Returns the layer index that contains the given z coordinate. Note that
        the layer includes the top surface boundary but not the bottom surface
        boundary, i.e. the z extent of the layer is [z_top, z_bottom), where
        z_bottom > z_top.

        Parameters
        ----------
        z: float
            Z coordinate of a point.

        Returns
        -------
        layer_index: int
            Index of the layer that contains the give point z.
        '''
        index = len(self._layers) - 1
        bottom = 0.0
        for pos, layer in enumerate(self._layers):
            if pos > 0:
                bottom += layer.d
            if z < bottom:
                index = pos
                break
        return index

    def thickness(self) -> float:
        '''
        Thickness of the layer stack excluding the topmost and bottommost
        layers that surround the sample.

        Returns
        -------
        thickness: float
            The sample thickness excluding the topmost and bottommost layers of
            the surrounding medium.
        '''
        d = 0.0
        for layer in self._layers[1:-1]:
            d += layer.d
        return d

    def cl_type(self, mc: mcobject.McObject) -> cltypes.Array:
        '''
        Returns an OpenCL array of ClLayer structures that is used to
        represent one instance of a layer stack.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.

        Returns
        ------- 
        clarray: cltypes.Structure*len(self)
            Array of ClLayers.
        '''
        return self._layers[0].fetch_cl_type(mc)*len(self._layers)

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Array = None) \
            -> cltypes.Array:
        '''
        Pack the layers into an OpenCL data type. The OpenCL data
        type is returned by the :py:meth:`Layers.cl_type` method.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: cltypes.Structure*len(self)
            A structure representing an array of Layers.

        Returns
        -------
        target: cltypes.Structure
            Filled structure received as an input argument or a new
            instance if the input argument target is None.

        Note
        ----
        This method only packs the fields that depend on the
        properties of the other/neighboring layers in the stack
        (PACKS top, bottom, cos_critical_top, cos_critical_bottom, but
        USES the LAYER instance to PACK all the remaining fields).
        '''

        num_layers = len(self._layers)

        if target is None or len(target) != num_layers:
            target_type = self.fetch_cl_type(mc)
            target = target_type()

        for index, layer in enumerate(self._layers):
            # pack the properties that do not depend on neighboring layers
            layer.cl_pack(mc, target[index])

            cc_top = cc_bottom = 0.0
            if index > 0:
                cc_top = boundary.cos_critical(
                    layer.n, self._layers[index - 1].n)
            if index + 1 < num_layers:
                cc_bottom = boundary.cos_critical(
                    layer.n, self._layers[index + 1].n)

            # target[index].thickness = layer.d
            if index == 0:
                target[index].top = -float('inf')
                target[index].bottom = 0.0
                target[index].thickness = float('inf')
            elif index == num_layers - 1:
                target[index].top = target[index - 1].bottom.value
                target[index].bottom = float('inf')
                target[index].thickness = float('inf')
            else:
                target[index].top = target[index - 1].bottom.value
                target[index].bottom = target[index - 1].bottom.value + layer.d

            target[index].cos_critical_top = cc_top
            target[index].cos_critical_bottom = cc_bottom

            layer.pf.cl_pack(mc, target[index].pf)

        return target

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'layers': [layer.todict() for layer in self._layers],
                'type': 'Layers'}

    def tolist(self) -> List[Layer or AnisotropicLayer]:
        '''
        Returns a weak copy of the list of managed layers.

        Returns
        -------
        layers: list[layers]
            List of managed layers
        '''
        return list(self._layers)

    @classmethod
    def fromdict(cls, data: dict) -> 'Layers':
        '''
        Create a new Layers object from a dict. The dict keys must match
        the parameter names defined by the constructor.
        '''
        data_ = dict(data)
        T = data_.pop('type')
        if T != 'Layers':
            raise ValueError(
                'Cannot create a Layers instance from the given data!')

        layers = []
        for item in data_.pop('layers'):
            T = {'Layer': Layer,
                 'AnisotropicLayer': AnisotropicLayer}.get(item.get('type'))

            layers.append(T.fromdict(item))

        return cls(layers=layers, **data_)

    def __getitem__(self, what):
        return self._layers[what]

    def __setitem__(self, what, value):
        self._layers[what] = value

    def __len__(self):
        return len(self._layers)

    def __str__(self):
        layers = ['    ' + str(layer) for layer in self._layers]
        layers[0] +=  ",  # medium above the sample"
        layers[-1] += "   # medium bellow the sample"
        layers_str = ',\n'.join(layers)
        return 'Layers([\n{}\n])'.format(layers_str)

    def __repr__(self):
        return  '{:s} # id 0x{:>08X}.'.format(self.__str__(), id(self))
