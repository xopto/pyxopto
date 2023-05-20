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

from typing import List, Tuple
from xopto.mcbase.mcpf.pfbase import PfBase

from xopto.mcml import mcobject
from xopto.mcml import mctypes
from xopto.mcml.mcutil import boundary
from xopto.mcml import cltypes
from xopto.mcml import mcpf


def ray_cylinder_intersection(r: float,
                              pos: Tuple[float, float, float],
                              dir: Tuple[float, float, float]) \
                                  -> Tuple[float, float] or None:
    '''
    Computes intersection between a cylinder :math:`x^2 + y^2 = r^2` that is
    centered on the :math:`z` axis and a ray
    with origin at :math:`pos=(x_0, y_0, z_0)` and direction
    :math:`dir=(p_x, p_y, p_z)`. The ray propagation satisfies
    a parametric equation
    :math:`(x, y, z) = (x_0, y_0, z_0) + t (p_x, p_y, p_z)`.
    The intersection is governed by a  quadratic equation for distance to
    intersection :math:`t`:

    .. math::

        t^2 (p_x^2 + p_y^2) + 2 t (x_0 p_x + y_0 p_y) + x_0^2 + y_0^2 - r^2.

    Solution is found as :math:`t = \\frac{-b \\pm\\sqrt{b^2 - 4 a c}}{2 a}`,
    where:

    .. math::

        a = p_x^2 + p_y^2,

        b = 2 (x_0 p_x + y_0 p_y),

        c = x_0^2 + y_0^2 - r^2.

    If :math:`a=0` or :math:`b^2 - 4 a c < 0` there is no intersection.


    Parameters
    ----------
    r: float
        Cylinder radius. The cylinder is aligned along the z axis.
    pos: Tuple[float, float, float]
        Ray origin as a tuple (x, y, z).
    dir: Tuple[float, float, float]
        Ray direction as a tuple (x, y, z). Note that the length of the
        direction vector must equal 1.

    Returns
    -------
    intersections: Tuple(float, float) or Tuple(None, None)
        Returns None if there is no intersection, else the two distances to the
        intersection as a tuple (d1, d2).
    '''
    a = dir[0]**2 + dir[1]**2
    b = 2.0*(pos[0]*dir[0] + pos[1]*dir[1])
    c = pos[0]**2 + pos[1]**2 - r**2
    D = b**2 - 4*a*c
    if D < 0.0 or a == 0.0:
        return None, None

    D = D**0.5
    d1 = (-b - D)/(2.0*a)
    d2 = (-b + D)/(2.0*a)

    return d1, d2

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
                ('r_inner', T.mc_fp_t),
                ('r_outer', T.mc_fp_t),
                ('n', T.mc_fp_t),
                ('cos_critical_inner', T.mc_fp_t),
                ('cos_critical_outer', T.mc_fp_t),
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
            '	mc_fp_t r_inner;					/**< Radius of the inner boundary. */',
            '	mc_fp_t r_outer;					/**< Radius of the outer boundary. */',
            '	mc_fp_t n;							/**< Layer index of refraction. */',
            '	mc_fp_t cos_critical_inner;			/**< Total internal reflection angle cosine for the inner boundary. */',
            '	mc_fp_t cos_critical_outer;			/**< Total internal reflection angle cosine for the outer boundary. */',
            '	mc_fp_t mus;						/**< Scattering coefficient. */',
            '	mc_fp_t mua;						/**< Absorption coefficient. */',
            '	mc_fp_t inv_mut;					/**< Reciprocal of the total attenuation coefficient. */',
            '	mc_fp_t mua_inv_mut;				/**< Absorption coefficient multiplied by the reciprocal of the total attenuation coefficient. */',
            '	McPf pf;							/**< Scattering phase function parameters. */',
            '};',
            '/**',
            ' * @}',
            ' */',
            '/** @brief Data type representing a sample layer. */',
            'typedef struct McLayer McLayer;',
            '',
            '/**',
            ' * @brief Evaluates to the radius of the inner boundary.',
            ' * @param[in] player Pointer to a layer object.',
            ' */',
            '#define mc_layer_r_inner(player) ((player)->r_inner)',
            '',
            '/**',
            ' * @brief Evaluates to the radius of the outer boundary.',
            ' * @param[in] player Pointer to a layer object.',
            ' */',
            '#define mc_layer_r_outer(player) ((player)->r_outer)',
            '',
            '/**',
            ' * @brief Evaluates to the z coordinate of the top layer surface.',
            ' * @param[in] player Pointer to a layer object.',
            ' */',
            '#define mc_layer_top(player) ((player)->top)',
            '',
            '/**',
            ' * @brief Evaluates to the z coordinate of the bottom layer surface.',
            ' * @param[in] player Pointer to a layer object.',
            ' */',
            '#define mc_layer_bottom(player) ((player)->bottom)',
            '',
            '/**',
            ' * @brief Evaluates to the layer refractive index.',
            ' * @param[in] player Pointer to a layer object.',
            ' */',
            ' #define mc_layer_n(player) ((player)->n)',
            '',
            '/**',
            ' * @brief Evaluates to the critical cosine (total internal reflection)',
            ' *        at the inner layer boundary.',
            ' * @details If the absolute cosine of the angle of incidence',
            ' *          (with respect to z axis) is less than the critical cosine,',
            ' @          the incident packet is reflected at the boundary.',
            ' * @param[in] player Pointer to a layer object.',
            ' */',
            ' #define mc_layer_cc_inner(player) ((player)->cos_critical_inner)',
            '',
            '/**',
            ' * @brief Evaluates to the critical cosine (total internal reflection)',
            ' *        at the outer layer boundary.',
            ' * @details If the absolute cosine of the angle of incidence',
            ' *          (with respect to z axis) is less than the critical cosine, the',
            ' *          incident packet is reflected from the boundary.',
            ' * @param[in] player Pointer to a layer object.',
            ' */',
            ' #define mc_layer_cc_outer(player) ((player)->cos_critical_outer)',
            '',
            '/**',
            ' * @brief Evaluates to the scattering coefficient of the layer.',
            ' * @param[in] player Pointer to a layer object.',
            ' * @param[in] pdir   Propagation direction vector.',
            ' */',
            '#define mc_layer_mus(player, pdir) ((player)->mus)',
            '',
            '/**',
            ' * @brief Evaluates to the absorption coefficient of the layer.',
            ' * @param[in] player Pointer to a layer object.',
            ' * @param[in] pdir   Propagation direction vector.',
            ' */',
            '#define mc_layer_mua(player, pdir) ((player)->mua)',
            '',
            '/**',
            ' * @brief Evaluates to the inverse of the absorption coefficient of the layer.',
            ' * @param[in] player Pointer to a layer object.',
            ' * @param[in] pdir   Propagation direction vector.',
            ' */',
            '#define mc_layer_inv_mua(player, pdir) \\',
            '	(((player)->mua != FP_0) ? mc_fdiv(FP_1, (player)->mua) : FP_INF)',
            '',
            '/**',
            ' * @brief Evaluates to the total attenuation coefficient, i.e. the sum of the',
            ' *        layer absorption and scattering coefficients.',
            ' * @param[in] player Pointer to a layer object.',
            ' * @param[in] pdir   Propagation direction vector.',
            ' */',
            '#define mc_layer_mut(player, pdir) ((player)->mua + (player)->mus)',
            '',
            '/**',
            ' * @brief Evaluates to the reciprocal of the total attenuation coefficient,',
            ' *        i.e. the reciprocal of the sum of the layer absorption and scattering',
            ' *        coefficients.',
            ' * @param[in] player Pointer to a layer object.',
            ' * @param[in] pdir   Propagation direction vector.',
            ' */',
            '#define mc_layer_inv_mut(player, pdir) ((player)->inv_mut)',
            '',
            '/**',
            ' * @brief Evaluates to the absorption coefficient of the layer multiplied',
            ' *        by the reciprocal of the total attenuation coefficient.',
            ' * @param[in] player Pointer to a layer object.',
            ' * @param[in] pdir   Propagation direction vector.',
            ' */',
            '#define mc_layer_mua_inv_mut(player, pdir) ((player)->mua_inv_mut)',
            '',
            '#if MC_ENABLE_DEBUG || defined(__DOXYGEN__)',
            '/**',
            ' * @brief Print one sample layer.',
            ' * param[in] prefix Can be used to pass indent string for the layer parameters."',
            ' * @param[in] player Pointer to a layer instance.',
            ' */',
            '#define dbg_print_layer(player, prefix) \\',
            '	printf( \\',
            '		prefix "r_inner: %.9f\\n" \\',
            '		prefix "r_outer: %.9f\\n" \\',
            '		prefix "n: %.9f\\n" \\',
            '		prefix "cos_critical_inner: %.9f\\n" \\',
            '		prefix "cos_critical_outer: %.9f\\n" \\',
            '		prefix "mua: %.9f\\n" \\',
            '		prefix "mus: %.9f\\n" \\',
            '		prefix "inv_mut: %.9f\\n" \\',
            '		prefix "mua_inv_mut: %.9f\\n", \\',
            '				(player)->r_inner, (player)->r_outer, (player)->n, \\',
            '				(player)->cos_critical_inner, (player)->cos_critical_outer, \\',
            '				(player)->mua, (player)->mus, \\',
            '				(player)->inv_mut, (player)->mua_inv_mut \\',
            '	); \\',
            '	{ McPf const _dbg_pf=(player)->pf; dbg_print_pf(&_dbg_pf);};',
            '',
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
            Layer diameter (m).
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
              Layer diameter (m).
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
        The inner layer boundary bellongs to the layer, but the outer does not.
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
    d = property(_get_d, _set_d, None, 'Layer diameter (m).')

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
        return cls(pf=pf_type.fromdict(pf_data), **data_)

    def __str__(self):
        return 'Layer(d={}, n={}, mua={}, mus={}, pf={})'.format(
            self._d, self._n, self._mua, self._mus, self._pf)

    def __repr__(self):
        return  '{:s} # id 0x{:>08X}.'.format(self.__str__(), id(self))


class Layers(mcobject.McObject):
    '''
    Class that represents a stack of concentric layers forming the sample.

    Note
    ----
    The first layer of the stack is used to describe the
    medium that surrounds the sample.
    Therefore, at least two layers must be always specified,
    namely the layer of the surrounding medium and one sample layer!
    The diameter of the outermost layer (first in the list of layers)
    will be automatically set to infinity regardless of the layer diameter
    set by the user.
    '''
    def __init__(self, layers: List[Layer] or 'Layers'):
        '''
        Constructs a managed sample layer stack from a list of sample layers.
        The outermost layer that sourrounds the sample must be specified first,
        followed by the sample layer from the outermost to the innermost.

        Note
        ----
        The first layer of the stack is used to describe the
        medium that surrounds the sample.
        Therefore, at least two layers must be always specified,
        namely one layer of the surrounding medium and one sample layer!

        The layers extend to infinity along both directions of the z axis.

        The diameter of the outermost layer will be automatically set
        to infinity when passed to the OpenCL kernel (regardless of the
        layer diameter set by the user).

        The layer diameters must be monotonically decreasing. Layers of
        zero thickness (diameters of the previous and next layer are the
        same) are not allowed.

        Note that all layers must use the same scattering phase function
        model.

        Parameters
        ----------
        layers: list or Layers
            A list of sample layers. Requires at least 2 items!
        '''
        self._pf_type = None

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
        if len(self._layers) < 2:
                ValueError('At least two layers are required, '
                           'but got only {:d}!'.format(len(self._layers)))

        if self._pf_type is None:
            self._pf_type = type(self._layers[1].pf)

        d_prev = float('inf')
        for layer_index, layer in enumerate(self._layers):
            if not isinstance(layer, Layer):
                raise TypeError('All layers must be instances of Layer '
                                'but found {:s}!'.format(type(layer).__name__))
            if type(layer.pf) != self._pf_type:
                raise TypeError(
                    'All the sample layer must use the same scattering phase '
                    'function model! Found {} and {}!'.format(
                        self._pf_type.__name__, type(layer.pf).__name__))
            if layer_index > 0 and layer.d > d_prev:
                raise ValueError('The diameters of layers must be '
                                 'monotonically decreasing!')
                d_prev = layer.d

    def layer(self, index: int) -> Layer:
        '''
        Returns layer at the specified index. Note that the first layer
        (index 0) represents the medium sourounding the sample.
        '''
        return self._layers[index]

    def layer_index(self, r: float) -> int:
        '''
        Returns the layer index that contains the given radius. Note that
        the layer includes the inner surface boundary but not the outer surface
        boundary, i.e. the r extent of the layer is [r_inner, r_outer), where
        r_outer > r_inner.

        Parameters
        ----------
        r: float
            R coordinate (radius) of a point.

        Returns
        -------
        layer_index: int
            Index of the layer that contains the give radius r.
        '''
        index = 0
        for pos, layer in enumerate(self._layers[::-1]):
            if 2.0*r < layer.d:
                index = len(self._layers) - 1 - pos
                break
        return index

    def diameter(self) -> float:
        '''
        Diameter of the layer stack excluding the outermost layer that
        sourrounds the sample.

        Returns
        -------
        diameter: float
            The sample diameter excluding the outermost layer of
            the surrounding medium.
        '''
        return self._layers[1].d

    def cl_options(self, mc: mcobject.McObject) -> str:
        '''
        OpenCL options of the scattering phase function.
        '''
        return self._layers[1].pf.fetch_cl_options(mc)

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        OpenCL declarations of the scattering phase function.
        '''
        return self._layers[1].pf.fetch_cl_declaration(mc)

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        '''
        OpenCL implementation of the scattering phase function.
        '''
        return self._layers[1].pf.fetch_cl_implementation(mc)

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
        '''
        self.check()

        num_layers = len(self._layers)

        if target is None or len(target) != num_layers:
            target_type = self.fetch_cl_type(mc)
            target = target_type()

        for index, layer in enumerate(self._layers):
            cc_outer = cc_inner = 0.0
            if index > 0:
                cc_outer = boundary.cos_critical(
                    layer.n, self._layers[index - 1].n)
            if index + 1 < num_layers:
                cc_inner = boundary.cos_critical(
                    layer.n, self._layers[index + 1].n)

            mut = layer.mua + layer.mus
            if mut > 0.0:
                inv_mut = 1.0/mut
            else:
                inv_mut = float('inf')
            if layer.mus == 0.0:
                mua_inv_mut = 1.0
            else:
                mua_inv_mut = layer.mua*inv_mut

            target[index].diameter = layer.d
            if index == 0:
                target[index].r_outer = float('inf')
                target[index].r_inner = self._layers[1].d*0.5
            else:
                target[index].r_outer = self._layers[index].d*0.5
                if index + 1 < num_layers:
                    target[index].r_inner = self._layers[index + 1].d*0.5
                else:
                    target[index].r_inner = 0.0
            target[index].n = layer.n
            target[index].cos_critical_outer = cc_outer
            target[index].cos_critical_inner = cc_inner
            target[index].mua = layer.mua
            target[index].mus = layer.mus
            target[index].inv_mut = inv_mut
            target[index].mua_inv_mut = mua_inv_mut

            layer.pf.cl_pack(mc, target[index].pf)

        return target

    def intersect(self, pos:Tuple[float, float, float],
                  dir:Tuple[float, float, float], entrance: bool = False) -> \
                      Tuple[Tuple[float, float, float] or None, \
                            Tuple[float, float, float] or None]:
        '''
        Intersect the sample with the ray and return the intersection.
        Intersection is considered to exist only if the distance to
        the intersection is positive or entrance is set to True

        Parameters
        ----------
        pos: (float, float, float)
            Ray origin.
        dir: (float, float, float)
            Ray direction.
        entrance: bool
            If True, compute the entrance point regardless of the current
            position of the ray (can be propagated backwards).

        Returns
        -------
        intersection: (float, float, float) or None
            Returns intersection if one exists, else None.
        normal: (float, float, float) or None
            Surface normal of the sample that points in the
            propagation direction.
        '''
        d1, d2 = ray_cylinder_intersection(self.layer(1).d*0.5, pos, dir)
        if d1 is None or d2 is None:
            return None, None

        if d1 < 0.0 and d2 < 0.0 and not entrance:
            return None, None

        if entrance:
            d = min(d1, d2)
        else:
            if d1 > 0.0 and d2 >= 0.0:
                d = min(d1, d2)
            else:
                d = max(d1, d2)

        intersection = pos[0] + dir[0]*d, pos[1] + dir[1]*d, pos[2] + dir[2]*d
        normal = (-intersection[0], -intersection[1], 0.0)

        return intersection, normal

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'layers': [layer.todict() for layer in self._layers],
                'type': 'Layers'}

    def tolist(self) -> List[Layer]:
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
        layers = [Layer.fromdict(item) for item in data_.pop('layers')]
        return cls(layers=layers, **data_)

    def __getitem__(self, what):
        return self._layers[what]

    def __setitem__(self, what, value):
        self._layers[what] = value

    def __len__(self):
        return len(self._layers)

    def __str__(self):
        layers = ['    ' + str(layer) for layer in self._layers]
        layers[0] +=  ",  # surrounding medium"
        layers_str = ',\n'.jtargetoin(layers)
        return 'Layers([\n{}\n])'.format(layers_str)

    def __repr__(self):
        return  '{:s} # id 0x{:>08X}.'.format(self.__str__(), id(self))
