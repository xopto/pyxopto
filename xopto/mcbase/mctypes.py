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

import math
from typing import List, Tuple

import numpy as np
from jinja2 import Template

from xopto.mcbase import cltypes
from xopto.mcbase.mcobject import McObject


class StructToListBasic:
    def tolist(self):
        return [getattr(self, name).value for name, _ in self._fields_]


class Matrix2Helper:
    def fromarray(self, array: np.ndarray):
        self.a_11, self.a_12 = array[0, 0], array[0, 1]
        self.a_21, self.a_22 = array[1, 0], array[1, 1]

    def toarray(self) -> np.ndarray:
        return np.array([
            [self.a_11.value, self.a_12.value],
            [self.a_21.value, self.a_22.value],
        ])

class Matrix3Helper:
    def fromarray(self, array: np.ndarray):
        self.a_11, self.a_12, self.a_13 = array[0, 0], array[0, 1], array[0, 2]
        self.a_21, self.a_22, self.a_23 = array[1, 0], array[1, 1], array[1, 2]
        self.a_31, self.a_32, self.a_33 = array[2, 0], array[2, 1], array[2, 2]

    def toarray(self) -> np.ndarray:
        return np.array([
            [self.a_11.value, self.a_12.value, self.a_13.value],
            [self.a_21.value, self.a_22.value, self.a_23.value],
            [self.a_31.value, self.a_32.value, self.a_33.value]
        ])


class Vector4Helper:
    def fromarray(self, array: np.ndarray):
        self.x, self.y, self.z, self.w = array[0], array[1], array[2], array[3]

    def toarray(self) -> np.ndarray:
        return np.array([self.x.value, self.y.value, self.z.value, self.w.value])

class Vector3Helper:
    def fromarray(self, array: np.ndarray):
        self.x, self.y, self.z = array[0], array[1], array[2]

    def toarray(self) -> np.ndarray:
        return np.array([self.x.value, self.y.value, self.z.value])


class Vector2Helper:
    def fromarray(self, array: np.ndarray):
        self.x, self.y = array[0], array[1]

    def toarray(self) -> np.ndarray:
        return np.array([self.x.value, self.y.value])

class McSizeT32:
    '''
    Class that represents 32-bit size type.
    The OpenCL data types are defined by the following attributes:

    - mc_size_t - Size type.
    - mc_size<n>_t for vectorized unsigned integers (<n> must be one of 2, 3, 4, 8, 16)
    - mc_point2s_t structure for a 2D point with mc_size_t coordinates.
    - mc_point3s_t structure for a 3D point with mc_size_t coordinates.

    The related numpy data types are defined by the following attributes:

    - np_size - For mc_size_t.

    Additional attributes of the class:

    - mc_maxsize - maximum unsigned integer value that can be represented by the size type.
    '''
    mc_size_t = cltypes.cl_uint
    mc_size2_t = cltypes.cl_uint2
    mc_size3_t = cltypes.cl_uint3
    mc_size4_t = cltypes.cl_uint4
    mc_size8_t = cltypes.cl_uint8
    mc_size16_t = cltypes.cl_uint16

    class mc_point2s_t(cltypes.Structure, StructToListBasic, Vector3Helper):
        ''' Structure representing a 2D point using mc_size_t. '''
        _fields_ = [
            ('x', cltypes.cl_uint),
            ('y', cltypes.cl_uint),
        ]

    class mc_point3s_t(cltypes.Structure, StructToListBasic, Vector3Helper):
        ''' Structure representing a 3D point using mc_size_t. '''
        _fields_ = [
            ('x', cltypes.cl_uint),
            ('y', cltypes.cl_uint),
            ('z', cltypes.cl_uint)
        ]

    class mc_point4s_t(cltypes.Structure, StructToListBasic, Vector4Helper):
        ''' Structure representing a 4D point using mc_size_t. '''
        _fields_ = [
            ('x', cltypes.cl_uint),
            ('y', cltypes.cl_uint),
            ('z', cltypes.cl_uint),
            ('w', cltypes.cl_uint)
        ]

    mc_maxsize = 4294967295
    np_size = 'uint32'

class McSizeT64:
    '''
    Class that represents 64-bit size type.
    The OpenCL data types are defined by the following attributes:

    - mc_size_t - Size type.
    - mc_size<n>_t - For vectorized unsigned integers (<n> must be one of 2, 3, 4, 8, 16)
    - mc_point2_t - For a structure representing a 2D point with mc_int_t coordinates.
    - mc_point3_t - For a structure representing a 3D point with mc_int_t coordinates.
    - mc_point4_t - For a structure representing a 4D point with mc_int_t coordinates.
    - mc_matrix3_t - For a structure representing a 3x3 matrix with mc_int_t components.

    The related numpy data types are defined by the following attributes:

    - np_size - For mc_size_t.

    Additional attributes of the class:

    - mc_maxsize - maximum unsigned integer value that can be represented by the size type.
    '''

    mc_size_t = cltypes.cl_ulong
    mc_size2_t = cltypes.cl_ulong2
    mc_size3_t = cltypes.cl_ulong3
    mc_size4_t = cltypes.cl_ulong4
    mc_size8_t = cltypes.cl_ulong8
    mc_size16_t = cltypes.cl_ulong16

    class mc_point2s_t(cltypes.Structure, StructToListBasic, Vector3Helper):
        ''' Structure representing a 2D point using mc_size_t. '''
        _fields_ = [
            ('x', cltypes.cl_ulong),
            ('y', cltypes.cl_ulong),
        ]

    class mc_point3s_t(cltypes.Structure, StructToListBasic, Vector3Helper):
        ''' Structure representing a 3D point using mc_size_t. '''
        _fields_ = [
            ('x', cltypes.cl_ulong),
            ('y', cltypes.cl_ulong),
            ('z', cltypes.cl_ulong)
        ]

    class mc_point4s_t(cltypes.Structure, StructToListBasic, Vector4Helper):
        ''' Structure representing a 4D point using mc_size_t. '''
        _fields_ = [
            ('x', cltypes.cl_ulong),
            ('y', cltypes.cl_ulong),
            ('z', cltypes.cl_ulong),
            ('w', cltypes.cl_ulong)
        ]

    mc_maxsize = 18446744073709551615
    np_size = 'uint64'


class McInt32:
    '''
    Class that represents a 32-bit default integer type.
    The OpenCL data types are defined by the following attributes:

    - mc_int_t - For signed integers.
    - mc_int<n>_t - For vectorized signed integers (<n> must be one of 2, 3, 4, 8, 16).
    - mc_uint_t - For unsigned integers.
    - mc_uint<n>_t - For vectorized unsigned integers (<n> must be one of 2, 3, 4, 8, 16).
    - mc_point2_t - For a structure representing a 2D point with mc_int_t coordinates.
    - mc_point3_t - For a structure representing a 3D point with mc_int_t coordinates.
    - mc_point4_t - For a structure representing a 4D point with mc_int_t coordinates.
    - mc_matrix3_t - For a structure representing a 3x3 matrix with mc_int_t components.

    The related numpy data types are defined by the following attributes:

    - np_int - For signed integers.
    - np_uint - For unsigned integers.

    Additional attributes of the class:

    - mc_uint_max - maximum unsigned integer value that can be represented by the data type.
    '''
    mc_int_t = cltypes.cl_int
    mc_int2_t = cltypes.cl_int2
    mc_int3_t = cltypes.cl_int3
    mc_int4_t = cltypes.cl_int4
    mc_int8_t = cltypes.cl_int8
    mc_int16_t = cltypes.cl_int16

    mc_uint_t = cltypes.cl_uint
    mc_uint2_t = cltypes.cl_uint2
    mc_uint3_t = cltypes.cl_uint3
    mc_uint4_t = cltypes.cl_uint4
    mc_uint8_t = cltypes.cl_uint8
    mc_uint16_t = cltypes.cl_uint16

    class mc_point4_t(cltypes.Structure, StructToListBasic, Vector4Helper):
        ''' Structure representing a 4D point using 32-bit integers. '''
        _fields_ = [
            ('x', cltypes.cl_int),
            ('y', cltypes.cl_int),
            ('z', cltypes.cl_int),
            ('w', cltypes.cl_int)
        ]

    class mc_point3_t(cltypes.Structure, StructToListBasic, Vector3Helper):
        ''' Structure representing a 3D point using 32-bit integers. '''
        _fields_ = [
            ('x', cltypes.cl_int),
            ('y', cltypes.cl_int),
            ('z', cltypes.cl_int)
        ]

    class mc_point2_t(cltypes.Structure, StructToListBasic, Vector2Helper):
        ''' Structure representing a 2D point using 32-bit integers. '''
        _fields_ = [
            ('x', cltypes.cl_int),
            ('y', cltypes.cl_int),
        ]

    class mc_matrix3_t(cltypes.Structure, StructToListBasic, Matrix3Helper):
        ''' Structure representing a 3x3 matrix using 32-bit integers. '''
        _fields_ = [
            ('a_11', cltypes.cl_int), ('a_12', cltypes.cl_int), ('a_13', cltypes.cl_int),
            ('a_21', cltypes.cl_int), ('a_22', cltypes.cl_int), ('a_23', cltypes.cl_int),
            ('a_31', cltypes.cl_int), ('a_32', cltypes.cl_int), ('a_33', cltypes.cl_int)
        ]

    mc_uint_max = 4294967296
    mc_int_max = 2147483647

    np_int = 'int32'
    np_uint = 'uint32'

class McInt64:
    '''
    Class that represents a 64-bit default integer type.
    The OpenCL data types are defined by the following attributes:

    - mc_int_t - For signed integers.
    - mc_int<n>_t - For vectorized signed integers (<n> must be one of 2, 3, 4, 8, 16).
    - mc_uint_t - For unsigned integers.
    - mc_uint<n>_t - For vectorized unsigned integers (<n> must be one of 2, 3, 4, 8, 16).
    - mc_point2_t - For a Structure representing a 2D point with mc_int_t coordinates.
    - mc_point3_t - For a Structure for representing a 3D point with mc_int_t coordinates.
    - mc_point4_t - For a Structure for representing a 4D point with mc_int_t coordinates.
    - mc_matrix3_t - For a Structure representing a 3x3 matrix with mc_int_t elements.

    The related numpy data types are defined by the following attributes:

    - mc_int_t - For signed integers.
    - mc_uint_t - For unsigned integers.

    Additional attributes of the class:

    - mc_uint_max - maximum unsigned integer value that can be represented by the data type.
    '''
    mc_int_t = cltypes.cl_long
    mc_int2_t = cltypes.cl_long2
    mc_int3_t = cltypes.cl_long3
    mc_int4_t = cltypes.cl_long4
    mc_int8_t = cltypes.cl_long8
    mc_int16_t = cltypes.cl_long16

    mc_uint_t = cltypes.cl_ulong
    mc_uint2_t = cltypes.cl_ulong2
    mc_uint3_t = cltypes.cl_ulong3
    mc_uint4_t = cltypes.cl_ulong4
    mc_uint8_t = cltypes.cl_ulong8
    mc_uint16_t = cltypes.cl_ulong16

    class mc_point4_t(cltypes.Structure, StructToListBasic, Vector4Helper):
        ''' Structure representing a 4D point using 64-bit integers. '''
        _fields_ = [
            ('x', cltypes.cl_long),
            ('y', cltypes.cl_long),
            ('z', cltypes.cl_long),
            ('w', cltypes.cl_long)
        ]

    class mc_point3_t(cltypes.Structure, StructToListBasic, Vector3Helper):
        ''' Structure representing a 3D point using 64-bit integers. '''
        _fields_ = [
            ('x', cltypes.cl_long),
            ('y', cltypes.cl_long),
            ('z', cltypes.cl_long)
        ]

    class mc_point2_t(cltypes.Structure, StructToListBasic, Vector2Helper):
        ''' Structure representing a 2D point using 64-bit integers. '''
        _fields_ = [
            ('x', cltypes.cl_long),
            ('y', cltypes.cl_long),
        ]

    class mc_matrix3_t(cltypes.Structure, StructToListBasic, Matrix3Helper):
        ''' Structure representing a 3x3 matrix using 64-bit integers. '''
        _fields_ = [
            ('a_11', cltypes.cl_long), ('a_12', cltypes.cl_long), ('a_13', cltypes.cl_long),
            ('a_21', cltypes.cl_long), ('a_22', cltypes.cl_long), ('a_23', cltypes.cl_long),
            ('a_31', cltypes.cl_long), ('a_32', cltypes.cl_long), ('a_33', cltypes.cl_long)
        ]

    mc_uint_max = 18446744073709551616
    mc_int_max = 9223372036854775807

    np_uint = 'int64'
    np_uint = 'uint64'

McInt = McInt32

class McAccu32:
    '''
    Class that represents a 32-bit detector accumulator type. Tends to quickly
    overflow, even with a relatively small number of simulated photon packets.
    This class exposes the following attributes:

    - mc_accu_t - Defines the related OpenCL data type.
    - np_accu - Defines the related numpy data type.
    - mc_accu_k - Defines the factor that converts floating-point photon packet weight to an integer.
    - mc_accu_max - Defines the maximum integer value that can be represented by the detector accumulator type.
    '''
    mc_accu_t = cltypes.cl_uint
    mc_accu_k = 0x7FFFFF # 23 bits
    mc_cnt_max = 4294967296
    np_accu = 'uint32'


class McAccu64:
    '''
    Class that represents a 64-bit detector accumulator type.
    This class exposes the following attributes:

    - mc_accu_t - Defines the related OpenCL data type.
    - np_accu - Defines the related numpy data type.
    - mc_accu_k - Defines the factor that converts floating-point photon packet weight to an integer.
    - mc_accu_max - Defines the maximum integer value that can be represented by the detector accumulator type.
    '''
    mc_accu_t = cltypes.cl_ulong
    mc_accu_k = 0x7FFFFF # 23 bits
    mc_accu_max = 18446744073709551616
    np_accu = 'uint64'


class McCnt32:
    '''
    Class that represents a 32-bit photon packet counter type.
    The maximum number of photon packets that can be processed in a single
    OpenCL call of the Monte Carlo kernel is limited to 4,294,967,296.
    This class exposes the following attributes:

    - mc_cnt_t - Defines the related OpenCL data type.
    - np_cnt - Defines the related numpy data type.
    - mc_cnt_max - Defines the maximum integer value that can be represented by the photon packet counter type.
    '''
    mc_cnt_t = cltypes.cl_uint
    mc_cnt_max = 4294967296
    np_cnt = 'uint32'


class McCnt64:
    '''
    Class that represents a 64-bit photon packet counter type.
    The maximum number of photon packets that can be processed in a single
    OpenCL call of the Monte Carlo kernel is for all practical purposes
    virtually unlimited (18,446,744,073,709,551,616).
    This class exposes the following attributes:

    - mc_cnt_t - Defines the related OpenCL data type.
    - np_cnt - Defines the related numpy data type.
    - mc_cnt_max - Defines the maximum integer value that can be represented by the photon packet counter type.
    '''
    mc_cnt_t = cltypes.cl_ulong
    mc_cnt_max = 18446744073709551616
    np_cnt = 'uint64'


McCnt = McCnt32


class McFloat:
    '''
    Class that represents a single precision floating-point data type/arithmetics.
    This class implements the following attributes that define the related
    OpenCL data types:

    - mc_fp_t - Scalar floating-point data type.
    - mc_fp<n>_t - Vectorized floating point data types (<n> must be one of 2, 3, 4, 8, 16).
    - mc_point2f_t - Structure representing a 2D point with mc_fp_t coordinates.
    - mc_point3f_t - Structure representing a 3D point with mc_fp_t coordinates.
    - mc_point4f_t - Structure representing a 4D point with mc_fp_t coordinates.
    - mc_matrix2f_t - Structure representing a 2x2 matrix with mc_fp_t components.
    - mc_matrix3f_t - Structure representing a 3x3 matrix with mc_fp_t components.

    Attribute that defined the related numpy data type:

    - np_float - Scalar floating-point data type.

    Additional attributes of the class:

    - eps - The difference between 1.0 and the next smallest representable float larger than 1.0.
    - mc_fp_maxint - Maximum unsigned integer value that can be represented by the float type.
    '''
    mc_fp_t = cltypes.cl_float
    mc_fp2_t = cltypes.cl_float2
    mc_fp3_t = cltypes.cl_float3
    mc_fp4_t = cltypes.cl_float4
    mc_fp8_t = cltypes.cl_float8
    mc_fp16_t = cltypes.cl_float16
    eps = 1.1920929e-07
    mc_fp_maxint = 0x7FFFFF # 23 bits
    np_float = 'float32'

    class mc_point2f_t(cltypes.Structure, StructToListBasic, Vector2Helper):
        ''' Structure representing a 2D point using single-precision floating-point. '''
        _fields_ = [
            ('x', cltypes.cl_float),
            ('y', cltypes.cl_float),
        ]

    class mc_point3f_t(cltypes.Structure, StructToListBasic, Vector3Helper):
        ''' Structure representing a 3D point using single-precision floating-point. '''
        _fields_ = [
            ('x', cltypes.cl_float),
            ('y', cltypes.cl_float),
            ('z', cltypes.cl_float)
        ]

    class mc_point4f_t(cltypes.Structure, StructToListBasic, Vector4Helper):
        ''' Structure representing a 4D point using single-precision floating-point. '''
        _fields_ = [
            ('x', cltypes.cl_float),
            ('y', cltypes.cl_float),
            ('z', cltypes.cl_float),
            ('w', cltypes.cl_float)
        ]

    class mc_matrix2f_t(cltypes.Structure, StructToListBasic, Matrix2Helper):
        ''' Structure representing a 2x2 matrix using single-precision floating-point. '''
        _fields_ = [
            ('a_11', cltypes.cl_float), ('a_12', cltypes.cl_float),
            ('a_21', cltypes.cl_float), ('a_22', cltypes.cl_float),
        ]

    class mc_matrix3f_t(cltypes.Structure, StructToListBasic, Matrix3Helper):
        ''' Structure representing a 3x3 matrix using single-precision floating-point. '''
        _fields_ = [
            ('a_11', cltypes.cl_float), ('a_12', cltypes.cl_float), ('a_13', cltypes.cl_float),
            ('a_21', cltypes.cl_float), ('a_22', cltypes.cl_float), ('a_23', cltypes.cl_float),
            ('a_31', cltypes.cl_float), ('a_32', cltypes.cl_float), ('a_33', cltypes.cl_float),
        ]

McFloat32 = McFloat
McSingle = McFloat


class McDouble:
    '''
    Class that represents a double precision floating-point data type/arithmetics.
    This class implements the following attributes that define the related
    OpenCL data types:

    - mc_fp_t - Scalar floating-point data type.
    - mc_fp<n>_t - Vectorized floating point data types (<n> must be one of 2, 3, 4, 8, 16).
    - mc_point2f_t - Structure representing a 2D point with mc_fp_t coordinates.
    - mc_point3f_t - Structure representing a 3D point with mc_fp_t coordinates.
    - mc_point4f_t - Structure representing a 4D point with mc_fp_t coordinates.
    - mc_matrix2f_t - Structure representing a 2x2 matrix with mc_fp_t components.
    - mc_matrix3f_t - Structure representing a 3x3 matrix with mc_fp_t components.

    Attribute that defined the related numpy data type:

    - np_float - Scalar floating-point data type.

    Additional attributes of the class:

    - eps - The difference between 1.0 and the next smallest representable float larger than 1.0. 
    - mc_fp_maxint - Maximum unsigned integer value that can be represented by the float type.
    '''
    mc_fp_t = cltypes.cl_double
    mc_fp2_t = cltypes.cl_double2
    mc_fp3_t = cltypes.cl_double3
    mc_fp4_t = cltypes.cl_double4
    mc_fp8_t = cltypes.cl_double8
    mc_fp16_t = cltypes.cl_double16
    eps = 2.220446049250313e-16
    mc_fp_maxint = 0xFFFFFFFFFFFFF # 52 bits
    np_float = 'float64'

    class mc_point2f_t(cltypes.Structure, StructToListBasic, Vector2Helper):
        ''' Structure representing a 2D point using double-precision floating-point. '''
        _fields_ = [
            ('x', cltypes.cl_double),
            ('y', cltypes.cl_double),
        ]

    class mc_point3f_t(cltypes.Structure, StructToListBasic, Vector3Helper):
        ''' Structure representing a 3D point using single-precision floating-point. '''
        _fields_ = [
            ('x', cltypes.cl_double),
            ('y', cltypes.cl_double),
            ('z', cltypes.cl_double)
        ]

    class mc_point4f_t(cltypes.Structure, StructToListBasic, Vector4Helper):
        ''' Structure representing a 4D point using single-precision floating-point. '''
        _fields_ = [
            ('x', cltypes.cl_double),
            ('y', cltypes.cl_double),
            ('z', cltypes.cl_double),
            ('w', cltypes.cl_double)
        ]

    class mc_matrix2f_t(cltypes.Structure, StructToListBasic, Matrix2Helper):
        ''' Structure representing a 2x2 matrix using double-precision floating-point. '''
        _fields_ = [
            ('a_11', cltypes.cl_double), ('a_12', cltypes.cl_double),
            ('a_21', cltypes.cl_double), ('a_22', cltypes.cl_double)
        ]

    class mc_matrix3f_t(cltypes.Structure, StructToListBasic, Matrix3Helper):
        ''' Structure representing a 3x3 matrix using double-precision floating-point. '''
        _fields_ = [
            ('a_11', cltypes.cl_double), ('a_12', cltypes.cl_double), ('a_13', cltypes.cl_double),
            ('a_21', cltypes.cl_double), ('a_22', cltypes.cl_double), ('a_23', cltypes.cl_double),
            ('a_31', cltypes.cl_double), ('a_32', cltypes.cl_double), ('a_33', cltypes.cl_double)
        ]

McFloat64 = McDouble


class _McDataTypesMeta(type):
    def __repr__(self):
        int_t = {False:McInt32, True:McInt64}.get(
            cltypes.sizeof(self.mc_int_t) == 8)
        accu_t = {False:McAccu32, True:McAccu64}.get(
            cltypes.sizeof(self.mc_accu_t) == 8)
        cnt_t = {False:McCnt32, True:McCnt64}.get(
            cltypes.sizeof(self.mc_cnt_t) == 8)
        fp_t = {False:McSingle, True:McDouble}.get(
            cltypes.sizeof(self.mc_fp_t) == 8)
        return 'class McDataTypes({}, {}, {}, {})'.format(
            int_t.__name__, accu_t.__name__, cnt_t.__name__, fp_t.__name__)

    def __str__(self):
        return self.__repr__()


class _McDataTypesHelper:
    @classmethod
    def cl_options(cls, *_):
        '''
        Returns a list of OpenCL options that set the selected default
        data types of the Monte Carlo kernel.
        '''
        return [
            ('MC_USE_DOUBLE_PRECISION', cltypes.sizeof(cls.mc_fp_t) == 8),
            ('MC_USE_64_BIT_SIZE_T', cltypes.sizeof(cls.mc_size_t) == 8),
            ('MC_USE_64_BIT_INTEGER', cltypes.sizeof(cls.mc_int_t) == 8),
            ('MC_USE_64_BIT_PACKET_COUNTER', cltypes.sizeof(cls.mc_cnt_t) == 8),
            ('MC_USE_64_BIT_ACCUMULATORS', cltypes.sizeof(cls.mc_accu_t) == 8),
            ('MC_INT_ACCUMULATOR_K', cls.mc_accu_k),
        ]

    @classmethod
    def double_precision(cls) -> 'McDataTypes':
        ''' Returns a double-precision variant of the data type. '''
        return cls.custom(McDouble)

    @classmethod
    def single_precision(cls) -> 'McDataTypes':
        ''' Returns a single-precision variant of the data type. '''
        return cls.custom(McFloat)

    @classmethod
    def custom(cls, *args) -> 'McDataTypes':
        '''
        Creates a custom data type configuration for the Monte Carlo kernel.

        Parameters
        ---------- 
        args: tuple
            Positional arguments that customize the default kernel data types:

            - size type

                - Use McSizeT32 for 32-bit unsigned integer size type.
                - Use McSizeT64 for 64-bit unsigned integer size type.
 
            - integer type

                - Use McInt32 for 32-bit integers.
                - Use McInt64 for 64-bit integers. 

            - photon packet counter type

                - Use McCnt32 for a 32-bit unsigned counter.
                - Use McCnt64 for a 64-bit unsigned counter.

            - floating point precision

                - Use McSingle for single precision floating-point arithmetics.
                - Use McDouble for double precision floating-point arithmetics. 

            - accumulator integer type

                - Use McAccu32 for 32-bit detector accumulators (not recommended). 
                - Use McAccu64 for 64-bit detector accumulators. 
        '''
        class McDataTypesCustom(*cls.make_base(*args),
                                metaclass =_McDataTypesMeta):
            pass
        return McDataTypesCustom

    @classmethod
    def make_base(cls, *args):
        '''
        Creates a customized list of base classes for the specified data
        type configuration for the OpenCL kernel.

        Parameters
        ---------- 
        args: tuple
            Positional arguments that customize the default kernel data types:

            - size type

                - Use McSizeT32 for 32-bit unsigned integer size type.
                - Use McSizeT64 for 64-bit unsigned integer size type. 

            - integer type

                - Use McInt32 for 32-bit integers.
                - Use McInt64 for 64-bit integers. 

            - photon packet counter type

                - Use McCnt32 for a 32-bit unsigned counter.
                - Use McCnt64 for a 64-bit unsigned counter.

            - floating point precision

                - Use McSingle for single precision floating-point arithmetics.
                - Use McDouble for double precision floating-point arithmetics. 

            - accumulator integer type

                - Use McAccu32 for 32-bit detector accumulators (not recommended). 
                - Use McAccu64 for 64-bit detector accumulators. 
        '''
        base = list(cls.__bases__)
        for item in args:
            if item == McDouble and McSingle in base:
                base[base.index(McSingle)] = McDouble
            elif item == McSingle and McDouble in base:
                base[base.index(McDouble)] = McSingle

            elif item == McInt32 and McInt64 in base:
                base[base.index(McInt64)] = McInt32
            elif item == McInt64 and McInt32 in base:
                base[base.index(McInt32)] = McInt64

            elif item == McCnt32 and McCnt64 in base:
                base[base.index(McCnt64)] = McCnt32
            elif item == McCnt64 and McCnt32 in base:
                base[base.index(McCnt32)] = McCnt64

            elif item == McSizeT32 and McSizeT64 in base:
                base[base.index(McSizeT64)] = McSizeT32
            elif item == McSizeT64 and McSizeT32 in base:
                base[base.index(McSizeT32)] = McSizeT64

            elif item == McAccu32 and McAccu64 in base:
                base[base.index(McAccu64)] = McAccu32
            elif item == McAccu64 and McAccu32 in base:
                base[base.index(McAccu32)] = McAccu64
        return base

    @classmethod
    def mc_fp_lut_t(cls):
        '''
        Returns a ctypes structure type that represents a linear lookup table.
        '''
        class ClFpLut(cltypes.Structure):
            ''' Structure representing a lookup table. '''
            _fields_ = [
                ('first', cls.mc_fpt_t),
                ('inv_span', cls.mc_fp_t),
                ('n', cls.mc_size_t),
                ('offset', cls.mc_size_t),
            ]
        return ClFpLut

class McDataTypesBase:
    pass

class McDataTypes(McFloat32, McCnt32, McInt32, McAccu64, McSizeT32,
                  McDataTypesBase, _McDataTypesHelper, McObject,
                  metaclass =_McDataTypesMeta):
    '''
    Class representing Monte Carlo kernel data types. Use the
    customize static method to customize the default integer,
    floating-point, detector accumulator and photon packet counter data types.

    This class sets up the default set of kernel data types:

    - size type
        Use McSizeT32 for 32-bit unsigned integer size type.
    - integer type
        Use McInt32 for 32-bit integers.
    - photon packet counter type
        Use McCnt32 for a 32-bit unsigned counter.
    - floating point precision
        Use McSingle for single precision floating-point arithmetics.
    - accumulator integer type
        Use McAccu64 for 64-bit detector accumulators.
    '''
    pass


class McDataTypesSingle(*McDataTypes.make_base(), metaclass =_McDataTypesMeta):
    '''
    Class will set up the following OpenCL data types:

    - size type
        Use McSizeT32 for 32-bit unsigned integer size type.
    - integer type
        Use McInt32 for 32-bit integers.
    - photon packet counter type
        Use McCnt32 for a 32-bit unsigned counter.
    - floating point precision
        Use McSingle for single precision floating-point arithmetics.
    - accumulator integer type
        Use McAccu64 for 64-bit detector accumulators.
    '''
    @classmethod
    def double_precision(cls) -> 'McDataTypesDouble':
        ''' Returns a double-precision floating-point variant of the data type. '''
        return McDataTypesDouble

    @classmethod
    def single_precision(cls) -> 'McDataTypesSingle':
        ''' Returns a single-precision floating-point variant of the data type. '''
        return cls


class McDataTypesSingleCnt64(*McDataTypes.make_base(McCnt64, McSizeT64),
                             metaclass =_McDataTypesMeta):
    '''
    Class will set up the following OpenCL data types:

    - size type
        Use McSizeT64 for 64-bit unsigned integer size type.
    - integer type
        Use McInt32 for 32-bit integers.
    - photon packet counter type
        Use McCnt64 for a 64-bit unsigned counter.
    - floating point precision
        Use McSingle for single precision floating-point arithmetics.
    - accumulator integer type
        Use McAccu64 for 64-bit detector accumulators.
    '''
    @classmethod
    def double_precision(cls) -> 'McDataTypesDoubleCnt64':
        ''' Returns a double-precision floating-point variant of the data type. '''
        return McDataTypesDoubleCnt64

    @classmethod
    def single_precision(cls) -> 'McDataTypesSingleCnt64':
        ''' Returns a single-precision floating-point variant of the data type. '''
        return cls


class McDataTypesDouble(*McDataTypes.make_base(McDouble),
                         metaclass =_McDataTypesMeta):
    '''
    Class will set up the following OpenCL data types:

    - size type
        Use McSizeT32 for 32-bit unsigned integer size type.
    - integer type
        Use McInt32 for 32-bit integers.
    - photon packet counter type
        Use McCnt32 for a 32-bit unsigned counter.
    - floating point precision
        Use McDouble for double precision floating-point arithmetics.
    - accumulator integer type
        Use McAccu64 for 64-bit detector accumulators.
    '''
    @classmethod
    def double_precsision(cls) -> 'McDataTypesDouble':
        ''' Returns a double-precision floating-point variant of the data type. '''
        return cls

    @classmethod
    def single_precision(cls) -> McDataTypesSingle:
        ''' Returns a single-precision floating-point variant of the data type. '''
        return McDataTypesSingle


class McDataTypesDoubleCnt64(
        *McDataTypes.make_base(McDouble, McCnt64, McSizeT64),
         metaclass =_McDataTypesMeta):
    '''
    Class will set up the following OpenCL data types:

    - size type
        Use McSizeT64 for 64-bit unsigned integer size type.
    - integer type
        Use McInt32 for 32-bit integers.
    - photon packet counter type
        Use McCnt32 for a 32-bit unsigned counter.
    - floating point precision
        Use McDouble for double precision floating-point arithmetics.
    - accumulator integer type
        Use McAccu64 for 64-bit detector accumulators.
    '''
    @classmethod
    def double_precsision(cls) -> 'McDataTypesDoubleCnt64':
        ''' Returns a double-precision floating-point variant of the data type. '''
        return cls

    @classmethod
    def single_precision(cls) -> McDataTypesSingleCnt64:
        ''' Returns a single-precision floating-point variant of the data type. '''
        return McDataTypesSingleCnt64


class OpenCLSource:
    def __init__(self, headers: List[str], sources: List[str]):
        '''
        Merges OpenCL header and source files in the listed order.

        Parameters
        ---------- 
        headers: list
            List of OpenCL headers.
        sources: list
            List of OpenCL headers.
        '''
        self._cl_sources = '\n'.join(headers) + '\n'.join(sources)
        self._cl_template = Template(self._cl_source)

    def render(self, data: dict):
        default = {
            'mc': {'options': ''},
            'source': {'declaration':'', 'implementation': ''},
            'pf': {'declaration':'', 'implementation': ''},
            'rt': {'declaration':'', 'implementation': ''},
            'trace': {'declaration':'', 'implementation': ''},
            'fluence': {'declaration':'', 'implementation': ''},
            'top_geometry': {'declaration':'', 'implementation': ''},
            'bottom_geometry': {'declaration':'', 'implementation': ''}
        }
        default.update(data)

        return self._cl_template.render(**default)
