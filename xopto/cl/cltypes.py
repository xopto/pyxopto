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

import ctypes

class _CLScalar:
    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.value)

    @classmethod
    def vector(cls, n: int):
        if n not in (2, 3, 4, 8, 16):
            raise ValueError(
                'OpenCL vector types must be of length 2, 3, 4, 8 or 16!')
        return cls._VECTORS[{2:0, 3:1, 4:2, 8:3, 16:4}.get(n)]

    def __str__(self):
        return self.__repr__()


class _CLVector2:
    def _get_x(self):
        return self[0]

    def _set_x(self, value):
        self[0] = value

    x = property(_get_x, _set_x)

    def _get_y(self):
        return self[1]

    def _set_y(self, value):
        self[1] = value

    y = property(_get_y, _set_y)

    def __repr__(self):
        return '{}({}, {})'.format(
            self.__class__.__name__, self[0], self[1])

    def __str__(self):
        return self.__repr__()


class _CLVector3(_CLVector2):
    def _get_z(self):
        return self[2]

    def _set_z(self, value):
        self[2] = value

    z = property(_get_z, _set_z)

    def __repr__(self):
        return '{}({}, {}, {})'.format(
            self.__class__.__name__, self[0], self[1], self[2])


class _CLVector4(_CLVector3):
    def _get_w(self):
        return self[3]

    def _set_w(self, value):
        self[3] = value

    w = property(_get_w, _set_w)

    def __repr__(self):
        return '{}{}'.format(self.__class__.__name__, tuple(self))


class _CLVector8(_CLVector4):
    pass

class _CLVector16(_CLVector8):
    pass

class cl_char2(_CLVector2, ctypes.c_char*2):
    TYPENAME='char2'
    dtype = 'int8'
class cl_char3(_CLVector3, ctypes.c_char*4):
    TYPENAME='char3'
    dtype = 'int8'
class cl_char4(_CLVector4, ctypes.c_char*4):
    TYPENAME='char4'
    dtype = 'int8'
class cl_char8(_CLVector8, ctypes.c_char*8):
    TYPENAME='char8'
    dtype = 'int8'
class cl_char16(_CLVector16, ctypes.c_char*16):
    TYPENAME='char16'
    dtype = 'int8'
class cl_char(_CLScalar, ctypes.c_char):
    TYPENAME='char'
    dtype = 'int8'
    _VECTORS = (cl_char2, cl_char3, cl_char4, cl_char8, cl_char16)

class cl_uchar2(_CLVector2, ctypes.c_ubyte*2):
    TYPENAME='uchar2'
    dtype = 'uint8'
class cl_uchar3(_CLVector3, ctypes.c_ubyte*4):
    TYPENAME='uchar3'
    dtype = 'uint8'
class cl_uchar4(_CLVector4, ctypes.c_ubyte*4):
    TYPENAME='uchar4'
    dtype = 'uint8'
class cl_uchar8(_CLVector8, ctypes.c_ubyte*8):
    TYPENAME='uchar8'
    dtype = 'uint8'
class cl_uchar16(_CLVector16, ctypes.c_ubyte*16):
    TYPENAME='uchar16'
    dtype = 'uint16'
class cl_uchar(_CLScalar, ctypes.c_ubyte):
    TYPENAME='uchar'
    dtype = 'uint8'
    _VECTORS = (cl_uchar2, cl_uchar3, cl_uchar4, cl_uchar8, cl_uchar16)

class cl_short2(_CLVector2, ctypes.c_int16*2):
    TYPENAME='short2'
    dtype = 'int16'
class cl_short3(_CLVector3, ctypes.c_int16*4):
    TYPENAME='short3'
    dtype = 'int16'
class cl_short4(_CLVector4, ctypes.c_int16*4):
    TYPENAME='short4'
    dtype = 'int16'
class cl_short8(_CLVector8, ctypes.c_int16*8):
    TYPENAME='short8'
    dtype = 'int16'
class cl_short16(_CLVector16, ctypes.c_int16*16):
    TYPENAME='short16'
    dtype = 'int16'
class cl_short(_CLScalar, ctypes.c_int16):
    TYPENAME='short'
    dtype = 'int16'
    _VECTORS = (cl_short2, cl_short3, cl_short4, cl_short8, cl_short16)

class cl_ushort2(_CLVector2, ctypes.c_uint16*2):
    TYPENAME='ushort2'
    dtype = 'uint16'
class cl_ushort3(_CLVector3, ctypes.c_uint16*4):
    TYPENAME='ushort3'
    dtype = 'uint16'
class cl_ushort4(_CLVector4, ctypes.c_uint16*4):
    TYPENAME='ushort4'
    dtype = 'uint16'
class cl_ushort8(_CLVector8, ctypes.c_uint16*8):
    TYPENAME='ushort8'
    dtype = 'int16'
class cl_ushort16(_CLVector16, ctypes.c_uint16*16):
    TYPENAME='ushort16'
    dtype = 'uint16'
class cl_ushort(_CLScalar, ctypes.c_uint16):
    TYPENAME='ushort'
    dtype = 'uint16'
    _VECTORS = (cl_ushort2, cl_ushort3, cl_ushort4, cl_ushort8, cl_ushort16)

class cl_int2(_CLVector2, ctypes.c_int32*2):
    TYPENAME='int2'
    dtype = 'int32'
class cl_int3(_CLVector3, ctypes.c_int32*4):
    TYPENAME='int3'
    dtype = 'int32'
class cl_int4(_CLVector4, ctypes.c_int32*4):
    TYPENAME='int4'
    dtype = 'int32'
class cl_int8(_CLVector8, ctypes.c_int32*8):
    TYPENAME='int8'
    dtype = 'int32'
class cl_int16(_CLVector16, ctypes.c_int32*16):
    TYPENAME='int16'
    dtype = 'int32'
class cl_int(_CLScalar, ctypes.c_int32):
    TYPENAME='int'
    dtype = 'int32'
    _VECTORS = (cl_int2, cl_int3, cl_int4, cl_int8, cl_int16)

class cl_uint2(_CLVector2, ctypes.c_uint32*2):
    TYPENAME='uint2'
    dtype = 'uint32'
class cl_uint3(_CLVector3, ctypes.c_uint32*4):
    TYPENAME='uint3'
    dtype = 'uint32'
class cl_uint4(_CLVector4, ctypes.c_uint32*4):
    TYPENAME='uint4'
    dtype = 'uint32'
class cl_uint8(_CLVector8, ctypes.c_uint32*8):
    TYPENAME='uint8'
    dtype = 'uint32'
class cl_uint16(_CLVector16, ctypes.c_uint32*16):
    TYPENAME='uint16'
    dtype = 'uint32'
class cl_uint(_CLScalar, ctypes.c_uint32):
    TYPENAME='uint'
    dtype = 'uint32'
    _VECTORS = (cl_uint2, cl_uint3, cl_uint4, cl_uint8, cl_uint16)

class cl_long2(_CLVector2, ctypes.c_int64*2):
    TYPENAME='long2'
    dtype = 'int64'
class cl_long3(_CLVector3, ctypes.c_int64*4):
    TYPENAME='long3'
    dtype = 'int64'
class cl_long4(_CLVector4, ctypes.c_int64*4):
    TYPENAME='long4'
    dtype = 'int64'
class cl_long8(_CLVector8, ctypes.c_int64*8):
    TYPENAME='long8'
    dtype = 'int64'
class cl_long16(_CLVector16, ctypes.c_int64*16):
    TYPENAME='long16'
    dtype = 'int64'
class cl_long(_CLScalar, ctypes.c_int64):
    TYPENAME='long'
    dtype = 'int64'
    _VECTORS = (cl_long2, cl_long3, cl_long4, cl_long8, cl_long16)

class cl_ulong2(_CLVector2, ctypes.c_uint64*2):
    TYPENAME='ulong2'
    dtype = 'uint64'
class cl_ulong3(_CLVector3, ctypes.c_uint64*4):
    TYPENAME='ulong3'
    dtype = 'uint64'
class cl_ulong4(_CLVector4, ctypes.c_uint64*4):
    TYPENAME='ulong4'
    dtype = 'uint64'
class cl_ulong8(_CLVector8, ctypes.c_uint64*8):
    TYPENAME='ulong8'
    dtype = 'uint64'
class cl_ulong16(_CLVector16, ctypes.c_uint64*16):
    TYPENAME='ulong16'
    dtype = 'uint64'
class cl_ulong(_CLScalar, ctypes.c_uint64):
    TYPENAME='ulong'
    dtype = 'uint64'
    _VECTORS = (cl_ulong2, cl_ulong3, cl_ulong4, cl_ulong8, cl_ulong16)

class cl_float2(_CLVector2, ctypes.c_float*2):
    TYPENAME='float2'
    dtype = 'float32'
class cl_float3(_CLVector3, ctypes.c_float*4):
    TYPENAME='float3'
    dtype = 'float32'
class cl_float4(_CLVector4, ctypes.c_float*4):
    TYPENAME='float4'
    dtype = 'float32'
class cl_float8(_CLVector8, ctypes.c_float*8):
    TYPENAME='float8'
    dtype = 'float32'
class cl_float16(_CLVector16, ctypes.c_float*16):
    TYPENAME='float16'
    dtype = 'float32'
class cl_float(_CLScalar, ctypes.c_float):
    TYPENAME='float'
    dtype = 'float32'
    _VECTORS = (cl_float2, cl_float3, cl_float4, cl_float8, cl_float16)

class cl_double2(_CLVector2, ctypes.c_double*2):
    TYPENAME='double2'
    dtype = 'float64'
class cl_double3(_CLVector3, ctypes.c_double*4):
    TYPENAME='double3'
    dtype = 'float64'
class cl_double4(_CLVector4, ctypes.c_double*4):
    TYPENAME='double4'
    dtype = 'float64'
class cl_double8(_CLVector8, ctypes.c_double*8):
    TYPENAME='double8'
    dtype = 'float64'
class cl_double16(_CLVector16, ctypes.c_double*16):
    TYPENAME='double16'
    dtype = 'float64'
class cl_double(_CLScalar, ctypes.c_double):
    TYPENAME='double'
    dtype = 'float64'
    _VECTORS = (cl_double2, cl_double3, cl_double4, cl_double8, cl_double16)

# aliases
cl_byte = cl_char
cl_ubyte = cl_uchar
cl_int8_t = cl_char
cl_uint8_t = cl_uchar
cl_int16_t = cl_short
cl_uint16_t = cl_ushort
cl_int32_t = cl_int
cl_uint32_t = cl_uint
cl_int64_t = cl_long
cl_uint64_t = cl_ulong
cl_float32_t = cl_float
cl_float64_t = cl_double

sizeof = ctypes.sizeof
ARRAY = ctypes.ARRAY
Array = ctypes.Array
byref = ctypes.byref
addressof = ctypes.addressof

class Structure(ctypes.Structure):
    def __str__(self):
        fields = []
        for field_name, field_type in self._fields_:
            field_value = getattr(self, field_name)
            value = getattr(field_value, 'value', field_value)

            fields.append('{:s}={:s}'.format(field_name, str(value)))

        return '{:s}({:s})'.format(self.__class__.__name__, ', '.join(fields))

    def __repr__(self):
        return '{:s} #id {}'.format(self.__str__(), id(self))


if __name__ == '__main__':
    class T(Structure):
        _pack_ = 1
        _fields_ = [
            ('a', cl_int3),
            ('b', ctypes.c_float*2)
        ]
    print(T((1,2,4), (2, 2)))
