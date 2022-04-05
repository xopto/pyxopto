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

from typing import List, Tuple, Type

import numpy as np

from xopto.mcvox import mcobject
from xopto.mcvox import mctypes
from xopto.mcvox import cltypes
from xopto.mcvox import mcpf
from xopto.mcvox.mcutil.axis import Axis

class Voxel(mcobject.McObject):
    '''
    Minimal implementation of a ctypes structure that represents
    a singe voxel in the Monte Carlo simulations.
    All the fields defined by this implementation are required.
    '''

    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClVoxel(cltypes.Structure):
            '''
            Ctypes structure that represents a singe voxel.

            Fields
            ------
            material: mc_int_t
                Voxel material as an index into the array of materials.
            '''
            _fields_ = [
                ('material_index', T.mc_int_t),
            ]
        return ClVoxel

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        '''
        Structure that represents a single voxel in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McVoxel {',
            '	mc_int_t material_index; /* < @brief Index of the voxel material. */',
            '};'
        ))

    @staticmethod
    def np_dtype(mc: mcobject.McObject) -> np.dtype:
        '''
        Numpy data type representing a single voxel. This implementation
        uses structured numpy array that can be directly/efficiently
        mapped to OpenCL data arrays.

        Parameters
        ----------
        mc: mcobject.McObject
            Target simulator instance.

        Returns
        -------
        dtype: np.dtype
            Numpy data type representing a single voxel.
        '''
        return np.dtype([('material_index', mc.types.np_int)])


class Voxels(mcobject.McObject):
    '''
    Minimal implementation of a voxelized sample.
    All the fields defined by this implementation are
    required by the Monte Carlo simulator.
    '''

    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClVoxels(cltypes.Structure):
            '''
            Ctypes structure that is passed to the Monte Carlo kernel.

            Fields
            ------
            top_left: mc_point3f_t
                Coordinates of the top left corner of the voxelized sample volume.
            bottom_right: mc_point3f_t
                Coordinates of the bottom right corner of the voxelized sample volume.
            direction: mc_point3f_t
                Size of the voxels along the x, y and z axis.
            reflectance: shape
                Shape of the voxelized sample volume along the x, y and z
                axis.
            '''
            _fields_ = [
                ('top_left', T.mc_point3f_t),
                ('bottom_right', T.mc_point3f_t),
                ('size', T.mc_point3f_t),
                ('shape', T.mc_point3_t)
            ]
        return ClVoxels

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Structure that defines the source in the Monte Carlo simulator and
        the data type that represents a single voxel.
        '''
        voxel_cl = self._voxel_instance.fetch_cl_declaration(mc)
        if voxel_cl is None:
            voxel_cl = ''
        return '\n'.join((
            voxel_cl,
            '',
            'struct MC_STRUCT_ATTRIBUTES McVoxelConfig {',
            '	mc_point3f_t top_left; /**< @brief Coordinates of the top left corner. */',
            '	mc_point3f_t bottom_right; /**< @brief Coordinates of the bottom right corner. */',
            '	mc_point3f_t size;     /**< @brief Size of a voxel. */',
            '	mc_point3_t shape;     /**< @brief Shape of the voxelized space. */',
            '};'
        ))

    def cl_options(self, mc: mcobject.McObject) -> str:
        '''
        OpenCL options of the scattering phase function.
        '''
        return self._voxel_instance.fetch_cl_options(mc)

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        '''
        OpenCL options of the scattering phase function.
        '''
        voxel_cl  = self._voxel_instance.fetch_cl_implementation(mc)
        if voxel_cl is None:
            voxel_cl = ''
        return '\n'.join((
            voxel_cl,
            '',
            'void dbg_print_voxel_cfg(__constant const McVoxelConfig *cfg){',
            '	dbg_print("McVoxelConfig:");',
            '	dbg_print_point3f(INDENT "top_left:", &cfg->top_left);',
            '	dbg_print_point3f(INDENT "bottom_right:", &cfg->bottom_right);',
            '	dbg_print_point3f(INDENT "size:", &cfg->size);',
            '	dbg_print_point3(INDENT "shape:", &cfg->shape);',
            '};',
        ))

    def __init__(self, xaxis: Axis, yaxis: Axis, zaxis: Axis,
                 voxel: Type[Voxel] = None):
        '''
        Object that manages the voxelized sample volume. Each voxel item
        is of integer type that represents the index of material filling
        the voxel volume.

        Parameters
        ----------
        xaxis, yaxis, zaxis: Axis
            Objects that specify the voxelized grid along the three orthogonal
            coordinate axis.
            Note that the Axis constructor takes the outer boundaries of the
            range and not of the voxel centers. Use Axis(-20e-3, 20e-3, 40) to
            have 40 voxels from -20 to 20 mm (voxel size is 1 mm).
            It is recommended to start the zaxis at z=0.

        voxel: type
            Data type representing a single voxel

        Note
        ----
        Initially the material_index of all voxels is set to 0, referencing the
        material of the surrounding medium. The data array indexing follows
        C scheme, i.e [z, y, x]. Use the data method to get the numpy voxel
        array.
        '''
        if xaxis.logscale or yaxis.logscale or zaxis.logscale:
            raise ValueError('Logarithmic scale is not supported!')

        self._grid = None
        self._x_axis = xaxis
        self._y_axis = yaxis
        self._z_axis = zaxis

        if voxel is None:
            voxel = Voxel

        self._voxel_size = (
            abs(self._x_axis.step),
            abs(self._y_axis.step),
            abs(self._z_axis.step)
        )
        self._top_left = (
            min(self._x_axis.span),
            min(self._y_axis.span),
            min(self._z_axis.span)
        )
        self._bottom_right = (
            max(self._x_axis.span),
            max(self._y_axis.span),
            max(self._z_axis.span)
        )

        self._data = None
        self._voxel_type = voxel
        self._voxel_instance = self._voxel_type()
        self._update = True

    def _create_data(self, dtype) -> np.ndarray:
        return np.zeros(
            (self._z_axis.n, self._y_axis.n, self._x_axis.n),
            dtype=dtype
        )

    def __getitem__(self, item):
        return self._data[item]

    def __setitem__(self, item, value):
        self._data[item] = value

    def index(self, position: Tuple[float, float, float]) -> \
            Tuple[int, int, int]:
        '''
        Returns index of the voxel at the given position.

        Parameters
        ----------
        position: (float, float, float)
            Position as a tuple of three coordinates (x, y, z).

        Returns
        -------
        index: (int, int, int)
            A tuple of indices into the voxel array.
        '''
        indx = int((position[0] - self.top_left[0])//self._voxel_size[0])
        indy = int((position[1] - self.top_left[1])//self._voxel_size[1])
        indz = int((position[2] - self.top_left[2])//self._voxel_size[2])

        return (indz, indy, indx)

    def center(self, index: Tuple[int, int, int]) -> \
            Tuple[float, float, float]:
        '''
        Returns the center of the voxel at the given index.

        Parameters
        ----------
        index: (int, int, int)
            Voxel index as a tuple of three indices (ind_x, ind_y, ind_z).

        Returns
        -------
        position: (float, float, float)
            A tuple of indices into the voxel array.
        '''
        x = self.top_left[0] + (index[0] + 0.5)*self._voxel_size[0]
        y = self.top_left[1] + (index[1] + 0.5)*self._voxel_size[1]
        z = self.top_left[2] + (index[2] + 0.5)*self._voxel_size[2]

        return (x, y, z)

    def isvalid(self, index: Tuple[int, int, int]) -> bool:
        '''
        Returns True if the given index is valid for the voxel array.

        Parameters
        ----------
        index: (int, int, int)
            Index into the voxel array.

        Returns
        -------
        valid: bool
            True if index is valid.
        '''
        return (self.shape[0] > index[0] >= 0) and \
               (self.shape[1] > index[1] >= 0) and \
               (self.shape[2] > index[2] >= 0)

    def contains(self, position: Tuple[float, float, float]) -> bool:
        '''
        Returns True if the position is within the voxel array boundaries. Note
        that the lower bound is included but the upper bound is not.

        Parameters
        ----------
        position: (float, float, float)
            Position as a tuple of three coordinates (x, y, z).

        Returns
        -------
        contains: bool
            Returns True if the given position lies within the voxel array
            boundaries.
        '''
        return self.xaxis.start <= position[0] < self.xaxis.stop and \
               self.yaxis.start <= position[1] < self.yaxis.stop and \
               self.zaxis.start <= position[2] < self.zaxis.stop

    def intersect(self, pos:Tuple[float, float, float],
                  dir:Tuple[float, float, float]) -> \
                      Tuple[Tuple[float, float, float] or None, \
                            Tuple[float, float, float] or None]:
        '''
        Intersect the voxelized box with the ray and return the intersection.
        Intersection is considered to exist only if the distance to
        the intersection is positive.

        Parameters
        ----------
        pos: (float, float, float)
            Ray origin.
        dir: (float, float, float)
            Ray direction.

        Returns
        -------
        intersection: (float, float, float) or None
            Returns intersection if one exists, else None.
        normal: (float, float, float) or None
            Surface normal of the voxelized sample box that points in the
            propagation direction.
        '''
        invdir = [np.inf, np.inf, np.inf]
        for ax in (0, 1, 2):
            if dir[ax] != 0.0:
                invdir[ax] = 1.0/dir[ax]

        xmin, xmax = self.xaxis.edges[[0, -1]]
        ymin, ymax = self.yaxis.edges[[0, -1]]
        zmin, zmax = self.zaxis.edges[[0, -1]]

        t1 = (xmin - pos[0])*invdir[0]
        t2 = (xmax - pos[0])*invdir[0]
        t3 = (ymin - pos[1])*invdir[1]
        t4 = (ymax - pos[1])*invdir[1]
        t5 = (zmin - pos[2])*invdir[2]
        t6 = (zmax - pos[2])*invdir[2]

        txmin, txmax = min(t1, t2), max(t1, t2)
        tymin, tymax = min(t3, t4), max(t3, t4)
        tzmin, tzmax = min(t5, t6), max(t5, t6)

        tmin = max(txmin, tymin, tzmin)
        tmax = min(txmax, tymax, tzmax)

        if tmax < tmin or tmax < 0.0:
            return None, None

        if tmin >= 0.0:
            t = tmin
        else:
            t = tmax

        intersection = (pos[0] + t*dir[0], pos[1] + t*dir[1], pos[2] + t*dir[2])
        normal = [0.0, 0.0, 0.0]
        normal[0] = int(txmin >= tymin and txmin >= tzmin)*np.sign(dir[0])
        normal[1] = int(not normal[0] and tymin >= tzmin)*np.sign(dir[1])
        normal[2] = int(not (normal[0] or normal[1]))*np.sign(dir[2])

        return intersection, normal

    def update_required(self, clear: bool = True) -> bool:
        '''
        Returns True if the voxel data have been changed and need to be
        uploaded to the OpenCL device before running the next simulation.
        The update flag is cleared (set to False) if the value of the clear
        argument is True.

        Parameters
        ----------
        clear: bool
            If True, the update flag is cleared (set to False) after returning
            its value.

        Returns
        -------
        update_required: bool
            True if the voxel data have changed and need to be uploaded to
            tho OpenCL device before running a new simulation.
        '''
        value = self._update
        if clear:
            self._update = False
        return value

    def update(self):
        '''
        Call this method if the voxel data needs to be updated/uploaded to
        the OpenCL device.
        This function should be called if the voxel data have been changed
        after the last simulation call. 
        '''
        self._update = True

    def _get_voxel_type(self) -> Type[Voxel]:
        return self._voxel_type
    voxel_type = property(_get_voxel_type, None, None,
                          'Voxel data type used by this voxelized sample.')

    def _get_voxel_instance(self) -> Voxel:
        return self._voxel_instance
    voxel = property(_get_voxel_instance, None, None,
                     'Instance of the Voxel data type.')

    def _get_x_axis(self) -> Axis:
        return self._x_axis
    xaxis = property(_get_x_axis, None, None, 'Axis x of the voxel array.')

    def _get_x(self) -> np.ndarray:
        return self._x_axis.centers
    x = property(_get_x, None, None,
                 'Coordinates of the voxel centers along the x axis.')

    def _get_dx(self) -> np.ndarray:
        return self._voxel_size[0]
    dx = property(_get_dx, None, None, 'The size of voxels along the x axis.')

    def _get_y_axis(self) -> Axis:
        return self._y_axis
    yaxis = property(_get_y_axis, None, None, 'Axis y of the voxel array.')

    def _get_y(self) -> np.ndarray:
        return self._y_axis.centers
    y = property(_get_y, None, None,
                 'Coordinates of the voxel centers along the y axis.')

    def _get_dy(self) -> np.ndarray:
        return self._voxel_size[1]
    dy = property(_get_dy, None, None, 'The size of voxels along the y axis.')

    def _get_z_axis(self) -> Axis:
        return self._z_axis
    zaxis = property(_get_z_axis, None, None, 'Axis z of the voxel array.')

    def _get_z(self) -> np.ndarray:
        return self._z_axis.centers
    z = property(_get_z, None, None,
                 'Coordinates of the voxel centers along the z axis.')

    def _get_dz(self) -> np.ndarray:
        return self._voxel_size[2]
    dz = property(_get_dz, None, None, 'The size of voxels along the z axis.')

    def _get_voxel_size(self) -> Tuple[float, float, float]:
        return self._voxel_size
    voxel_size = property(_get_voxel_size, None, None,
                          'Size of a voxel (m) as (dx, dy, dz).')

    def _get_shape(self) -> Tuple[int, int, int]:
        return self._data.shape
    shape = property(_get_shape, None, None, 'Voxel array shape (z, y, x).')

    def _get_size(self) -> int:
        return self._data.size
    size = property(_get_size, None, None, 'Total number of voxels.')

    def _get_grid(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        return self.meshgrid()
    grid = property(_get_grid, None, None,
                    'Coordinate matrices Z, Y and X of voxel centers as '
                    'returned by the :py:meth:`Voxels.meshgrid` method.')

    def _get_top_left(self) -> Tuple[float, float, float]:
        return tuple(self._top_left)
    top_left = property(_get_top_left, None, None,
                        'Coordinates (x, y, z) of the top left corner of '
                        'the voxel array.')

    def meshgrid(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        '''
        Coordinate matrices with x, y and z coordinates of all the voxels.
        The shape of the matrices equals the shape of the voxel array.

        Returns
        -------
        Z: np.ndarray
            Z coordinates of all the voxels as a 3D array.
        Y: np.ndarray
            Y coordinates of all the voxels as a 3D array.
        X: np.ndarray
            X coordinates of all the voxels as a 3D array.

        Note
        ----
        The coordinate matrices are created on the first call.
        '''
        if self._grid is None:
            self._grid = np.meshgrid(
                self._z_axis.centers,
                self._y_axis.centers,
                self._x_axis.centers, indexing='ij')

        return self._grid

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure) \
            -> cltypes.Structure:
        '''
        Fills the ctypes structure (target) with the data required by the
        Monte Carlo simulator. See :py:meth:`Voxels.cl_type` for a detailed
        list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: mcypes.Structure
            Ctypes structure that is filled with the voxelization configuration.

        Returns
        -------
        target: mctypes.Structures
            Filled ctypes structure received as an input argument or a new
            instance if the input argument target is None.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        target.top_left.fromarray(self._top_left)
        target.bottom_right.fromarray(self._bottom_right)

        target.size.fromarray(self._voxel_size)

        target.shape.x = self._x_axis.n
        target.shape.y = self._y_axis.n
        target.shape.z = self._z_axis.n

        return target

    def data(self, mc: mcobject.McObject = None) -> np.ndarray:
        '''
        Return the voxel data as a 3D array representing the voxelized sample
        space.

        The data array should be addressed as data[z, y, x], wher x, y and z
        are the indices along the corresponding spatial axes passed to the
        constructor. The individual data fields of the selected voxels
        are addressed as data[x, y, z]['field_name'],
        e.g. data[0, :, :]['material_index'] = 1 sets the material of the
        first layer of voxels along the z axis to 1 (the second material in
        the materials list). 

        Parameters
        ----------
        mc: mcobject.McObject
            A simulator instance that is required only during the first call
            to derive the proper voxel data type and allocate the a 3D
            structured array that represents the voxels.

        Returns
        -------
        data: np.ndarray
            Voxel data as a 3D numpy array.
        '''
        if self._data is None:
            if mc is None:
                raise TypeError(
                    'Cannot allocate voxel data array without a '
                    'simulator instance!'
                )
            self._data = self._create_data(self._voxel_instance.np_dtype(mc))
        else:
            if mc is not None and self._voxel_instance.np_dtype(mc) != self._data.dtype:
                raise TypeError('The data type of voxels cannot be changed!')

        return self._data

    def _get_material_index(self):
        return self.data()['material_index']
    material = property(_get_material_index, None, None, 
                        '3D array of material indices.')


    def plot(self, axis: str ='z', autoscale: bool = True,
             show: bool = True):
        '''
        Show slices through the materials array.

        Parameters
        ----------
        scale: str
            Data scaling can be "log" for logarithmic or "lin" for linear.
        axis: str
            The axis of slicing ("x", "y" or "z").
        autoscale: bool
            Scale the color coding of individual slices to the corresponding
            range of weights. If True, the color coding changes from slice
            to slice.
        show: bool 
        '''
        from xopto.util import sliceview

        data = self.material

        ax = {'z':0, 'y':1, 'x':2}.get(axis[0], 0)
        title = 'Slice {{slice}}/{} @ {} = {{pos:.6f}}'.format(
            data.shape[ax], axis)

        fig = None

        if ax == 0:
            extent = [self._x_axis.start, self._x_axis.stop,
                      self._y_axis.start, self._y_axis.stop]
            slices = self._z_axis.centers
            xlabel, ylabel = 'x', 'y'
        elif ax == 1:
            extent = [self._x_axis.start, self._x_axis.stop,
                      self._z_axis.start, self._z_axis.stop]
            slices = self._y_axis.centers
            xlabel, ylabel = 'x', 'z'
        elif ax == 2:
            extent = [self._y_axis.start, self._y_axis.stop,
                      self._z_axis.start, self._z_axis.stop]
            slices = self._x_axis.centers
            xlabel, ylabel = 'y', 'z'

        sv = sliceview.SliceView(
            data, axis=ax, slices=slices, title=title, logscale=False,
            extent=extent, xlabel=xlabel, ylabel=ylabel, origin='lower',
            autoscale=autoscale)
        sv.fig.canvas.manager.set_window_title('Voxel Material SliceView')
        if show:
            sv.show()

    def render(self, show: bool = True):
        '''
        Show the voxelized volume in a 3D viewer.

        Parameters
        ----------
        show: bool
            If True, shows the window and starts the Qt event loop that
            will block until the window is closed.

        Returns
        -------
        viewer: slicer3d.Slicer3D
            Use the :py:meth:`~xopto.util.widgets.visualization.slicer3d.Slicer3D.exec`
            method to show the viewer and start the Qt event loop that will
            block until the window is closed.

        Note
        ----
        The 3D viewer requires PySide6 and PyQtGraph.
        '''
        from xopto.util.widgets import common
        from xopto.util.widgets.visualization import slicer3d

        app = common.prepareqt()

        data = self.material
        span = (data.min(), data.max())

        slicer = slicer3d.Slicer3D()
        slicer.setData(data, range_=span,
                       x=self.x*1e3, y=self.y*1e3, z=self.z*1e3)
        slicer.setXLabel('x (mm)')
        slicer.setYLabel('y (mm)')
        slicer.setZLabel('z (mm)')
        slicer.view3D().setStandardCameraView('isometric')
        slicer.setWindowTitle('Voxelized medium')

        if show:
            slicer.show()
            app.exec()

        return slicer

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {
            'type': 'Voxels',
            'xaxis': self._x_axis.todict(),
            'yaxis': self._y_axis.todict(),
            'zaxis': self._z_axis.todict(),
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'Voxels':
        '''
        Create a new instance of :py:class:`Voxels` from a dict.

        Parameters
        ----------
        data: dict
            A dict with instance data.

        Returns
        -------
        instance: Voxels
            A new instance of :py:class:`Voxels`.
        '''
        data_ = dict(data)
        type_name = data_.pop('type')
        if type_name != cls.__name__:
            raise TypeError('Expected data of instance Voxels '
                            'but got "{}"!'.format(type_name))
        return Voxels(
            Axis.fromdict(data_.pop('xaxis')),
            Axis.fromdict(data_.pop('yaxis')),
            Axis.fromdict(data_.pop('zaxis')),
            **data_
        )

    def __str__(self):
        return 'Voxels(xaxis={}, yaxis={}, zaxis={}, voxel={})'.format(
            self._x_axis, self._y_axis, self._z_axis, self._voxel_type)

    def __repr__(self):
        return '{:s} # id 0x{:>08X}.'.format(self.__str__(), id(self))
