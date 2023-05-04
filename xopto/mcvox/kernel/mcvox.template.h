/******************************** Begin license ********************************
 * Copyright (C) Laboratory of Imaging technologies,
 *               Faculty of Electrical Engineering,
 *               University of Ljubljana.
 *
 * This file is part of PyXOpto.
 *
 * PyXOpto is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PyXOpto is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PyXOpto. If not, see <https://www.gnu.org/licenses/>.
 ******************************** End license *********************************/

#ifndef __MCVOX_H
#define __MCVOX_H

#include "mcbase.template.h"

/*################## Start simulator configuration options ###################*/
/**
 * @addtogroup mc_simulator_options Simulator options
 * @{
 */

#if !defined(MC_USE_TOP_SURFACE_LAYOUT) || defined(__DOXYGEN__)
	/** @brief Turn on/off the use of advanced layout at the top sample surface. */
	#define MC_USE_TOP_SURFACE_LAYOUT			FALSE
#endif

#if !defined(MC_USE_BOTTOM_SURFACE_LAYOUT) || defined(__DOXYGEN__)
	/** @brief Turn on/off the use of advanced layout at the bottom sample surface. */
	#define MC_USE_BOTTOM_SURFACE_LAYOUT		FALSE
#endif

#if !defined(MC_MATERIAL_ARRAY_MEMORY) || defined(__DOXYGEN__)
	/** @brief Memory space of the materials array: must be __global or __constant. */
	#define MC_MATERIAL_ARRAY_MEMORY			__global
#endif

/**
 * @} // end @addtogroup mc_simulator_options
 */
/*################### End simulator configuration options ####################*/


/*############## Start GPU memory type used for different data ###############*/
/**
* @addtogroup mc_memory_types Memory types simulator data
* @{
*/
/** @brief Keep the voxelization configuration in the constant memory. */
#define __mc_voxelcfg_mem		__constant
/** @brief Keep the photon packet source configuration data in constant memory. */
#define __mc_source_mem			__constant
/** @brief Keep the scattering phase function parameters in constant memory. */
#define __mc_pf_mem				__constant
/** @brief Keep the surface reflectance / transmittance detector configuration data in constant memory. */
#define __mc_detector_mem		__constant
/** @brief Configure the floating-point lookup table memory type as set by configuration options. */
#define __mc_fp_lut_mem			MC_FP_LUT_ARRAY_MEMORY
/** @brief Configure the integer lookup table memory type as set by configuration options. */
#define __mc_int_lut_mem		MC_INT_LUT_ARRAY_MEMORY
/** @brief Keep fluence configuration data in constant memory. */
#define __mc_fluence_mem		__constant
/** @brief Keep the configuration data of the advanced layout at top surface of the sample in constant memory. */
#define __mc_surface_mem		__constant
/** @brief Keep the configuration data of the advanced layout at bottom surface of the sample in constant memory. */
#define __mc_trace_mem			__constant
/** @brief Keep the photon packet source configuration data in constant memory. */
#define __mc_material_mem		MC_MATERIAL_ARRAY_MEMORY
/**
 * @} // end @addtogroup mc_memory_types
 */
/*############### End GPU memory type used for different data ################*/


/*######################## Start forward declarations ########################*/
/**
* @addtogroup mc_simulator_state Simulator state
* @{
*/

/* Forward declaration of the simulator state and data object. */
struct McSim;
/** @brief Simulator state and data object. */
typedef struct McSim McSim;

/**
 * @} // end @addtogroup mc_simulator_state
 */
/*######################### End forward declarations #########################*/


/*############## Start scattering phase function declarations ################*/
/**
* @addtogroup mc_scattering_phase_functions Scattering phase function
* @{
*/

/**
* @brief Data type describing a scattering phase function.
* @{
*/
struct McPf;
/**
 * @}
 */
/** @brief Scattering phase function. */
typedef struct McPf McPf;

/* Scattering phase function structure is defined as part of materials! */

/**
 * @} // end @addtogroup mc_scattering_phase_functions
 */
/*############### End scattering phase function declarations #################*/


/*################# Start Photon packet source declarations ##################*/
/**
* @addtogroup mc_photon_packet_source Photon packet source
* @{
*/

/**
* @brief Data type describing a photon packet source.
* @{
*/
struct McSource;
/**
 * @}
 */
/** @brief Photon packet source. */
typedef struct McSource McSource;

/* Photon packet source structure definition goes here - DO NOT EDIT! */
/* START_SOURCE_DECLARATION_BLOCK */
{{ source.declaration or '' }}
/* END_SOURCE_DECLARATION_BLOCK */

/**
 * @brief Launches a single photon packet.
 *
 * @param mcsim A simulator instance.
 */
inline void mcsim_launch(McSim *mcsim);

/**
 * @} // end @addtogroup mc_photon_packet_source
 */
/*################## End Photon packet source declarations ###################*/


/*##################### Start declarations of detectors ######################*/
/**
* @addtogroup mc_surface_reflectance_transmittance Photon packet accumulator
* @{
*/

/**
* @brief Data type describing the surface reflectance / transmittance detector.
* @{
*/
struct McDetectors;
/**
 * @}
 */
/** @brief Surface reflectance / transmittance detector. */
typedef struct McDetectors McDetectors;

/* Surface reflectance / transmittance detector structure definition goes here - DO NOT EDIT! */
/* START_DETECTORS_DECLARATION_BLOCK */
{{ detectors.declaration or 'struct McDetectors{mc_int_t dummy;};' }}
/* END_DETECTORS_DECLARATION_BLOCK */

/**
 * @} // end @addtogroup mc_detectors
 */
/*###################### End declarations of detectors #######################*/


/*################# Start photon packet trace declarations ###################*/
/**
* @addtogroup mc_photon_packet_trace Photon packet trace
* @{
*/

/**
 * @addtogroup mcsim_multilayer_option_trace flags Trace flags
 * @{
 */

 /** @brief If no other flags are defined, trace functionality is disabled. */
#define MC_USE_TRACE_NONE	0
/** @brief Trace the initial state of photon packets. */
#define MC_USE_TRACE_START	1
/** @brief Trace the final state of the photon packets. */
#define MC_USE_TRACE_END	2
/** @brief Trace all states of the photon packets. */
#define MC_USE_TRACE_ALL	7

/**
 * @} // end @addtogroup mcsim_multilayer_option_trace
 */

/**
* @brief Data type for configuring event tracing.
* @{
*/
struct McTrace;
/**
 * @}
 */
/** @brief Event tracing type. */
typedef struct McTrace McTrace;

/* Photon packet trace structure definition goes here - DO NOT EDIT! */
/* START_TRACE_DECLARATION_BLOCK */
{{ trace.declaration or '' }}
/* END_TRACE_DECLARATION_BLOCK */

/**
 * @} // end @addtogroup mc_photon_packet_trace
 */
/*################## End photon packet trace declarations ####################*/


/*################# Start fluence accumulator declarations ###################*/
/**
* @addtogroup mc_fluence_accumulator Fluence accumulator
* @{
*/

/**
* @brief Data type describing fluence.
* @{
*/
struct McFluence;
/**
 * @}
 */
/** @brief Fluence type. */
typedef struct McFluence McFluence;

/* Fluence structure definition goes here - DO NOT EDIT! */
/* START_FLUENCE_DECLARATION_BLOCK */
{{ fluence.declaration or '' }}
/* END_FLUENCE_DECLARATION_BLOCK */

/**
 * @} // end @addtogroup mc_fluence_accumulator
 */
/*################## End fluence accumulator declarations ####################*/


/*############## Start advanced sample surface layout declaration ############*/
/**
* @addtogroup mc_advanced_surface_layout Advanced sample surface layout
* @{
*/

/** @brief This value should be returned by ::mcsim_top_surface_layout_handler
 *         and ::mcsim_bottom_surface_layout_handler, if the boundary requires
 *         standard handling with the provided parameters. */
#define MC_SURFACE_LAYOUT_CONTINUE -1

/** @brief Data type describing advanced layout at the sample top surface. */
/* Definition of the structure that describes the surface layouts
 * goes here - DO NOT EDIT! */
/* START_SURFACE_LAYUTS_DECLARATION_BLOCK */
{{ surface_layouts.declaration or 'typedef mc_int_t McSurfaceLayouts;' }}
/* END_SURFACE_LAYOUTS_DECLARATION_BLOCK */

/**
 * @brief Handles layout at the top (z = 0) sample surface.
 *
 * @param[in] psim Pointer to a simulator instance.
 * @param[in, out] n2 Initialized with the refractive index of the surrounding
 *                    medium. Update with a custom value of the refractive
 *                    index on return.
 * @return Must return ::MC_SURFACE_LAYOUT_CONTINUE if the surface
 *         boundary should be processed using the returned values of n2 and cc.
 *         Must return MC_REFLECTED or MC_REFRACTED if the boundary was
 *         fully processed by the call (including updating of the photon
 *         packet propagation direction and weight, and current layer index
 *         if required), the value will be returned by mcsim_boundary.
 */
inline int mcsim_top_surface_layout_handler(McSim *psim, mc_fp_t *n2);

/**
 * @brief Handles advanced layout at the bottom sample surface (z = sample_thickness).
 *
 * @param[in] psim Pointer to a simulator instance.
 * @param[in, out] n2 Initialized with the refractive index of the surrounding
 *                    medium. Update with a custom value of the refractive
 *                    index on return.
 * @return Must return ::MC_SURFACE_LAYOUT_CONTINUE if the surface
 *         boundary should be processed using the returned values of n2 and cc.
 *         Must return MC_REFLECTED or MC_REFRACTED if the boundary was
 *         fully processed by the call (including updating of the photon
 *         packet propagation direction and weight, and current layer index
 *         if required), the value will be returned by mcsim_boundary.
 */
inline int mcsim_bottom_surface_layout_handler(McSim *psim, mc_fp_t *n2);

/**
 * @} // end @addtogroup mc_advanced_surface_layout
 */
/*############### End advanced sample surface layout declaration #############*/


/*######################### Start simulator materials ########################*/
/**
* @addtogroup mc_materials Materials
* @{
*/

/**
 * @brief Data type describing a single sample material.
 * @note The members of this object are constant and do not change during the simulation.
 */
struct McMaterial;

/** @brief A material data type. */
typedef struct McMaterial McMaterial;

/* Material and scattering phase function structure and API definitions go here
   - DO NOT EDIT! */
/* START_MATERIALS_DECLARATION_BLOCK */
{{ materials.declaration }}
/* END_MATERIALS_DECLARATION_BLOCK */

/**
* @} // mc_materials Materials
*/
/*########################## End simulator materials #########################*/


/*#################### Start simulation state and data #######################*/
/**
* @addtogroup mcsim_voxel Sample volume voxelization
* @{
*/

/**
 * @brief Voxel array configuration
 * @{
*/

/**
* @brief Data type describing a single sample voxel.
* @{
*/
struct McVoxel;
/**
 * @}
 */
/** @brief Sample voxel. */
typedef struct McVoxel McVoxel;

/**
* @brief Data type describing the sample voxelization configuration.
* @{
*/
struct McVoxelConfig;
/**
 * @}
 */
/** @brief Sample voxelization configuration. */
typedef struct McVoxelConfig McVoxelConfig;

/* Voxel structure definition goes here - DO NOT EDIT! */
/* START_VOXEL_DECLARATION_BLOCK*/
{{ voxels.declaration }}
/* END_VOXEL_DECLARATION_BLOCK*/

/**
* @} // Sample volume voxelization
*/
/*##################### End simulation state and data ########################*/


/*#################### Start simulation state and data #######################*/
/**
* @addtogroup mcvox_simulator_core Simulator core
* @{
*/

/**
 * @brief Simulation state
 * @{
*/
struct McSimState{
	mc_point3f_t position;		/**< @brief Photon packet position. */
	mc_point3f_t direction;		/**< @brief Photon packet propagation direction. */
	mc_point3_t voxel_index;		/**< @brief Current voxel index. */
	mc_size_t voxel_material_index;	/**< @brief Index of the current voxel material. */
	uint64_t rngStateX;			/**< Random generator state X (changes on each mcsim_random call). */
	uint32_t rngStateA;			/**< Random generator state A (not changed by mcsim_random call). */
	mc_fp_t weight;				/**< @brief Photon packet weight. */
	mc_cnt_t photon_index;		/**< @brief Absolute photon packet index. */
	#if MC_USE_EVENTS || defined(__DOXYGEN__)
	mc_uint_t event_flags;		/**< @brief All events that the packet underwent during this step. */
	#endif
	#if MC_TRACK_OPTICAL_PATHLENGTH || defined(__DOXYGEN__)
	mc_fp_t optical_pathlength;		/**< Optical pathlength traveled by the photon packet. */
	#endif
	#if MC_USE_TRACE || defined(__DOXYGEN__)
	mc_uint_t trace_count;		/**< @brief Number of logged trace events since packet launch. */
	#endif
	#if MC_USE_FLUENCE && MC_USE_FLUENCE_CACHE
	mc_accucache_t fluence_cache;	/**< @brief Fluence cache object. */
	#endif
};
 /** @} */
typedef struct McSimState McSimState;

/**
* @brief Data type holding the simulator state and all required data
* @details  The McSim::remainingStep is normalized and must be multiplied by the
*			inverse of the sum of material absorption and scattering coefficient
*			to obtain the remaining step size for a particular material.
* @{
*/
struct McSim{
	/** @brief Simulation state. */
	McSimState state;

	/** @brief Configuration data of the voxel array. */
	__mc_voxelcfg_mem McVoxelConfig const *voxel_cfg;
	/** @brief Voxels array - default type represents indices into the materials array. */
	__global McVoxel const *voxels;

	/** @brief Number of materials used in the sample. */
	mc_size_t num_materials;
	/** @brief A list of McMaterial objects. */
	__mc_material_mem McMaterial const *materials;

	/** @brief Photon packet source object. */
	__mc_source_mem McSource const *source;

	#if MC_USE_SURFACE_LAYOUTS || defined(__DOXYGEN__)
		__mc_surface_mem McSurfaceLayouts const *surface_layouts;	/**< Advanced layout at the sample top and bottom surfaces. */
	#endif

	#if MC_USE_FP_LUT || defined(__DOXYGEN__)
		__mc_fp_lut_mem mc_fp_t const *fp_lut_array;	/**< Floating-point lookup table(s) data. */
	#endif

	#if MC_USE_TRACE || defined(__DOXYGEN__)
		__mc_trace_mem McTrace const *trace;	/**< @brief Trace configuration data. */;
	#endif

	#if MC_USE_FLUENCE || defined(__DOXYGEN__)
		__mc_fluence_mem McFluence const *fluence;	/**< @brief Fluence array strides (dimensions) as [nx, ny, nz]. */
	#endif

	#if MC_USE_DETECTORS || defined(__DOXYGEN__)
		__mc_detector_mem McDetectors const *detectors;	/**< @brief Reflectance/transmittance detector configuration data. */
	#endif

	__global mc_int_t *integer_buffer; 		/**< @brief Common integer buffer. */
	__global mc_fp_t *float_buffer; 		/**< @brief Common floating-point buffer. */
	__global mc_accu_t *accumulator_buffer; /**< @brief Common accumulator buffer. */
};
/**
 * @}
 */

/**
 * @brief Evaluates to the photon packet index.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_packet_index(psim) ((psim)->state.photon_index)

/**
 * @brief Evaluates to the current position (::mc_point_t type) of the photon packet.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_position(psim) (&(psim)->state.position)

/**
 * @brief Evaluates to the x coordinate of the current photon packet position.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_position_x(psim) ((psim)->state.position.x)
/**
 * @brief Evaluates to the y coordinate of the current photon packet position.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_position_y(psim) ((psim)->state.position.y)
/**
 * @brief Evaluates to the z coordinate of the current photon packet position.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_position_z(psim) ((psim)->state.position.z)
/**
 * @brief Evaluates to the r squared (polar radius) coordinate of the
 *			photon packet position.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_position_r2(psim) \
	((psim)->state.position.x*(psim)->state.position.x + \
	(psim)->state.position.y*(psim)->state.position.y)
/**
 * @brief Evaluates to the r (polar radius) coordinate of the
 *			photon packet position.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_position_r(psim) mc_sqrt(mcsim_position_r2(psim))

/**
 * @brief Returns the r squared (polar radius) coordinate of the
 *        photon packet position with respect to the given origin.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in]  Coordinate x of the origin.
 * @param[in] y Coordinate y of the origin.
 * @returns r2 Squared distance from the origin.
 */
inline mc_fp_t mcsim_position_r2_ex(
	McSim const *psim, mc_fp_t x_center, mc_fp_t y_center);

/**
 * @brief Evaluates to the r (polar radius) coordinate of the
 *        photon packet position with  respect to the given origin.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in]  Coordinate x of the origin.
 * @param[in] y Coordinate y of the origin.
 */
#define mcsim_position_r_ex(psim, x, y) \
	mc_sqrt(mcsim_position_r2_ex((psim), (x), (y)))

/**
 * @brief Sets the z coordinante of the photon packet position to the specified value.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] zpos New packet position z coordinate.
 */
#define mcsim_set_position_z(psim, zpos) ((psim)->state.position.z = (zpos))

/**
 * @brief Sets the current position of the photon packet to the specified value (::mc_point_t type).
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] ppoint Pointer to the new photon packet position (::mc_point_t type).
 */
#define mcsim_set_position(psim, ppoint) ((psim)->state.position = *(ppoint))

/**
 * @brief Sets the current position of the photon packet to the specified value (::mc_point_t type).
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] posx X coordinate of the new photon packet position.
 * @param[in] posy Y coordinate of the new photon packet position.
 * @param[in] posz Z coordinate of the new photon packet position.
 */
#define mcsim_set_position_coordinates(psim, posx, posy, posz) \
	{(psim)->state.position.x = (posx); \
	(psim)->state.position.y = (posy); \
	(psim)->state.position.z = (posz);}

/**
 * @brief Adjusts the current photon packet position by moving along
 * 		the current photon packet propagation direction by d.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] d Step length along the current photon packet propagation direction.
 */
#define mcsim_move(psim, d) \
	{(psim)->state.position.x += (d)*(psim)->state.direction.x; \
	(psim)->state.position.y += (d)*(psim)->state.direction.y; \
	(psim)->state.position.z += (d)*(psim)->state.direction.z;}

/**
 * @brief Evaluates to the current propagation direction (::mc_point_t type) of the
 * photon packet.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_direction(psim) (&(psim)->state.direction)

/**
 * @brief Sets the propagation direction of the photon packet.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] pdir Pointer to the new propagation direction vector.
 */
#define mcsim_set_direction(psim, pdir) \
	{ \
		(psim)->state.direction.x = (pdir)->x; \
		(psim)->state.direction.y = (pdir)->y; \
		(psim)->state.direction.z = (pdir)->z; \
	}

/**
 * @brief Reverse the z component of the photon packet propagation direction.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_reverse_direction_z(psim) \
	{(psim)->state.direction.z = -(psim)->state.direction.z;}

/**
 * @brief Sets the propagation direction of the photon packet.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] px X component of the new propagation direction vector.
 * @param[in] py Y component of the new propagation direction vector.
 * @param[in] pz Z component of the new propagation direction vector.
 */
#define mcsim_set_direction_coordinates(psim, px, py, pz) \
	{(psim)->state.direction.x = px; \
		(psim)->state.direction.y = py; \
		(psim)->state.direction.z = pz;}

/**
 * @brief Evaluates to the x component of the photon packet propagation direction.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_direction_x(psim) ((psim)->state.direction.x)
/**
 * @brief Evaluates to the y component of the photon packet propagation direction.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_direction_y(psim) ((psim)->state.direction.y)
/**
 * @brief Evaluates to the z component of the photon packet propagation direction.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_direction_z(psim) ((psim)->state.direction.z)

/**
 * @brief Evaluates to the x component of the top left corner of the voxel array.
 */
#define mcsim_top_left_x(psim) ((psim)->voxel_cfg->top_left.x)
/**
 * @brief Evaluates to the y component of the top left corner of the voxel array.
 */
#define mcsim_top_left_y(psim) ((psim)->voxel_cfg->top_left.y)
/**
 * @brief Evaluates to the z component of the top left corner of the voxel array.
 */
#define mcsim_top_left_z(psim) ((psim)->voxel_cfg->top_left.z)
/**
 * @brief Evaluates to the top left corner of the voxel array (pointer to mc_point3f_t).
 */
#define mcsim_top_left(psim) (&((psim)->voxel_cfg->top_left))

/**
 * @brief Evaluates to the size of the voxel array.
 */
#define mcsim_size(psim) \
	( \
		((mc_size_t)(psim)->voxel_cfg->shape.x)* \
		((mc_size_t)(psim)->voxel_cfg->shape.y)* \
		((mc_size_t)(psim)->voxel_cfg->shape.z) \
	)

/**
 * @brief Evaluates to the size of the voxel array along the x axis.
 */
#define mcsim_shape_x(psim) ((psim)->voxel_cfg->shape.x)
/**
 * @brief Evaluates to the size of the voxel array along the y axis.
 */
#define mcsim_shape_y(psim) ((psim)->voxel_cfg->shape.y)
/**
 * @brief Evaluates to the size of the voxel array along the z axis.
 */
#define mcsim_shape_z(psim) ((psim)->voxel_cfg->shape.z)
/**
 * @brief Evaluates to the shape of the voxel array (pointer to mc_point3_t).
 */
#define mcsim_shape(psim) (&((psim)->voxel_cfg->shape))

/**
 * @brief Evaluates to the x component of the bottom right corner of the voxel array.
 */
#define mcsim_bottom_right_x(psim) ((psim)->voxel_cfg->bottom_right.x)
/**
 * @brief Evaluates to the y component of the bottom right corner of the voxel array.
 */
#define mcsim_bottom_right_y(psim) ((psim)->voxel_cfg->bottom_right.y)
/**
 * @brief Evaluates to the z component of the bottom right corner of the voxel array.
 */
#define mcsim_bottom_right_z(psim) ((psim)->voxel_cfg->bottom_right.z)
/**
 * @brief Evaluates to the bottom right corner of the voxel array (pointer to mc_point3f_t).
 */
#define mcsim_bottom_right(psim) (&((psim)->voxel_cfg->bottom_right))

/**
 * @brief Set the current voxel index by components and update the index of material.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] x Component x of the index.
 * @param[in] y Component y of the index.
 * @param[in] z Component z of the index.
 */
#define mcsim_set_voxel_index_components_and_material(psim, x, y, z) \
	{\
		(psim)->state.voxel_index.x = x;\
		(psim)->state.voxel_index.y = y;\
		(psim)->state.voxel_index.z = z;\
		(psim)->state.voxel_material_index = mcsim_voxel_material_index( \
			psim, &(psim)->state.voxel_index); \
	}

/**
 * @brief Set the current voxel index by components and keep the current index of material.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] x Component x of the index.
 * @param[in] y Component y of the index.
 * @param[in] z Component z of the index.
 */
#define mcsim_set_voxel_index_components(psim, x, y, z) \
	{\
		(psim)->state.voxel_index.x = x;\
		(psim)->state.voxel_index.y = y;\
		(psim)->state.voxel_index.z = z;\
	}

/**
 * @brief Set the current voxel index and update the current material index.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] pindex Pointer to an instance mc_point_t.
 */
#define mcsim_set_voxel_index_and_material(psim, pindex) \
	{\
		(psim)->state.voxel_index.x = (pindex)->x;\
		(psim)->state.voxel_index.y = (pindex)->y;\
		(psim)->state.voxel_index.z = (pindex)->z;\
		(psim)->state.voxel_material_index = mcsim_voxel_material_index( \
			psim, &(psim)->state.voxel_index); \
	}

/**
 * @brief Set the current voxel material index.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] material_index Index of the material.
 */
#define mcsim_set_voxel_material_index(psim, material_index) \
	{\
		(psim)->state.voxel_material_index = (material_index); \
	}

/**
 * @brief Set the current voxel index and keep the current material index.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] pindex Pointer to an instance mc_point_t.
 */
#define mcsim_set_voxel_index(psim, pindex) \
	{\
		(psim)->state.voxel_index.x = (pindex)->x;\
		(psim)->state.voxel_index.y = (pindex)->y;\
		(psim)->state.voxel_index.z = (pindex)->z;\
	}

/**
 * @brief Evaluates to a pointer to McMaterial at the specified material index.
 *
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] index Material index.
 */
#define mcsim_material(psim, index) (&((psim)->materials[index]))

/**
 * @brief Evaluates to the material of the surrounding medium.
 *
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_surrounding_material(psim) (mcsim_material(psim, 0))

/**
 * @brief Evaluates to the index of the surrounding medium material.
 *
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_surrounding_material_index(psim) 0

/**
 * @brief Evaluates to a flat voxel index.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] pindex Pointer to a voxel index (mc_point3_t value).
 */
#define mcsim_flat_voxel_index(psim, pindex) \
	(((pindex)->z*mcsim_shape_y(psim) + (pindex)->y)*mcsim_shape_x(psim) + (pindex)->x)

/**
 * @brief Evaluates to voxel array.
 *
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_voxels(psim) ((psim)->voxels)

/**
 * @brief Evaluates to the index of the material for the given voxel index.
 * @note Index must be valid, no check is performed!
 *
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] pindex Pointer to a voxel index (mc_point3_t).
 */
#define mcsim_voxel_material_index(psim, pindex) \
	((psim)->voxels[mcsim_flat_voxel_index(psim, pindex)].material_index)

/**
 * @brief evaluates to the material at the specified voxel index.
 * @note The specified index is not checked and must be valid!
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] pindex Pointer to a voxel index.
 */
#define mcsim_voxel_material(psim, pindex) \
	mcsim_material(psim, mcsim_voxel_material_index(psim, pindex))

/**
 * @brief Evaluates to the material index of the current voxel.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_current_voxel_material_index(psim) \
	((psim)->state.voxel_material_index)

/**
 * @brief Evaluates to the material of the current voxel.
 * @note The current voxel index is not checked and must be valid!
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] flat_index Flat voxel index.
 */
#define mcsim_current_voxel_material(psim) \
	mcsim_material(psim, (psim)->state.voxel_material_index)

/**
 * @brief Evaluates to the current photon packet weight (from [0.0, 1.0]).
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_weight(psim) ((psim)->state.weight)

/**
 * @brief Convert floating-point weight to integer weight.
 * @param[in] psim Pointer to a simulator instance.
 */
#define weight_to_int(weight) ((weight)*MC_INT_ACCUMULATOR_K + FP_0p5)

/**
 * @brief Evaluates to the current photon packet weight converted to integer.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_int_weight(psim) weight_to_int((psim)->state.weight)

/**
 * @brief Evaluates to the current photon packet weight converted to an integer.
 * @param[in] scale Additional multiplicative factor applied to the photon packet weight.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_int_weight_ex(psim, scale) \
	((psim)->state.weight*scale*MC_INT_ACCUMULATOR_K + FP_0p5)

/**
 * @brief Sets the current photon packet weight (from [0.0, 1.0]) to the specified value.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] w the new photon packet weight from [0.0, 1.0].
 */
#define mcsim_set_weight(psim, w) ((psim)->state.weight = (w))

/**
 * @brief Adjusts the current photon packet weight (from [0.0, 1.0]) by
 *			subtracting the specified value.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] delta Subtracted from the current photon packet weight.
 * @note The resulting photon packet weight is NOT adjusted to valid range [0.0, 1.0].
 */
#define mcsim_adjust_weight(psim, delta) ((psim)->state.weight -= (delta))

/**
* @brief Evaluates to the x component of the current voxel index.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_voxel_index_x(psim) ((psim)->state.voxel_index.x)

/**
* @brief Evaluates to the y component of the current voxel index.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_voxel_index_y(psim) ((psim)->state.voxel_index.y)

/**
* @brief Evaluates to the z component of the current voxel index.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_voxel_index_z(psim) ((psim)->state.voxel_index.z)

/**
* @brief Evaluates to a pointer to the current voxel index.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_voxel_index(psim) (&((psim)->state.voxel_index))

/**
* @brief Evaluates to the x component of the current voxel index.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_voxel_size_x(psim) ((psim)->voxel_cfg->size.x)

/**
* @brief Evaluates to the y component of the current voxel index.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_voxel_size_y(psim) ((psim)->voxel_cfg->size.y)

/**
* @brief Evaluates to the z component of the current voxel index.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_voxel_size_z(psim) ((psim)->voxel_cfg->size.z)

/**
* @brief Evaluates to the x coordinate of the voxel center.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_voxel_center_x(psim) \
	(mcsim_top_left_x(psim) + \
	(mcsim_voxel_index_x(psim) + FP_0p5)*mcsim_voxel_size_x(psim))

/**
* @brief Evaluates to the y coordinate of the voxel center.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_voxel_center_y(psim) \
	(mcsim_top_left_y(psim) + \
	(mcsim_voxel_index_y(psim) + FP_0p5)*mcsim_voxel_size_y(psim))

/**
* @brief Evaluates to the z coordinate of the voxel center.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_voxel_center_z(psim) \
	(mcsim_top_left_z(psim) + \
	(mcsim_voxel_index_z(psim) + FP_0p5)*mcsim_voxel_size_z(psim))

/**
* @brief Evaluates to the coordinate of the voxel center.
* @param[in] psim Pointer to a simulator instance.
* @param[out] point3 Pointer to a mc_point3_t instance that will be filled with coordinates.
*/
#define mcsim_voxel_center(psim, pcenter) \
	(pcenter->x = mcsim_voxel_center_x(psim) \
		pcenter->y = mcsim_voxel_center_y(psim) \
		pcenter->z = mcsim_voxel_center_z(psim))

/**
* @brief Evaluates to the total number of voxels.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_voxel_count(psim) \
	( \
		(psim)->voxel_cfg->shape.x* \
		(psim)->voxel_cfg->shape.y* \
		(psim)->voxel_cfg->shape.z \
	)

/**
 * @brief evaluates to nonzero if the voxel index is valid
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] pindex Pointer to a mc_point_t value.
 */
#define mcsim_is_valid_voxel_index(psim, pindex) \
	( \
		((pindex)->x >= 0 && (pindex)->x < (psim)->voxel_cfg->shape.x) && \
		((pindex)->y >= 0 && (pindex)->y < (psim)->voxel_cfg->shape.y) && \
		((pindex)->z >= 0 && (pindex)->z < (psim)->voxel_cfg->shape.z) \
	)

/**
 * @brief Evaluates to nonzero if the photon packet has escaped the sample.
 *
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_packet_escaped_sample(psim) \
	(!mcsim_is_valid_voxel_index(psim, mcsim_voxel_index(psim)))

/**
* @brief Evaluates to a pointer to the current voxel material scattering phase
*		function.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_current_voxel_pf(psim) \
	(&(mcsim_current_voxel_material(psim)->pf))

/**
 * Portable alias for mcsim_current_voxel_pf
 */
#define mcsim_current_pf(psim) mcsim_current_voxel_pf(psim)

/**
 * @brief Evaluates to a pointer to the photon packet source object (::McSource).
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_source(psim) ((psim)->source)

/**
 * @brief Evaluates to the total optical pathlength of the photon packet.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_optical_pathlength(psim) ((psim)->state.optical_pathlength)

/**
 * @brief Adds a value to the total optical pathlength of the photon packet.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] opl  Optical pathlength to add to the total pathlength.
 */
#define mcsim_optical_pathlength_add(psim, opl) \
	((psim)->state.optical_pathlength += (opl))

/**
* @brief Evaluates to the number of trace events since packet launch.
* @param[in] psim Pointer to a simulator instance.
* @note Depending on the available trace length, not all events might be logged.
*/
#define mcSimTraceCount(psim) ((psim)->state.traceCount)

#if MC_USE_FP_LUT || defined(__DOXYGEN)
    /**
    * @brief Evaluates to the main array of floating-point lookup table data.
    * @param[in] psim Pointer to a simulator instance.
    */
    #define mcsim_fp_lut_array(psim) ((psim)->fp_lut_array)

    /**
     * @brief In this implementation, the scattering phase function uses the main
     *        lookup table array.
     */
    #define mcsim_pf_lut_array mcsim_fp_lut_array

#endif

/**
* @brief Evaluates to the scattering phase function lookup table array
*        that starts at the specified offset.
* @param[in] psim Pointer to a simulator instance.
* @param[in] offset Offset from the start of the lookup table array.
*/
#define mcsim_fp_lut_array_ex(psim, offset) ((psim)->fp_lut_array + (offset))

/**
 * @brief Evaluates to a pointer to the detectors configuration structure.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_detectors(psim) ((psim)->detectors)

/**
 * @brief Evaluates to a pointer to the common accumulator buffer.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_accumulator_buffer(psim) ((psim)->accumulator_buffer)

/**
 * @brief Evaluates to a pointer to the common accumulator buffer
 *        at the specified offset.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] offset Offset from the beginning of the buffer expressed
 *                   in number of buffer items.
 */
#define mcsim_accumulator_buffer_ex(psim, offset) \
	((psim)->accumulator_buffer + (offset))

/**
 * @brief Evaluates to a pointer to the common floating-point buffer.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_float_buffer(psim) ((psim)->float_buffer)

/**
 * @brief Evaluates to a pointer to the common floating-point buffer
 *        at the specified offset.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] offset Offset from the beginning of the buffer expressed
 *                   in number of buffer items.
 */
#define mcsim_float_buffer_ex(psim, offset) \
	((psim)->float_buffer + offset)

/**
 * @brief Evaluates to a pointer to the common integer buffer.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_integer_buffer(psim) ((psim)->integer_buffer)

/**
 * @brief Evaluates to a pointer to the common integer buffer at the specified
 *        offset.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] offset Offset from the beginning of the buffer expressed
 *                   in number of buffer items.
 */
#define mcsim_integer_buffer_ex(psim, offset) \
	((psim)->integer_buffer + offset)

/**
* @brief Evaluates to a pointer to the surface layouts.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_surface_layouts(psim) ((psim)->surface_layouts)

/**
* @brief Evaluates to a pointer to the fluence configuration structure.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_fluence(psim) ((psim)->fluence)

#if MC_USE_FLUENCE || defined(__DOXYGEN)
/**
 * @brief Deposit the given weight to the fluence accumulator. Takes care
 *        of the different function signatures in case of the weight deposition
 *        and fluence rate mode.
 *
 * @param[in] psim     Simulator instance.
 * @param[in] pos      Deposit the weight at this position.
 * @param[in] deposit  Weight to deposit.
 */
inline void mcsim_fluence_deposit_weight(
	McSim *psim, mc_point3f_t const *pos, mc_fp_t deposit);

#if MC_FLUENCE_MODE_RATE
/**
 * @brief Deposit the given weight to the fluence accumulator that operates
 *        in fluence rate mode.
 *
 * @param[in] mcsim   Simulator instance.
 * @param[in] pos     Deposit the weight at this position.
 * @param[in] weight  Weight to deposit.
 * @param[in] mua     Absorption coefficient of the medium that is required to
 *                    to compute the fluence rate.
 *
 * @note The source code of this function is implemented in related python modules.
 */
inline void mcsim_fluence_deposit_at(
	McSim *mcsim, mc_point3f_t const *pos, mc_fp_t weight, mc_fp_t mua);
#else
/**
 * @brief Deposit the given weight to the fluence accumulator that operates
 *        in energy deposition mode.
 *
 * @param[in] mcsim   Simulator instance.
 * @param[in] pos     Deposit the weight at this position.
 * @param[in] weight  Weight to deposit.
 *
 * @note The source code of this function is implemented in related python modules.
 */
inline void mcsim_fluence_deposit_at(
	McSim *mcsim, mc_point3f_t const *pos, mc_fp_t weight);
#endif

#if MC_USE_FLUENCE_CACHE
	/**
	 * @brief Evaluates to a pointer to the fluence cache.
	 *
	 * @return   Pointer to a fluence cache instance.
	 */
	#define mcsim_fluence_cache(psim) (&((psim)->state.fluence_cache))
#endif

/**
 * @brief Low-level deposition function that can use an intermediate cache if
 *        configured so through the ::MC_USE_FLUENCE_CACHE option.
 *
 * @param psim     Simulator instance.
 * @param offset   Deposition address/offset.
 * @param weight   Weight to deposit.
 */
inline void mcsim_fluence_weight_deposit_ll(
	McSim *psim, size_t offset, uint32_t weight);

#endif /* MC_USE_FLUENCE */

/**
* @brief Evaluates to a pointer to the trace configuration structure.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_trace(psim) ((psim)->trace)

/**
* @brief Evaluates to the number of trace events since the launch of the photon packet.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_trace_count(psim) ((psim)->state.trace_count)

/**
 * @brief Logs one trace event.
 *
 * @param[in] psim        Pointer to a simulator instance.
 * @param[in] event_count Number of already traced events for this photon packet.
 *
 * @return                Returns nonzero if the event counter needs to be incremented.
 *
 * @note The source code of this function is implemented in related python modules.
 */
inline int mcsim_trace_event(McSim *psim, mc_uint_t event_count);

/**
 * @brief Finalizes the trace.
 *
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] event_count The number of traced events for this photon packet.
 *
 * @note The source code of this function is implemented in related python modules.
 */
inline void mcsim_trace_complete(McSim *psim, mc_uint_t event_count);

/**
 * @} // end @addtogroup mcvox_simulator_core
 */
/*############## End Monte Carlo simulator state declarations ################*/


/*######################## Start geometry handling ###########################*/
/**
 * @addtogroup mc_voxel_tools Tools for the voxelized geometry
 * @{
 */

/**
 * @brief Checks if the sample box contains the given point.
 *
 * @param psim Pointer to a simulator instance.
 * @param pos Point to check.
 * @return Returns nonzero if the given point lies inside the sample volume.
 */
#define mcsim_sample_box_contains(psim, pt) \
	( \
		(pt)->x >= mcsim_top_left_x(psim) && \
		(pt)->x < mcsim_bottom_right_x(psim) && \
		(pt)->y >= mcsim_top_left_y(psim) && \
		(pt)->y < mcsim_bottom_right_y(psim) && \
		(pt)->z >= mcsim_top_left_z(psim) && \
		(pt)->z < mcsim_bottom_right_z(psim) \
	)

/**
 * @brief Intersect the sample box and return the intersection and surface
 *        normal that points in the propagation direction.
 *
 * @param psim Pointer to a simulator instance.
 * @param pos Origin of the intersecting ray.
 * @param dir Propagation direction of the intersecting ray.
 * @param intersection Computed intersection
 * @param normal Computed surface normal at the intersection.
 *
 * @return Nonzero if an intersection exists.
 */
inline int mcsim_sample_box_intersect(
		McSim *psim, mc_point3f_t const *pos, mc_point3f_t const *dir,
		mc_point3f_t *intersection, mc_point3f_t *normal);

/**
 * @brief Computes distance to intersection along the propagation direction
 * 			with all the boundaries of the current voxel.
 *
 * @param[in] psim Pointer to a simulator instance.
 * @param[out] distance Filled with the computed distances on return;
 *
 * @return Distance to the nearest boundary along the propagation direction.
 */
inline mc_fp_t mcsim_intersect(McSim *psim, mc_point3f_t *distance);

/**
 * @brief Compute the voxel index from the given position. No check are made
 *        if the position is within the sample box.
 *
 * param[in]  psim Pointer to a simulator instance.
 * param[in]  pos  Position.
 * param[out] ind  Voxel index.
 */
inline void mcsim_position_to_voxel_index(
	McSim *psim, mc_point3f_t const *pos, mc_point3_t *ind);

/**
 * @brief Compute the voxel index from the given position. Clip the indices
 *        to the valid range.
 *
 * param[in]  psim Pointer to a simulator instance.
 * param[in]  pos  Position.
 * param[out] ind  Voxel index.
 */
inline void mcsim_position_to_voxel_index_safe(
	McSim *psim, mc_point3f_t const *pos, mc_point3_t *ind);

/**
 * @brief Update the current voxel index from the current position of the
 *			photon packet. Note that the lower bounds of the voxel boundaries
 *			are included in the voxel, bu the upper bounds are not.
 *
 * param[in] psim Pointer to a simulator instance.
 */
inline void mcsim_set_voxel_index_from_position(McSim *psim);

/**
 * @} // Tools for the voxelized geometry
 */
/*######################### End geometry handling ############################*/


/*############### Start layer boundary handler declarations ##################*/
/**
 * @addtogroup mc_boundary_crossing Boundary crossing
 * @{
 */

/**
 * @addtogroup mc_boundary_interaction_outcome Boundary interaction outcomes
 *             returned by ::mcsim_boundary function.
 * @{
 */
 /**@brief Photon packet has been reflected from the boundary. */
#define MC_REFLECTED	0
/**@brief Photon packet has been refracted into the next layer. */
#define MC_REFRACTED	1

/**
 * @} // end @addtogroup mc_boundary_interaction_outcome
 */

/**
 * @brief Handles voxel boundary interactions (refraction/reflection).
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] distances Distances to the voxel boundaries.
 * @return Returns MC_REFLECTED if the photon packet is reflected from the
 *         boundary or MC_REFRECTED if the photon packet is refracted across
 *         the voxel boundary.
 */
inline mc_int_t mcsim_boundary(McSim *psim, const mc_point3f_t *distances);

/**
 * @} // end @addtogroup mc_boundary_crossing
 */
/*################ End layer boundary handler declarations ###################*/


/*############### Start random number generator declarations #################*/
/**
* @addtogroup mc_random Random number generator
* @{
*/

/**
 * @brief Generates a single precision random number from (0.0, 1.0) and
 *        update the generator state.
 * @note Due to precision issues
 *       the open nature of the interval is not guaranteed.
 *
 * @param[in]	psim Pointer to a simulator instance.
 * @return		A random number from (0.0, 1.0). Due to precision issues
 *				the open nature of the interval is not guaranteed.
 *
 @details George Marsaglia's Random Number Generator. A single precision
 *        floating-point number has a 23-bit mantisa, hence only integers from
 *        [0, 2^23 - 1] can be represented without loss of information.
 */
inline mc_fp_t mcsim_random_single(McSim *psim);

#if MC_USE_DOUBLE_PRECISION || defined(__DOXYGEN__)
	/**
	 * @brief Generates a double precision random number from (0.0, 1.0) and
	 *        update the generator state.
	 * @note Due to precision issues
	 *       the open nature of the interval is not guaranteed.
	 * @param[in]	psim Pointer to a simulator instance.
	 * @return		A random number from (0.0, 1.0).
	 *
	 * @details George Marsaglia's Random Number Generator. A double precision
	 *          floating-point number has a 52-bit mantisa, hence only integers from
	 *          0, 2^52 - 1] can be represented without loss of information.
	 */
	inline mc_fp_t mcsim_random_double(McSim *psim);

	/** @brief Using double precision random generator when ::MC_USE_DOUBLE_PRECISION is defined. **/
	#define mcsim_random	mcsim_random_double
#else
	/** @brief Using single precision random generator when ::MC_USE_DOUBLE_PRECISION is NOT TRUE. **/
	#define mcsim_random	mcsim_random_single
#endif

/**
 * @} // end @addtogroup mc_random
 */
/*################ End random number generator declarations ##################*/


/*################# Start scattering handling declarations ###################*/
/**
* @addtogroup mc_scattering Scattering
* @{
*/

/**
* @brief Called when the photon packet needs to be scattered.
* @param[in, out] psim Simulator satate.
* @details Function computes and updates the packet propagation direction
*          and associated states photon packet states.
*/
inline void mcsim_scatter(McSim *psim);

/**
 * @brief Implementation of the scattering phase function.
 *
 * @param[in, out] psim    Simulator satate.
 * @param[out] azimuth     Azimuth scattering angle. Is usually a uniformly
 *                         distributed number from interval [0.0, 2*pi]
 * @return Deflection angle cosine.
 *
 * @note Scattering phase functions are implemented in python modules.
 */
inline mc_fp_t mcsim_sample_pf(McSim *psim, mc_fp_t *azimuth);

/**
 * @} // end @addtogroup mc_scattering
 */
/*################## End scattering handling declarations ####################*/


/*#################### Start debug support declarations ######################*/
/**
* @addtogroup mc_debug_support Debug support
* @{
*/

#if MC_ENABLE_DEBUG || defined(__DOXYGEN__)
/**
 * @brief Print simulator status.
 * @param[in] psim Pointer to a simulator instance.
 * param[in] label Label that will be printed before the simulator status"
 */
#define dbg_print_status(psim, label) \
	dbg_print("Status: " label); \
	dbg_print_cnt(INDENT     "Packet index/id:", mcsim_packet_index(psim)); \
	dbg_print_point3f(INDENT "Position       :", mcsim_position(psim)); \
	dbg_print_point3f(INDENT "Direction      :", mcsim_direction(psim)); \
	dbg_print_point3(INDENT  "Voxel index    :", mcsim_voxel_index(psim)); \
	dbg_print_size_t(INDENT  "Material index :", mcsim_current_voxel_material_index(psim)); \
	dbg_print_float(INDENT   "Weight         :", mcsim_weight(psim));

#else

#define dbg_print_status(psim, label) ;

#endif

/**
 * @} // end @addtogroup mc_debug_support
 */
/*##################### End debug support declarations #######################*/


/*#################### Start events support declarations #####################*/

/**
* @addtogroup mc_events
* @{
*/
#if MC_USE_EVENTS || defined(__DOXYGEN__)
    /**
     * @brief Clear all the event flags.
     *
     * @param[in] psim      Simulator instance.
     */
    static inline void mcsim_event_flags_clear(McSim *psim) {
        psim->state.event_flags = 0U;
    };
    
    /**
     * @brief Add packet event flags to the event mask.
     *
     * @param[in] psim      Simulator instance.
     * @param[in] flags     Flags that will be added to the event mask. 
     */
    static inline void mcsim_event_flags_add(McSim *psim, mc_uint_t flags) {
        psim->state.event_flags |= flags;
    };
    
    /**
     * @brief Returns the event flags.
     */
    static inline unsigned int mcsim_event_flags(McSim *psim) {
        return psim->state.event_flags;
    };
#else
	/** @brief Default implementation is void. */
	#define mcsim_event_flags_clear(psim)			((void)(psim))
	/** @brief Default implementation is void. */
	#define mcsim_event_flags_add(psim, flags)		((void)(psim), (void)(flags))
	/** @brief Default implementation returns all flags set. */
	#define mcsim_event_flags(psim)                  (0xFFFFFFFFU)
#endif

/**
 * @} // end @addtogroup mc_events
 */

/*##################### End events support declarations ######################*/

#endif /* #define __MCVOX_H */
