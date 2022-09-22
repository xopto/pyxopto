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

#ifndef __MCML_H
#define __MCML_H

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

/**
 * @} // end @addtogroup mc_simulator_options
 */
/*################### End simulator configuration options ####################*/


/*############## Start GPU memory type used for different data ###############*/
/**
* @addtogroup mc_memory_types Memory types simulator data
* @{
*/

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

/* Scattering phase function structure definition goes here - DO NOT EDIT! */
/* START_PF_DECLARATION_BLOCK*/
{{ pf.declaration or '' }}
/* END_PF_DECLARATION_BLOCK */

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
 * @param[out] cc     Update with the cosine of the critical incidence angle
 *                    for the boundary transition.
 * @return Must return ::MC_SURFACE_LAYOUT_CONTINUE if the surface
 *         boundary should be processed using the returned values of n2 and cc.
 *         Must return MC_REFLECTED or MC_REFRACTED if the boundary was
 *         fully processed by the call (including updating of the photon
 *         packet propagation direction and weight, and current layer index
 *         if required), the value will be returned by mcsim_boundary.
 */
inline int mcsim_top_surface_layout_handler(
	McSim *psim, mc_fp_t *n2, mc_fp_t *cc);

/**
 * @brief Handles advanced layout at the bottom sample surface (z = sample_thickness).
 *
 * @param[in] psim Pointer to a simulator instance.
 * @param[in, out] n2 Initialized with the refractive index of the surrounding
 *                    medium. Update with a custom value of the refractive
 *                    index on return.
 * @param[out] cc     Update with the cosine of the critical incidence angle
 *                    for the boundary transition.
 * @return Must return ::MC_SURFACE_LAYOUT_CONTINUE if the surface
 *         boundary should be processed using the returned values of n2 and cc.
 *         Must return MC_REFLECTED or MC_REFRACTED if the boundary was
 *         fully processed by the call (REFLincluding updating of the photon
 *         packet propagation direction and weight, and current layer index
 *         if required), the value will be returned by mcsim_boundary.
 */
inline int mcsim_bottom_surface_layout_handler(
	McSim *psim, mc_fp_t *n2, mc_fp_t *cc);

/**
 * @} // end @addtogroup mc_advanced_surface_layout
 */
/*############### End advanced sample surface layout declaration #############*/


/*######################### Start layer declarations #########################*/
/**
* @addtogroup mc_layer Sample Layer
* @{
*/

 /**
 * @brief Data type describing a single sample layer.
 * @note The members of this object are constant and do not change during the simulation.
 * @{
 */
struct MC_STRUCT_ATTRIBUTES McLayer {
	mc_fp_t thickness;					/**< Layer thickness. */
	mc_fp_t top;						/**< Z coordinate of the layer top surface (z coordinate increases with the layer index). */
	mc_fp_t bottom;						/**< Z coordinate of the layer bottom surface (z coordinate increases with the layer index). */
	mc_fp_t n;							/**< Layer index of refraction. */
	mc_fp_t cos_critical_top;			/**< Total internal reflection angle cosine for the above layer. */
	mc_fp_t cos_critical_bottom;		/**< Total internalreflection angle cosine for the bellow layer. */
	mc_fp_t mus;						/**< Scattering coefficient. */
	mc_fp_t mua;						/**< Absorption coefficient. */
	mc_fp_t inv_mut;					/**< Reciprocal of the total attenuation coefficient. */
	mc_fp_t mua_inv_mut;				/**< Absorption coefficient multiplied by the reciprocal of the total attenuation coefficient. */
	McPf pf;							/**< Scattering phase function parameters. */
};
/**
 * @}
 */
/** @brief Data type representing a sample layer. */
typedef struct McLayer McLayer;

/**
 * @brief Evaluates to the layer thickness.
 * @param[in] player Pointer to a layer object.
 */
#define mc_layer_thickness(player) ((player)->thickness)

/**
 * @brief Evaluates to the z coordinate of the top layer surface.
 * @param[in] player Pointer to a layer object.
 */
#define mc_layer_top(player) ((player)->top)

/**
 * @brief Evaluates to the z coordinate of the bottom layer surface.
 * @param[in] player Pointer to a layer object.
 */
#define mc_layer_bottom(player) ((player)->bottom)

/**
 * @brief Evaluates to the layer refractive index.
 * @param[in] player Pointer to a layer object.
 */
 #define mc_layer_n(player) ((player)->n)

/**
 * @brief Evaluates to the critical cosine (total internal reflection)
 *        at the top layer boundary.
 * @details If the absolute cosine of the angle of incidence
 *          (with respect to z axis) is less than the critical cosine,
 @          the incident packet is reflected at the boundary.
 * @param[in] player Pointer to a layer object.
 */
 #define mc_layer_cc_top(player) ((player)->cos_critical_top)

/**
 * @brief Evaluates to the critical cosine (total internal reflection)
 *        at the bottom layer boundary.
 * @details If the absolute cosine of the angle of incidence
 *          (with respect to z axis) is less than the critical cosine, the
 *          incident packet is reflected from the boundary.
 * @param[in] player Pointer to a layer object.
 */
 #define mc_layer_cc_bottom(player) ((player)->cos_critical_bottom)

/**
 * @brief Evaluates to the scattering coefficient of the layer.
 * @param[in] player Pointer to a layer object.
 */
#define mc_layer_mus(player) ((player)->mus)

/**
 * @brief Evaluates to the absorption coefficient of the layer.
 * @param[in] player Pointer to a layer object.
 */
#define mc_layer_mua(player) ((player)->mua)

/**
 * @brief Evaluates to the inverse of the absorption coefficient of the layer.
 * @param[in] player Pointer to a layer object.
 */
#define mc_layer_inv_mua(player) \
	(((player)->mua != FP_0) ? mc_fdiv(FP_1, (player)->mua) : FP_INF)

/**
 * @brief Evaluates to the total attenuation coefficient, i.e. the sum of the
 *        layer absorption and scattering coefficients.
 * @param[in] player Pointer to a layer object.
 */
#define mc_layer_mut(player) ((player)->mua + (player)->mus)

/**
 * @brief Evaluates to the reciprocal of the total attenuation coefficient,
 *        i.e. the reciprocal of the sum of the layer absorption and scattering
 *        coefficients.
 * @param[in] player Pointer to a layer object.
 */
#define mc_layer_inv_mut(player) ((player)->inv_mut)

/**
 * @brief Evaluates to the absorption coefficient of the layer multiplied
 *        by the reciprocal of the total attenuation coefficient.
 * @param[in] player Pointer to a layer object.
 */
#define mc_layer_mua_inv_mut(player) ((player)->mua_inv_mut)

/**
 * @} // end @addtogroup mc_layer
 */
/*########################## End layer declarations ##########################*/


/*############# Start Monte Carlo simulator state declarations ###############*/
/**
* @addtogroup mc_simulator_core Simulator core
* @{
*/

/**
 * @brief Simulation core state
 * @{
*/
struct McSimState{
	mc_point3f_t position;		/**< @brief Photon packet position. */
	mc_point3f_t direction;		/**< @brief Photon packet propagation direction. */
	uint64_t rngStateX;			/**< Random generator state X (changes on each mcsim_random call). */
	uint32_t rngStateA;			/**< Random generator state A (not changed by mcsim_random call). */
	mc_fp_t weight;				/**< @brief Photon packet weight. */
	mc_cnt_t photon_index;		/**< @brief Absolute photon packet index. */
	mc_int_t layer_index;		/**< Current layer index. */
	#if MC_TRACK_OPTICAL_PATHLENGTH || defined(__DOXYGEN__)
	mc_fp_t optical_pathlength;		/**< Optical pathlength traveled by the photon packet. */
	#endif
	#if MC_USE_TRACE || defined(__DOXYGEN__)
	mc_uint_t trace_count;		/**< @brief Number of logged trace events since the packet launch. */
	#endif
	#if MC_USE_FLUENCE && MC_USE_FLUENCE_CACHE
	mc_accucache_t fluence_cache;	/**< @brief Fluence cache object. */
	#endif
};
 /** @} */
 /** @brief Data type representing the Monte Carlo simulator core state. */
typedef struct McSimState McSimState;

/**
* @brief Data type holding the simulator state and all required data
* @details  The McSim::remainingStep is normalized and must be multiplied by the
*			inverse of the sum of layer absorption and scattering coefficient
*			to obtain the remaining step size for a particular layer.
* @{
*/
struct McSim{
	McSimState state;		/**< Simulation state. */

	mc_int_t num_layers;	/**< Number of layers including the two outermost layers. */
	__constant McLayer const *layers;	/**< Layer objects. */

	__mc_source_mem McSource const *source; 	/**< Photon packet source object. */

	#if MC_USE_SURFACE_LAYOUTS || defined(__DOXYGEN__)
		__mc_surface_mem McSurfaceLayouts const *surface_layouts;	/**< Advanced layout at the sample top and bottom surfaces. */
	#endif

	#if MC_USE_FP_LUT || defined(__DOXYGEN__)
		__mc_fp_lut_mem mc_fp_t const *fp_lut_array;	/**< Floating-point lookup table(s) data. */
	#endif

	#if MC_USE_TRACE || defined(__DOXYGEN__)
		__mc_trace_mem const McTrace *trace;	/**< @brief Trace configuration data. */;
	#endif

	#if MC_USE_FLUENCE || defined(__DOXYGEN__)
		__mc_fluence_mem const McFluence *fluence;	/**< @brief Fluence array strides (dimensions) as [nx, ny, nz]. */
	#endif

	#if MC_USE_DETECTORS || defined(__DOXYGEN__)
		__mc_detector_mem const McDetectors *detectors;	/**< @brief Reflectance/transmittance detector configuration data. */
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
 * @brief Evaluates to the current position (::mc_point3f_t type) of the photon packet.
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
 *        photon packet position.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_position_r2(psim) \
	((psim)->state.position.x*(psim)->state.position.x + \
	(psim)->state.position.y*(psim)->state.position.y)

/**
 * @brief Evaluates to the r (polar radius) coordinate of the
 *        photon packet position.
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
 * @brief Sets the z coordinate of the photon packet position to the specified value.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] zpos New packet position z coordinate.
 */
#define mcsim_set_position_z(psim, zpos) ((psim)->state.position.z = (zpos))

/**
 * @brief Sets the current position of the photon packet to the specified value (::mc_point3f_t type).
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] ppoint Pointer to the new photon packet position (::mc_point3f_t type).
 */
#define mcsim_set_position(psim, ppoint) ((psim)->state.position = *(ppoint))

/**
 * @brief Sets the current position of the photon packet to the specified value (::mc_point3f_t type).
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
 *        the current photon packet propagation direction by d.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] d Step length along the current photon packet propagation direction.
 */
#define mcsim_move(psim, d) \
	{(psim)->state.position.x += (d)*(psim)->state.direction.x; \
	(psim)->state.position.y += (d)*(psim)->state.direction.y; \
	(psim)->state.position.z += (d)*(psim)->state.direction.z;}

/**
 * @brief Evaluates to the current propagation direction (::mc_point3f_t type) of the
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
	{(psim)->state.direction.x = (pdir)->x; \
	(psim)->state.direction.y = (pdir)->y; \
	(psim)->state.direction.z = (pdir)->z;}

/**
 * @brief Reverse the x component of the photon packet propagation direction.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_reverse_direction_x(psim, pdir) \
	{(psim)->state.direction.x = -(paim)->state.direction.x;}

/**
 * @brief Reverse the y component of the photon packet propagation direction.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_reverse_direction_y(psim, pdir) \
	{(psim)->state.direction.y = -(paim)->state.direction.y;}

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
 * @brief Adjusts the current photon packet weight (from [0.0, 1.0])by
 *        subtracting the specified value.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] delta Subtracted from the current photon packet weight.
 * @note The resulting photon packet weight is NOT adjusted to valid range [0.0, 1.0].
 */
#define mcsim_adjust_weight(psim, delta) ((psim)->state.weight -= (delta))

/**
* @brief Evaluates to the current layer index.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_current_layer_index(psim) ((psim)->state.layer_index)

/**
* @brief Sets the current layer index.
* @param[in] psim Pointer to a simulator instance.
* @param[in] index Current layer index to be set.
* @note	The supplied index is NOT checked for valid range!
*/
#define mcsim_set_current_layer_index(psim, index) ((psim)->state.layer_index = (index))

/**
* @brief Evaluates to the next layer index.
* @param[in] psim Pointer to a simulator instance.
* @note	The calculated index might be out of bounds
*/
#define mcsim_next_layer_index(psim) \
	((psim)->state.layer_index + mc_fsign((psim)->state.direction.z))

/**
* @brief Evaluates to the next layer.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_next_layer(psim) \
	(&(psim)->layers[(psim)->state.layer_index + mc_fsign((psim)->state.direction.z)])

/**
* @brief Evaluates to the total number of layers including the two outer layers.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_layer_count(psim) ((psim)->num_layers)

/**
* @brief Evaluates to nonzero if the layer index is valid. Note that the first
*        and last layer of the stack represent the surrounding medium and are
*        and are not considered valid.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_layer_index_is_sample(psim, index) \
	((index) > 0 && (index) < (psim)->num_layers - 1)

/**
 * @brief Evaluates to nonzero if the photon packet has escaped the sample.
 *        Note that the first and last layer of the stack represent the
 *        surrounding medium and the photon packets located inside the
 *        two outermost layers are considered escaped.
 *
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_packet_escaped_sample(psim) \
	(!mcsim_layer_index_is_sample(psim, mcsim_current_layer_index(psim)))

/**
* @brief Evaluates to a pointer to the current layer.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_current_layer(psim) (&(psim)->layers[(psim)->state.layer_index])

/**
* @brief Evaluates to a pointer to the current layer scattering phase function.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_current_layer_pf(psim) (&mcsim_current_layer(psim)->pf)

/**
 * Portable alias for mcsim_current_layer_pf
 */
#define mcsim_current_pf(psim) mcsim_current_layer_pf(psim)

/**
* @brief Evaluates to a pointer to the specified layer.
* @param[in] psim Pointer to a simulator instance.
* @param[in] index The requested layer index.
* @note The specified layer index is NOT checked for valid range.
*/
#define mcsim_layer(psim, index) (&(psim)->layers[index])

/**
* @brief Evaluates to the index of the top layer.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_top_layer_index(psim) (0)

/**
* @brief Evaluates to a pointer to the top layer.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_top_layer(psim) (&(psim)->layers[0])

/**
* @brief Evaluates to the index of the top layer of the sample.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_top_sample_layer_index(psim) (1)

/**
* @brief Evaluates to a pointer to the top layer of the sample.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_top_sample_layer(psim) (&(psim)->layers[1])

/**
* @brief Evaluates to a pointer to the bottom layer.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_bottom_layer(psim) (&(psim)->layers[mcsim_layer_count(psim) - 1])

/**
* @brief Evaluates to the index of the bottom layer.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_bottom_layer_index(psim) (mcsim_layer_count(psim) - 1)

/**
* @brief Evaluates to a pointer to the bottom layer of the sample.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_bottom_sample_layer(psim) (&(psim)->layers[mcsim_layer_count(psim) - 2])

/**
* @brief Evaluates to the index of the bottom sample layer.
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_bottom_sample_layer_index(psim) (mcsim_layer_count(psim) - 2)

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
* @brief Evaluates to the number of trace events since the photon packet launch.
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
 * @param[in] pos      Position at which to deposit the weight.
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
 *                to compute the fluence rate.
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

#if MC_USE_FLUENCE_CACHE
	/**
	 * @brief Evaluates to a pointer to the fluence cache.
	 *
	 * @return   Pointer to a fluence cache instance.
	 */
	#define mcsim_fluence_cache(psim) (&((psim)->state.fluence_cache))
#endif

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
 * @} // end @addtogroup mc_simulator_core
 */
/*############## End Monte Carlo simulator state declarations ################*/


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
 * @brief Handles layer boundary interactions (refraction/reflection).
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] Next layer nexLayerIndex.
 * @return Returns MC_REFLECTED if the photon packet is reflected from the
 *         boundary or MC_REFRECTED if the photon packet is refracted across
 *         the layer boundary.
 */
inline mc_int_t mcsim_boundary(McSim *psim, mc_int_t nexLayerIndex);

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
	printf("Status: " label "\n"); \
	printf("  Packet index/id: %" FMT_CNT "\n", \
		mcsim_packet_index(psim)); \
	printf("  Position: %.6f, %.6f, %.6f\n", mcsim_position_x(psim), \
			mcsim_position_y(psim), mcsim_position_z(psim)); \
	printf("  Direction: %.6f, %.6f, %.6f\n", mcsim_direction_x(psim), \
			mcsim_direction_y(psim), mcsim_direction_z(psim)); \
	printf("  Weight: %.6f\n", mcsim_weight(psim));

/**
 * @brief Print one sample layer.
 * param[in] prefix Can be used to pass indent string for the layer parameters."
 * @param[in] player Pointer to a layer instance.
 */
#define dbg_print_layer(player, prefix) \
	printf(prefix "d: %.9f\n" \
			prefix "top: %.9f\n" \
			prefix "bottom: %.9f\n" \
			prefix "n: %.9f\n" \
			prefix "cctop: %.9f\n" \
			prefix "ccbottom: %.9f\n" \
			prefix "mua: %.9f\n" \
			prefix "mus: %.9f\n" \
			prefix "inv_mut: %.9f\n" \
			prefix "mua_inv_mut: %.9f\n", \
					(player)->thickness, (player)->top, (player)->bottom, (player)->n, \
					(player)->cos_critical_top, (player)->cos_critical_bottom, \
					(player)->mua, (player)->mus, \
					(player)->inv_mut, (player)->mua_inv_mut); \
			{ McPf const _dbg_pf=(player)->pf; dbg_print_pf(&_dbg_pf); };

#else

#define dbg_print_status(psim, label) ;
#define dbg_print_layer(player, label) ;

#endif

/**
 * @} // end @addtogroup mc_simulator_core
 */
/*##################### End debug support declarations #######################*/

#endif /* #define __MCML_H */
