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

#include "mcvox.template.h"


/**
 * @brief Returns the r squared (polar radius) coordinate of the
 *        photon packet position with respec to the given origin.
 * @param[in] psim Pointer to a simulator instance.
 * @param[in]  Coordinate x of the origin.
 * @param[in] y Coordinate y of the origin.
 * @returns r2 Squared distance from the origin.
 */
inline mc_fp_t mcsim_position_r2_ex(
		McSim const *psim, mc_fp_t x_center, mc_fp_t y_center){
	mc_fp_t dx = mcsim_position_x(psim) - x_center;
	mc_fp_t dy = mcsim_position_y(psim) - y_center;

	return dx*dx + dy*dy;
};

/*############## Start random number generator implementation ################*/
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
	 inline mc_fp_t mcsim_random_double(McSim *psim){
		return fp_random_double(&psim->state.rngStateX, psim->state.rngStateA);
	};
#endif

/**
 * @brief Generates a single precision random number from (0.0, 1.0) and
 *        update the generator state.
 * @note Due to precision issues
 *       the open nature of the interval is not guaranteed.
 *
 * @param[in]	psim Pointer to a simulator instance.
 * @return		A random number from (0.0, 1.0). Due to precision issues
 * 				the open nature of the interval is not guaranteed.
 *
 @details George Marsaglia's Random Number Generator. A single precision
 *        floating-point number has a 23-bit mantisa, hence only integers from
 *        [0, 2^23 - 1] can be represented without loss of information.
 */
 inline mc_fp_t mcsim_random_single(McSim *psim){
	return fp_random_single(&psim->state.rngStateX, psim->state.rngStateA);
};
/*############### End random number generator implementation #################*/


/*######################## Start geometry handling ###########################*/

/**
 * @brief Intersect the sample box and return the intersection and surface
 *        normal that points in the propagation direction.
 *
 * @param psim Pointer to a simulator instance.
 * @param pos Origin of the intersecting ray.
 * @param dir Propagation direction of the intersecting ray
 * @param intersection Computed intersection
 * @param normal Computed surface normal at the intersection.
 *
 * @return Nonzero if an intersection exists.
 */
inline int mcsim_sample_box_intersect(
		McSim *psim, mc_point3f_t const *pos, mc_point3f_t const *dir,
		mc_point3f_t *intersection, mc_point3f_t *normal){
	mc_fp3_t invdir;

	invdir.x = (dir->x != FP_0) ? mc_fdiv(FP_1, dir->x) : FP_INF;
	invdir.y = (dir->y != FP_0) ? mc_fdiv(FP_1, dir->y) : FP_INF;
	invdir.z = (dir->z != FP_0) ? mc_fdiv(FP_1, dir->z) : FP_INF;

	mc_fp_t t1 = (mcsim_top_left_x(psim) - pos->x)*invdir.x;
	mc_fp_t t2 = (mcsim_bottom_right_x(psim) - pos->x)*invdir.x;
	mc_fp_t t3 = (mcsim_top_left_y(psim) - pos->y)*invdir.y;
	mc_fp_t t4 = (mcsim_bottom_right_y(psim) - pos->y)*invdir.y;
	mc_fp_t t5 = (mcsim_top_left_z(psim) - pos->z)*invdir.z;
	mc_fp_t t6 = (mcsim_bottom_right_z(psim) - pos->z)*invdir.z;

	mc_fp_t txmin = mc_fmin(t1, t2);
	mc_fp_t txmax = mc_fmax(t1, t2);
	mc_fp_t tymin = mc_fmin(t3, t4);
	mc_fp_t tymax = mc_fmax(t3, t4);
	mc_fp_t tzmin = mc_fmin(t5, t6);
	mc_fp_t tzmax = mc_fmax(t5, t6);

	dbg_print_float("txmin:", txmin);
	dbg_print_float("txmax:", txmax);
	dbg_print_float("tymin:", tymin);
	dbg_print_float("tymax:", tymax);
	dbg_print_float("tzmin:", tzmin);
	dbg_print_float("tzmax:", tzmax);

	mc_fp_t tmin = mc_fmax(mc_fmax(txmin, tymin), tzmin);
	mc_fp_t tmax = mc_fmin(mc_fmin(txmax, tymax), tzmax);

	dbg_print_float("tmin 1:", tmin);
	dbg_print_float("tmax 1:", tmax);

	if (tmin > tmax)
		return 0;

	normal->x = (txmin >= tymin && txmin >= tzmin);
	normal->y = ((normal->x == FP_0) && tymin >= tzmin);
	normal->z = normal->x == FP_0 && normal->y == FP_0;

	normal->x *= mc_fsign(dir->x);
	normal->y *= mc_fsign(dir->y);
	normal->z *= mc_fsign(dir->z);

	mc_fp_t t = (tmin > FP_0) ? tmin : tmax;

	dbg_print_float("t:", t);

	intersection->x = pos->x + t*dir->x;
	intersection->y = pos->y + t*dir->y;
	intersection->z = pos->z + t*dir->z;

	dbg_print_point3f("Intersection:", intersection);

	intersection->x = mc_fclip(
		intersection->x, mcsim_top_left_x(psim),
		mcsim_bottom_right_x(psim) - FP_EPS);
	intersection->y = mc_fclip(
		intersection->y, mcsim_top_left_y(psim),
		mcsim_bottom_right_y(psim) - FP_EPS);
	intersection->z = mc_fclip(
		intersection->z, mcsim_top_left_z(psim),
		mcsim_bottom_right_z(psim) - FP_EPS);

	return tmax >= tmin && tmax >= FP_0;
};

/**
 * @brief Computes distance to intersection along the propagation direction
 * 			with all the boundaries of the current voxel.
 *
 * @param[in] psim Pointer to a simulator instance.
 * @param[out] distance Filled with the computed distances on return;
 *
 * @return Distance to the nearest boundary along the propagation direction.
 */
inline mc_fp_t mcsim_intersect(McSim *psim, mc_point3f_t *distance){
	distance->x = mcsim_top_left_x(psim) +
		((mcsim_direction_x(psim) >= FP_0) + mcsim_voxel_index_x(psim))*
		mcsim_voxel_size_x(psim) - mcsim_position_x(psim);

	distance->y = mcsim_top_left_y(psim) +
		((mcsim_direction_y(psim) >= FP_0) + mcsim_voxel_index_y(psim))*
		mcsim_voxel_size_y(psim) - mcsim_position_y(psim);

	distance->z = mcsim_top_left_z(psim) +
		((mcsim_direction_z(psim) >= FP_0) + mcsim_voxel_index_z(psim))*
		mcsim_voxel_size_z(psim) - mcsim_position_z(psim);

	distance->x = (mcsim_direction_x(psim) != FP_0) ?
		mc_fdiv(distance->x, mcsim_direction_x(psim)) : FP_INF;

	distance->y = (mcsim_direction_y(psim) != FP_0) ?
		mc_fdiv(distance->y, mcsim_direction_y(psim)) : FP_INF;

	distance->z = (mcsim_direction_z(psim) != FP_0) ?
		mc_fdiv(distance->z, mcsim_direction_z(psim)) : FP_INF;

	return mc_fmin(distance->x, mc_fmin(distance->y, distance->z));
};

/**
 * @brief Compute the voxel index from the given position. No check are made
 *        if the position is within the sample box.
 *
 * param[in]  psim Pointer to a simulator instance.
 * param[in]  pos  Position.
 * param[out] ind  Voxel index.
 */
inline void mcsim_position_to_voxel_index(
		McSim *psim, mc_point3f_t const *pos, mc_point3_t *ind){
	ind->x = mc_int(
		mc_fdiv(mcsim_position_x(psim) - mcsim_top_left_x(psim),
			mcsim_voxel_size_x(psim))
	);
	ind->y = mc_int(
		mc_fdiv(mcsim_position_y(psim) - mcsim_top_left_y(psim),
			mcsim_voxel_size_y(psim))
	);
	ind->z = mc_int(
		mc_fdiv(mcsim_position_z(psim) - mcsim_top_left_z(psim),
			mcsim_voxel_size_z(psim))
	);
};

/**
 * @brief Compute the voxel index from the given position. Clip the indices
 *        to the valid range.
 *
 * param[in]  psim Pointer to a simulator instance.
 * param[in]  pos  Position.
 * param[out] ind  Voxel index.
 */
inline void mcsim_position_to_voxel_index_safe(
		McSim *psim, mc_point3f_t const *pos, mc_point3_t *ind){
	mcsim_position_to_voxel_index(psim, pos, ind);
	ind->x = mc_clip(ind->x, 0, mcsim_shape_x(psim) - 1);
	ind->y = mc_clip(ind->y, 0, mcsim_shape_y(psim) - 1);
	ind->z = mc_clip(ind->z, 0, mcsim_shape_z(psim) - 1);
};

/**
 * @brief Update the current voxel index from the current position of the
 *			photon packet. Note that the lower bounds of the voxel boundaries
 *			are included in the voxel, but the upper bounds are not.
 *
 * param[in] psim Pointer to a simulator instance.
 */
inline void mcsim_set_voxel_index_from_position(McSim *psim){
	mcsim_position_to_voxel_index(
		psim, &psim->state.position, &psim->state.voxel_index);

	psim->state.voxel_material_index = mcsim_voxel_material_index(
		psim, &psim->state.voxel_index);
};

/**
 * @brief Call the user defined photon packet launcher and update the
 *		  current voxel index.
 *
 * param[in] psim Pointer to a simulator instance.
*/
inline void mcsim_launch_one(McSim *psim){
	mcsim_launch(psim);
	mcsim_set_voxel_index_from_position(psim);
};
/*######################### End geometry handling ############################*/


/*############## Start layer boundary handler implementation #################*/
/**
 * @brief Handles voxel boundary interactions (refraction/reflection).
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] distances Distances to the voxel boundaries.
 * @return Returns MC_REFLECTED if the photon packet is reflected from the
 *         boundary or MC_REFRECTED if the photon packet is refracted across
 *         the voxel boundary.
 */
inline mc_int_t mcsim_boundary(McSim *psim, const mc_point3f_t *distances){

	mc_point3_t normal;
	mc_point3_t next_voxel_index;

	mc_fp_t cc_n1_n2, cos1, R;
	mc_fp_t n1, n2;

	#if MC_GEOMETRY_EX
		mc_int_t exNextMaterialIndex;
	#endif

	/* find the boundary normal that points outwards */
	normal.x = (distances->x <= distances->y) && (distances->x <= distances->z);
	normal.y = normal.x ?
		0 : (distances->y <= distances->x) && (distances->y <= distances->z);
	normal.z = !(normal.x + normal.y);
	normal.x = mcsim_direction_x(psim) < FP_0 ? -normal.x : normal.x;
	normal.y = mcsim_direction_y(psim) < FP_0 ? -normal.y : normal.y;
	normal.z = mcsim_direction_z(psim) < FP_0 ? -normal.z : normal.z;

	/* compute the next voxel index */
	next_voxel_index.x = mcsim_voxel_index_x(psim) + normal.x;
	next_voxel_index.y = mcsim_voxel_index_y(psim) + normal.y;
	next_voxel_index.z = mcsim_voxel_index_z(psim) + normal.z;

	/* check if the photon packet might escape the sample */
	int escaping_sample = !mcsim_is_valid_voxel_index(psim, &next_voxel_index);

	/* fetch the materials of this and the next voxel */
	__mc_material_mem const McMaterial *current_material =
		mcsim_current_voxel_material(psim);

	mc_int_t next_voxel_material_index = !escaping_sample ?
		mcsim_voxel_material_index(psim, &next_voxel_index) :
		mcsim_surrounding_material_index(psim);

	__mc_material_mem const McMaterial *next_material =
		mcsim_material(psim, next_voxel_material_index);

	#if MC_ENABLE_DEBUG_
		dbg_print("Boundary:");
		dbg_print_float(INDENT "n1:", current_material->n);
		dbg_print_float(INDENT "n2:", next_material->n);
		dbg_print_point3(INDENT "Next voxel index:", &next_voxel_index);
		dbg_print_int(INDENT "Photon packet escaping sample:", escaping_sample);
	#endif

	/* the same material and not escaping the sample */
	if (current_material == next_material && !escaping_sample){
		/* can keep the current material */
		mcsim_set_voxel_index(psim, &next_voxel_index);
		return MC_REFRACTED;
	};

	n1 = current_material->n;
	n2 = next_material->n;

	#if MC_USE_TOP_SURFACE_LAYOUT
	/* advance layout at the top sample surface */
	if(next_voxel_index.z < 0)
	{
        int top_res = mcsim_top_surface_layout_handler(psim, &n2);
		if(top_res != MC_SURFACE_LAYOUT_CONTINUE){
			/* boundary handling done/completed by the top layout handler */
			return top_res;
		}
	}else
	#endif
	#if MC_USE_BOTTOM_SURFACE_LAYOUT
	if(next_voxel_index.z >= mcsim_shape_z(psim))
	{
		int bottom_res = mcsim_bottom_surface_layout_handler(psim, &n2);
		if (bottom_res != MC_SURFACE_LAYOUT_CONTINUE){
			/* boundary handling done/completed by the top layout handler */
			return bottom_res;
		}
	}
	#endif

	/* refractive indices match */
	if(n1 == n2){
		mcsim_set_voxel_index(psim, &next_voxel_index);
		mcsim_set_voxel_material_index(psim, next_voxel_material_index);
		return MC_REFRACTED;
	};

	cc_n1_n2 = cos_critical(n1, n2);

	/* the incidence cosine - positive since the normal is pointing outwards */
	cos1 = normal.x*mcsim_direction_x(psim) +
		normal.y*mcsim_direction_y(psim) +
		normal.z*mcsim_direction_z(psim);

	mc_point3f_t fpnormal =
		{(mc_fp_t)normal.x, (mc_fp_t)normal.y, (mc_fp_t)normal.z};

	/* check if above the critical angle - otherwise reflection holds */
	if(cos1 > cc_n1_n2){
		R = reflectance(n1, n2, cos1, cc_n1_n2);
		if (R < mcsim_random(psim)){ /* [1.0, 0.0) */
			/* we have a refraction */
			refract(mcsim_direction(psim), &fpnormal, n1, n2, mcsim_direction(psim));
			mcsim_set_voxel_index(psim, &next_voxel_index);
			mcsim_set_voxel_material_index(psim, next_voxel_material_index);
			return MC_REFRACTED;
		};
	};

	/* this is a reflection */
	reflect(mcsim_direction(psim), &fpnormal, mcsim_direction(psim));

	return MC_REFLECTED;
};
/*############### End layer boundary handler implementation ##################*/


/*############### Start photon packet fluence implementation #################*/
#if MC_USE_FLUENCE || defined(__DOXYGEN__)
/**
 * @brief Deposit the given weight to the fluence accumulator. Takes care
 *        of the different function signatures in case of the weight deposition
 *        and fluence rate mode.
 *
 * @param[in] sim      Simulator instance.
 * @param[in] pos      Position at which to deposit the weight.
 * @param[in] deposit  Weight to deposit.
 */
inline void mcsim_fluence_deposit_weight(
		McSim *psim, mc_point3f_t const *pos, mc_fp_t deposit){
	#if MC_FLUENCE_MODE_RATE
		mcsim_fluence_deposit_at(psim, pos, deposit,
			mc_material_mua(
				mcsim_current_voxel_material(psim),  mcsim_direction(psim))
		);
	#else
		mcsim_fluence_deposit_at(psim, pos, deposit);
	#endif
};

/**
 * @brief Low-level deposition function that can use an intermediate cache if
 *        configured so through the ::MC_USE_FLUENCE_CACHE option.
 *
 * @param psim     Simulator instance.
 * @param offset   Deposition address/offset.
 * @param weight   Weight to deposit.
 */
inline void mcsim_fluence_weight_deposit_ll(
		McSim *psim, size_t offset, uint32_t weight){
	#if MC_USE_FLUENCE_CACHE
		mc_accucache_weight_add(
			mcsim_fluence_cache(psim), offset, weight,
			mcsim_accumulator_buffer(psim)
		);
	#else
		accumulator_deposit(
			(__global void *)mcsim_accumulator_buffer_ex(psim, offset),
			weight
		);
	#endif
};

#endif /* MC_USE_FLUENCE */
/*################ End photon packet fluence implementation ##################*/


/*################ Start photon packet trace implementation ##################*/
#if MC_USE_TRACE || defined(__DOXYGEN__)
/**
 * @brief Call the user-defined trace function and increment the number of
 * logged trace events if function returns nonzero.
 * @param[in] psim Pointer to a simulator instance.
 */
static inline void mcsim_trace_this_event(McSim *psim){
	if (mcsim_trace_event(psim, psim->state.trace_count))
		psim->state.trace_count++;
};

/**
 * Finalizes the trace - saves the number of trace events.
 */
static inline void mcsim_trace_finalize(McSim *psim){
	mcsim_trace_complete(psim, mcsim_trace_count(psim));
};
#endif
/*################# End photon packet trace implementation ###################*/


/*################ Start scattering handling implementations ##################*/
/**
* @brief Called when the photon packet needs to be scattered.
* @param[in, out] psim Simulator satate.
* @details Function computes and updates the packet propagation direction
*          and associated states photon packet states.
*/
inline void mcsim_scatter(McSim *psim){
	#if defined(MC_PF_SAMPLE_DIRECTION)
	mc_point3f_t new_dir;
	mcsim_pf_sample_dir(psim, &new_dir);
	mcsim_set_direction(psim, &new_dir);
	#else
	mc_fp_t fi, cos_theta;

	/* sample the scattering phase functions */
	cos_theta = mcsim_pf_sample_angles(psim, &fi);

	scatter_direction(mcsim_direction(psim), cos_theta, fi);
	#endif
};
/*################# End scattering handling implementation ###################*/


/*#################### Start sample voxel implementation #####################*/
/* User-defined implementation of sample voxels - DO NOT EDIT! */
/* START_VOXELS_IMPLEMENTATION_BLOCK */
{{ voxels.implementation or '' }}
/* END_VOXELS_IMPLEMENTATION_BLOCK */
/*##################### End sample voxel implementation ######################*/


/*##################### Start materials implementation #######################*/
/* User-defined implementation of materials including
   the scattering phase function and related API- DO NOT EDIT! */
/* START_MATERIALS_IMPLEMENTATION_BLOCK */
{{ materials.implementation or '' }}
/* END_MATERIALS_IMPLEMENTATION_BLOCK */
/*###################### End materials implementation ########################*/


/*##################### Start materials implementation #######################*/
/* User-defined implementation of the photon packet source - DO NOT EDIT! */
/* START_SOURCE_IMPLEMENTATION_BLOCK */
{{ source.implementation or '' }}
/* END_SOURCE_IMPLEMENTATION_BLOCK */
/*###################### End materials implementation ########################*/


/*##################### Start detectors implementation #######################*/
/* User-defined surface reflectance/transmittance detector implementation
   goes here - DO NOT EDIT! */
/* START_DETECTORS_IMPLEMENTATION_BLOCK */
{{ detectors.implementation or '' }}
/* END_DETECTORS_IMPLEMENTATION_BLOCK */
/*####################### End detectors implementation #######################*/


/*################ Start photon packet trace implementation ##################*/
/* User-defined photon packet trace implementation goes here - DO NOT EDIT! */
/* START_TRACE_IMPLEMENTATION_BLOCK */
{{ trace.implementation or '' }}
/* END_TRACE_IMPLEMENTATION_BLOCK */
/*################# End photon packet trace implementation ###################*/


/*############### Start photon packet fluence implementation #################*/
/* User-defined fluence implementation goes here - DO NOT EDIT! */
/* START_FLUENCE_IMPLEMENTATION_BLOCK */
{{ fluence.implementation or ''}}
/* END_FLUENCE_IMPLEMENTATION_BLOCK */
/*################ End photon packet fluence implementation ##################*/


/*############# Start top sample surface layout implementation ###############*/
/* User-defined sample surface layout implementation
   goes here - DO NOT EDIT! */
/* START_SURFACE_LAYOUTS_IMPLEMENTATION_BLOCK */
{{ surface_layouts.implementation or '' }}
/* END_SURFACE_LAYOUTS_IMPLEMENTATION_BLOCK */
/*############## End top sample surface layout implementation ################*/


/*#################### Start Monte Carlo OpenCL kernel #######################*/
__kernel void McKernel(
	mc_cnt_t num_packets,
	__global mc_cnt_t *num_packets_done,

	__global uint32_t *num_kernels,

	mc_fp_t mc_rmax,

	__global uint64_t *rng_state_x,
	__global uint32_t const *rng_state_a,	// must nut be constant - insufficient memory on GPU

	__mc_voxelcfg_mem McVoxelConfig const *voxel_cfg,
	__global McVoxel const *voxels,

	mc_size_t num_materials,
	__mc_material_mem McMaterial const *materials,

	__mc_source_mem McSource const *source,
	__mc_surface_mem McSurfaceLayouts const *surface_layouts,

	__mc_trace_mem McTrace const *trace,

	__mc_fluence_mem McFluence const *fluence,

	__mc_detector_mem McDetectors const *detectors,

	__mc_fp_lut_mem mc_fp_t const *fp_lut_array,

	__global mc_int_t *integer_buffer,
	__global mc_fp_t *float_buffer,
	__global mc_accu_t *accumulator_buffer
)
{
	mc_fp_t d, d_ok, step, deposit;
	mc_point3f_t distances;
	bool done = false;

	mc_point3f_t src_pos = source->position;

	/* create and initialize a simulation structure */
	McSim sim = {
		{				/* McSimState state: simulator state */
			{FP_0, FP_0, FP_0}		/* mc_point3f_t position: Current photon packet position. */
			,{FP_0, FP_0, FP_0}		/* mc_point3f_t direction: Current photon packet propagation direction. */
			,{0, 0, 0}				/* mc_point_t voxel_index: Index of the current voxel. */
			,0						/* mc_size_t voxel_material_index: Index of the current voxel material. */
			,rng_state_x[get_global_id(0)]	/* Random generator state X (changes on each mcsim_random call). */
			,rng_state_a[get_global_id(0)]	/* Random generator state A (not changed by mcsim_random call). */
			,FP_1							/* mc_fp_t weight: Current photon packet weight from [0.0, 1.0]. */
			,0								/* mc_cnt_t: Absolute photon packet index. */
			#if MC_USE_EVENTS || defined(__DOXYGEN__)
			,0						/* mc_uint_t event_flags: All events that the packet underwent during this step. */
			#endif
			#if MC_TRACK_OPTICAL_PATHLENGTH || defined(__DOXYGEN__)
				,FP_0				/* mc_fp_t optical_pathlength: Optical pathlength traveled by the photon packet. */
			#endif
			#if MC_USE_TRACE || defined(__DOXYGEN__)
			,0						/* mc_uint_t trace_count: Number of traced events. */
			#endif
			#if MC_USE_FLUENCE && MC_USE_FLUENCE_CACHE
			,mc_accucache_initializer	/* mc_accucache_t fluence_cache: Fluence cache object. */
			#endif
		},

		voxel_cfg,		/* __mc_voxelcfg_mem McVoxelConfig const *voxel_cfg: Voxel array configuration. */
		voxels,			/* __global McVoxel const *voxels: Voxel array. */

		num_materials,	/* mc_size_t num_materials: Number of materials. */
		materials, 		/* __mc_material_mem McMaterial const *materials: Array of all materials. */

		source			/* __mc_source_mem McSource const *source: Photon packet source object. */

		#if MC_USE_SURFACE_LAYOUTS
			,surface_layouts		/* __mc_surface_mem McSurfaceLayouts const *surface_layout: Advanced layout at the sample top and bottom surfaces. */
		#endif

		#if MC_USE_FP_LUT
			,fp_lut_array			/* __mc_fp_lut_mem mc_fp_t const *fp_lut_array: Lookup table(s) data. */
		#endif

		#if MC_USE_TRACE
			,trace					/* __mc_trace_mem McTrace const *trace: Trace configuration object. */
		#endif

		#if MC_USE_FLUENCE
			,fluence				/* __mc_fluence_mem McFluence const *fluence: Fluence configuration struct */
		#endif

		#if MC_USE_DETECTORS
			,detectors				/* __mc_detector_mem McDetectors const *detectors: Surface and specular detectors. */
		#endif

		,integer_buffer		/* __global mc_int_t *integer_buffer: Common integer buffer. */
		,float_buffer		/* __global mc_fp_t *float_buffer: Common floating-point buffer. */
		,accumulator_buffer	/* __global mc_accu_t *accumulator_buffer: Common accumulator buffer. */
	};

	/* make unused variables void to surpresss compiler warnings */
	#if !MC_USE_SURFACE_LAYOUTS
		(void)surface_layouts;
	#endif

	#if !MC_USE_FP_LUT
		(void)fp_lut_array;
	#endif

	#if !MC_USE_TRACE
		(void)trace;
	#endif

	#if !MC_USE_FLUENCE
		(void)fluence;
	#endif

	#if !MC_USE_DETECTORS
		(void)detectors;
	#endif

	/* check if photon packets have to be launched by this thread */
	if((sim.state.photon_index = pkt_cnt_atomic_inc(num_packets_done)) < num_packets){

		atomic_inc(num_kernels);

		#if MC_ENABLE_DEBUG
		/* print debug information ony for the first photon packet */
		if (sim.state.photon_index == 0) {
			dbg_print("");
			dbg_print("OpenCL runtime:");
			dbg_printf(INDENT "supported version: %.1f\n",
				__OPENCL_VERSION__/FP_LITERAL(100.0));
			dbg_printf(INDENT "built for version: %.1f\n",
				__OPENCL_C_VERSION__/FP_LITERAL(100.0));

			dbg_print("");
			dbg_print("Simulation:");
			dbg_print_cnt(INDENT   "num_packets     :", num_packets);
			dbg_print_cnt(INDENT   "num_packets_done:", *num_packets_done);
			dbg_print_cnt(INDENT   "num_materials   :", num_materials);
			dbg_print_float(INDENT "mc_rmax (mm)    :", mc_rmax*FP_LITERAL(1e3));

			/* Print the voxel array for small arrays */
			dbg_print("");
			if (mcsim_size(&sim) <= 512) {
				dbg_print("Voxel array:");
				dbg_printf("[\n");
				for (int k=0; k < mcsim_shape_x(&sim); ++k) {
					dbg_printf(" [\n");
					for (int j=0; j < mcsim_shape_x(&sim); ++j) {
						dbg_printf("  [");
						for (int i=0; i < mcsim_shape_x(&sim); ++i) {
							mc_point3_t index = {i, j, k};
							mc_size_t flat_index = mcsim_flat_voxel_index(&sim, &index);
							mc_int_t voxel = mcsim_voxels(&sim)[flat_index].material_index;
							if (i < mcsim_shape_x(&sim) - 1)
								dbg_printf("%4d, ", voxel);
							else
								dbg_printf("%4d", voxel);
						}
						if (j < mcsim_shape_y(&sim) - 1)
							dbg_printf("],\n");
						else
							dbg_printf("]\n");
					}
					if (k < mcsim_shape_z(&sim) - 1)
						dbg_printf(" ],\n");
					else
						dbg_printf(" ]\n");
				}
				dbg_printf("]\n");
			};

			dbg_print("");
			for (mc_size_t i=0; i < num_materials; ++i) {
				dbg_print_size_t("Material:", i);
				dbg_print_material(&materials[i]);
			};

			dbg_print(""); dbg_print_voxel_cfg(voxel_cfg);

			dbg_print(""); dbg_print_source(source);

			#if MC_USE_SURFACE_LAYOUTS
				dbg_print(""); dbg_print_surface_layouts(surface_layouts);
			#endif

			#if MC_USE_TRACE
				dbg_print(""); dbg_print_trace(trace);
			#endif

			#if MC_USE_FLUENCE
				dbg_print_fluence(fluence);
			#endif

			#if MC_USE_DETECTORS
				dbg_print(""); dbg_print_detectors(detectors);
			#endif

			dbg_print("#########################################");
		};
		#endif

		/* launch a new photon packet */
		mcsim_launch_one(&sim);
		mcsim_event_flags_add(&sim, MC_EVENT_PACKET_LAUNCH);

		#if MC_USE_TRACE & MC_USE_TRACE_START
			/* initial photon packet state */
			mcsim_trace_this_event(&sim);
		#endif

		/* loop through the simulation steps until all the photon packets 
			have been processed */
		while (!done) {

			/* generate a new step */
			#if MC_METHOD == MICROSCOPIC_BEER_LAMBERT
				step = mc_fdiv(-mc_log(mcsim_random(&sim)),
					mc_material_mus(mcsim_current_voxel_material(&sim), mcsim_direction(&sim)));
			#else
				step = -mc_log(mcsim_random(&sim))*
					mc_material_inv_mut(mcsim_current_voxel_material(&sim), mcsim_direction(&sim));
			#endif
			step = mc_fmin(step, FP_MAX);

			/* compute distance to the voxel boundaries */
			d = mcsim_intersect(&sim, &distances);

			#if MC_ENABLE_DEBUG
				dbg_print("Step to boundary:\n");
				dbg_print_float(INDENT "step", step);
				dbg_print_point3f(INDENT "Distances:", &distances);
				dbg_print_float(INDENT "d:", d);
				dbg_print_status(&sim, "Status:");
			#endif

			/* propagate the photon packet */
			d_ok = mc_fmin(d, step);
			mcsim_set_position_coordinates(
				&sim,
				mcsim_position_x(&sim) + d_ok*mcsim_direction_x(&sim),
				mcsim_position_y(&sim) + d_ok*mcsim_direction_y(&sim),
				mcsim_position_z(&sim) + d_ok*mcsim_direction_z(&sim)
			);

			/* update total optical pathlength of the photon packet*/
			#if MC_TRACK_OPTICAL_PATHLENGTH
				mcsim_optical_pathlength_add(
					&sim,
					mc_material_n(mcsim_current_voxel_material(&sim))*d_ok
				);
			#endif

			#if MC_METHOD == MICROSCOPIC_BEER_LAMBERT
				mc_fp_t mua = mc_material_mua(
					mcsim_current_voxel_material(&sim), mcsim_direction(&sim));
				mc_fp_t deposit_fraction = FP_1 - mc_exp(-mua*d_ok);
				deposit = deposit_fraction*mcsim_weight(&sim);
				mcsim_adjust_weight(&sim, deposit);
				mcsim_event_flags_add(&sim, MC_EVENT_PACKET_ABSORPTION);

				#if MC_USE_FLUENCE
					mc_fp_t step_back = (mua != FP_0) ?
						d_ok - mc_fdiv(-mc_log(FP_1 - mcsim_random(&sim)*deposit_fraction), mua) :
						FP_0;
					//mc_fp_t w_c = FP_0p5*step;
					mc_point3f_t deposit_pos = {
						mcsim_position_x(&sim) - step_back*mcsim_direction_x(&sim),
						mcsim_position_y(&sim) - step_back*mcsim_direction_y(&sim),
						mcsim_position_z(&sim) - step_back*mcsim_direction_z(&sim)
					};
					mcsim_fluence_deposit_weight(&sim, &deposit_pos, deposit);
				#endif

				/* process boundary hit or handle absorption */
				if (d < step){
					/* deal with the boundary hit ... returns MC_REFRACTED or MC_REFLECTED */
					mc_uint_t _boundary_flags = mcsim_boundary(&sim, &distances);
					mcsim_event_flags_add(&sim, _boundary_flags | MC_EVENT_BOUNDARY_HIT);
					(void)_boundary_flags;

					#if MC_ENABLE_DEBUG
						dbg_print_status(&sim, "Boundary hit:");
					#endif

					if (mcsim_packet_escaped_sample(&sim)) {
						#if MC_USE_TOP_DETECTOR || MC_USE_BOTTOM_DETECTOR
							if ( (mcsim_voxel_index_z(&sim) < 0) ){
								#if MC_USE_TOP_DETECTOR
									mcsim_top_detector_deposit(
										&sim,
										mcsim_position(&sim),
										mcsim_direction(&sim),
										mcsim_weight(&sim)
									);
								#endif
							} else if (mcsim_voxel_index_z(&sim) >=
									mcsim_shape_z(&sim)){
								#if MC_USE_BOTTOM_DETECTOR
									mcsim_bottom_detector_deposit(
										&sim,
										mcsim_position(&sim),
										mcsim_direction(&sim),
										mcsim_weight(&sim)
									);
								#endif
							}
						#endif

						done = true; /* photon packet escaped the sample volume */
					}

				} else {
					/* Scatter the photon packet, */
					mcsim_scatter(&sim);
					mcsim_event_flags_add(&sim, MC_EVENT_PACKET_SCATTERING);
				}

				/* Perform survival lottery if required (packet not done).*/
				if(mcsim_weight(&sim) < MC_PACKET_WEIGHT_MIN){
					#if MC_USE_LOTTERY
						/* perform lottery */
						if(mcsim_random(&sim) > MC_PACKET_LOTTERY_CHANCE){
							/* should leave the weight as is */
							/* mcsim_set_weight(&sim, FP_0); */
							done = true;
						}else{
							mcsim_set_weight(&sim,
								mc_fdiv(mcsim_weight(&sim), MC_PACKET_LOTTERY_CHANCE));
						};
					#else
						/* no lottery - terminate on weight threshold */
						done = true;
					#endif

				};

			#else /* Albedo Weight or Albedo Rejection */
				/* process boundary hit or handle absorption */
				if (d < step){
					/* deal with the boundary hit ... returns MC_REFRACTED or MC_REFLECTED */
					mc_uint_t _boundary_flags = mcsim_boundary(&sim, &distances);
					mcsim_event_flags_add(&sim, _boundary_flags | MC_EVENT_BOUNDARY_HIT);
					(void)_boundary_flags;
					/* check if photon escaped the sample */

					#if MC_ENABLE_DEBUG
						dbg_print_status(&sim, "Boundary hit:");
					#endif

					if (mcsim_packet_escaped_sample(&sim)) {
						#if MC_USE_TOP_DETECTOR || MC_USE_BOTTOM_DETECTOR
							if ( (mcsim_voxel_index_z(&sim) < 0) ){
								#if MC_USE_TOP_DETECTOR
									mcsim_top_detector_deposit(
										&sim,
										mcsim_position(&sim),
										mcsim_direction(&sim),
										mcsim_weight(&sim)
									);
								#endif
							} else if (mcsim_voxel_index_z(&sim) >=
									mcsim_shape_z(&sim)){
								#if MC_USE_BOTTOM_DETECTOR
									mcsim_bottom_detector_deposit(
										&sim,
										mcsim_position(&sim),
										mcsim_direction(&sim),
										mcsim_weight(&sim)
									);
								#endif
							}
						#endif

						done = true; /* photon packet escaped the sample volume */
					}

				} else {
					#if MC_METHOD == ALBEDO_REJECTION
						/* Do absorption or scattering only when no layer boundary has been hit.*/
						if (mcsim_random(&sim) < mc_material_mua_inv_mut(
								mcsim_current_voxel_material(&sim), mcsim_direction(&sim))) {
							/* Deposit the entire weight of the packet. */
							deposit = mcsim_weight(&sim);
							mcsim_adjust_weight(&sim, deposit);
							mcsim_event_flags_add(&sim, MC_EVENT_PACKET_ABSORPTION);
							done = true;
							#if MC_USE_FLUENCE
								/* update the fluence data in fluence mode */
								mcsim_fluence_deposit_weight(
									&sim, mcsim_position(&sim), deposit);
							#endif
						} else {
							/* Scatter the photon packet, */
							mcsim_scatter(&sim);
							mcsim_event_flags_add(&sim, MC_EVENT_PACKET_SCATTERING);
						}

					#else	/* Albedo Weight */
						/* Do absorption only when no voxel boundary has been hit.*/
						deposit = mcsim_weight(&sim)*
							mc_material_mua_inv_mut(mcsim_current_voxel_material(&sim), mcsim_direction(&sim));
						mcsim_adjust_weight(&sim, deposit);
						mcsim_event_flags_add(&sim, MC_EVENT_PACKET_ABSORPTION);
						#if MC_USE_FLUENCE
							/* update the fluence data in fluence mode */
							mcsim_fluence_deposit_weight(
								&sim, mcsim_position(&sim), deposit);
						#endif

						/* Scatter the photon packet, */
						mcsim_scatter(&sim);
						mcsim_event_flags_add(&sim, MC_EVENT_PACKET_SCATTERING);

						/* Perform survival lottery if required (packet not done).*/
						if(mcsim_weight(&sim) < MC_PACKET_WEIGHT_MIN){
							#if MC_USE_LOTTERY
								/* perform lottery */
								if(mcsim_random(&sim) > MC_PACKET_LOTTERY_CHANCE){
									/* should leave the weight as is */
									/* mcsim_set_weight(&sim, FP_0); */
									done = true;
								}else{
									mcsim_set_weight(&sim,
										mc_fdiv(mcsim_weight(&sim), MC_PACKET_LOTTERY_CHANCE));
								};
							#else
								/* no lottery - terminate on weight threshold */
								done = true;
							#endif
						};
					#endif
				}
			#endif

			/* check if photon escaped the predefined simulation domain */
			if (mc_distance2_point3f(mcsim_position(&sim), &src_pos) >
					mc_rmax*mc_rmax || mcsim_weight(&sim) <= FP_0){
				done = true;
				mcsim_event_flags_add(&sim, MC_EVENT_PACKET_ESCAPED);
			};

			mcsim_event_flags_add(&sim, (done) ? MC_EVENT_PACKET_TERMINATED : 0);

			#if MC_USE_TRACE == MC_USE_TRACE_ALL
				mcsim_trace_this_event(&sim);
			#elif MC_USE_TRACE & MC_USE_TRACE_END
				if (done)
					mcsim_trace_this_event(&sim);
			#endif

			/* clear the event flags for this step */
			mcsim_event_flags_clear(&sim);

			if (done) {
				/* Finalize the trace - only saves the number of trace events. */
				#if MC_USE_TRACE
					mcsim_trace_finalize(&sim);
				#endif

				/* call user defined termination */
				#if defined(MC_TERMINAL_HOOK)
					/* Photon packet has escaped the sample voxels -
						call the provided user defined hook. */
					MC_TERMINAL_HOOK(&sim);
				#endif

				if ((sim.state.photon_index = pkt_cnt_atomic_inc(num_packets_done)) < num_packets) {
					/* initialize a new photon packet trace if required */
					#if MC_USE_TRACE
						/* initialize the trace */
						sim.state.trace_count = 0;
					#endif

					/* initialize the optical pathlength */
					#if MC_TRACK_OPTICAL_PATHLENGTH
						/* initialize the trace */
						sim.state.optical_pathlength = FP_0;
					#endif

					/* launch a new photon packet */
					mcsim_launch_one(&sim);
					mcsim_event_flags_add(&sim, MC_EVENT_PACKET_LAUNCH);

					#if MC_USE_TRACE & MC_USE_TRACE_START
						/* initial photon packet state */
						mcsim_trace_this_event(&sim);
					#endif
					done = false; /* have a new packet ... not done yet */
				}
			}
		};

		#if MC_USE_FLUENCE && MC_USE_FLUENCE_CACHE
			mc_accucache_flush(
				mcsim_fluence_cache(&sim), mcsim_accumulator_buffer(&sim))
		#endif

		/* save/update the random number generator state to the global memory */
		rng_state_x[get_global_id(0)] = sim.state.rngStateX;
	};
};
/*##################### End Monte Carlo OpenCL kernel ########################*/

__kernel void sizeof_datatypes(__global uint *n){
	if (get_global_id(0) == 0){
		n[0] = sizeof(McVoxel);
		n[1] = sizeof(McMaterial);
		n[2] = sizeof(McSource);
		#if MC_USE_SURFACE_LAYOUTS
		n[3] = sizeof(McSurfaceLayouts);
		#endif
		#if MC_USE_DETECTORS
		n[4] = sizeof(McDetectors);
		#endif
		#if MC_USE_TRACE
		n[5] = sizeof(McTrace);
		#endif
		#if MC_USE_FLUENCE
		n[6] = sizeof(McFluence);
		#endif
	}
};
