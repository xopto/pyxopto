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

#include "mccyl.template.h"


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

/**
 * @brief Returns the r squared (polar radius) coordinate of the
 *        photon packet position with respect to the given origin.
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

/*########### Start geometry boundary intersection implementation ############*/

/**
* @brief Calculates distance to the intersections of a cylinder with
*        radius r (centered on the z axis) and ray described by initial
*        position and propagation direction.
*
* @param[in] r Radius of the cylinder.
* @param[in] pos Initial position of the ray.
* @param[in] dir Propagation direction of the ray (norm(dir) = 1).
* @param[in] d1 Distance to the first intersection.
* @param[in] d2 Distance to the second intersection.
*
* @return Returns nonzero if intersection exists (either positive or
*         negative distance to intersection).
*
* @details Intersection between a cylinder x^2 + y^2 = r^2 and a ray
*			T = pos + d*dir = (x0, y0, z0) + d*(px, py, pz) yields a
*			quadratic equation for d:
*			t^2*(px^2 + py^2) + 2*t*(x0*px + y0*py) + x0^2 + y0^2 - r^2.
*			Solution is found as d = (-b * sqrt(b^2 - 4*a*c))/2a, where:
*				a = px^2 + py^2,
*				b = 2*(x0*px + y0*py),
*				c = x0^2 + y0^2 - r^2.
*			If a=0 or b^2 - 4*a*c < 0 there is no intersection.
*/
inline int ray_cylinder_intersections(
		mc_fp_t r, mc_point3f_t const *pos, mc_point3f_t const *dir,
		mc_fp_t *d1, mc_fp_t *d2){

	mc_fp_t a = dir->x*dir->x + dir->y*dir->y;
	mc_fp_t b = FP_2*(pos->x*dir->x + pos->y*dir->y);
	/* mc_fp_t c = pos->x*pos->x + pos->y*pos->y - r*r; */
	mc_fp_t D = b*b - FP_4*a*(pos->x*pos->x + pos->y*pos->y - r*r);

	if(a != FP_0 && D > FP_0){
		D = mc_sqrt(D);
		*d1 = mc_fdiv(-b - D, FP_2*a);
		*d2 = mc_fdiv(-b + D, FP_2*a);
		return 1;
	}
	return 0;
};

/**
* @brief Calculates distance to the closest boundary of the current layer.
*
* @param[in] psim Pointer to a simulator instance.
* @param[out] d_inner Distance to the inner boundary intersection.
* @param[out] d_outer Distance to the outer boundary intersection.
*
* @return Return distance to the closest intersection.
*
* @details Intersection between a cylinder x^2 + y^2 = r^2 and a ray
*			T = pos + d*dir = (x0, y0, z0) + d*(px, py, pz) yields a
*			quadratic equation for d:
*			t^2*(px^2 + py^2) + 2*t*(x0*px + y0*py) + x0^2 + y0^2 - r^2.
*			Solution is found as d = (-b * sqrt(b^2 - 4*a*c))/2a, where:
*				a = px^2 + py^2,
*				b = 2*(x0*px + y0*py),
*				c = x0^2 + y0^2 - r^2.
*			If a=0 or b^2 - 4*a*c < 0 there is no intersection.
*/
inline mc_fp_t mcsim_distance_to_boundary(
		McSim *psim, mc_fp_t *d_inner, mc_fp_t *d_outer){
	mc_point3f_t const *pos = mcsim_position(psim);
	mc_point3f_t const *dir = mcsim_direction(psim);

	/* quadratic equation for distance to intersection: a*d^2 + b*d + c = 0 */
	mc_fp_t inv_2a;
	mc_fp_t a = dir->x*dir->x + dir->y*dir->y;
	mc_fp_t b = FP_2*(pos->x*dir->x + pos->y*dir->y);
	mc_fp_t c = pos->x*pos->x + pos->y*pos->y; /* - r*r; */
	mc_fp_t D;
	mc_fp_t d1, d2;
	*d_outer = FP_INF; /* intersection at infinity - sign is irrelevant */
	*d_inner = FP_INF; /* intersection at infinity - sign is irrelevant */

	/*
	 * Note that b is also the scalar/dot product of the radial normal vector
	 * (length is not normalized to 1) with the propagation direction vector.
	 * The inner boundary can be intersected ony if b is negative!
	*/
	inv_2a = mc_fdiv(FP_1, FP_2*a);

	/* Inner boundary intersection:
	 * If a <= 0 then x and y components of the propagation vector are 0 and
	 * there is no intersection with the boundary. If D < 0, the photon
	 * packet does not hit the boundary. If D == 0 then the photon packet
	 * touches the boundary - considered as intersection.
	*/
	/* if( a > FP_0){ This is already checked in the kernel ... */
	D = b*b - FP_4*a*
		(c - mc_fsquare(mc_layer_r_inner(mcsim_current_layer(psim))));
	#if MC_ENABLE_DEBUG
		dbg_print("Computing boundary intersections:");
		dbg_print_float(INDENT "Discriminant of the inner boundary:", D);
	#endif

	/* intersection with the inner boundary is only possible if the
	 *  boundary radius is nonzero (not the innermost layer ) */
	if (mc_layer_r_inner(mcsim_current_layer(psim)) > FP_0 && D > FP_0){
		D = mc_sqrt(D);
		d1 = (-b - D)*inv_2a;
		d2 = (-b + D)*inv_2a;
		dbg_print_float(INDENT "Inner boundary d1:", d1);
		dbg_print_float(INDENT "Inner boundary d2:", d2);
		*d_inner = (d2 > FP_2*FP_EPS) ? mc_fmax(d1, FP_0) : FP_INF;
	}

	/* Outer boundary intersection */
	D = b*b - FP_4*a*(c -
		mc_fsquare(mc_layer_r_outer(mcsim_current_layer(psim))));
	#if MC_ENABLE_DEBUG
		dbg_print_float(INDENT "Discriminant of the outer boundary:", D);
	#endif

	if (D >= FP_0){
		D = mc_sqrt(D);
		d1 = (-b - D)*inv_2a;
		d2 = (-b + D)*inv_2a;
		dbg_print_float(INDENT "Outer boundary d1:", d1);
		dbg_print_float(INDENT "Outer boundary d2:", d2);
		*d_outer = mc_fmax(d2, FP_0);
	};

	return mc_fmin(*d_outer, *d_inner);
};

/**
 * @brief Computes radial normal vector at the given position. The normal
 *        points outwards.
 *
 * @param[in] pos      Position.
 * @param[out] normal  Filled with components of the radial normal on return.
 *
 * @return Pointer to the initialized radial normal vector.
 */
inline mc_point3f_t *radial_normal(const mc_point3f_t *pos, mc_point3f_t *normal){
	mc_fp_t k = mc_sqrt(pos->x*pos->x + pos->y*pos->y);
	k = (k > FP_0) ? mc_fdiv(FP_1, k) : FP_0;
	normal->x = pos->x*k;
	normal->y = pos->y*k;
	normal->z = FP_0;

	return normal;
};

/**
 * @brief Computes radial normal vector at the given position. The normal
 *        points inwards.
 *
 * @param[in] pos      Position.
 * @param[out] normal  Filled with components of the radial normal on return.
 *
 * @return Pointer to the initialized radial normal vector.
 */
inline mc_point3f_t *radial_normal_inv(const mc_point3f_t *pos, mc_point3f_t *normal){
	mc_fp_t k = mc_sqrt(pos->x*pos->x + pos->y*pos->y);
	k = (k > FP_0) ? mc_fdiv(FP_1, k) : FP_0;
	normal->x = -pos->x*k;
	normal->y = -pos->y*k;
	normal->z = FP_0;

	return normal;
};
/*############ End geometry boundary intersection implementation #############*/

/*############## Start layer boundary handler implementation #################*/
inline mc_int_t mcsim_boundary(McSim *psim, mc_int_t nextLayerIndex){

	mc_fp_t cos_critical, cos1, sin1, cos2, sin2;
	mc_fp_t n_cos1, n_cos2;
	mc_fp_t n2, n1_d_n2, R, Rs, Rp;
	mc_point3f_t normal;
	mc_point3f_t *dir = mcsim_direction(psim);

	__constant const McLayer *currentLayer = mcsim_current_layer(psim);
	__constant const McLayer *nextLayer = mcsim_layer(psim, nextLayerIndex);

	/* normal must point in the propagation direction; normal*dir = cos1 is positive */
	/* get the critical cosine an normal for the boundary */
	if (mcsim_current_layer_index(psim) < nextLayerIndex){
		/* intersecting inner layer boundary */
		radial_normal_inv(mcsim_position(psim), &normal);
		cos_critical = mc_layer_cc_inner(currentLayer);
	}else{
		/* intersecting outer layer boundary */
		radial_normal(mcsim_position(psim), &normal);
		cos_critical = mc_layer_cc_outer(currentLayer);
	}
	dbg_print_float("\nNormal x:", normal.x);
	dbg_print_float("Normal y:", normal.y);
	dbg_print_int("Current layer index:", mcsim_current_layer_index(psim));
	dbg_print_int("Next layer index   :", nextLayerIndex);

	n2 = mc_layer_n(nextLayer);
	mc_fp_t n1 = mc_layer_n(currentLayer);

	#if MC_USE_OUTER_SURFACE_LAYOUT
	/* advance layout at the outer sample surface */
	if(nextLayerIndex == 0)
	{
		int outer_res = mcsim_outer_surface_layout_handler(
			psim, 0, &n2, &cos_critical);
		if(outer_res != MC_SURFACE_LAYOUT_CONTINUE){
			/* boundary handling done/completed by the outer layout handler */
			return outer_res;
		}
	}else
	#endif
	#if MC_USE_INTERNAL_SURFACE_LAYOUT
	if(nextLayerIndex >= 1)
	{
		int internal_res = mcsim_internal_surface_layout_handler(
			psim, nextLayerIndex, &n2, &cos_critical);
		if (internal_res != MC_SURFACE_LAYOUT_CONTINUE){
			/* boundary handling done/completed by the inner layout handler */
			return internal_res;
		}
	}
	else
	#endif

	dbg_print_float("n1", n1);
	dbg_print_float("n2", n2);

	/* refractive indices match*/
	if(n1 == n2){
		mcsim_set_current_layer_index(psim, nextLayerIndex);
		return MC_REFRACTED;
	}

	cos1 = mc_fmin(
		mc_fabs(normal.x*mcsim_direction_x(psim) +
				normal.y*mcsim_direction_y(psim)),
		FP_1
	);
	dbg_print_float("cos1", cos1);
	dbg_print_float("cos_critical", cos_critical);

	/* check if under the critical angle - otherwise reflection holds */
	if(cos1 > cos_critical){
		/* under the critical cosine ... reflect or refract */

		/* under the critical cosine ... reflect or refract */

		/* calculate reflectance */
		n1_d_n2 = mc_fdiv(mc_layer_n(currentLayer), n2);

		sin1 = mc_sqrt(FP_1 - cos1*cos1);
		if(cos1 >= FP_COS_0)
			sin1 = FP_0;

		sin2 = mc_fmin(FP_1, n1_d_n2*sin1);
		cos2 = mc_sqrt(FP_1 - sin2*sin2);

		#if MC_ENABLE_DEBUG
			printf("At the boundary: sin1 %.9f; sin2 %.9f, cos2 %.9f\n",
				sin1, sin2, cos2);
		#endif

		/*
		Fresnel's equations
		Rs = (n*cos1 - cos2)^2/(n*cos1 + cos2)^2
		Rp = (n*cos2 - cos1)^2/(n*cos2 + cos1)^2
		n = n1/n2
		R = 1/2*(Rs + Rp)
		*/

		n_cos1 = n1_d_n2*cos1;
		n_cos2 = n1_d_n2*cos2;

		Rs = mc_fdiv(n_cos1 - cos2, n_cos1 + cos2);
		Rs *= Rs;
		Rp = mc_fdiv(n_cos2 - cos1, n_cos2 + cos1);
		Rp *= Rp;

		R = FP_0p5*(Rp + Rs);

		/* Handle special case cos2 == 0 (90 deg) */
		if (cos1 <= FP_COS_90 || sin2 == FP_1)
			R = FP_1;

		#if MC_ENABLE_DEBUG
			printf("At the boundary: reflectance: %.9f\n", R);
		#endif

		if (R < mcsim_random(psim)){
			/* we have a refraction */
			mcsim_set_current_layer_index(psim, nextLayerIndex);

			dbg_print("\nRefracted");
			dbg_print_float("n1_d_n2", n1_d_n2);

			/* dir_f = n1/n2*dir + (n1/n2*cos1 - cos2)*normal */
			mc_fp_t k = n1_d_n2*cos1 - cos2;
			dir->x = n1_d_n2*dir->x - k*normal.x;
			dir->y = n1_d_n2*dir->y - k*normal.y;
			dir->z = n1_d_n2*dir->z; /* nz = 0 */

			return MC_REFRACTED;
		}else{
			/* reflect */
			dbg_print("\nReflected");

			/* reflection dir_r = dir - 2*(dir*normal)*normal */
			dir->x = dir->x - FP_2*cos1*normal.x;
			dir->y = dir->y - FP_2*cos1*normal.y;

		}
	} else {
		/* exceeding critical cosine ... reflect */
		/* reflection dir_r = dir - 2*(dir*normal)*normal */
		dir->x = dir->x - FP_2*cos1*normal.x;
		dir->y = dir->y - FP_2*cos1*normal.y;

	}

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
 * @param sim      Simulator instance.
 * @param deposit  Weight to deposit.
 */
inline void mcsim_fluence_deposit_weight(
		McSim *psim, mc_point3f_t const *pos, mc_fp_t deposit){
	#if MC_FLUENCE_MODE_RATE
		mcsim_fluence_deposit(psim, pos, deposit,
			mc_layer_mua(mcsim_current_layer(psim), mcsim_direction(psim))
		);
	#else
		mcsim_fluence_deposit(psim, pos, deposit);
	#endif
};

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


/*################ Start photon packet source implementation #################*/
/* User-defined photon packet source implementation goes here - DO NOT EDIT! */
/* START_SOURCE_IMPLEMENTATION_BLOCK */
{{ source.implementation  or '' }}
/* END_SOURCE_IMPLEMENTATION_BLOCK */
/*################# End photon packet source implementation ##################*/


/*############# Start scattering phase function implementation ###############*/
/* User-defined scattering phase function implementation goes here -
   DO NOT EDIT! */
/* START_PF_IMPLEMENTATION_BLOCK */
{{ pf.implementation or '' }}
/* END_PF_IMPLEMENTATION_BLOCK */
/*############## End scattering phase function implementation ###############*/


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

	uint32_t num_layers,
	__constant McLayer const *layers,

	__mc_source_mem McSource const *source,
	__mc_surface_mem McSurfaceLayouts const *surface_layouts,

	__mc_trace_mem const McTrace *trace,

	__mc_fluence_mem McFluence const *fluence,

	__mc_detector_mem McDetectors const *detectors,

	__mc_fp_lut_mem mc_fp_t const *fp_lut_array,

	__global mc_int_t *integer_buffer,
	__global mc_fp_t *float_buffer,
	__global mc_accu_t *accumulator_buffer
)
{
	mc_fp_t step, deposit;
	mc_int_t nexLayerIndex;
	bool done = false;

	mc_point3f_t src_pos = source->position;

	/* create and initialize the simulator structure */
	McSim sim = {
		{						/* McSimState state: simulator state */
			{FP_0, FP_0, FP_0}			/* mc_point3f_t position: Current photon packet position. */
			,{FP_0, FP_0, FP_0}			/* mc_point3f_t direction: Current photon packet propagation direction. */
			,rng_state_x[get_global_id(0)]		/* Random generator state X (changes on each mcsim_random call). */
			,rng_state_a[get_global_id(0)]		/* Random generator state A (not changed by mcsim_random call). */
			,FP_1								/* mc_fp_t weight: Current photon packet weight from [0.0, 1.0]. */
			,0									/* McUint: Absolute photon packet index. */
			,0 									/* mc_int_t layer_index: Current layer index. */
			#if MC_USE_EVENTS || defined(__DOXYGEN__)
			,0		/* mc_uint_t event_flags: All events that the packet underwent during this step. */
			#endif
			#if MC_TRACK_OPTICAL_PATHLENGTH || defined(__DOXYGEN__)
			,FP_0	/* mc_fp_t optical_pathlength: Optical pathlength traveled by the photon packet. */
			#endif
			#if MC_USE_TRACE || defined(__DOXYGEN__)
				,0	/* mc_uint_t trace_count: Number of traced events. */
			#endif
			#if MC_USE_FLUENCE && MC_USE_FLUENCE_CACHE
			,mc_accucache_initializer	/* mc_accucache_t fluence_cache: Fluence cache object. */
			#endif
		},

		num_layers,					/* mc_int_t num_layers: Number of layers including the two outermost layers. */
		layers, 					/* __constant McLayer *layers: Layer objects. */

		source						/* __constant McSource *source: Photon packet source object. */

		#if MC_USE_SURFACE_LAYOUTS
			,surface_layouts		/* __mc_surface_mem McSurfaceLayouts - Advanced layout at the sample top and bottom surfaces. */
		#endif

		#if MC_USE_FP_LUT
			,fp_lut_array				/* __constant mc_fp_t *fp_lut_array: Lookup table(s) data. */
		#endif

		#if MC_USE_TRACE
			,trace					/* __mc_trace_mem McTrace: Trace configuration object. */
		#endif

		#if MC_USE_FLUENCE
			,fluence				/* __mc_fluence_mem const McFluence *fluence; Fluence configuration struct */
		#endif

		#if MC_USE_DETECTORS
			,detectors				/* McDetectors : Surface and specular detectors. */
		#endif

		,integer_buffer		/* __global mc_int_t *fp_buffer: Common integer buffer. */
		,float_buffer		/* __global mc_fp_t *fp_buffer: Common floating-point buffer. */
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
	if ((sim.state.photon_index = pkt_cnt_atomic_inc(num_packets_done)) < num_packets) {

		atomic_inc(num_kernels);

		#if MC_ENABLE_DEBUG
		/* print debug information ony for the first photon packet */
		if (sim.state.photon_index == 0){
			dbg_print("");
			dbg_print("OpenCL runtime:");
			dbg_printf(INDENT     "supported version:  %.1f\n",
				__OPENCL_VERSION__/FP_LITERAL(100.0));
			dbg_printf(INDENT     "built for version:  %.1f\n",
				__OPENCL_C_VERSION__/FP_LITERAL(100.0));

			dbg_print("");
			dbg_print("Simulation:");
			dbg_print_cnt(INDENT   "num_packets     :", num_packets);
			dbg_print_cnt(INDENT   "num_packets_done:", *num_packets_done);
			dbg_print_uint(INDENT  "num_layers      :", num_layers);
			dbg_print_float(INDENT "mc_rmax (mm)    :", mc_rmax*FP_LITERAL(1e3));

			dbg_print("");
			for (int i=0; i < sim.num_layers; ++i){
				dbg_print_int("Layer:", i);
				dbg_print_layer(&sim.layers[i], INDENT);
			};

			dbg_print(""); dbg_print_source(source);

			#if MC_USE_SURFACE_LAYOUTS
				dbg_print(""); dbg_print_surface_layouts(surface_layouts);
			#endif
			#if MC_USE_TRACE
				dbg_print(""); dbg_print_trace(trace);
			#endif
			#if MC_USE_FLUENCE
				dbg_print(""); dbg_print_fluence(fluence);
			#endif
			#if MC_USE_DETECTORS
				dbg_print(""); dbg_print_detectors(detectors);
			#endif
			printf("Test:\n\t(float)0x100000000UL: %.9f\n"
					"\n\tInitial rng states X: %lu A: %u"
					"\n\tTest mcsim_random() %.9f"
					"\n\tNew rng states X: %lu A: %u\n",
					(float)0x100000000,
					sim.state.rngStateX, sim.state.rngStateA,
					mcsim_random(&sim),
					sim.state.rngStateX, sim.state.rngStateA);

			dbg_print("#########################################");
		};
		#endif

		/* launch a new photon packet */
		mcsim_launch(&sim);
		mcsim_event_flags_add(&sim, MC_EVENT_PACKET_LAUNCH);

		#if MC_USE_TRACE & MC_USE_TRACE_START
			/* initial photon packet state */
			mcsim_trace_this_event(&sim);
		#endif

		/* loop through the simulation steps until all the photon packets
			have been processed */
		int num_steps = 0;
		while (!done && num_steps++ < 1000000) {

			/* generate a new step */
			//step = -mc_log(mc_fclip(mcsim_random(&sim), FP_EPS, FP_1 - FP_EPS))*
			#if MC_METHOD == MICROSCOPIC_BEER_LAMBERT
				step = mc_fdiv(-mc_log(mcsim_random(&sim)),
					mc_layer_mus(mcsim_current_layer(&sim), mcsim_direction(&sim)));
			#else
				step = -mc_log(mcsim_random(&sim))*
					mc_layer_inv_mut(mcsim_current_layer(&sim), mcsim_direction(&sim));
			#endif
			step = mc_fmin(step, FP_MAX);

			/* initialize the next layer index with the current layer index */
			nexLayerIndex = mcsim_current_layer_index(&sim);

			dbg_print_status(&sim, "Packet start step");
			dbg_print_float(INDENT "Step taken:", step);
			dbg_print_int(INDENT "Current layer index:", mcsim_current_layer_index(&sim));
			dbg_print("\n");

			/* intersection with boundaries possible only if photon packet is
				not moving in the direction of the z axis */
			if  ((mcsim_direction_x(&sim) != FP_0) ||
					(mcsim_direction_y(&sim) != FP_0)){

				mc_fp_t d, d_inner, d_outer;
				d = mcsim_distance_to_boundary(&sim, &d_inner, &d_outer);
				dbg_print_float("Distance to inner boundary:", d_inner);
				dbg_print_float("Distance to outer boundary:", d_outer);
				dbg_print_float("Free distance d:", d);

				nexLayerIndex = (step > d) ?
					(nexLayerIndex + ((d_inner <= d_outer) ? 1 : -1)) :
					nexLayerIndex;

				/* limit the step size if required */
				step = mc_fmin(d, step);
			}

			/* move the photon packet by the estimated step */
			mcsim_set_position_coordinates(&sim,
				mcsim_position_x(&sim) + mcsim_direction_x(&sim)*step,
				mcsim_position_y(&sim) + mcsim_direction_y(&sim)*step,
				mcsim_position_z(&sim) + mcsim_direction_z(&sim)*step
			);

			/* update total optical pathlength of the photon packet*/
			#if MC_TRACK_OPTICAL_PATHLENGTH
				mcim_optical_pathlength_add(
					&sim, mc_layer_n(mcsim_current_layer(&sim))*step);
			#endif

			#if MC_METHOD == MICROSCOPIC_BEER_LAMBERT
				mc_fp_t mua = mc_layer_mua(mcsim_current_layer(&sim), mcsim_direction(&sim));
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
				if (nexLayerIndex != mcsim_current_layer_index(&sim)) {
					/* deal with the boundary hit ... returns MC_REFRACTED or MC_REFLECTED */
					dbg_print_status(&sim, "\nEntering mcsim_boundary:");
					mc_uint_t _boundary_flags = mcsim_boundary(&sim, nexLayerIndex);
					mcsim_event_flags_add(&sim, _boundary_flags | MC_EVENT_BOUNDARY_HIT);
					(void)_boundary_flags;
					dbg_print_status(&sim, "Leaving mcsim_boundary:");

					/* handle the sample surface detector */
					if (mcsim_packet_escaped_sample(&sim)) {
						dbg_print("Packet escaped sample.");
						#if MC_USE_DETECTORS
							if (mcsim_current_layer_index(&sim) <= 0){
								#if MC_USE_OUTER_DETECTOR
									mcsim_outer_detector_deposit(
										&sim,
										mcsim_position(&sim),
										mcsim_direction(&sim),
										mcsim_weight(&sim)
									);
								#endif
							}
						#endif

						/* check if the packet escaped the sample */
						done = true;
					}
				} else {
					mcsim_scatter(&sim);
					mcsim_event_flags_add(&sim, MC_EVENT_PACKET_SCATTERING);
				}

				/* Perform survival lottery if required (packet not done). */
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

			#else /* Albedo weight or Albedo Rejection */
				/* process boundary hit or handle absorption */
				if (nexLayerIndex != mcsim_current_layer_index(&sim)) {
					/* deal with the boundary hit ... returns MC_REFRACTED or MC_REFLECTED */
					dbg_print_status(&sim, "\nEntering mcsim_boundary:");
					mc_uint_t _boundary_flags = mcsim_boundary(&sim, nexLayerIndex);
					mcsim_event_flags_add(&sim, _boundary_flags | MC_EVENT_BOUNDARY_HIT);
					(void)_boundary_flags;
					dbg_print_status(&sim, "Leaving mcsim_boundary:");

					/* handle the sample surface detector */
					if (mcsim_packet_escaped_sample(&sim)) {
						dbg_print("Packet escaped sample.");
						#if MC_USE_DETECTORS
							if (mcsim_current_layer_index(&sim) <= 0){
								#if MC_USE_OUTER_DETECTOR
									mcsim_outer_detector_deposit(
										&sim,
										mcsim_position(&sim),
										mcsim_direction(&sim),
										mcsim_weight(&sim)
									);
								#endif
							}
						#endif

						/* check if the packet escaped the sample */
						done = true;
					}

				} else {
					#if MC_METHOD == ALBEDO_REJECTION
						/* Do absorption or scattering only when no layer boundary has been hit.*/
						if (mcsim_random(&sim) < mc_layer_mua_inv_mut(mcsim_current_layer(&sim)), mcsim_direction(&sim)){
							/* Deposit the entire weight of the packet. */
							deposit = mcsim_weight(&sim);
							done = true;
							mcsim_adjust_weight(&sim, deposit);
							mcsim_event_flags_add(&sim, MC_EVENT_PACKET_ABSORPTION);
							#if MC_USE_FLUENCE
								mcsim_fluence_deposit_weight(
									&sim, mcsim_position(&sim), deposit);
							#endif
						} else {
							/* Scatter the photon packet. */
							mcsim_scatter(&sim);
							mcsim_event_flags_add(&sim, MC_EVENT_PACKET_SCATTERING);
						}
					#else /* Albedo Weight */
						/* Do absorption only when no layer boundary has been hit.*/
						deposit = mcsim_weight(&sim)*
							mc_layer_mua_inv_mut(mcsim_current_layer(&sim), mcsim_direction(&sim));
						mcsim_adjust_weight(&sim, deposit);
						mcsim_event_flags_add(&sim, MC_EVENT_PACKET_ABSORPTION);
						#if MC_USE_FLUENCE
							mcsim_fluence_deposit_weight(
								&sim, mcsim_position(&sim), deposit);
						#endif

						dbg_print_float("\nDepositing", deposit);

						/* Scatter the photon packet. */
						mcsim_scatter(&sim);
						mcsim_event_flags_add(&sim, MC_EVENT_PACKET_SCATTERING);
						dbg_print_status(&sim, "\nScattered");

						/* Perform survival lottery if required (packet not done). */
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

			mc_point3f_t *dir = mcsim_direction(&sim);
			mc_fp_t dir_len = mc_sqrt(dir->x*dir->x + dir->y*dir->y + dir->z*dir->z);
			if ( mc_fabs(dir_len - FP_1) > 10*FP_EPS){
				printf("\nDIRECTION ERROR len=%.8f: (%.8f, %.8f, %.8f)",
					dir_len, dir->x, dir->y, dir->z);
				done = true;
			}
			/* check if photon escaped the predefined simulation domain */
			if (mc_distance2_point3f(mcsim_position(&sim), &src_pos) >
					mc_rmax*mc_rmax || mcsim_weight(&sim) <= FP_0) {
				mcsim_event_flags_add(&sim, MC_EVENT_PACKET_ESCAPED);
				done = true;
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
					/* Photon packet has escaped the sample layers -
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
					mcsim_launch(&sim);
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

		rng_state_x[get_global_id(0)] = sim.state.rngStateX;
	};

};
/*##################### End Monte Carlo OpenCL kernel ########################*/

__kernel void sizeof_datatypes(__global uint *n){
	if (get_global_id(0) == 0){
		n[0] = sizeof(McLayer);
		n[1] = sizeof(McSource);
		#if MC_USE_SURFACE_LAYOUTS
		n[2] = sizeof(McSurfaceLayouts);
		#endif
		#if MC_USE_DETECTORS
		n[3] = sizeof(McDetectors);
		#endif
		#if MC_USE_TRACE
		n[4] = sizeof(McTrace);
		#endif
		#if MC_USE_FLUENCE
		n[5] = sizeof(McFluence);
		#endif
	}
};
