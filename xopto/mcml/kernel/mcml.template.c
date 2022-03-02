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

#include "mcml.template.h"


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

/*############## Start layer boundary handler implementation #################*/
inline mc_int_t mcsim_boundary(McSim *psim, mc_int_t nextLayerIndex){

	mc_fp_t cos_critical, cos1, sin1, cos2, sin2;
	mc_fp_t n_cos1, n_cos2;
	mc_fp_t n2, n1_d_n2, R, Rs, Rp;

	#if MC_GEOMETRY_EX
		mc_int_t exNextMaterialIndex;
	#endif

	__constant const McLayer *currentLayer = mcsim_current_layer(psim);
	__constant const McLayer *nextLayer = mcsim_next_layer(psim);

	mc_point3f_t *dir = mcsim_direction(psim);


	/* get the critical cosine for the boundary */
	if (dir->z < FP_0)
		cos_critical = mc_layer_cc_top(currentLayer);
	else
		cos_critical = mc_layer_cc_bottom(currentLayer);

	n2 = mc_layer_n(nextLayer);


	#if MC_USE_TOP_SURFACE_LAYOUT
	/* advance layout at the top sample surface */
	if(nextLayerIndex == 0)
	{
		int top_res = mcsim_top_surface_layout_handler(psim, &n2, &cos_critical);
		if(top_res != MC_SURFACE_LAYOUT_CONTINUE){
			/* boundary handling done/completed by the top layout handler */
			return top_res;
		}
	}else
	#endif
	#if MC_USE_BOTTOM_SURFACE_LAYOUT
	if(nextLayerIndex == mcsim_layer_count(psim) - 1)
	{
		int bottom_res = mcsim_bottom_surface_layout_handler(psim, &n2, &cos_critical);
		if (bottom_res != MC_SURFACE_LAYOUT_CONTINUE){
			/* boundary handling done/completed by the top layout handler */
			return bottom_res;
		}
	}
	else
	#endif

	/* refractive indices match*/
	if(mc_layer_n(currentLayer) == n2){
		mcsim_set_current_layer_index(psim, nextLayerIndex);
		return MC_REFRACTED;
	}

	/* reflect or refract */
	/* default reflect */
	dir->z = -dir->z;

	cos1 = mc_fabs(dir->z);

	/* check if under the critical angle - otherwise reflection holds */	
	if(cos1 > cos_critical){

		/* normal incidence */
		if (cos1 > FP_COS_0){
			R = mc_fdiv((mc_layer_n(currentLayer) - n2),
				(mc_layer_n(currentLayer) + n2));
			if (R*R < mcsim_random(psim)){
				/* refraction */
				mcsim_set_current_layer_index(psim, nextLayerIndex);
				/* change back from default reflect */
				dir->z = -dir->z;
				return MC_REFRACTED;
			}else{
				return MC_REFLECTED;
			};
		};

		n1_d_n2 = mc_fdiv(mc_layer_n(currentLayer), n2);

		sin1 = mc_sqrt(FP_1 - cos1*cos1);
		if(cos1 >= FP_COS_0) 
			sin1 = FP_0;

		sin2 = mc_fmin(FP_1, n1_d_n2*sin1);
		cos2 = mc_sqrt(FP_1 - sin2*sin2);

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

		R = FP_0p5*(Rs + Rp);

		// Handle special case cos2 == 0 (90 deg) - special if statement used 
		// in previous implementation
		if (cos1 <= FP_COS_90 || sin2 == FP_1) 
			R = FP_1;

		if (R < mcsim_random(psim)){ /* [1.0, 0.0) */
			/* we have a refraction */
			mcsim_set_current_layer_index(psim, nextLayerIndex);

			dir->x *= n1_d_n2;
			dir->y *= n1_d_n2;
			dir->z = -mc_fcopysign(cos2, dir->z);

			return MC_REFRACTED;
		};
	};

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
		mcsim_fluence_deposit(psim, deposit, pos,
			mc_layer_mua(mcsim_current_layer(psim))
		);
	#else
		mcsim_fluence_deposit(&psim, pos, deposit);
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
	mc_fp_t fi, cos_theta;

	/* sample the scattering phase functions */
	cos_theta = mcsim_sample_pf(psim, &fi);

	scatter_direction(mcsim_direction(psim), cos_theta, fi);
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

	/* prepare a simulation structure - the part that does not change between simulations */
	McSim sim = {
		{
			{FP_0, FP_0, FP_0}			/* mc_point3f_t position: Current photon packet position. */
			,{FP_0, FP_0, FP_0}			/* mc_point3f_t direction: Current photon packet propagation direction. */
			,rng_state_x[get_global_id(0)]		/* Random generator state X (changes on each mcsim_random call). */
			,rng_state_a[get_global_id(0)]		/* Random generator state A (not changed by mcsim_random call). */
			,FP_1								/* mc_fp_t weight: Current photon packet weight from [0.0, 1.0]. */
			,0									/* McUint: Absolute photon packet index. */
			,0 									/* mc_int_t layer_index: Current layer index. */
			#if MC_TRACK_OPTICAL_PATHLENGTH || defined(__DOXYGEN__)
			,FP_0	/* mc_fp_t optical_pathlength: Optical pathlength traveled by the photon packet. */
			#endif
			#if MC_USE_TRACE || defined(__DOXYGEN__)
				,0	/* mc_uint_t trace_count: Number of traced events. */
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

		/* initialize a new photon packet trace */
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

		#if MC_USE_TRACE & MC_USE_TRACE_START
			/* initial photon packet state */
			mcsim_trace_this_event(&sim);
		#endif

		/* loop through the simulation steps until all the photon packets 
			have been processed */
		while (!done) {

			/* generate a new step */
			//step = -mc_log(mc_fclip(mcsim_random(&sim), FP_EPS, FP_1 - FP_EPS))*
			#if MC_METHOD == MICROSCOPIC_BEER_LAMBERT
				step = mc_fdiv(-mc_log(mcsim_random(&sim)),
					mc_layer_mus(mcsim_current_layer(&sim)));
			#else
				step = -mc_log(mcsim_random(&sim))*
					mc_layer_inv_mut(mcsim_current_layer(&sim));
			#endif
			step = mc_fmin(step, FP_INV_EPS);
			//step = -mc_log(mcsim_random(&sim));
			//step = mc_fclip(step, FP_EPS, FP_1 - FP_EPS)*mc_layer_inv_mut(mcsim_current_layer(&sim));
			

			/* initialize the next layer index with the current layer index */ 
			nexLayerIndex = mcsim_current_layer_index(&sim);

			/* check if the photon packet hits the top layer boundary and if so,
				adjust the step size - go only to the layer boundary */
			if (mcsim_position_z(&sim) + step*mcsim_direction_z(&sim) <
					mc_layer_top(mcsim_current_layer(&sim))) {
				--nexLayerIndex;
				if (mc_fabs(mcsim_direction_z(&sim)) != FP_0) {
					step = mc_fdiv( mc_layer_top(mcsim_current_layer(&sim)) - 
						mcsim_position_z(&sim), mcsim_direction_z(&sim) );
				};
			};
			/* check if the photon packet hits the bottom layer boundary and if so,
				adjust the step size - go only to the layer boundary */
			if (mcsim_position_z(&sim) + 
					step*mcsim_direction_z(&sim) >=
					mc_layer_bottom(mcsim_current_layer(&sim))) {
				++nexLayerIndex;
				if (mc_fabs(mcsim_direction_z(&sim)) != FP_0){
					step = mc_fdiv( mc_layer_bottom(mcsim_current_layer(&sim)) - 
						mcsim_position_z(&sim), mcsim_direction_z(&sim) );
				};
			};
			
			/* move the photon packet by the estimated step */
			mcsim_set_position_coordinates(&sim,
				mcsim_position_x(&sim) + mcsim_direction_x(&sim)*step,
				mcsim_position_y(&sim) + mcsim_direction_y(&sim)*step,
				mcsim_position_z(&sim) + mcsim_direction_z(&sim)*step
			);

			/* update total optical pathlength of the photon packet*/
			#if MC_TRACK_OPTICAL_PATHLENGTH
				mcsim_optical_pathlength_add(
					&sim, mc_layer_n(mcsim_current_layer(&sim))*step);
			#endif

			/* limit the movement to the layer boundaries - FPU precision related */
			if (mcsim_current_layer_index(&sim) < nexLayerIndex)
				mcsim_set_position_z(&sim, 
					mc_layer_bottom(mcsim_current_layer(&sim)));
			if (mcsim_current_layer_index(&sim) > nexLayerIndex)
				mcsim_set_position_z(&sim, 
					mc_layer_top(mcsim_current_layer(&sim)));

			#if MC_METHOD == MICROSCOPIC_BEER_LAMBERT
				mc_fp_t mua = mc_layer_mua(mcsim_current_layer(&sim));
				mc_fp_t deposit_fraction = FP_1 - mc_exp(-mua*step);
				deposit = deposit_fraction*mcsim_weight(&sim);
				mcsim_adjust_weight(&sim, deposit);

				#if MC_USE_FLUENCE
					mc_fp_t step_back = (mua != FP_0) ?
						step - mc_fdiv(-mc_log(FP_1 - mcsim_random(&sim)*deposit_fraction), mua) :
						FP_0;
					//mc_fp_t w_c = FP_0p5*step;
					mc_point3f_t deposit_pos = {
						mcsim_position_x(&sim) - step_back*mcsim_direction_x(&sim),
						mcsim_position_y(&sim) - step_back*mcsim_direction_y(&sim),
						mcsim_position_z(&sim) - step_back*mcsim_direction_z(&sim)
					};

					/* update the fluence data in fluence mode */
					mcsim_fluence_deposit_weight(&sim, &deposit_pos, deposit),
				#endif

				/* process boundary hit or handle absorption */
				if (nexLayerIndex != mcsim_current_layer_index(&sim)) {
					/* deal with the boundary hit ... returns MC_REFRACTED or MC_REFLECTED */
					mcsim_boundary(&sim, nexLayerIndex);

					/* handle the sample surface detector */
					if (mcsim_packet_escaped_sample(&sim)) {
						#if MC_USE_TOP_DETECTOR || MC_USE_BOTTOM_DETECTOR
							if ( (mcsim_current_layer_index(&sim) <= 0) ){
								#if MC_USE_TOP_DETECTOR
									mcsim_top_detector_deposit(
										&sim,
										mcsim_position(&sim),
										mcsim_direction(&sim),
										mcsim_weight(&sim)
									);
								#endif
							} else if (mcsim_current_layer_index(&sim) >=
									mcsim_layer_count(&sim) - 1){
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

						/* check if the packet escaped the sample */
						done = true;
					}

				} else {
					/* Scatter the photon packet. */
						mcsim_scatter(&sim);
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

				}
			#else /* Albedo Weight or Albedo Rejection */
				/* process boundary hit or handle absorption */
				if (nexLayerIndex != mcsim_current_layer_index(&sim)) {
					/* deal with the boundary hit ... returns MC_REFRACTED or MC_REFLECTED */
					mcsim_boundary(&sim, nexLayerIndex);

					/* handle the sample surface detector */
					if (mcsim_packet_escaped_sample(&sim)) {
						#if MC_USE_TOP_DETECTOR || MC_USE_BOTTOM_DETECTOR
							if ( (mcsim_current_layer_index(&sim) <= 0) ){
								#if MC_USE_TOP_DETECTOR
									mcsim_top_detector_deposit(
										&sim,
										mcsim_position(&sim),
										mcsim_direction(&sim),
										mcsim_weight(&sim)
									);
								#endif
							} else if (mcsim_current_layer_index(&sim) >=
									mcsim_layer_count(&sim) - 1){
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

						/* check if the packet escaped the sample */
						done = true;
					}

				} else {
					#if MC_METHOD == ALBEDO_REJECTION
						/* Do absorption or scattering only when no layer boundary has been hit.*/
						if (mcsim_random(&sim) <= mc_layer_mua_inv_mut(mcsim_current_layer(&sim))){
							/* Deposit the entire weight of the packet. */
							deposit = mcsim_weight(&sim);
							done = true;
							mcsim_adjust_weight(&sim, deposit);
							#if MC_USE_FLUENCE
								mcsim_fluence_deposit_weight(
									&sim, mcsim_position(&sim), deposit);
							#endif
						} else {
							/* Scatter the photon packet. */
							mcsim_scatter(&sim);
						}
					#else	/* Albedo Weight */.
						/* Do absorption only when no layer boundary has been hit.*/
						deposit = mcsim_weight(&sim)*
							mc_layer_mua_inv_mut(mcsim_current_layer(&sim));
						mcsim_adjust_weight(&sim, deposit);
						#if MC_USE_FLUENCE
							mcsim_fluence_deposit_weight(
								&sim, mcsim_position(&sim), deposit);
						#endif

						/* Scatter the photon packet. */
						mcsim_scatter(&sim);

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

			/* check if photon escaped the predefined simulation domain */
			if (point3f_distance_squared(mcsim_position(&sim), &src_pos) >
					mc_rmax*mc_rmax) {
				done = true;
			};

			#if MC_USE_TRACE == MC_USE_TRACE_ALL
				mcsim_trace_this_event(&sim);
			#elif MC_USE_TRACE & MC_USE_TRACE_END
				if (done)
					mcsim_trace_this_event(&sim);
			#endif

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

					#if MC_USE_TRACE & MC_USE_TRACE_START
						/* initial photon packet state */
						mcsim_trace_this_event(&sim);
					#endif
					done = false; /* have a new packet ... not done yet */
				}
			}
		};
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
