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

#include "mcsv.template.h"

/*#################### Start Sampling Volume Extension #######################*/

#if MC_USE_SAMPLING_VOLUME || defined(__DOXYGEN__)

/**
 * @brief Fetches a trace event from the trace data buffer.
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @param[out] ev Pointer to an event instance that fill be filled with the event data.
 * @param[in]  index Event index (0-based).
 * @returns Pointer to the event instance (input argument) filled with event data.
 */
inline TraceEvent *sva_get_event(McSamplingVolumeAnalyzer const *sva,
	TraceEvent *ev, mc_int_t ev_index){

	mc_int_t offset = TRACE_ENTRY_LEN*ev_index;

	ev->pos.x = sva->trace_data[offset++];
	ev->pos.y = sva->trace_data[offset++];
	ev->pos.z = sva->trace_data[offset++];
	ev->dir.x = sva->trace_data[offset++];
	ev->dir.y = sva->trace_data[offset++];
	ev->dir.z = sva->trace_data[offset++];
	ev->weight = sva->trace_data[offset];

	return ev;
};

/**
 * @brief Computes the voxel coordinates at the position of the event.
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @param[in] ev Pointer to an event instance.
 * @param[out] voxel Pointer to a 3D integer point that will be filled with the voxel indices.
 * @returns Input argument voxel.
 */ 
inline mc_point3_t *sva_voxel(McSamplingVolumeAnalyzer const *sva,
		TraceEvent const *ev, mc_point3_t *voxel){
	voxel->x = mc_fdiv((ev->pos.x - sva->sv->top_left.x), sva->sv->voxel_size.x);
	voxel->y = mc_fdiv((ev->pos.y - sva->sv->top_left.y), sva->sv->voxel_size.y);
	voxel->z = mc_fdiv((ev->pos.z - sva->sv->top_left.z), sva->sv->voxel_size.z);

	return voxel;
};

/**
 * @brief Return nonzero if the given voxel is valid / within the voxelized
 *        sampling volume domain.
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @param[out] voxel Voxel to test.
 * @returns Nonzero if the voxel lies within the sampling volume domain.
 */ 
inline int sva_is_valid_voxel(McSamplingVolumeAnalyzer const *sva,
		mc_point3_t const *voxel){

	return voxel->x >= 0 && voxel->y >= 0 && voxel->z >= 0 &&
			voxel->x < sva->sv->shape.x &&
			voxel->y < sva->sv->shape.y &&
			voxel->z < sva->sv->shape.z;
};

/**
 * @brief Intersect the voxel at the current position.
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @param[in] ev Event at the current position.
 * @param[in] voxel Voxel at the current poosition.
 * @param[out] distances Distances to intersections with the voxel in x, y and z direction.
 * @returns Distance from ev1 to the voxel intersection.
 */ 
inline mc_fp_t sva_intersect(
		McSamplingVolumeAnalyzer const *sva, 
		TraceEvent const *ev, mc_point3_t const *voxel,
		mc_point3f_t *distances){

	distances->x = sva->sv->top_left.x +
		((ev->dir.x >= FP_0) + voxel->x)*sva->sv->voxel_size.x - ev->pos.x;

	distances->y = sva->sv->top_left.y +
		((ev->dir.y >= FP_0) + voxel->y)*sva->sv->voxel_size.y - ev->pos.y;

	distances->z = sva->sv->top_left.z +
		((ev->dir.z >= FP_0) + voxel->z)*sva->sv->voxel_size.z - ev->pos.z;

	distances->x = (ev->dir.x != FP_0) ? mc_fdiv(distances->x, ev->dir.x) : FP_INF;
	distances->y = (ev->dir.y != FP_0) ? mc_fdiv(distances->y, ev->dir.y) : FP_INF;
	distances->z = (ev->dir.z != FP_0) ? mc_fdiv(distances->z, ev->dir.z) : FP_INF;

	return mc_fmin(distances->z, mc_fmin(distances->x, distances->y));
};

/**
 * @brief Compute the next voxel index.
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @param[in] ev  Pointer to the current event instance.
 * @param[in] distances Distances to the current voxel boundaries along the x, y and z axis.
 * @param[in, out] voxel Pointer to the current voxel indices.
 * @returns Input argument voxel.
 */
inline mc_point3_t *sva_next_voxel(McSamplingVolumeAnalyzer const *sva,
		TraceEvent const *ev, mc_point3f_t const *distances,
		mc_point3_t *voxel){

	mc_point3_t normal;

	/* find the boundary normal that points outwards */
	normal.x = (distances->x <= distances->y) && (distances->x <= distances->z);
	normal.y = normal.x ?
		0 : (distances->y <= distances->x) && (distances->y <= distances->z);
	normal.z = !(normal.x + normal.y);
	normal.x = ev->dir.x < FP_0 ? -normal.x : normal.x;
	normal.y = ev->dir.y < FP_0 ? -normal.y : normal.y;
	normal.z = ev->dir.z < FP_0 ? -normal.z : normal.z;

	/* compute the next voxel index */
	voxel->x += normal.x;
	voxel->y += normal.y;
	voxel->z += normal.z;

	return voxel;
};

/**
 * @brief Compute the total path length of the photon packet.
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @returns Total path length of the photon packet.
 */
inline mc_fp_t sva_total_path(McSamplingVolumeAnalyzer const *sva){
	mc_int_t i;
	mc_point3f_t pos1, pos2;
	mc_fp_t total_path = FP_0;
	mc_int_t pos = 0;

	pos1.x = sva->trace_data[pos++];
	pos1.y = sva->trace_data[pos++];
	pos1.z = sva->trace_data[pos];

	for(i=1; i < sva->num_events; ++i){
		pos = i*TRACE_ENTRY_LEN;
		pos2.x = sva->trace_data[pos++];
		pos2.y = sva->trace_data[pos++];
		pos2.z = sva->trace_data[pos];
		total_path += point3f_distance(&pos1, &pos2);
		pos1 = pos2;
	};
	return total_path;
};

/**
 * @brief Compute unsigned 32-bit integer weight that will be deposited to the
 *        sampling volume voxel..
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @param[in] terminal_weight Terminal / final weight of the photon packet.
 * @param[in] l_voxel Length of the path traveled in this voxel.
 * @returns Weight to deposit in this voxel.
 */
inline uint32_t sva_deposit_weight(McSamplingVolumeAnalyzer const *sva,
		mc_fp_t terminal_weight, mc_fp_t l_voxel){

	return (uint32_t)(terminal_weight*l_voxel*sva->sv->multiplier*
		sva->sv->k + FP_0p5);
};

/**
 * @brief Returnflat voxel index into a flat array.
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @param[in] voxel X, y and z indices of the voxel.
 * @return Flat index.
 */
inline mc_int_t sva_flat_voxel_index(McSamplingVolumeAnalyzer const *sva,
		mc_point3_t const *voxel){
	return (voxel->z*sva->sv->shape.y + voxel->y)*sva->sv->shape.x + voxel->x;
};

__kernel void SamplingVolume(
	mc_cnt_t npackets,
	__global mc_cnt_t *npackets_processed, 

	__global uint32_t *num_kernels,

	__constant McTrace *trace,
	__constant McSamplingVolume *sampling_volume,

	__global mc_accu_t *total_weight,

	__global mc_int_t *int_buffer,
	__global mc_fp_t *fp_buffer,	
	__global mc_accu_t *accu_buffer
){
	mc_int_t trace_data_offset, trace_length_offset;
	mc_cnt_t packet_index;

	if((packet_index = pkt_cnt_atomic_inc(npackets_processed)) < npackets){

		atomic_inc(num_kernels);

		#if MC_ENABLE_DEBUG
			/* print debug information ony for the first photon packet */
			if (packet_index == 0){
				dbg_print("Sampling volume:");
				dbg_print_cnt("Number of photon packets (traces):", npackets);
				dbg_print_trace(trace);
				dbg_print_sampling_volume(sampling_volume);
			}
		#endif

		/* Initialize the sampling volume satate */
		TraceEvent ev1, ev2;
		mc_int_t event_index, n_valid;
		mc_point3_t voxel;
		mc_point3f_t distances;
		mc_fp_t d_voxel, d_ev1_ev2, d_move, end_weight, l_voxel;
		bool done = false;

		/* offset of the first trace entry for the processed photon packet */
		trace_data_offset = trace->data_buffer_offset +
			packet_index*trace->max_events*TRACE_ENTRY_LEN;
		trace_length_offset = trace->count_buffer_offset + packet_index;

		n_valid = mc_min(trace->max_events, int_buffer[trace_length_offset]);
	
		/* initializing the samping volume state data type */
		McSamplingVolumeAnalyzer sva = {
			sampling_volume,			/* __const McSamplingVolume *sv - sampling volume configuration */
			trace,						/* __const McTrace *trace - trace configuration */
			packet_index,				/* mc_cnt_t packet_index - photon packet index */
			n_valid,					/* mc_int_t num_events - number of events recorded for the photon packet */
			fp_buffer + trace_data_offset,		/* __global mc_fp_t *trace_data - trace event data buffer */
		};

		/* starting with the first event */
		event_index = 0;

		/* compute the total path length of the photon packet */
		// l_total = sva_total_path(&sva);

		/* terminal/end weight of the photon packet */
		end_weight = sva.trace_data[sva.num_events*TRACE_ENTRY_LEN - 1];
		/* save the end weight */
		accumulator_deposit(total_weight,
			(uint32_t)(end_weight*sva.sv->k + FP_0p5));

		/* get the first event data */
		sva_get_event(&sva, &ev1, event_index);

		/* get the second event data */
		sva_get_event(&sva, &ev2, mc_min(event_index + 1, sva.num_events - 1));

		/* compute distance between the events */
		d_ev1_ev2 = point3f_distance(&ev1.pos, &ev2.pos);

		/* get the voxel index at the location of the event */
		sva_voxel(&sva, &ev1, &voxel);

		/* the distance traveled in this voxel */
		l_voxel = FP_0;

		/* start processing the trace events */
		while (!done){
	
			/* get distances to the voxel intersection */
			d_voxel = sva_intersect(&sva, &ev1, &voxel, &distances);

			/* move to the voxel boundaries or to the next event position - whichever is closer */
			d_move = mc_fmin(d_voxel, d_ev1_ev2);

			l_voxel += d_move; /* the distance traveled in this voxel */

			if (d_voxel < d_ev1_ev2){
				/* voxel intersection is closer than the position of the next event */
				
				/* leaving this voxel - update the voxel weight */
				if (sva_is_valid_voxel(&sva, &voxel)){
					mc_int_t index = sva_flat_voxel_index(&sva, &voxel);
					__global mc_accu_t *address = accu_buffer + sva.sv->offset + index;
					uint32_t ui32w = sva_deposit_weight(&sva, end_weight, l_voxel);
					accumulator_deposit(address, ui32w);
				}

				/* update ev1 since it is partially processed */
				ev1.pos.x += ev1.dir.x*d_move;
				ev1.pos.y += ev1.dir.y*d_move;
				ev1.pos.z += ev1.dir.z*d_move;

				/* update the distance to ev2 */
				d_ev1_ev2 -= d_move;

				/* moving to the next voxel - updating the voxel index */
				sva_next_voxel(&sva, &ev1, &distances, &voxel);

				/* reset the distance traveled in this voxel */
				l_voxel = FP_0;
			}else{
				/* voxel intersection is beyond the next event */
				/* ev1 is fully processed */
				ev1 = ev2;
				event_index += 1; /* processing the next event */

				if (event_index < sva.num_events - 1){
					/* fetch the next event */
					sva_get_event(&sva, &ev2, event_index + 1);

					/* compute the distance between the two events */
					d_ev1_ev2 = point3f_distance(&ev1.pos, &ev2.pos);
				} else {
					if (sva_is_valid_voxel(&sva, &voxel)){
						mc_int_t index = sva_flat_voxel_index(&sva, &voxel);
						__global mc_accu_t *address = accu_buffer + sva.sv->offset + index;
						uint32_t ui32w = sva_deposit_weight(&sva, end_weight, l_voxel);
						accumulator_deposit(address, ui32w);
					}

					/* no more events to process - fetch the next trace */
					if((packet_index = pkt_cnt_atomic_inc(npackets_processed)) < npackets){

						trace_data_offset = trace->data_buffer_offset +
							packet_index*trace->max_events*TRACE_ENTRY_LEN;
						trace_length_offset = trace->count_buffer_offset + packet_index;

						n_valid = mc_min(trace->max_events, int_buffer[trace_length_offset]);

						sva.trace_data = fp_buffer + trace_data_offset;
						sva.packet_index = packet_index;
						sva.num_events = n_valid;

						/* starting with the first event */
						event_index = 0;

						/* compute the total path length of the photon packet */
						// l_total = sva_total_path(&sva);

						/* terminal/end weight of the photon packet */
						end_weight = sva.trace_data[sva.num_events*TRACE_ENTRY_LEN - 1];
						/* save the end weight */
						accumulator_deposit(total_weight,
							(uint32_t)(end_weight*sva.sv->k + FP_0p5));

						/* get the first event data */
						sva_get_event(&sva, &ev1, event_index);

						/* get the second event data */
						sva_get_event(&sva, &ev2, mc_min(event_index + 1, sva.num_events - 1));

						/* compute distance between the events */
						d_ev1_ev2 = point3f_distance(&ev1.pos, &ev2.pos);

						/* get the voxel index at the location of the event */
						sva_voxel(&sva, &ev1, &voxel);

						/* the distance traveled in this voxel */
						l_voxel = FP_0;
					}else{
						/* no more packets to process */
						done = true;
					};
				};
			};
		};
	};
};

#endif /* #if MC_USE_SAMPLING_VOLUME */

/*##################### End Sampling Volume Extension ########################*/
