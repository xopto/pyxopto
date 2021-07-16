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

#ifndef __MC_SV_H
#define __MC_SV_H

#include "mcbase.template.h"

/*#################### Start Sampling Volume Extension #######################*/
/**
*@brief Structure that represents a voxelized sampling volume.
*@{
*/
struct McSamplingVolume {
	/** @brief Coordinates of the top left corner of the voxelized sampling volume. */
	mc_point3f_t top_left;
	/** @brief Size of a voxel. */
	mc_point3f_t voxel_size;
	/** @brief Shape of the voxelized sampling volume data buffer. */
	mc_point3_t shape;
	/** @brief Multiplier that is applied to the product of the packet weight
	 *         and path length traveled in a voxel before converting to int */
	mc_fp_t multiplier;
	/** @brief Offset of the top-left (first) voxel in the accumulator buffer. */
	mc_size_t offset;
	/** @brief Photon packet weight conversion from floating-poiny to integer accumulator type **/
	mc_int_t k;
};
/**
 * @}
 */
/** @brief Data type that represents a voxelized sampling volume. */
typedef struct McSamplingVolume McSamplingVolume;

/**
*@brief Structure that represents the state of the sampling volume analyzer.
*@{
*/
struct McSamplingVolumeAnalyzer{
	__constant McSamplingVolume *sv;		/** < @brief Sampling volume configuration. */
	__constant McTrace *trace;	/** < @brief Trace configuration. */
	mc_cnt_t packet_index;			/** < @brief Current photon packet index. */
	mc_int_t num_events;			/** < @brief Number of events recorded for this photon packet. */
	__global mc_fp_t *trace_data;	/** < @brief Event data buffer for this photon packet. */
};
/**
 * @}
 */
/** @brief Data type that represents the state of the sampling volume analyzer. */
typedef struct McSamplingVolumeAnalyzer McSamplingVolumeAnalyzer;

/**
*@brief Structure that represents a single trace event.
*@{
*/
struct TraceEvent{
	mc_point3f_t	pos;	/** < @brief Photon packet position. */
	mc_point3f_t	dir;	/** < @brief Photon packet propagation direction. */
	mc_fp_t	weight;	/** < @brief Photon packet weight. */
};
/**
 * @}
 */
/** @brief Data type that represents a single trace event. */
typedef struct TraceEvent TraceEvent;

/**
 * @brief Fetches a trace event from the trace data buffer.
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @param[out] ev Pointer to an event instance that fill be filled with the event data.
 * @param[in]  index Event index (0-based).
 * @returns Pointer to the event instance (input argument) filled with event data.
 */
inline TraceEvent *sva_get_event(McSamplingVolumeAnalyzer const *sva,
	TraceEvent *ev, mc_int_t ev_index);

/**
 * @brief Computes the voxel coordinates at the position of the event.
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @param[in] ev Pointer to an event instance.
 * @param[out] voxel Pointer to a 3D integer point that will be filled with the voxel indices.
 * @returns voxel Input argument voxel.
 */ 
inline mc_point3_t *sva_voxel(McSamplingVolumeAnalyzer const *sva,
		TraceEvent const *ev, mc_point3_t *voxel);

/**
 * @brief Computes the voxel coordinates at the position of the event.
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @param[in] ev Pointer to an event instance.
 * @param[out] voxel Pointer to a 3D integer point that will be filled with the voxel indices.
 * @returns voxel Input argument voxel.
 */ 
inline mc_point3_t *sva_voxel(McSamplingVolumeAnalyzer const *sva,
		TraceEvent const *ev, mc_point3_t *voxel);

/**
 * @brief Return nonzero if the given voxel is valid / within the voxelized
 *        sampling volume domain.
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @param[out] voxel Voxel to test.
 * @returns Nonzero if the voxel lies within the sampling volume domain.
 */ 
inline int sva_is_valid_voxel(McSamplingVolumeAnalyzer const *sva,
		mc_point3_t const *voxel);

/**
 * @brief Intersect the voxel at the current position.
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @param[in] ev Event at the current position.
 * @param[in] voxel Voxel at the current poosition.
 * @param[out] distances Distances to intersections with the voxel in x, y and z direcyion.
 * @returns Distance from ev1 to the voxel intersection.
 */ 
inline mc_fp_t sva_intersect(
		McSamplingVolumeAnalyzer const *sva, 
		TraceEvent const *ev, mc_point3_t const *voxel,
		mc_point3f_t *distances);

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
		mc_point3_t *voxel);

/**
 * @brief Compute the total path length of the photon packet.
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @returns Total path length of the photon packet.
 */
inline mc_fp_t sva_total_path(McSamplingVolumeAnalyzer const *sva);

/**
 * @brief Compute unsigned 32-bit integer weight that will be deposited to the
 *        sampling volume voxel..
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @param[in] terminal_weight Terminal / final weight of the photon packet.
 * @param[in] l_voxel Length of the path traveled in this voxel.
 * @returns Weight to deposit in this voxel.
 */
inline uint32_t sva_deposit_weight(McSamplingVolumeAnalyzer const *sva,
		mc_fp_t terminal_weight, mc_fp_t l_voxel);

/**
 * @brief Returnflat voxel index into a flat array.
 * @param[in] sva Pointer to a sampling volume analyzer instance.
 * @param[in] voxel X, y and z indices of the voxel.
 * @return Flat index.
 */
inline mc_int_t sva_flat_voxel_index(McSamplingVolumeAnalyzer const *sva,
	mc_point3_t const *voxel);

#if MC_ENABLE_DEBUG
	#define dbg_print_sampling_volume(psv) \
		dbg_print("Sampling volume:"); \
		dbg_print_point3f(INDENT "top_left:", &(psv)->top_left); \
		dbg_print_point3f(INDENT "voxel_size:", &(psv)->voxel_size); \
		dbg_print_point3(INDENT "shape:", &(psv)->shape); \
		dbg_print_float(INDENT "multiplier:", (psv)->multiplier); \
		dbg_print_size_t(INDENT "offset:", (psv)->offset); \
		dbg_print_int(INDENT "k:", (psv)->k);

	#define dbg_print_trace_event(label, pev) \
		dbg_print("Trace event: " label); \
		dbg_print_point3f(INDENT "Position:", &(pev)->pos); \
		dbg_print_point3f(INDENT "Direction:", &(pev)->dir); \
		dbg_print_float(INDENT "Weight:", (pev)->weight);

#else
	#define dbg_print_sampling_volume(psv) ;
	#define dbg_print_trace_event(label, pev) ;
#endif
/*##################### End Sampling Volume Extension ########################*/

#endif /* #define __MC_SV_H */
