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

#ifndef __MCBASE_H
#define __MCBASE_H

#include "mcbase.template.h"

/*############## Start 64-bit atomic increment implementation ################*/

#if MC_USE_64_BIT_ACCUMULATORS || defined(__DOXYGEN__)
/** 
 * @brief Atomic deposit a 32-bit integer weight to 64-bit accumulator
 * @param[in] address Accumulator addres.
 * @param[in] weight Weight to deposit / add to the accumulator.
 */
inline void accu_64_deposit_32(volatile __global ulong *address, uint weight){
	if (atomic_add((volatile __global uint *)address, weight) + weight < weight)
		atomic_add((volatile __global uint *)address + 1, 1);
};
#endif

#if MC_USE_64_BIT_PACKET_COUNTER || defined(__DOXYGEN__)
/**
 * @brief Software implemantation of 64-bit atomic increment operation using
 *        OpenCL 32-bit atomic operations.
 * @param[in] ptr	Pointer to a 64-bit unsigned integer to be incremented.
 * @return			Returns old value at the pointer location.
 */
inline uint64_t atomic_inc_uint64(__global uint64_t *ptr){
	__global uint32_t *ui32_ptr = (__global uint32_t *)ptr;
	uint32_t low = atomic_inc(ui32_ptr);
	uint32_t high = 0;
	if (low + 1 < low)
		high = atomic_inc(ui32_ptr + 1);
	else
		high = atomic_cmpxchg(ui32_ptr + 1, 0, 0);
	return low + ((uint64_t)(high) << 32);
};
#endif
/*############### End 64-bit atomic increment implementation #################*/


/*################### Start vector/matrix implementation #####################*/
/**
 * @brief Normalizes vector length to unity. 
 * @param[in, out] pv Pointer to a vector normalized on return.
 */	
inline void point3f_normalize(mc_point3f_t *v){
	mc_fp_t k=mc_rsqrt(v->x*v->x + v->y*v->y + v->z*v->z);
	v->x *= k;
	v->y *= k;
	v->z *= k;
};

/**
 * @brief Calculates square of the Euclidean distance between two points.
 * @param[in] pT1 Pointer to the first point.
 * @param[in] pT2 Pointer to the second point.
 * @return Returns the square of Euclidean distance between points T1 and T2.
 */
inline mc_fp_t point2f_distance_squared(
		const mc_point2f_t * pT1, const mc_point2f_t * pT2){
	mc_fp_t d, tmp;
	tmp = pT1->x - pT2->x; d = tmp*tmp;
	tmp = pT1->y - pT2->y; d += tmp*tmp;
	return d;
};

/**
 * @brief Calculates square of the Euclidean distance between two points.
 * @param[in] pT1 Pointer to the first point.
 * @param[in] pT2 Pointer to the second point.
 * @return Returns the square of Euclidean distance between points T1 and T2.
 */
inline mc_fp_t point3f_distance_squared(
		const mc_point3f_t * pT1, const mc_point3f_t * pT2){
	mc_fp_t d, tmp;
	tmp = pT1->x - pT2->x; d = tmp*tmp;
	tmp = pT1->y - pT2->y; d += tmp*tmp;
	tmp = pT1->z - pT2->z; d += tmp*tmp;
	return d;
};

/**
 * @brief Calculates square of the Euclidean distance between two 2D points.
 * @param[in] x1 Coordinate x of the first point.
 * @param[in] y1 Coordinate y of the first point.
 * @param[in] x2 Coordinate x the second point.
 * @param[in] x2 Coordinate y the second point.
 * @return Returns the square of Euclidean distance between (x1, y1) and (x2, y2).
 */
inline mc_fp_t xyf_distance_squared(const mc_fp_t x1, mc_fp_t y1,
		const mc_fp_t x2, const mc_fp_t y2){
	mc_fp_t tmp;
	mc_fp_t d = FP_0;

	tmp = x1 - x2;
	d += tmp*tmp;
	tmp = y1 -y2;
	d += tmp;

	return d;
}

/**
 * @brief Checks if a rectangle contains the given point.
 * @param[in] top_left_x Coordinate x of the top-left edge of the rectangle (mc_fp_t).
 * @param[in] top_left_y Coordinate y of the top-left edge of the rectangle (mc_fp_t).
 * @param[in] width Rectangle width (mc_fp_t).
 * @param[in] height Rectangle height (mc_fp_t).
 * @param[in] x Coordinate x of the point (mc_fp_t).
 * @param[in] y Coordinate y of the point (mc_fp_t).
 * @return Returns nonzero if rectangle contains the given point.
 */
inline int rectf_contains_ex(
		mc_fp_t top_left_x, mc_fp_t top_left_y, mc_fp_t width, mc_fp_t height,
		mc_fp_t x, mc_fp_t y){
	return	(top_left_x <= x) && 
			(top_left_x + width >= x) &&
			(top_left_y <= y) &&
			(top_left_y + height) >= y;
};

/**
 * @brief Checks if a circle contains the given point.
 * @param[in] center_x Coordinate x of the circle center (mc_fp_t).
 * @param[in] center_y Coordinate y of the circle center (mc_fp_t).
 * @param[in] r Circle radius (mc_fp_t).
 * @param[in] x Coordinate x of the point (mc_fp_t).
 * @param[in] y Coordinate y of the point (mc_fp_t).
 * @return Returns nonzero if rectangle contains the given point.
 */
inline int circf_contains_ex(
		mc_fp_t center_x, mc_fp_t center_y, mc_fp_t r, mc_fp_t x, mc_fp_t y){
	mc_fp_t dx = ((center_x) - (x));
	mc_fp_t dy = ((center_y) - (y));

	return dx*dx + dy*dy <= r*r;
}

/**
 * @brief Checks if slot contains a point.
 * @param[in] cx Slot center coordinate x.
 * @param[in] cy Slot center coordinate y.
 * @param[in] width Slot width.
 * @param[in] height Slot height.
 * @param[in] x Coordinate x of the point.
 * @param[in] y Coordinate y of the point.
 * @returns Returns nonzero if the point is within the slot.
 */
inline int slotf_contains_ex(mc_fp_t cx, mc_fp_t cy,
		mc_fp_t width, mc_fp_t height, mc_fp_t x, mc_fp_t y){
	mc_fp_t dx = x - cx;
	mc_fp_t dy = y - cy;
	mc_fp_t c_offset = FP_0p5*mc_fabs(width - height);
	
	mc_fp_t dc = (width >= height) ? 
		(mc_fabs(dx) - c_offset) : (mc_fabs(dy) - c_offset);
	
	return (width >= height) ?
		(
			(mc_fabs(dy) < height*FP_0p5) && (
				(mc_fabs(dx) <= c_offset) || 
				((dy*dy + dc*dc) < height*height*FP_0p25)
			) 
		) :
		(
			(mc_fabs(dx) < width*FP_0p5) && ( 
				(mc_fabs(dy) <= c_offset) ||
				((dx*dx + dc*dc) < width*width*FP_0p25)
			)
		);
}

/*#################### End vector/matrix implementation ######################*/


/*################## Start boundary physics implementation ###################*/
/**
 * @brief Computes reflectance for the given boundary conditions.
 * @param[in] n1 Refractive index of the incident layer.
 * @param[in] n2 Refractive index of the layer across the boundary.
 * @param[in] cosIncidence Incidence angle (with respect to the z axis) cosine.
 * @param[in] cosCritical Critical angle cosine for the interface n1 => n2.
 * @return Returns the calculated reflectance probability from [0.0, 1.0].
 */
inline mc_fp_t reflectance(
		mc_fp_t n1, mc_fp_t n2, mc_fp_t cos1, mc_fp_t cosCritical){
	mc_fp_t Rp, Rs, R = FP_1;
	mc_fp_t n1_d_n2;
	mc_fp_t sin1, sin2, cos2;
	mc_fp_t n_cos1, n_cos2;

	cos1 = mc_fabs(cos1);

	if (n1 == n2)
		return FP_0;

	if(cos1 > cosCritical){
		n1_d_n2 = mc_fdiv(n1, n2);

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

		R = FP_0p5*(Rp + Rs);

		/* Handle special case cos2 == 0 (90 deg) */
		if (cos1 <= FP_COS_90 || sin2 == FP_1) 
			return FP_1;
	}
	return R;
};

/**
 * @brief Computes reflectance for the given boundary conditions by
 * using data from the secondary side.
 * @param[in] n1 Refractive index of the incident layer.
 * @param[in] n2 Refractive index of the layer across the boundary.
 * @param[in] cosIncidence2 Incidence angle (with respect to z axis) cosine.
 * @return Returns the calculated reflectance probability from [0.0, 1.0].
 */
inline mc_fp_t reflectance_cos2(
		mc_fp_t n1, mc_fp_t n2, mc_fp_t cos2){
	mc_fp_t Rp, Rs, R = FP_1;
	mc_fp_t n1_d_n2;
	mc_fp_t sin1, sin2, cos1;
	mc_fp_t n_cos1, n_cos2;

	cos2 = mc_fabs(cos2);

	if (n1 == n2)
		return FP_0;

	sin2 = mc_sqrt(FP_1 - cos2*cos2);
	if(cos2 >= FP_COS_0)
		sin2 = FP_0;

	sin1 = mc_fdiv(n2, n1)*sin2;

	if (sin1 < FP_1){
		cos1 = mc_sqrt(FP_1 - sin1*sin1);

		/*
		Fresnel's equations
		Rs = (n*cos1 - cos2)^2/(n*cos1 + cos2)^2
		Rp = (n*cos2 - cos1)^2/(n*cos2 + cos1)^2
		n = n1/n2
		R = 1/2*(Rs + Rp)
		*/
		n1_d_n2 = mc_fdiv(n1, n2);
		n_cos1 = n1_d_n2*cos1;
		n_cos2 = n1_d_n2*cos2;

		Rs = mc_fdiv(n_cos1 - cos2, n_cos1 + cos2);
		Rs *= Rs;
		Rp = mc_fdiv(n_cos2 - cos1, n_cos2 + cos1);
		Rp *= Rp;

		R = FP_0p5*(Rp + Rs);

		/* Handle special case cos2 == 0 (90 deg) */
		if (cos1 <= FP_COS_90 || sin2 == FP_1) 
			return FP_1;
	}
	return R;
};


/**
 * Compute cosine of the critical incidence angle beyond which the incident
 * beam is reflected at the boundary of materials with refractive indices
 * n1 and n2.
 * @param[in] n1 Refractive index of the incident medium.
 * @param[in] n2 Refractive index of the material across the boundary.
 * @return Cosine of the incident angle beyond which the incident beam is
 *         reflected at the given boundary.
 */
inline mc_fp_t cos_critical(mc_fp_t n1, mc_fp_t n2){
    return (n1 > n2) ? mc_sqrt(FP_1 - mc_fdiv(n2*n2, n1*n1)): FP_0;
};

/**
 * Computes propagation direction of the reflected beam for the given
 * incident propagation direction and boundary normal.
 * @param[in] p Incident propagation direction.
 * @param[in] n Boundary normal (can be pointing outwards or inwards).
 * @param[out] r Reflected beam propagation direction on return.
 *
 * @return Input parameter r.
 * @note Reflected beam propagation direction is computed as p - 2*n*(p*n)
 */
inline mc_point3f_t *reflect(mc_point3f_t const *p, mc_point3f_t const *n, mc_point3f_t *r){
	mc_fp_t p_n_2 = FP_2*dot3f(p, n);
	r->x = p->x - n->x*p_n_2;
	r->y = p->y - n->y*p_n_2;
	r->z = p->z - n->z*p_n_2;

	return r;
}

/**
 * Computes propagation direction of the refracted beam for the given
 * incident propagation direction, boundary normal and refractive indices
 * of the materials on each side of the boundary. Requires signed incident
 * angle cosine computed as n*p.
 * @param[in] p Incident propagation direction.
 * @param[in] n Boundary normal (pointing outwards or inwards, requires signed
 *              cos1 to resolve the surface normal direction!).
 * @param[in] n1 Refractive index of the material on the incident side of the boundary.
 * @param[in] n2 Refractive index of the material across the boundary.
 * @param[in] cos1 SIGNED incident angle that must be calculated as n*p.
 * @param[out] r Refrected beam propagation direction filled in on return.
 *
 * @note Refracted beam propagation direction is computed as p - 2*n*(p*n).
 * @details Refraction for a normal that points inward (n2 => n1) is computed
 *          as:
 *            kn = n1/n2
 *            cos1 = n*p
 *            sin1 = (1 - cos1^2)^0.5
 *            cos2 = (1 - kn^2*sin1^2)^0.5
 *            r = kn*p + (kn*cos1 - cos2)*n
 */

inline mc_point3f_t *refract_cos1(
	mc_point3f_t const *p, mc_point3f_t const *n,
		mc_fp_t n1, mc_fp_t n2, mc_fp_t cos1, mc_point3f_t *r){

	mc_fp_t n1_d_n2 = mc_fdiv(n1, n2);

	mc_fp_t sin2_squared = n1_d_n2*n1_d_n2*(FP_1 - cos1*cos1);

	mc_fp_t k = mc_fsign(cos1)*(n1_d_n2*mc_fabs(cos1) - mc_sqrt(FP_1 - sin2_squared));

	r->x = n1_d_n2*p->x + k*n->x;
	r->y = n1_d_n2*p->y + k*n->y;
	r->z = n1_d_n2*p->z + k*n->z;

	return r;
}

/**
 * @brief Computes propagation direction of the refracted beam for the given
 * incident propagation direction, boundary normal and refractive indices
 * of the materials on each side of the boundary. If the beam is actually
 * reflected, this call can produce unexpected results!
 * 
 * @param[in] p Incident propagation direction.
 * @param[in] n Boundary normal (pointing outwards or inwards).
 * @param[in] n1 Refractive index of the material on the incident side of the boundary.
 * @param[in] n2 Refractive index of the material across the boundary.
 * @param[out] r Refrected beam propagation direction filled in on return.
 *
 * @return Returns the refracted propagation direction as a pointer (input argument r).
 *
 * @note Refracted beam propagation direction is computed as p - 2*n*(p*n).
 * @note No checks are made if the beam is actually reflected. This will lead
 *       to NaN components of the refracted vector or unexpected results!
 * @details Refraction for a normal that points inward (n2 => n1) is computed
 *          as:
 *            kn = n1/n2
 *            cos1 = n*p
 *            sin1 = (1 - cos1^2)^0.5
 *            cos2 = (1 - kn^2*sin1^2)^0.5
 *            r = kn*p + (kn*cos1 - cos2)*n
 */
inline mc_point3f_t *refract(
	mc_point3f_t const *p, mc_point3f_t const *n,
		mc_fp_t n1, mc_fp_t n2, mc_point3f_t *r){

	/* For outwards pointing normal, cos1 is negative. */
	mc_fp_t cos1 = dot3f(p, n);

	mc_fp_t n1_d_n2 = mc_fdiv(n1, n2);

	mc_fp_t sin2_squared = n1_d_n2*n1_d_n2*(FP_1 - cos1*cos1);

	mc_fp_t k = mc_fsign(cos1)*(n1_d_n2*mc_fabs(cos1) - mc_sqrt(FP_1 - sin2_squared));

	r->x = n1_d_n2*p->x - k*n->x;
	r->y = n1_d_n2*p->y - k*n->y;
	r->z = n1_d_n2*p->z - k*n->z;

	return r;
}

/**
 * @brief Computes propagation direction of the refracted beam for the given
 * incident propagation direction, boundary normal and refractive indices
 * of the materials on each side of the boundary. If the beam
 * is actually reflected, the call returns nonzero.
 * 
 * @param[in] p Incident propagation direction.
 * @param[in] n Boundary normal (pointing outwards or inwards).
 * @param[in] n1 Refractive index of the material on the incident side of the boundary.
 * @param[in] n2 Refractive index of the material across the boundary.
 * @param[out] r Refrected beam propagation direction filled in on return.
 *
 * @return Returns zero if the beam is refracted else nonzero
 *         (the parameter r ios left unchanged).
 *
 * @note Refracted beam propagation direction is computed as p - 2*n*(p*n).
 * @note No checks are made if the beam is actually reflected. This will lead
 *       to NaN components of the refracted vector or unexpected results!
 * @details Refraction for a normal that points inward (n2 => n1) is computed
 *          as:
 *            kn = n1/n2
 *            cos1 = n*p
 *            sin1 = (1 - cos1^2)^0.5
 *            cos2 = (1 - kn^2*sin1^2)^0.5
 *            r = kn*p + (kn*cos1 - cos2)*n
 */
inline int refract_safe(mc_point3f_t const *p, mc_point3f_t const *n,
		mc_fp_t n1, mc_fp_t n2, mc_point3f_t *r){

	/* For outwards pointing normal, cos1 is negative. */
	mc_fp_t cos1 = dot3f(p, n);

	mc_fp_t n1_d_n2 = mc_fdiv(n1, n2);

	mc_fp_t sin2_squared = n1_d_n2*n1_d_n2*(FP_1 - cos1*cos1);
	
	if (sin2_squared > FP_1)
		return 1;

	mc_fp_t k = mc_fsign(cos1)*(n1_d_n2*mc_fabs(cos1) - mc_sqrt(FP_1 - sin2_squared));

	r->x = n1_d_n2*p->x - k*n->x;
	r->y = n1_d_n2*p->y - k*n->y;
	r->z = n1_d_n2*p->z - k*n->z;

	return 0;
}

/*#################### End boundary physics implementation ###################*/


/*#################### Begin buffer initialization kernels ###################*/
__kernel void fill_int8(
		__global char *buffer, char const value,
		mc_size_t nb, mc_size_t offset){

	if (get_global_id(0) < nb)
		buffer[offset + get_global_id(0)] = value;
};

__kernel void fill_uint8(
		__global uchar *buffer, uchar const value,
		mc_size_t nb, mc_size_t offset){

	if (get_global_id(0) < nb)
		buffer[offset + get_global_id(0)] = value;
};

__kernel void fill_int16(
		__global short *buffer, short const value,
		mc_size_t nb, mc_size_t offset){

	if (get_global_id(0) < nb)
		buffer[offset + get_global_id(0)] = value;
};

__kernel void fill_uint16(
		__global ushort *buffer, ushort const value,
		mc_size_t nb, mc_size_t offset){

	if (get_global_id(0) < nb)
		buffer[offset + get_global_id(0)] = value;
};

__kernel void fill_int32(
		__global int *buffer, int const value,
		mc_size_t nb, mc_size_t offset){

	if (get_global_id(0) < nb)
		buffer[offset + get_global_id(0)] = value;

};

__kernel void fill_uint32(
		__global uint *buffer, uint const value,
		mc_size_t nb, mc_size_t offset){

	if (get_global_id(0) < nb)
		buffer[offset + get_global_id(0)] = value;
};

__kernel void fill_int64(
		__global long *buffer, long const value,
		mc_size_t nb, mc_size_t offset){

	if (get_global_id(0) < nb)
		buffer[offset + get_global_id(0)] = value;

};

__kernel void fill_uint64(
		__global ulong *buffer, ulong const value,
		mc_size_t nb, mc_size_t offset){

	if (get_global_id(0) < nb)
		buffer[offset + get_global_id(0)] = value;

};

__kernel void fill_float(
		__global float *buffer, float const value,
		mc_size_t nb, mc_size_t offset){

	if (get_global_id(0) < nb)
		buffer[offset + get_global_id(0)] = value;
};

#if MC_USE_DOUBLE_PRECISION || defined(__DOXYGEN__)
__kernel void fill_double(
		__global double *buffer, double const value,
		mc_size_t nb, mc_size_t offset){

	if (get_global_id(0) < nb)
		buffer[offset + get_global_id(0)] = value;

};
#endif

__kernel void ScalarFillMcAccu(
		__global mc_accu_t *buffer, mc_accu_t const value,
		mc_size_t nb, mc_size_t offset){

	if (get_global_id(0) < nb)
		buffer[offset + get_global_id(0)] = value;

};

/*##################### End buffer initialization kernels ####################*/


/*############## Start random number generator implementation ################*/

/**
 * @brief Generates a single precision random number from [0.0, 1.0] and
 *        updates the generator state.
 * @note Due to precision issues
 *       the open nature of the interval is not guaranteed.
 *
 * @param[in,out]	x Mutable state of the random number generator.
 * @param[in]		a Immutable state of the random number generator.
 * @return			A random number from [0.0, 1.0].
 * 
 @details George Marsaglia's Random Number Generator. A single precision 
 *        floating-point number has a 23-bit mantisa, hence only integers from 
 *        [0, 2^23 - 1] can be represented without loss of information.
 */
inline mc_fp_t fp_random_single(unsigned long *x, unsigned a){
	*x = (*x & (uint64_t)0xFFFFFFFFUL)*(a) + (*x >> 32);

	return mc_fdiv(
		((float)( ((unsigned int)(*x) & 0x7FFFFFU) )),
		(float)0x7FFFFFU
	);
};

#if MC_USE_DOUBLE_PRECISION || defined(__DOXYGEN__)
	/**
	 * @brief Generates a double precision random number from [0.0, 1.0] and
	 *        update the generator state.
	 * @note Due to precision issues
	 *       the open nature of the interval is not guaranteed.
	 * @param[in,out]	x Mutable state of the random number generator.
	 * @param[in]		a Immutable state of the random number generator.
	 * @return			A random number from [0.0, 1.0].
	 *
	 * @details George Marsaglia's Random Number Generator. A double precision 
	 *          floating-point number has a 52-bit mantisa, hence only integers from 
	 *          0, 2^52 - 1] can be represented without loss of information.
	 */
	inline mc_fp_t fp_random_double(unsigned long *x, unsigned a){
		*x = (*x & (uint64_t)0xFFFFFFFFUL)*(a) + (*x >> 32);

		/* Generate a random number (0,1). */
		return mc_fdiv(
			((double)( (*x & (uint64_t)0xFFFFFFFFFFFFFUL)) ),
			(double)0xFFFFFFFFFFFFFUL
		);
	};

#endif

__kernel void RngKernel(
	unsigned long x, unsigned a, mc_size_t n, __global mc_fp_t *buffer){
		if (get_global_id(0) == 0){
			for(mc_size_t i=0; i < n; ++i)
				buffer[i] = fp_random(&x, a);
		}
};

/*############### End random number generator implementation #################*/


#endif /* #define __MCBASE_H */
