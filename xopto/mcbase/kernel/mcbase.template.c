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

/*################## Start simulator configuration options ###################*/
/*################### End simulator configuration options ####################*/


/*########## Start basic data types, constants and math functions  ###########*/
/*########### End basic data types, constants and math functions  ############*/


/*############## Start 64-bit atomic increment implementation ################*/

#if MC_USE_64_BIT_ACCUMULATORS || defined(__DOXYGEN__)
/**
 * @brief Atomic deposit of a 32-bit integer weight into 64-bit accumulator.
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
 * @brief Transforms a 2D vector using a 2D matrix with elements
 *        of type ::mc_int_t.
 *
 * @param[in]  m   Pointer to a 2D transformation matrix.
 * @param[in]  v   Pointer to a 2D input vector.
 * @param[out] r   Pointer to a 2D output vector.
 *                 Must NOT be the input vector.
 *
 * @returns        Pointer to the transformed input vector r.
 */
inline mc_intv2_t *mc_transform_intv2(
        mc_matrix2_int_t const *m, mc_intv2_t const *v,  mc_intv2_t *r){
    r->x = m->a_11*v->x + m->a_12*v->y;
    r->y = m->a_21*v->x + m->a_22*v->y;
    return r;
};

/**
 * @brief Transforms a 3D vector using a 3D matrix with elements
 *        of type ::mc_int_t.
 *
 * @param[in]  m   Pointer to a 3D transformation matrix.
 * @param[in]  v   Pointer to a 3D input vector.
 * @param[out] r   Pointer to a 3D output vector.
 *                 Must NOT be the input vector.
 *
 * @returns        Pointer to the transformed input vector r.
 */
inline mc_intv3_t *mc_transform_intv3(
        mc_matrix3_int_t const *m, mc_intv3_t const *v,  mc_intv3_t *r){
    r->x = m->a_11*v->x + m->a_12*v->y + m->a_13*v->z;
    r->y = m->a_21*v->x + m->a_22*v->y + m->a_23*v->z;
    r->z = m->a_31*v->x + m->a_32*v->y + m->a_33*v->z;

    return r;
};

/**
 * @brief Transforms a 4D vector using a 3D matrix with elements
 *        of type ::mc_int_t.
 *
 * @param[in]  m   Pointer to a 4D transformation matrix.
 * @param[in]  v   Pointer to a 4D input vector.
 * @param[out] r   Pointer to a 4D output vector.
 *                 Must NOT be the input vector.
 *
 * @returns        Pointer to the transformed input vector r.
 */
inline mc_intv4_t *mc_transform_intv4(
        mc_matrix4_int_t const *m, mc_intv4_t const *v,  mc_intv4_t *r){
    r->x = m->a_11*v->x + m->a_12*v->y + m->a_13*v->z + m->a_14*v->w;
    r->y = m->a_21*v->x + m->a_22*v->y + m->a_23*v->z + m->a_24*v->w;
    r->z = m->a_31*v->x + m->a_32*v->y + m->a_33*v->z + m->a_34*v->w;
    r->w = m->a_41*v->x + m->a_42*v->y + m->a_43*v->z + m->a_44*v->w;

    return r;
};

/**
 * @brief Multiply two 2D matrices with elements
 *        of type ::mc_int_t.
 *
 * @param[in]  m1  Pointer to the first input 2D matrix.
 * @param[in]  m2  Pointer to the second input 2D matrix.
 * @param[out] r   Pointer to the output 2D matrix.
 *                 Must NOT be any of the two input matrices.
 *
 * @returns        Pointer to the resulting matrix r.
 */
inline mc_matrix2_int_t *mc_matrix2_mul_int(
        mc_matrix2_int_t const *m1, mc_matrix2_int_t const *m2,
        mc_matrix2_int_t *r){
    r->a_11 = m1->a_11*m2->a_11 + m1->a_12*m2->a_21;
    r->a_12 = m1->a_11*m2->a_12 + m1->a_12*m2->a_22;

    r->a_21 = m1->a_21*m2->a_11 + m1->a_22*m2->a_21;
    r->a_22 = m1->a_21*m2->a_12 + m1->a_22*m2->a_22;

    return r;
};

/**
 * @brief Multiply two 3D matrices with elements
 *        of type ::mc_int_t.
 *
 * @param[in]  m1  Pointer to the first input 3D matrix.
 * @param[in]  m2  Pointer to the second input 3D matrix.
 * @param[out] r   Pointer to the output 3D matrix.
 *                 Must NOT be any of the two input matrices.
 *
 * @returns        Pointer to the resulting matrix r.
 */
inline mc_matrix3_int_t *mc_matrix3_mul_int(
        mc_matrix3_int_t const *m1, mc_matrix3_int_t const *m2,
        mc_matrix3_int_t *r){
    r->a_11 = m1->a_11*m2->a_11 + m1->a_12*m2->a_21 + m1->a_13*m2->a_31;
    r->a_12 = m1->a_11*m2->a_12 + m1->a_12*m2->a_22 + m1->a_13*m2->a_32;
    r->a_13 = m1->a_11*m2->a_13 + m1->a_12*m2->a_23 + m1->a_13*m2->a_33;

    r->a_21 = m1->a_21*m2->a_11 + m1->a_22*m2->a_21 + m1->a_23*m2->a_31;
    r->a_22 = m1->a_21*m2->a_12 + m1->a_22*m2->a_22 + m1->a_23*m2->a_32;
    r->a_23 = m1->a_21*m2->a_13 + m1->a_22*m2->a_23 + m1->a_23*m2->a_33;

    r->a_31 = m1->a_31*m2->a_11 + m1->a_32*m2->a_21 + m1->a_33*m2->a_31;
    r->a_32 = m1->a_31*m2->a_12 + m1->a_32*m2->a_22 + m1->a_33*m2->a_32;
    r->a_33 = m1->a_31*m2->a_13 + m1->a_32*m2->a_23 + m1->a_33*m2->a_33;

    return r;
};

/**
 * @brief Multiply two 3D matrices with elements
 *        of type ::mc_int_t.
 *
 * @param[in]  m1  Pointer to the first input 3D matrix.
 * @param[in]  m2  Pointer to the second input 3D matrix.
 * @param[out] r   Pointer to the output 3D matrix.
 *                 Must NOT be any of the two input matrices.
 *
 * @returns        Pointer to the resulting matrix r.
 */
inline mc_matrix4_int_t *mc_matrix4_mul_int(
        mc_matrix4_int_t const *m1, mc_matrix4_int_t const *m2,
        mc_matrix4_int_t *r){
    r->a_11 = m1->a_11*m2->a_11 + m1->a_12*m2->a_21 + m1->a_13*m2->a_31 + m1->a_14*m2->a_41;
    r->a_12 = m1->a_11*m2->a_12 + m1->a_12*m2->a_22 + m1->a_13*m2->a_32 + m1->a_14*m2->a_42;
    r->a_13 = m1->a_11*m2->a_13 + m1->a_12*m2->a_23 + m1->a_13*m2->a_33 + m1->a_14*m2->a_43;
    r->a_14 = m1->a_11*m2->a_14 + m1->a_12*m2->a_24 + m1->a_13*m2->a_34 + m1->a_14*m2->a_44;

    r->a_21 = m1->a_21*m2->a_11 + m1->a_22*m2->a_21 + m1->a_23*m2->a_31 + m1->a_24*m2->a_41;
    r->a_22 = m1->a_21*m2->a_12 + m1->a_22*m2->a_22 + m1->a_23*m2->a_32 + m1->a_24*m2->a_42;
    r->a_23 = m1->a_21*m2->a_13 + m1->a_22*m2->a_23 + m1->a_23*m2->a_33 + m1->a_24*m2->a_43;
    r->a_24 = m1->a_21*m2->a_14 + m1->a_22*m2->a_24 + m1->a_23*m2->a_34 + m1->a_24*m2->a_44;

    r->a_31 = m1->a_31*m2->a_11 + m1->a_32*m2->a_21 + m1->a_33*m2->a_31 + m1->a_34*m2->a_41;
    r->a_32 = m1->a_31*m2->a_12 + m1->a_32*m2->a_22 + m1->a_33*m2->a_32 + m1->a_34*m2->a_42;
    r->a_33 = m1->a_31*m2->a_13 + m1->a_32*m2->a_23 + m1->a_33*m2->a_33 + m1->a_34*m2->a_43;
    r->a_34 = m1->a_31*m2->a_14 + m1->a_32*m2->a_24 + m1->a_33*m2->a_34 + m1->a_34*m2->a_44;

    r->a_41 = m1->a_41*m2->a_11 + m1->a_42*m2->a_21 + m1->a_43*m2->a_31 + m1->a_44*m2->a_41;
    r->a_42 = m1->a_41*m2->a_12 + m1->a_42*m2->a_22 + m1->a_43*m2->a_32 + m1->a_44*m2->a_42;
    r->a_43 = m1->a_41*m2->a_13 + m1->a_42*m2->a_23 + m1->a_43*m2->a_33 + m1->a_44*m2->a_43;
    r->a_44 = m1->a_41*m2->a_14 + m1->a_42*m2->a_24 + m1->a_43*m2->a_34 + m1->a_44*m2->a_44;

    return r;
};

/**
 * @brief Reverses the direction of a 2D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_intv2_t).
 * @param[out] r   Pointer to the output vector (::mc_intv2_t).
 *                 Can be the input vector a.
 *
 * @returns        Pointer to the reversed input vector r.
 */
inline mc_intv2_t *mc_reverse_intv2(mc_intv2_t const *a, mc_intv2_t *r){
    r->x = -a->x;
    r->y = -a->y;
    return r;
};

/**
 * @brief Reverses the direction of a 3D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_intv3_t).
 * @param[out] r   Pointer to the output vector (::mc_intv3_t).
 *                 Can be the input vector a.
 *
 * @returns        Pointer to the reversed input vector r.
 */
inline mc_intv3_t *mc_reverse_intv3(mc_intv3_t const *a, mc_intv3_t *r){
    r->x = -a->x;
    r->y = -a->y;
    r->z = -a->z;
    return r;
};

/**
 * @brief Reverses the direction of a 4D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_intv4_t).
 * @param[out] r   Pointer to the output vector (::mc_intv4_t).
 *                 Can be the input vector a.
 *
 * @returns        Pointer to the reversed input vector r.
 */
inline mc_intv4_t *mc_reverse_intv4(mc_intv4_t const *a, mc_intv4_t *r){
    r->x = -a->x;
    r->y = -a->y;
    r->z = -a->z;
    r->w = -a->w;
    return r;
};

/**
 * @brief Computes the dot product of two 2D vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_intv2_t).
 * @param[in]  b   Pointer to the first input vector (::mc_intv2_t).
 *
 * @returns        Dot product of the two vectors.
 */
inline mc_int_t mc_dot_intv2(mc_intv2_t const *a, mc_intv2_t const *b){
    return a->x*b->x + a->y*b->y;
};

/**
 * @brief Computes the dot product of two 3D vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_intv3_t).
 * @param[in]  b   Pointer to the first input vector (::mc_intv3_t).
 *
 * @returns        Dot product of the two vectors.
 */
inline mc_int_t mc_dot_intv3(mc_intv3_t const *a, mc_intv3_t const *b){
    return a->x*b->x + a->y*b->y + a->z*b->z;
};

/**
 * @brief Computes the dot product of two 4D vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_intv4_t).
 * @param[in]  b   Pointer to the first input vector (::mc_intv4_t).
 *
 * @returns        Dot product of the two vectors.
 */
inline mc_int_t mc_dot_intv4(mc_intv4_t const *a, mc_intv4_t const *b){
    return a->x*b->x + a->y*b->y + a->z*b->z + a->w*b->w;
};

/**
 * @brief Computes the length of a 2D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_intv2_t).
 *
 * @returns        Length of the input vector.
 */
inline mc_fp_t mc_length_intv2(mc_intv2_t const *a){
    return mc_sqrt(a->x*a->x + a->y*a->y);
};

/**
 * @brief Computes the length of a 3D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_intv3_t).
 *
 * @returns        Length of the input vector.
 */
inline mc_fp_t mc_length_intv3(mc_intv3_t const *a){
    return mc_sqrt(a->x*a->x + a->y*a->y + a->z*a->z);
};

/**
 * @brief Computes the length of a 4D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_intv4_t).
 *
 * @returns        Length of the input vector.
 */
inline mc_fp_t mc_length_intv4(mc_intv4_t const *a){
    return mc_sqrt(a->x*a->x + a->y*a->y + a->z*a->z + a->w*a->w);
};

/**
 * @brief Computes the cross product of two vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_intv3_t).
 * @param[in]  b   Pointer to the second input vector (::mc_intv3_t).
 *
 * @param[out] r   Pointer to the output vector (::mc_intv3_t)
 *                 filled with the cross product a x b.
 *                 Must NOT be any of the two input vectors.
 *
 * @returns        Pointer to the cross product vector r.
 */
inline mc_intv3_t *mc_cross_intv3(
        mc_intv3_t const *a, mc_intv3_t const *b, mc_intv3_t *r){
    r->x = a->y*b->z - a->z*b->y;
    r->y = a->z*b->x - a->x*b->z;
    r->z = a->x*b->y - a->y*b->x;

    return r;
};



/**
 * @brief Transforms a 2D vector using a 2D matrix with elements
 *        of type ::mc_size_t.
 *
 * @param[in]  m   Pointer to a 2D transformation matrix.
 * @param[in]  v   Pointer to a 2D input vector.
 * @param[out] r   Pointer to a 2D output vector.
 *                 Must NOT be the input vector.
 *
 * @returns        Pointer to the transformed input vector r.
 */
inline mc_sizev2_t *mc_transform_sizev2(
        mc_matrix2_size_t const *m, mc_sizev2_t const *v,  mc_sizev2_t *r){
    r->x = m->a_11*v->x + m->a_12*v->y;
    r->y = m->a_21*v->x + m->a_22*v->y;
    return r;
};

/**
 * @brief Transforms a 3D vector using a 3D matrix with elements
 *        of type ::mc_size_t.
 *
 * @param[in]  m   Pointer to a 3D transformation matrix.
 * @param[in]  v   Pointer to a 3D input vector.
 * @param[out] r   Pointer to a 3D output vector.
 *                 Must NOT be the input vector.
 *
 * @returns        Pointer to the transformed input vector r.
 */
inline mc_sizev3_t *mc_transform_sizev3(
        mc_matrix3_size_t const *m, mc_sizev3_t const *v,  mc_sizev3_t *r){
    r->x = m->a_11*v->x + m->a_12*v->y + m->a_13*v->z;
    r->y = m->a_21*v->x + m->a_22*v->y + m->a_23*v->z;
    r->z = m->a_31*v->x + m->a_32*v->y + m->a_33*v->z;

    return r;
};

/**
 * @brief Transforms a 4D vector using a 3D matrix with elements
 *        of type ::mc_size_t.
 *
 * @param[in]  m   Pointer to a 4D transformation matrix.
 * @param[in]  v   Pointer to a 4D input vector.
 * @param[out] r   Pointer to a 4D output vector.
 *                 Must NOT be the input vector.
 *
 * @returns        Pointer to the transformed input vector r.
 */
inline mc_sizev4_t *mc_transform_sizev4(
        mc_matrix4_size_t const *m, mc_sizev4_t const *v,  mc_sizev4_t *r){
    r->x = m->a_11*v->x + m->a_12*v->y + m->a_13*v->z + m->a_14*v->w;
    r->y = m->a_21*v->x + m->a_22*v->y + m->a_23*v->z + m->a_24*v->w;
    r->z = m->a_31*v->x + m->a_32*v->y + m->a_33*v->z + m->a_34*v->w;
    r->w = m->a_41*v->x + m->a_42*v->y + m->a_43*v->z + m->a_44*v->w;

    return r;
};

/**
 * @brief Multiply two 2D matrices with elements
 *        of type ::mc_size_t.
 *
 * @param[in]  m1  Pointer to the first input 2D matrix.
 * @param[in]  m2  Pointer to the second input 2D matrix.
 * @param[out] r   Pointer to the output 2D matrix.
 *                 Must NOT be any of the two input matrices.
 *
 * @returns        Pointer to the resulting matrix r.
 */
inline mc_matrix2_size_t *mc_matrix2_mul_size(
        mc_matrix2_size_t const *m1, mc_matrix2_size_t const *m2,
        mc_matrix2_size_t *r){
    r->a_11 = m1->a_11*m2->a_11 + m1->a_12*m2->a_21;
    r->a_12 = m1->a_11*m2->a_12 + m1->a_12*m2->a_22;

    r->a_21 = m1->a_21*m2->a_11 + m1->a_22*m2->a_21;
    r->a_22 = m1->a_21*m2->a_12 + m1->a_22*m2->a_22;

    return r;
};

/**
 * @brief Multiply two 3D matrices with elements
 *        of type ::mc_size_t.
 *
 * @param[in]  m1  Pointer to the first input 3D matrix.
 * @param[in]  m2  Pointer to the second input 3D matrix.
 * @param[out] r   Pointer to the output 3D matrix.
 *                 Must NOT be any of the two input matrices.
 *
 * @returns        Pointer to the resulting matrix r.
 */
inline mc_matrix3_size_t *mc_matrix3_mul_size(
        mc_matrix3_size_t const *m1, mc_matrix3_size_t const *m2,
        mc_matrix3_size_t *r){
    r->a_11 = m1->a_11*m2->a_11 + m1->a_12*m2->a_21 + m1->a_13*m2->a_31;
    r->a_12 = m1->a_11*m2->a_12 + m1->a_12*m2->a_22 + m1->a_13*m2->a_32;
    r->a_13 = m1->a_11*m2->a_13 + m1->a_12*m2->a_23 + m1->a_13*m2->a_33;

    r->a_21 = m1->a_21*m2->a_11 + m1->a_22*m2->a_21 + m1->a_23*m2->a_31;
    r->a_22 = m1->a_21*m2->a_12 + m1->a_22*m2->a_22 + m1->a_23*m2->a_32;
    r->a_23 = m1->a_21*m2->a_13 + m1->a_22*m2->a_23 + m1->a_23*m2->a_33;

    r->a_31 = m1->a_31*m2->a_11 + m1->a_32*m2->a_21 + m1->a_33*m2->a_31;
    r->a_32 = m1->a_31*m2->a_12 + m1->a_32*m2->a_22 + m1->a_33*m2->a_32;
    r->a_33 = m1->a_31*m2->a_13 + m1->a_32*m2->a_23 + m1->a_33*m2->a_33;

    return r;
};

/**
 * @brief Multiply two 3D matrices with elements
 *        of type ::mc_size_t.
 *
 * @param[in]  m1  Pointer to the first input 3D matrix.
 * @param[in]  m2  Pointer to the second input 3D matrix.
 * @param[out] r   Pointer to the output 3D matrix.
 *                 Must NOT be any of the two input matrices.
 *
 * @returns        Pointer to the resulting matrix r.
 */
inline mc_matrix4_size_t *mc_matrix4_mul_size(
        mc_matrix4_size_t const *m1, mc_matrix4_size_t const *m2,
        mc_matrix4_size_t *r){
    r->a_11 = m1->a_11*m2->a_11 + m1->a_12*m2->a_21 + m1->a_13*m2->a_31 + m1->a_14*m2->a_41;
    r->a_12 = m1->a_11*m2->a_12 + m1->a_12*m2->a_22 + m1->a_13*m2->a_32 + m1->a_14*m2->a_42;
    r->a_13 = m1->a_11*m2->a_13 + m1->a_12*m2->a_23 + m1->a_13*m2->a_33 + m1->a_14*m2->a_43;
    r->a_14 = m1->a_11*m2->a_14 + m1->a_12*m2->a_24 + m1->a_13*m2->a_34 + m1->a_14*m2->a_44;

    r->a_21 = m1->a_21*m2->a_11 + m1->a_22*m2->a_21 + m1->a_23*m2->a_31 + m1->a_24*m2->a_41;
    r->a_22 = m1->a_21*m2->a_12 + m1->a_22*m2->a_22 + m1->a_23*m2->a_32 + m1->a_24*m2->a_42;
    r->a_23 = m1->a_21*m2->a_13 + m1->a_22*m2->a_23 + m1->a_23*m2->a_33 + m1->a_24*m2->a_43;
    r->a_24 = m1->a_21*m2->a_14 + m1->a_22*m2->a_24 + m1->a_23*m2->a_34 + m1->a_24*m2->a_44;

    r->a_31 = m1->a_31*m2->a_11 + m1->a_32*m2->a_21 + m1->a_33*m2->a_31 + m1->a_34*m2->a_41;
    r->a_32 = m1->a_31*m2->a_12 + m1->a_32*m2->a_22 + m1->a_33*m2->a_32 + m1->a_34*m2->a_42;
    r->a_33 = m1->a_31*m2->a_13 + m1->a_32*m2->a_23 + m1->a_33*m2->a_33 + m1->a_34*m2->a_43;
    r->a_34 = m1->a_31*m2->a_14 + m1->a_32*m2->a_24 + m1->a_33*m2->a_34 + m1->a_34*m2->a_44;

    r->a_41 = m1->a_41*m2->a_11 + m1->a_42*m2->a_21 + m1->a_43*m2->a_31 + m1->a_44*m2->a_41;
    r->a_42 = m1->a_41*m2->a_12 + m1->a_42*m2->a_22 + m1->a_43*m2->a_32 + m1->a_44*m2->a_42;
    r->a_43 = m1->a_41*m2->a_13 + m1->a_42*m2->a_23 + m1->a_43*m2->a_33 + m1->a_44*m2->a_43;
    r->a_44 = m1->a_41*m2->a_14 + m1->a_42*m2->a_24 + m1->a_43*m2->a_34 + m1->a_44*m2->a_44;

    return r;
};

/**
 * @brief Reverses the direction of a 2D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_sizev2_t).
 * @param[out] r   Pointer to the output vector (::mc_sizev2_t).
 *                 Can be the input vector a.
 *
 * @returns        Pointer to the reversed input vector r.
 */
inline mc_sizev2_t *mc_reverse_sizev2(mc_sizev2_t const *a, mc_sizev2_t *r){
    r->x = -a->x;
    r->y = -a->y;
    return r;
};

/**
 * @brief Reverses the direction of a 3D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_sizev3_t).
 * @param[out] r   Pointer to the output vector (::mc_sizev3_t).
 *                 Can be the input vector a.
 *
 * @returns        Pointer to the reversed input vector r.
 */
inline mc_sizev3_t *mc_reverse_sizev3(mc_sizev3_t const *a, mc_sizev3_t *r){
    r->x = -a->x;
    r->y = -a->y;
    r->z = -a->z;
    return r;
};

/**
 * @brief Reverses the direction of a 4D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_sizev4_t).
 * @param[out] r   Pointer to the output vector (::mc_sizev4_t).
 *                 Can be the input vector a.
 *
 * @returns        Pointer to the reversed input vector r.
 */
inline mc_sizev4_t *mc_reverse_sizev4(mc_sizev4_t const *a, mc_sizev4_t *r){
    r->x = -a->x;
    r->y = -a->y;
    r->z = -a->z;
    r->w = -a->w;
    return r;
};

/**
 * @brief Computes the dot product of two 2D vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_sizev2_t).
 * @param[in]  b   Pointer to the first input vector (::mc_sizev2_t).
 *
 * @returns        Dot product of the two vectors.
 */
inline mc_size_t mc_dot_sizev2(mc_sizev2_t const *a, mc_sizev2_t const *b){
    return a->x*b->x + a->y*b->y;
};

/**
 * @brief Computes the dot product of two 3D vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_sizev3_t).
 * @param[in]  b   Pointer to the first input vector (::mc_sizev3_t).
 *
 * @returns        Dot product of the two vectors.
 */
inline mc_size_t mc_dot_sizev3(mc_sizev3_t const *a, mc_sizev3_t const *b){
    return a->x*b->x + a->y*b->y + a->z*b->z;
};

/**
 * @brief Computes the dot product of two 4D vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_sizev4_t).
 * @param[in]  b   Pointer to the first input vector (::mc_sizev4_t).
 *
 * @returns        Dot product of the two vectors.
 */
inline mc_size_t mc_dot_sizev4(mc_sizev4_t const *a, mc_sizev4_t const *b){
    return a->x*b->x + a->y*b->y + a->z*b->z + a->w*b->w;
};

/**
 * @brief Computes the length of a 2D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_sizev2_t).
 *
 * @returns        Length of the input vector.
 */
inline mc_fp_t mc_length_sizev2(mc_sizev2_t const *a){
    return mc_sqrt(a->x*a->x + a->y*a->y);
};

/**
 * @brief Computes the length of a 3D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_sizev3_t).
 *
 * @returns        Length of the input vector.
 */
inline mc_fp_t mc_length_sizev3(mc_sizev3_t const *a){
    return mc_sqrt(a->x*a->x + a->y*a->y + a->z*a->z);
};

/**
 * @brief Computes the length of a 4D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_sizev4_t).
 *
 * @returns        Length of the input vector.
 */
inline mc_fp_t mc_length_sizev4(mc_sizev4_t const *a){
    return mc_sqrt(a->x*a->x + a->y*a->y + a->z*a->z + a->w*a->w);
};

/**
 * @brief Computes the cross product of two vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_sizev3_t).
 * @param[in]  b   Pointer to the second input vector (::mc_sizev3_t).
 *
 * @param[out] r   Pointer to the output vector (::mc_sizev3_t)
 *                 filled with the cross product a x b.
 *                 Must NOT be any of the two input vectors.
 *
 * @returns        Pointer to the cross product vector r.
 */
inline mc_sizev3_t *mc_cross_sizev3(
        mc_sizev3_t const *a, mc_sizev3_t const *b, mc_sizev3_t *r){
    r->x = a->y*b->z - a->z*b->y;
    r->y = a->z*b->x - a->x*b->z;
    r->z = a->x*b->y - a->y*b->x;

    return r;
};



/**
 * @brief Transforms a 2D vector using a 2D matrix with elements
 *        of type ::mc_fp_t.
 *
 * @param[in]  m   Pointer to a 2D transformation matrix.
 * @param[in]  v   Pointer to a 2D input vector.
 * @param[out] r   Pointer to a 2D output vector.
 *                 Must NOT be the input vector.
 *
 * @returns        Pointer to the transformed input vector r.
 */
inline mc_fpv2_t *mc_transform_fpv2(
        mc_matrix2_fp_t const *m, mc_fpv2_t const *v,  mc_fpv2_t *r){
    r->x = m->a_11*v->x + m->a_12*v->y;
    r->y = m->a_21*v->x + m->a_22*v->y;
    return r;
};

/**
 * @brief Transforms a 3D vector using a 3D matrix with elements
 *        of type ::mc_fp_t.
 *
 * @param[in]  m   Pointer to a 3D transformation matrix.
 * @param[in]  v   Pointer to a 3D input vector.
 * @param[out] r   Pointer to a 3D output vector.
 *                 Must NOT be the input vector.
 *
 * @returns        Pointer to the transformed input vector r.
 */
inline mc_fpv3_t *mc_transform_fpv3(
        mc_matrix3_fp_t const *m, mc_fpv3_t const *v,  mc_fpv3_t *r){
    r->x = m->a_11*v->x + m->a_12*v->y + m->a_13*v->z;
    r->y = m->a_21*v->x + m->a_22*v->y + m->a_23*v->z;
    r->z = m->a_31*v->x + m->a_32*v->y + m->a_33*v->z;

    return r;
};

/**
 * @brief Transforms a 4D vector using a 3D matrix with elements
 *        of type ::mc_fp_t.
 *
 * @param[in]  m   Pointer to a 4D transformation matrix.
 * @param[in]  v   Pointer to a 4D input vector.
 * @param[out] r   Pointer to a 4D output vector.
 *                 Must NOT be the input vector.
 *
 * @returns        Pointer to the transformed input vector r.
 */
inline mc_fpv4_t *mc_transform_fpv4(
        mc_matrix4_fp_t const *m, mc_fpv4_t const *v,  mc_fpv4_t *r){
    r->x = m->a_11*v->x + m->a_12*v->y + m->a_13*v->z + m->a_14*v->w;
    r->y = m->a_21*v->x + m->a_22*v->y + m->a_23*v->z + m->a_24*v->w;
    r->z = m->a_31*v->x + m->a_32*v->y + m->a_33*v->z + m->a_34*v->w;
    r->w = m->a_41*v->x + m->a_42*v->y + m->a_43*v->z + m->a_44*v->w;

    return r;
};

/**
 * @brief Multiply two 2D matrices with elements
 *        of type ::mc_fp_t.
 *
 * @param[in]  m1  Pointer to the first input 2D matrix.
 * @param[in]  m2  Pointer to the second input 2D matrix.
 * @param[out] r   Pointer to the output 2D matrix.
 *                 Must NOT be any of the two input matrices.
 *
 * @returns        Pointer to the resulting matrix r.
 */
inline mc_matrix2_fp_t *mc_matrix2_mul_fp(
        mc_matrix2_fp_t const *m1, mc_matrix2_fp_t const *m2,
        mc_matrix2_fp_t *r){
    r->a_11 = m1->a_11*m2->a_11 + m1->a_12*m2->a_21;
    r->a_12 = m1->a_11*m2->a_12 + m1->a_12*m2->a_22;

    r->a_21 = m1->a_21*m2->a_11 + m1->a_22*m2->a_21;
    r->a_22 = m1->a_21*m2->a_12 + m1->a_22*m2->a_22;

    return r;
};

/**
 * @brief Multiply two 3D matrices with elements
 *        of type ::mc_fp_t.
 *
 * @param[in]  m1  Pointer to the first input 3D matrix.
 * @param[in]  m2  Pointer to the second input 3D matrix.
 * @param[out] r   Pointer to the output 3D matrix.
 *                 Must NOT be any of the two input matrices.
 *
 * @returns        Pointer to the resulting matrix r.
 */
inline mc_matrix3_fp_t *mc_matrix3_mul_fp(
        mc_matrix3_fp_t const *m1, mc_matrix3_fp_t const *m2,
        mc_matrix3_fp_t *r){
    r->a_11 = m1->a_11*m2->a_11 + m1->a_12*m2->a_21 + m1->a_13*m2->a_31;
    r->a_12 = m1->a_11*m2->a_12 + m1->a_12*m2->a_22 + m1->a_13*m2->a_32;
    r->a_13 = m1->a_11*m2->a_13 + m1->a_12*m2->a_23 + m1->a_13*m2->a_33;

    r->a_21 = m1->a_21*m2->a_11 + m1->a_22*m2->a_21 + m1->a_23*m2->a_31;
    r->a_22 = m1->a_21*m2->a_12 + m1->a_22*m2->a_22 + m1->a_23*m2->a_32;
    r->a_23 = m1->a_21*m2->a_13 + m1->a_22*m2->a_23 + m1->a_23*m2->a_33;

    r->a_31 = m1->a_31*m2->a_11 + m1->a_32*m2->a_21 + m1->a_33*m2->a_31;
    r->a_32 = m1->a_31*m2->a_12 + m1->a_32*m2->a_22 + m1->a_33*m2->a_32;
    r->a_33 = m1->a_31*m2->a_13 + m1->a_32*m2->a_23 + m1->a_33*m2->a_33;

    return r;
};

/**
 * @brief Multiply two 3D matrices with elements
 *        of type ::mc_fp_t.
 *
 * @param[in]  m1  Pointer to the first input 3D matrix.
 * @param[in]  m2  Pointer to the second input 3D matrix.
 * @param[out] r   Pointer to the output 3D matrix.
 *                 Must NOT be any of the two input matrices.
 *
 * @returns        Pointer to the resulting matrix r.
 */
inline mc_matrix4_fp_t *mc_matrix4_mul_fp(
        mc_matrix4_fp_t const *m1, mc_matrix4_fp_t const *m2,
        mc_matrix4_fp_t *r){
    r->a_11 = m1->a_11*m2->a_11 + m1->a_12*m2->a_21 + m1->a_13*m2->a_31 + m1->a_14*m2->a_41;
    r->a_12 = m1->a_11*m2->a_12 + m1->a_12*m2->a_22 + m1->a_13*m2->a_32 + m1->a_14*m2->a_42;
    r->a_13 = m1->a_11*m2->a_13 + m1->a_12*m2->a_23 + m1->a_13*m2->a_33 + m1->a_14*m2->a_43;
    r->a_14 = m1->a_11*m2->a_14 + m1->a_12*m2->a_24 + m1->a_13*m2->a_34 + m1->a_14*m2->a_44;

    r->a_21 = m1->a_21*m2->a_11 + m1->a_22*m2->a_21 + m1->a_23*m2->a_31 + m1->a_24*m2->a_41;
    r->a_22 = m1->a_21*m2->a_12 + m1->a_22*m2->a_22 + m1->a_23*m2->a_32 + m1->a_24*m2->a_42;
    r->a_23 = m1->a_21*m2->a_13 + m1->a_22*m2->a_23 + m1->a_23*m2->a_33 + m1->a_24*m2->a_43;
    r->a_24 = m1->a_21*m2->a_14 + m1->a_22*m2->a_24 + m1->a_23*m2->a_34 + m1->a_24*m2->a_44;

    r->a_31 = m1->a_31*m2->a_11 + m1->a_32*m2->a_21 + m1->a_33*m2->a_31 + m1->a_34*m2->a_41;
    r->a_32 = m1->a_31*m2->a_12 + m1->a_32*m2->a_22 + m1->a_33*m2->a_32 + m1->a_34*m2->a_42;
    r->a_33 = m1->a_31*m2->a_13 + m1->a_32*m2->a_23 + m1->a_33*m2->a_33 + m1->a_34*m2->a_43;
    r->a_34 = m1->a_31*m2->a_14 + m1->a_32*m2->a_24 + m1->a_33*m2->a_34 + m1->a_34*m2->a_44;

    r->a_41 = m1->a_41*m2->a_11 + m1->a_42*m2->a_21 + m1->a_43*m2->a_31 + m1->a_44*m2->a_41;
    r->a_42 = m1->a_41*m2->a_12 + m1->a_42*m2->a_22 + m1->a_43*m2->a_32 + m1->a_44*m2->a_42;
    r->a_43 = m1->a_41*m2->a_13 + m1->a_42*m2->a_23 + m1->a_43*m2->a_33 + m1->a_44*m2->a_43;
    r->a_44 = m1->a_41*m2->a_14 + m1->a_42*m2->a_24 + m1->a_43*m2->a_34 + m1->a_44*m2->a_44;

    return r;
};

/**
 * @brief Reverses the direction of a 2D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_fpv2_t).
 * @param[out] r   Pointer to the output vector (::mc_fpv2_t).
 *                 Can be the input vector a.
 *
 * @returns        Pointer to the reversed input vector r.
 */
inline mc_fpv2_t *mc_reverse_fpv2(mc_fpv2_t const *a, mc_fpv2_t *r){
    r->x = -a->x;
    r->y = -a->y;
    return r;
};

/**
 * @brief Reverses the direction of a 3D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_fpv3_t).
 * @param[out] r   Pointer to the output vector (::mc_fpv3_t).
 *                 Can be the input vector a.
 *
 * @returns        Pointer to the reversed input vector r.
 */
inline mc_fpv3_t *mc_reverse_fpv3(mc_fpv3_t const *a, mc_fpv3_t *r){
    r->x = -a->x;
    r->y = -a->y;
    r->z = -a->z;
    return r;
};

/**
 * @brief Reverses the direction of a 4D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_fpv4_t).
 * @param[out] r   Pointer to the output vector (::mc_fpv4_t).
 *                 Can be the input vector a.
 *
 * @returns        Pointer to the reversed input vector r.
 */
inline mc_fpv4_t *mc_reverse_fpv4(mc_fpv4_t const *a, mc_fpv4_t *r){
    r->x = -a->x;
    r->y = -a->y;
    r->z = -a->z;
    r->w = -a->w;
    return r;
};

/**
 * @brief Computes the dot product of two 2D vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv2_t).
 * @param[in]  b   Pointer to the first input vector (::mc_fpv2_t).
 *
 * @returns        Dot product of the two vectors.
 */
inline mc_fp_t mc_dot_fpv2(mc_fpv2_t const *a, mc_fpv2_t const *b){
    return a->x*b->x + a->y*b->y;
};

/**
 * @brief Computes the dot product of two 3D vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv3_t).
 * @param[in]  b   Pointer to the first input vector (::mc_fpv3_t).
 *
 * @returns        Dot product of the two vectors.
 */
inline mc_fp_t mc_dot_fpv3(mc_fpv3_t const *a, mc_fpv3_t const *b){
    return a->x*b->x + a->y*b->y + a->z*b->z;
};

/**
 * @brief Computes the dot product of two 4D vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv4_t).
 * @param[in]  b   Pointer to the first input vector (::mc_fpv4_t).
 *
 * @returns        Dot product of the two vectors.
 */
inline mc_fp_t mc_dot_fpv4(mc_fpv4_t const *a, mc_fpv4_t const *b){
    return a->x*b->x + a->y*b->y + a->z*b->z + a->w*b->w;
};

/**
 * @brief Computes the length of a 2D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_fpv2_t).
 *
 * @returns        Length of the input vector.
 */
inline mc_fp_t mc_length_fpv2(mc_fpv2_t const *a){
    return mc_sqrt(a->x*a->x + a->y*a->y);
};

/**
 * @brief Computes the length of a 3D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_fpv3_t).
 *
 * @returns        Length of the input vector.
 */
inline mc_fp_t mc_length_fpv3(mc_fpv3_t const *a){
    return mc_sqrt(a->x*a->x + a->y*a->y + a->z*a->z);
};

/**
 * @brief Computes the length of a 4D vector.
 *
 * @param[in]  a   Pointer to the input vector (::mc_fpv4_t).
 *
 * @returns        Length of the input vector.
 */
inline mc_fp_t mc_length_fpv4(mc_fpv4_t const *a){
    return mc_sqrt(a->x*a->x + a->y*a->y + a->z*a->z + a->w*a->w);
};

/**
 * @brief Computes the cross product of two vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv3_t).
 * @param[in]  b   Pointer to the second input vector (::mc_fpv3_t).
 *
 * @param[out] r   Pointer to the output vector (::mc_fpv3_t)
 *                 filled with the cross product a x b.
 *                 Must NOT be any of the two input vectors.
 *
 * @returns        Pointer to the cross product vector r.
 */
inline mc_fpv3_t *mc_cross_fpv3(
        mc_fpv3_t const *a, mc_fpv3_t const *b, mc_fpv3_t *r){
    r->x = a->y*b->z - a->z*b->y;
    r->y = a->z*b->x - a->x*b->z;
    r->z = a->x*b->y - a->y*b->x;

    return r;
};

/**
 * @brief Normalizes the length of the input 2D vector to unity.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv2_t).
 * @param[out] r   Pointer to output vector (::mc_fpv2_t).
 *                 Can be the input vector a.
 *
 * @returns        Pointer to the unit length vector r.
 */
inline mc_fpv2_t *mc_normalize_fpv2(mc_fpv2_t const *a, mc_fpv2_t *r){
    mc_fp_t k = mc_fdiv(FP_1, mc_length_fpv2(a));
    r->x = a->x*k;
    r->y = a->y*k;
    return r;
};

/**
 * @brief Normalizes the length of the input 3D vector to unity.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv3_t).
 * @param[out] r   Pointer to output vector (::mc_fpv3_t).
 *                 Can be the input vector a.
 *
 * @returns        Pointer to the unit length vector r.
 */
inline mc_fpv3_t *mc_normalize_fpv3(mc_fpv3_t const *a, mc_fpv3_t *r){
    mc_fp_t k = mc_fdiv(FP_1, mc_length_fpv3(a));
    r->x = a->x*k;
    r->y = a->y*k;
    r->z = a->z*k;
    return r;
};

/**
 * @brief Normalizes the length of the input 4D vector to unity.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv4_t).
 * @param[out] r   Pointer to output vector (::mc_fpv4_t).
 *                 Can be the input vector a.
 *
 * @returns        Pointer to the unit length vector r.
 */
inline mc_fpv4_t *mc_normalize_fpv4(mc_fpv4_t const *a, mc_fpv4_t *r){
    mc_fp_t k = mc_fdiv(FP_1, mc_length_fpv4(a));
    r->x = a->x*k;
    r->y = a->y*k;
    r->z = a->z*k;
    r->w = a->w*k;
    return r;
};

/**
 * @brief Computes the squared distance between two 2D vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv2_t).
 * @param[in]  b   Pointer to the second input vector (::mc_fpv2_t).
 *
 * @returns        Squared distance between the two vectors.
 */
inline mc_fp_t mc_distance2_fpv2(mc_fpv2_t const *a, mc_fpv2_t const *b){
    mc_fpv2_t diff = {a->x - b->x, a->y - b->y};
    return mc_dot_fpv2(&diff, &diff);
};

/**
 * @brief Computes the squared distance between two 3D vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv3_t).
 * @param[in]  b   Pointer to the second input vector (::mc_fpv3_t).
 *
 * @returns        Squared distance between the two vectors.
 */
inline mc_fp_t mc_distance2_fpv3(mc_fpv3_t const *a, mc_fpv3_t const *b){
    mc_fpv3_t diff = {a->x - b->x, a->y - b->y, a->z - b->z};
    return mc_dot_fpv3(&diff, &diff);
};

/**
 * @brief Computes the squared distance between two 4D vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv4_t).
 * @param[in]  b   Pointer to the second input vector (::mc_fpv4_t).
 *
 * @returns        Squared distance between the two vectors.
 */
inline mc_fp_t mc_distance2_fpv4(mc_fpv4_t const *a, mc_fpv4_t const *b){
    mc_fpv4_t diff = {a->x - b->x, a->y - b->y, a->z - b->z, a->w - b->w};
    return mc_dot_fpv4(&diff, &diff);
};

/**
 * @brief Computes the distance between two 2D vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv2_t).
 * @param[in]  b   Pointer to the second input vector (::mc_fpv2_t).
 *
 * @returns        Distance between the two vectors.
 */
inline mc_fp_t mc_distance_fpv2(mc_fpv2_t const *a, mc_fpv2_t const *b){
    mc_fpv2_t diff = {a->x - b->x, a->y - b->y};
    return mc_length_fpv2(&diff);
};

/**
 * @brief Computes the distance between two 3D vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv3_t).
 * @param[in]  b   Pointer to the second input vector (::mc_fpv3_t).
 *
 * @returns        Distance between the two vectors.
 */
inline mc_fp_t mc_distance_fpv3(mc_fpv3_t const *a, mc_fpv3_t const *b){
    mc_fpv3_t diff = {a->x - b->x, a->y - b->y, a->z - b->z};
    return mc_length_fpv3(&diff);
};

/**
 * @brief Computes the distance between two 4D vectors.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv4_t).
 * @param[in]  b   Pointer to the second input vector (::mc_fpv4_t).
 *
 * @returns        Distance between the two vectors.
 */
inline mc_fp_t mc_distance_fpv4(mc_fpv4_t const *a, mc_fpv4_t const *b){
    mc_fpv4_t diff = {a->x - b->x, a->y - b->y, a->z - b->z, a->w - b->w};
    return mc_length_fpv4(&diff);
};

/**
 * @brief Computes a + b*c - multiply-add for 2D vectors of type ::mc_fpv2_t.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv2_t).
 * @param[in]  b   Pointer to the second input vector (::mc_fpv2_t).
 * @param[in]  c   Multiplication factor
 * @param[out] r   Pointer to the output/result vector.
 *                 Can be any of the two input vectors.
 *
 * @returns        Pointer to the result vector r.
 */
inline mc_fpv2_t *mc_mad_fpv2(
        mc_fpv2_t const *a, mc_fpv2_t const *b, mc_fp_t c, mc_fpv2_t *r){
    r->x = a->x + b->x*c;
    r->y = a->y + b->y*c;
    return r;
};

/**
 * @brief Computes a + b*c - multiply-add for 3D vectors of type ::mc_fpv3_t.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv2_t).
 * @param[in]  b   Pointer to the second input vector (::mc_fpv2_t).
 * @param[in]  c   Multiplication factor
 * @param[out] r   Pointer to the output/result vector.
 *                 Can be any of the two input vectors.
 *
 * @returns        Pointer to the result vector r.
 */
inline mc_fpv3_t *mc_mad_fpv3(
        mc_fpv3_t const *a, mc_fpv3_t const *b, mc_fp_t c, mc_fpv3_t *r){
    r->x = a->x + b->x*c;
    r->y = a->y + b->y*c;
    r->z = a->z + b->z*c;
    return r;
};

/**
 * @brief Computes a + b*c - multiply-add for 4D vectors of type ::mc_fpv4_t.
 *
 * @param[in]  a   Pointer to the first input vector (::mc_fpv4_t).
 * @param[in]  b   Pointer to the second input vector (::mc_fpv4_t).
 * @param[in]  c   Multiplication factor
 * @param[out] r   Pointer to the output/result vector.
 *                 Can be any of the two input vectors.
 *
 * @returns        Pointer to the result vector r.
 */
inline mc_fpv4_t *mc_mad_fpv4(
        mc_fpv4_t const *a, mc_fpv4_t const *b, mc_fp_t c, mc_fpv4_t *r){
    r->x = a->x + b->x*c;
    r->y = a->y + b->y*c;
    r->z = a->z + b->z*c;
    r->w = a->w + b->w*c;
    return r;
};

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
inline int mc_rectf_contains_ex(
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
inline int mc_circf_contains_ex(
		mc_fp_t center_x, mc_fp_t center_y, mc_fp_t r, mc_fp_t x, mc_fp_t y){
	mc_fp_t dx = ((center_x) - (x));
	mc_fp_t dy = ((center_y) - (y));

	return dx*dx + dy*dy <= r*r;
};

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
inline int mc_slotf_contains_ex(mc_fp_t cx, mc_fp_t cy,
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
};
/*#################### End vector/matrix implementation ######################*/


/*################### Start debug support implementation #####################*/
/*#################### End debug support implementation ######################*/


/*################# Start linear lookup table implementation #################*/
/*################## End linear lookup table implementation ##################*/


/*################## Start boundary physics implementation ###################*/

/**
 * @brief Compute cosine of the critical incidence angle beyond which
 *        the incident beam is reflected at the boundary n1 => n2 of materials
 *        with refractive indices n1 and n2.
 *
 * @param[in] n1 Refractive index of the incident medium.
 * @param[in] n2 Refractive index of the material across the boundary.
 *
 * @return       Cosine of the incident angle beyond which the incident beam
 *               is reflected at the given boundary.
 */
inline mc_fp_t cos_critical(mc_fp_t n1, mc_fp_t n2){
    return (n1 > n2) ? mc_sqrt(FP_1 - mc_fdiv(n2*n2, n1*n1)): FP_0;
};

/**
 * @brief Computes reflectance for the given boundary conditions.
 *
 * @param[in] n1           Refractive index of the incident layer.
 * @param[in] n2           Refractive index of the layer across the boundary.
 * @param[in] cos1         Incidence angle (with respect to the z axis) cosine.
 * @param[in] cos_critical Critical angle cosine for the interface n1 => n2.
 *
 * @return                 Returns the calculated reflectance probability
 *                         from [0.0, 1.0].
 */
inline mc_fp_t reflectance(
		mc_fp_t n1, mc_fp_t n2, mc_fp_t cos1, mc_fp_t cos_critical){
	mc_fp_t Rp, Rs, R = FP_1;
	mc_fp_t n1_d_n2;
	mc_fp_t sin1, sin2, cos2;
	mc_fp_t n_cos1, n_cos2;

	cos1 = mc_fabs(cos1);

	if (n1 == n2)
		return FP_0;

	if(cos1 > cos_critical){
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
 *        using data from the secondary side.
 *
 * @param[in] n1            Refractive index of the incident layer.
 * @param[in] n2            Refractive index of the layer across the boundary.
 * @param[in] cos2          Incidence angle cosine (with respect to the
 *                          boundary normal).
 *
 * @return                  Returns the calculated reflectance probability
 *                          from [0.0, 1.0].
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
 * @brief Computes propagation direction of the reflected beam for the given
 *        incident propagation direction and boundary normal.
 *
 * @param[in] p  Incident propagation direction.
 * @param[in] n  Boundary normal (can be pointing outwards or inwards).
 * @param[out] r Reflected beam propagation direction on return.
 *
 * @return Pointer to the propagation direction r of the reflected beam.
 *
 * @note Reflected beam propagation direction is computed as p - 2*n*dot(p, n)
 */
inline mc_point3f_t *reflect(
        mc_point3f_t const *p, mc_point3f_t const *n, mc_point3f_t *r){
	mc_fp_t p_n_2 = FP_2*mc_dot_point3f(p, n);
	r->x = p->x - n->x*p_n_2;
	r->y = p->y - n->y*p_n_2;
	r->z = p->z - n->z*p_n_2;

	return r;
};

/**
 * @brief Computes propagation direction of the refracted beam for the given
 *        incident propagation direction, boundary normal and refractive
 *        indices of the materials on both sides of the boundary.
 *        Requires signed incident angle cosine computed as dot(n, p).
 *
 * @param[in] p    Incident propagation direction.
 * @param[in] n    Boundary normal (pointing outwards or inwards, requires
 *                 signed cos1 to resolve the surface normal direction!).
 * @param[in] n1   Refractive index of the material on the incident side
 *                 of the boundary.
 * @param[in] n2   Refractive index of the material across the boundary.
 * @param[in] cos1 SIGNED incident angle that must be calculated as dot(n,p).
 * @param[out] r   Refrected beam propagation direction filled in on return.
 *
 * @return         Pointer to the propagation direction r of the refracted beam.
 *
 * @note Refracted beam propagation direction is computed as p - 2*n*(p*n).
 *
 * @details Refraction for a normal that points inward (n2 => n1) is computed
 *          as:
 *            kn = n1/n2
 *            cos1 = dot(n, p)
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
};

/**
 * @brief Computes propagation direction of the refracted beam for the given
 *        incident propagation direction, boundary normal and refractive
 *        indices of the materials on both sides of the boundary.
 *
 * @note  If the beam is actually reflected, this call can produce unexpected
 *        results!
 *
 * @param[in] p  Incident propagation direction.
 * @param[in] n  Boundary normal (pointing outwards or inwards).
 * @param[in] n1 Refractive index of the material on the incident side
 *               of the boundary.
 * @param[in] n2 Refractive index of the material across the boundary.
 * @param[out] r Refrected beam propagation direction filled in on return.
 *
 * @return       Pointer to the propagation direction r of the refracted beam.
 *
 * @note    No checks are made if the beam is actually reflected. This can lead
 *          to NaN components int the direction vector of the refracted beam
 *          or other unexpected results!
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
	mc_fp_t cos1 = mc_dot_point3f(p, n);

	mc_fp_t n1_d_n2 = mc_fdiv(n1, n2);

	mc_fp_t sin2_squared = n1_d_n2*n1_d_n2*(FP_1 - cos1*cos1);

	mc_fp_t k = mc_fsign(cos1)*(n1_d_n2*mc_fabs(cos1) - mc_sqrt(FP_1 - sin2_squared));

	r->x = n1_d_n2*p->x - k*n->x;
	r->y = n1_d_n2*p->y - k*n->y;
	r->z = n1_d_n2*p->z - k*n->z;

	return r;
};

/**
 * @brief Computes propagation direction of the refracted beam for the given
 *        incident propagation direction, boundary normal and refractive
 *        indices of the materials on each side of the boundary.
 *        If the beam is actually reflected, the call returns nonzero.
 *
 * @param[in] p  Incident propagation direction.
 * @param[in] n  Boundary normal (pointing outwards or inwards).
 * @param[in] n1 Refractive index of the material on the incident side
 *               of the boundary.
 * @param[in] n2 Refractive index of the material across the boundary.
 * @param[out] r Refrected beam propagation direction filled in on return.
 *
 * @return       Returns zero if the beam is refracted else nonzero
 *               (parameter r is left unchanged if returning nonzero).
 *
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
	mc_fp_t cos1 = mc_dot_point3f(p, n);

	mc_fp_t n1_d_n2 = mc_fdiv(n1, n2);

	mc_fp_t sin2_squared = n1_d_n2*n1_d_n2*(FP_1 - cos1*cos1);

	if (sin2_squared > FP_1)
		return 1;

	mc_fp_t k = mc_fsign(cos1)*(n1_d_n2*mc_fabs(cos1) - mc_sqrt(FP_1 - sin2_squared));

	r->x = n1_d_n2*p->x - k*n->x;
	r->y = n1_d_n2*p->y - k*n->y;
	r->z = n1_d_n2*p->z - k*n->z;

	return 0;
};

/*#################### End boundary physics implementation ###################*/


/**
 * @brief Called when the photon packet needs to be scattered.
 *
 * @param[in, out] dir       Current propagation direction updated on return.
 * @param[in]      cos_theta Scattering angle cosine.
 * @param[in]      fi        Azimuth angle (radians).
 *
 * @details                  Function computes and updates the packet
 *                           propagation direction.
 */
inline void scatter_direction(mc_point3f_t *dir, mc_fp_t cos_theta, mc_fp_t fi){
	mc_fp_t sin_fi, cos_fi;
	mc_fp_t sin_theta;
	mc_fp_t px, k;
	mc_fp_t sin_theta_cos_fi, sin_theta_sin_fi;

	sin_theta = mc_sqrt(FP_1 - cos_theta*cos_theta);

	mc_sincos(fi, &sin_fi, &cos_fi);

	sin_theta_cos_fi = sin_theta*cos_fi;
	sin_theta_sin_fi = sin_theta*sin_fi;

	px = dir->x;

	if(mc_fabs(dir->z) >= FP_COS_0){
		dir->x = sin_theta_cos_fi;
		dir->y = sin_theta_sin_fi;
		dir->z = mc_fcopysign(cos_theta, dir->z*cos_theta);
	}else{
		k = mc_sqrt(FP_1 - dir->z*dir->z);

		dir->x = mc_fdiv(sin_theta_cos_fi*px*dir->z - sin_theta_sin_fi*dir->y, k) +
			px*cos_theta;
		dir->y = mc_fdiv(sin_theta_cos_fi*dir->y*dir->z + sin_theta_sin_fi*px, k) +
			dir->y*cos_theta;
		dir->z = (-sin_theta_cos_fi)*k + dir->z*cos_theta;
	};

	/* Single precision can lose unity vector length. */
	#if !defined(MC_DOUBLE_PRECISION)
		mc_normalize_point3f(dir);
	#endif
};


/*################## Start accumulator cache implementation ##################*/
/*################### End accumulator cache implementation ###################*/


/*############## Start random number generator implementation ################*/

/**
 * @brief Generates a single precision random number from [0.0, 1.0] and
 *        updates the generator state.
 *
 * @param[in,out]	x Mutable state of the random number generator.
 * @param[in]		a Immutable state of the random number generator.
 * @return			A random number from [0.0, 1.0].
 *
 @details George Marsaglia's Random Number Generator. A single precision
 *        floating-point number has a 23-bit mantissa, hence only integers
 *        from [0, 2^23 - 1] can be represented without loss of information.
 */
inline mc_fp_t fp_random_single(unsigned long *x, unsigned a){
	#if MC_USE_ENHANCED_RNG
		*x = (*x & (uint64_t)0xFFFFFFFFUL)*(a) + (*x >> 32);
		uint32_t high = *x;
		*x = (*x & (uint64_t)0xFFFFFFFFUL)*(a) + (*x >> 32);
		uint32_t low = *x;

		return mc_fdiv(
			(float) ((((uint64_t)high) << 32U) +low),
			(float) 0xFFFFFFFFFFFFFFFFUL
		);
	#else
		*x = (*x & (uint64_t)0xFFFFFFFFUL)*(a) + (*x >> 32);
		return mc_fdiv(
			(float)(uint32_t)*x,
			(float)0xFFFFFFFFU
		);
	#endif
};

#if MC_USE_DOUBLE_PRECISION || defined(__DOXYGEN__)
	/**
	 * @brief Generates a double precision random number from [0.0, 1.0] and
	 *        update the generator state.
	 *
	 * @param[in,out]	x Mutable state of the random number generator.
	 * @param[in]		a Immutable state of the random number generator.
	 * @return			A random number from [0.0, 1.0].
	 *
	 * @details George Marsaglia's Random Number Generator. A double precision
	 *          floating-point number has a 52-bit mantissa, hence only
	 *          integers from  [0, 2^52 - 1] can be represented without
	 *          loss of information.
	 */
	inline mc_fp_t fp_random_double(unsigned long *x, unsigned a){
		#if MC_USE_ENHANCED_RNG
			*x = (*x & (uint64_t)0xFFFFFFFFUL)*(a) + (*x >> 32);
			uint32_t high = *x;
			*x = (*x & (uint64_t)0xFFFFFFFFUL)*(a) + (*x >> 32);
			uint32_t low = *x;

			return mc_fdiv(
				(double) ((((uint64_t)high) << 32U) + low),
				(double) 0xFFFFFFFFFFFFFFFFUL
			);
		#else
			*x = (*x & (uint64_t)0xFFFFFFFFUL)*(a) + (*x >> 32);
			return mc_fdiv(
				(double)(uint32_t)*x,
				(double)0xFFFFFFFFU
			);
		#endif
	};
#endif

/**
 * @brief Default floating-point type random number generator kernel.
 *
 * @param[in, out] x    Random number generator seed - updated on each call.
 * @param[in] a         Random number generator seed - fixed.
 * @param[in] n         Number of random numbers to generate
 * @param buffer        Buffer for the n random numbers.
 *
 */
__kernel void RngKernel(
	unsigned long x, unsigned a, mc_size_t n, __global mc_fp_t *buffer){
		if (get_global_id(0) == 0){
			for(mc_size_t i=0; i < n; ++i)
				buffer[i] = fp_random(&x, a);
		}
};

/*############### End random number generator implementation #################*/


#ifndef __MCCLKERNELS_H
#define __MCCLKERNELS_H

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

#endif /* #define __MCCLKERNELS_H */