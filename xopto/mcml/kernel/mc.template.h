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

#ifndef __MC_BASE_H
#define __MC_BASE_H

#if !defined(TRUE) || defined(__DOXYGEN__)
	/** @brief Logical type true value. */
	#define TRUE 1
#endif

#if !defined(FALSE) || defined(__DOXYGEN__)
	/** @brief Logical type false value. */
	#define FALSE (!(TRUE))
#endif


/*########### Start user-defined simulator options - DO NOT EDIT! ############*/
/* START_MONTE_CARLO_SIMULATOR_OPTIONS_BLOCK */
{{ mc.options or '' }}
/* END_MONTE_CARLO_SIMULATOR_OPTIONS_BLOCK */
/*########### End user-defined simulator options - DO NOT EDIT! ##############*/


/*################## Start simulator configuration options ###################*/
/**
 * @addtogroup mc_simulator_options Simulator options
 * @{
 */

#if !defined(MC_ENABLE_DEBUG) || defined(__DOXYGEN__)
	/** @brief Enables stdout debug information on CPU defices. */
	#define MC_ENABLE_DEBUG						FALSE
#endif
 
/**
 * @addtogroup mc_option Simulation options
 * @{
 */
#if !defined(MC_USE_DOUBLE_PRECISION) || defined(__DOXYGEN__)
	/** @brief Enable double floating-point precision arithmetics. 
	 *  @note Double precision arithmetic is significantly slower than
	 *        the default single precision arithmetics */
	#define MC_USE_DOUBLE_PRECISION				FALSE
#endif

#if !defined(MC_USE_64_BIT_PACKET_COUNTER) || defined(__DOXYGEN__)
	/** @brief Enable 64-bit photon packet counter that allows simulation of
	 *         up to 2^64 photon packets in a kernel single call. 
	 *  @note  This requires 64-bit atomic increment operation that is part of
	 *         the cl_khr_int64_base_atomics extension. To enable the extension
	 *         define cl_khr_int64_base_atomics. If this define is not found,
	 *         the a soft implementation of the 64-bit atomic increment
	 *         operation is used (can slightly degrade the performance).
	 * the default single precision arithmetics */
	#define MC_USE_64_BIT_PACKET_COUNTER		FALSE
#endif

#if !defined(MC_USE_SOFT_64_ATOMICS) || defined(__DOXYGEN__)
	/** @brief Force the use of software implementation of the 64-bit atomic
	 *         operations regardless of the availability of the
	 *         cl_khr_int64_base_atomics OpenCL extension.
	 *  @note It seems that some platforms advertise cl_khr_int64_base_atomics
	 *        but the implementation fails. Setting this flag to TRUE will
	 *        force the use of software-implemented atomic operations even
	 *        if cl_khr_int64_base_atomics is supported by the OpenCL device.
	 */
	#define MC_USE_SOFT_64_ATOMICS				FALSE
#endif

#if !defined(MC_USE_64_BIT_INTEGER) || defined(__DOXYGEN__)
	/** @brief Enable 64-bit default integers. 
	 *  @note  This option is experimental. */
	#define MC_USE_64_BIT_INTEGER					FALSE
#endif 

#if !defined(MC_USE_NATIVE_MATH) || defined(__DOXYGEN__)
	/** @brief Enable native math function calls (divide, sqrt, rsqrt). 
	 *  #note Native math function are generally faster, however, depending on 
	 *	the underlying hardware implementation might be less accurate. */
	#define MC_USE_NATIVE_MATH					FALSE
#endif

#if !defined(MC_USE_HALF_MATH) || defined(__DOXYGEN__)
	/** @brief Enable half math function calls (divide, sqrt, rsqrt). 
	 *  #note Half math function are generally faster, however, depending on 
	 *	the underlying hardware implementation might be insufficiently accurate! */
	#define MC_USE_HALF_MATH 					FALSE
#endif

#if !defined(MC_USE_LOTTERY) || defined(__DOXYGEN__)
	/** @brief Terminate photon packet by lottery. */
	#define MC_USE_LOTTERY						TRUE
#endif 

#if !defined(MC_PACKET_WEIGHT_MIN) || defined(__DOXYGEN__)
	/** @brief Minimum photon packet weight before termination/lottery. */
	#define MC_PACKET_WEIGHT_MIN				(1.0e-4f)
#endif

#if !defined(MC_PACKET_LOTTERY_CHANCE) || defined(__DOXYGEN__)
	/** @brief Terminate photon packet by lottery. If the value of a 
		uniform random number exceeds ::MC_PACKET_LOTTERY_CHANCE, the packet is
		terminated. */
	#define MC_PACKET_LOTTERY_CHANCE			(1e-1f)
#endif

#if !defined(MC_USE_TRACE) || defined(__DOXYGEN__)
	/** @brief Define to a combination of flags ::MC_USE_TRACE_ALL, 
		::MC_USE_TRACE_START and/or MC_USE_TRACE_END. Use ::MC_USE_TRACE_NONE 
		to disable trace functionality.*/
	#define MC_USE_TRACE						FALSE
#endif

#if !defined(MC_USE_FP_LUT) || defined(__DOXYGEN__)
	/** @brief Define to TRUE if floating-point lookup table is used. */
	#define MC_USE_FP_LUT						FALSE
#endif

#if !defined(MC_FP_LUT_ARRAY_MEMORY) || defined(__DOXYGEN__)
	/** @brief Floating-point lookup table memory space: must be __global or __constant. */
	#define MC_FP_LUT_ARRAY_MEMORY				__global
#endif

#if !defined(MC_USE_INT_LUT) || defined(__DOXYGEN__)
	/** @brief Define to TRUE if integer lookup table is used. */
	#define MC_USE_INT_LUT						FALSE
#endif

#if !defined(MC_INT_LUT_ARRAY_MEMORY) || defined(__DOXYGEN__)
	/** @brief Integer lookup table memory space: must be __global or __constant. */
	#define MC_INT_LUT_ARRAY_MEMORY				__global
#endif

#if !defined(MC_LUT_ARRAY_MEMORY) || defined(__DOXYGEN__)
	/** @brief Integer lookup table memory space: must be __global or __constant. */
	#define MC_INT_LUT_ARRAY_MEMORY				__global
#endif

#if !defined(MC_USE_PACKED_STRUCTS) || defined(__DOXYGEN__)
	/** @brief Define to TRUE to force packed structures. */
	#define MC_USE_PACKED_STRUCTS				FALSE
#endif

#if !defined(MC_TRACK_OPTICAL_PATHLENGTH) || defined(__DOXYGEN__)
	/** @brief Define to TRUE to track packet optical pathlength. */
	#define MC_TRACK_OPTICAL_PATHLENGTH			FALSE
#endif

#if !defined(MC_USE_64_BIT_ACCUMULATORS)  || defined(__DOXYGEN__)
	/** @brief Define to TRUE if 64-bit detector accumulators 
		are to be used. */
	#define MC_USE_64_BIT_ACCUMULATORS			TRUE
#endif

#if !defined(MC_INT_ACCUMULATOR_K)  || defined(__DOXYGEN__)
	/** @brief Factor that is used to convert a floating-point photon packet
	 *         weight (value from 0.0 to 1.0) to an integer. */
	#define MC_INT_ACCUMULATOR_K				0x7FFFFF
#endif


#if !defined(MC_USE_FLUENCE) || defined(__DOXYGEN__)
	/** @brief Define to non-zero if fluence data should be collected by the simulation. */
	#define MC_USE_FLUENCE						FALSE
#endif

#if !defined(MC_USE_USER_DATA) || defined(__DOXYGEN__)
	/** @brief Define this macro if user-defined parameters and buffer are passed and used by
		the kernel hooks. */
	#define MC_USE_USER_DATA					FALSE
#endif

#if !defined(MC_USER_DATA_BUFFER_TYPE) || defined(__DOXYGEN__)
	/** @brief User-defined buffer data type. */
	#define MC_USER_DATA_BUFFER_TYPE			mc_fp_t
#endif

#if !defined(MC_USER_DATA_PARAMETERS_TYPE) || defined(__DOXYGEN__)
	/** @brief User-defined parameters data type. */
	#define MC_USER_DATA_PARAMETERS_TYPE		mc_fp_t
#endif

#if !defined(MC_N_USER_PARAMETERS)  || defined(__DOXYGEN__)
	/** @brief Maximum number of user-defined parameters data type. */
	#define MC_N_USER_PARAMETERS				16
#endif

#if MC_USE_PACKED_STRUCTS
	/** @brief Use packed struct definitions if MC_USE_PACKED_STRUCTS is set to TRUE. */
	#define MC_STRUCT_ATTRIBUTES 				__attribute__ ((packed))
#else
	/** @brief Use unpacked struct definitions if MC_USE_PACKED_STRUCTS is set to FALSE. */
	#define MC_STRUCT_ATTRIBUTES
#endif

/**
 * @} // end @addtogroup mc_simulator_options
 */
/*################### End simulator configuration options ####################*/


/*########## Start basic data types, constants and math functions  ###########*/
/**
 * @addtogroup mc_types_constants_math Math, Types and Constants
 * @{
 */

/** @brief Standard integer types. */
typedef ulong uint64_t;
typedef long int64_t;
typedef uint uint32_t;
typedef int int32_t;
typedef ushort uint16_t;
typedef short int16_t;
typedef uchar uint8_t;
typedef char int8_t;


#if MC_USE_DOUBLE_PRECISION || defined(__DOXYGEN__)
	#if defined(cl_khr_fp64)
		/** @brief Enable OpenCL double precision extension. */
		#pragma OPENCL EXTENSION cl_khr_fp64 : enable
	#else
		#error "Double precision floating point not supported by OpenCL implementation."
	#endif
	/** @brief  Define a double precision floating-point literal when ::MC_USE_DOUBLE_PRECISION is TRUE. */
	#define FP_LITERAL(value) value
	/** @brief Using double precision floating-point. */
	typedef double mc_fp_t;
	typedef double2 mc_fp2_t;
	typedef double3 mc_fp3_t;
	typedef double4 mc_fp4_t;
	typedef double8 mc_fp8_t;
	typedef double16 mc_fp16_t;
#else
	/** @brief  Define a single precision floating-point literal when ::MC_USE_DOUBLE_PRECISION is NOT TRUE. */
	#define FP_LITERAL(value) value##f
	/** @brief Using single precision floating-point. */
	typedef float mc_fp_t;
	typedef float2 mc_fp2_t;
	typedef float3 mc_fp3_t;
	typedef float4 mc_fp4_t;
	typedef float8 mc_fp8_t;
	typedef float16 mc_fp16_t;
#endif

/** @brief A short way to define a floating-point literal using the default precision. */
#define FPL(value) FP_LITERAL(value)

#if MC_USE_64_BIT_ACCUMULATORS || defined(__DOXYGEN__)
	typedef ulong mc_accu_t;
	#define MC_ACCU_MAX ((ulong)0xFFFFFFFFFFFFFFFFul)
#else
	typedef uint mc_accu_t;
	#define MC_ACCU_MAX ((uint)0xFFFFFFFF)
#endif

#if MC_USE_64_BIT_PACKET_COUNTER || defined(__DOXYGEN__)
	/** @brief Using a 64-bit unsigned integer as a packet counter when ::MC_USE_64_BIT_PACKET_COUNTER is TRUE. */
	typedef ulong mc_cnt_t;
	/** @brief enable 64-bit atomic extension if available. */
	#if MC_USE_SOFT_64_ATOMICS || !defined(cl_khr_int64_base_atomics)
		/** @brief Using software implementation for atomic incrementation of 64-bit unsigned integers. */
		#define pkt_cnt_atomic_inc(pcnt) atomic_inc_uint64(pcnt)
	#else
		/** @brief Using 64-bit OpenCL atomic extension. */
		#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
		#define pkt_cnt_atomic_inc(pcnt) atom_inc(pcnt)
	#endif
	#define MC_CNT_MAX ((ulong)0xFFFFFFFFFFFFFFFFul)
#else
	/** @brief Using a 32-bit unsigned integer as a packet counter when ::MC_USE_64_BIT_PACKET_COUNTER is NOT TRUE. */
	typedef uint mc_cnt_t;
	/** @brief Using OpenCL implementation of atomic incrementation of 32-bit unsigned integers. */
	#define pkt_cnt_atomic_inc(pcnt) atomic_inc(pcnt)
	#define MC_CNT_MAX ((uint)0xFFFFFFFF)
#endif

#if MC_USE_64_BIT_INTEGER || defined(__DOXYGEN__)
	/** @brief Using 64-bit signed integers when ::MC_USE_64_BIT_INTEGER is TRUE. */
	typedef long mc_int_t;
	/** @brief Using 64-bit signed integer vector[2] when ::MC_USE_64_BIT_INTEGER is TRUE. */
	typedef long2 mc_int2_t;
	/** @brief Using 64-bit signed integer vector[3] when ::MC_USE_64_BIT_INTEGER is TRUE. */
	typedef long3 mc_int3_t;
	/** @brief Using 64-bit signed integer vector[4] when ::MC_USE_64_BIT_INTEGER is TRUE. */
	typedef long4 mc_int4_t;
	/** @brief Using 64-bit signed integer vector[8] when ::MC_USE_64_BIT_INTEGER is TRUE. */
	typedef long8 mc_int8_t;
	/** @brief Using 64-bit signed integer vector[16] when ::MC_USE_64_BIT_INTEGER is TRUE. */
	typedef long16 mc_int16_t;
	
	/** @brief Using 64-bit unsigned integers when ::MC_USE_64_BIT_INTEGER is TRUE. */
	typedef ulong mc_uint_t;
	/** @brief Using 64-bit unsigned integer vector[2] when ::MC_USE_64_BIT_INTEGER is TRUE. */
	typedef ulong2 mc_uint2_t;
	/** @brief Using 64-bit unsigned integer vector[3] when ::MC_USE_64_BIT_INTEGER is TRUE. */
	typedef ulong3 mc_uint3_t;
	/** @brief Using 64-bit unsigned integer vector[4] when ::MC_USE_64_BIT_INTEGER is TRUE. */
	typedef ulong4 mc_uint4_t;
	/** @brief Using 64-bit unsigned integer vector[8] when ::MC_USE_64_BIT_INTEGER is TRUE. */
	typedef ulong8 mc_uint8_t;
	/** @brief Using 64-bit unsigned integer vector[16] when ::MC_USE_64_BIT_INTEGER is TRUE. */
	typedef ulong16 mc_uint16_t;

	#define MC_INT_MAX ((long)0x7FFFFFFFFFFFFFFF)
	#define MC_UINT_MAX ((ulong)0xFFFFFFFFFFFFFFFF)
#else
	/** @brief Using 32-bit signed integers when ::MC_USE_64_BIT_INTEGER is NOT TRUE. */
	typedef int mc_int_t;
	/** @brief Using 32-bit signed integer vector[2] when ::MC_USE_64_BIT_INTEGER is NOT TRUE. */
	typedef int2 mc_int2_t;
	/** @brief Using 32-bit signed integer vector[3] when ::MC_USE_64_BIT_INTEGER is NOT TRUE. */
	typedef int3 mc_int3_t;
	/** @brief Using 32-bit signed integer vector[4] when ::MC_USE_64_BIT_INTEGER is NOT TRUE. */
	typedef int4 mc_int4_t;
	/** @brief Using 32-bit signed integer vector[8] when ::MC_USE_64_BIT_INTEGER is NOT TRUE. */
	typedef int8 mc_int8_t;
	/** @brief Using 32-bit signed integer vector[16] when ::MC_USE_64_BIT_INTEGER is NOT TRUE. */
	typedef int16 mc_int16_t;

	/** @brief Using 32-bit unsigned integers when ::MC_USE_64_BIT_INTEGER is NOT TRUE. */
	typedef uint mc_uint_t;
	/** @brief Using 32-bit unsigned integer vector[2] when ::MC_USE_64_BIT_INTEGER is NOT TRUE. */
	typedef uint2 mc_uint2_t;
	/** @brief Using 32-bit unsigned integer vector[3] when ::MC_USE_64_BIT_INTEGER is NOT TRUE. */
	typedef uint3 mc_uint3_t;
	/** @brief Using 32-bit unsigned integer vector[4] when ::MC_USE_64_BIT_INTEGER is NOT TRUE. */
	typedef uint4 mc_uint4_t;
	/** @brief Using 32-bit unsigned integer vector[8] when ::MC_USE_64_BIT_INTEGER is NOT TRUE. */
	typedef uint8 mc_uint8_t;
	/** @brief Using 32-bit unsigned integer vector[16] when ::MC_USE_64_BIT_INTEGER is NOT TRUE. */
	typedef uint16 mc_uint16_t;

	#define MC_INT_MAX ((long)0x7FFFFFFF)
	#define MC_UINT_MAX ((ulong)0xFFFFFFFF)
#endif

#if MC_USE_DOUBLE_PRECISION || defined(__DOXYGEN__)
	/** @brief Double precision EPS. */
	#define FP_EPS			FP_LITERAL(2.220446049250313e-16)
	/** @brief Maximum integer (4503599627370495) that can be represented by a double precision floating-point number. */
	#define FP_MAX_INT		((uint64_t)0xFFFFFFFFFFFFFul)
#else
	/** @brief Single precision EPS. */
	#define FP_EPS			FP_LITERAL(1.1920929e-07)
	/** @brief Maximum integer (8388607) that can be represented by a single precision floating-point number. */
	#define FP_MAX_INT		(0x7FFFFF)
#endif

/** @brief Minimum radius for logscale radial accumulators. */
#define FP_RMIN			FP_LITERAL(1e-12)
/** @brief Minimum path length for path length accumulators. */
#define FP_PLMIN		FP_LITERAL(1e-12)
/** @brief A floating-point constant: zero. */
#define FP_0			FP_LITERAL(0.0)
/** @brief A floating-point constant: 1/27. */
#define FP_1d27			FP_LITERAL(0.037037037037037035f)
/** @brief A floating-point constant: quarter. */
#define FP_0p25			FP_LITERAL(0.25)
/** @brief A floating-point constant: zero point five. */
#define FP_0p5			FP_LITERAL(0.5)
/** @brief A floating-point constant: one point zero. */
#define FP_1			FP_LITERAL(1.0)
/** @brief A floating-point constant: one point five. */
#define FP_1p5			FP_LITERAL(1.5)
/** @brief A floating-point constant: two point zero. */
#define FP_2			FP_LITERAL(2.0)
/** @brief A floating-point constant: two point five. */
#define FP_2p5			FP_LITERAL(2.5)
/** @brief A floating-point constant: half pi. */
#define FP_HALF_PI		FP_LITERAL(1.5707963267948966)
/** @brief A floating-point constant: pi. */
#define FP_PI			FP_LITERAL(3.141592653589793)
/** @brief A floating-point constant: 2 times pi. */
#define FP_2PI			FP_LITERAL(6.283185307179586)
/** @brief A floating-point constant: cos(30deg). */
#define FP_COS_30		FP_LITERAL(0.8660254037844386f)
/** @brief Cosine of 90 deg. */
#define FP_COS_90		FP_0
/** @brief Cosine of 90 deg. */
#define FP_COS_0		(FP_1 - FP_COS_90)
/** @brief floating-point infinity constant. */
#define FP_INF			INFINITY

#if USE_NATIVE_MATH || defined(__DOXYGEN__)
	/** @brief Native sine function. */
	#define mc_sin(x)					native_sin(x)
	/** @brief Native cosine function. */
	#define mc_cos(x)					native_cos(x)
	/** @brief Native division function. */
	#define mc_fdiv(a,b)				native_divide((mc_fp_t)(a), (mc_fp_t)(b))
	/** @brief Native reciprocal. */
	#define mc_reciprocal(x)			native_recip((mc_fp_t)(x))
	/** @brief Native natural logarithm function. */
	#define mc_log(x)					native_log(x)
	/** @brief Native square root function. */
	#define mc_sqrt(x)					native_sqrt((mc_fp_t)(x))
	/** @brief Native reciprocal square root function. */
	#define mc_rsqrt(x)					native_rsqrt((mc_fp_t)(x))
	/** @brief Native power function, where x >= 0. */
	#define mc_pow(x, y)				native_powr(x, y)
	/** @brief Standard exponential function, where x >= 0. */
	#define mc_exp(x)					native_exp(x)
	/** @brief Native Simultaneous sine and cosine computation. */
	#define mc_sincos(fi, ps, pc)		(*(ps)=sincos(fi, (pc)))
	/* Native sin and cos functions seem to be highly inacurate on many platforms. */
	/*#define mc_sincos(fi, ps, pc)		{*(ps) = native_sin(fi); *(pc) = native_cos(fi);} */
#elif USE_HALF_MATH
	/** @brief Native half precision sine function. */
	#define mc_sin(x)					half_sin(x)
	/** @brief Native half precision cosine function. */
	#define mc_cos(x)					half_cos(x)
	/** @brief Native half precision division function. */
	#define mc_fdiv(a,b)				half_divide((mc_fp_t)(a), (mc_fp_t)(b))
	/** @brief Native half precision reciprocal. */
	#define mc_reciprocal(x)			half_recip((mc_fp_t)(x))
	/** @brief Native half precision natural logarithm function. */
	#define mc_log(x)					half_log(x)
	/** @brief Native half precision square root function. */
	#define mc_sqrt(x)					half_sqrt((mc_fp_t)(x))
	/** @brief Native half precision reciprocal square root function. */
	#define mc_rsqrt(x)					half_rsqrt((mc_fp_t)(x))
	/** @brief Native half precision power function, where x >= 0. */
	#define mc_pow(x, y)				half_powr(x, y)
	/** @brief Half precision exponential function, where x >= 0. */
	#define mc_exp(x)					half_exp(x)
	/** @brief Half precision simultaneous sine and cosine computation. */
	#define mc_sincos(fi, ps, pc)		(*(ps)=sincos(fi, (pc)))
	/* Half sin and cos functions are far too inacurate. */
	/* #define mc_sincos(fi, ps, pc)	{*(ps) = half_sin(fi); *(pc) = half_cos(fi);} */
#else
	/** @brief Standard sine function. */
	#define mc_sin(x)					sin(x)
	/** @brief Standard cosine function. */
	#define mc_cos(x)					cos(x)
	/** @brief Standard single precision division function. */
	#define mc_fdiv(a,b)				((mc_fp_t)(a)/(mc_fp_t)(b))
	/** @brief Native single precision reciprocal. */
	#define mc_reciprocal(x)			mc_fdiv((mc_fp_t)(1), (mc_fp_t)(x))
	/** @brief Standard single precision reciprocal square root function. */
	#define mc_rsqrt(x)					rsqrt((mc_fp_t)(x))
	/** @brief Standard single precision square root function. */
	#define mc_sqrt(x)					sqrt((mc_fp_t)(x))
	/** @brief Standard natural logarithm function. */
	#define mc_log(x)					log(x)
	/** @brief Standard power function, where x >= 0. */
	#define mc_pow(x, y)				powr(x, y)
	/** @brief Standard exponential function, where x >= 0. */
	#define mc_exp(x)					exp(x)
	/** @brief Simultaneous sine and cosine computation. */
	#define mc_sincos(fi, ps, pc)		(*(ps)=sincos(fi, (pc)))
#endif

/** @brief Copy the sign of a floating-point number. */
#define mc_fcopysign(to, from)		(copysign(to, from))
/** @brief Evaluate to integer sign of a floating-point number. */
#define mc_fsign(x)					(((x) >= FP_0) ? 1 : -1)
/** @brief Evaluate to integer sign of an integer number. */
#define mc_sign(x)					(((x) >= 0) ? 1 : -1)
/** @brief Absolute value of a floating-point number. */
#define mc_fabs(x)					fabs((mc_fp_t)(x))
/** @brief Minimum of two floating-point numbers. */
#define mc_fmin(x,y)				fmin((mc_fp_t)(x), (mc_fp_t)(y))
/** @brief Maximum of two floating-point numbers. */
#define mc_fmax(x,y)				fmax((mc_fp_t)(x), (mc_fp_t)(y))
/** @brief Minimum of two integer numbers. */
#define mc_min(x,y)					min((mc_int_t)(x), (mc_int_t)(y))
/** @brief Maximum of two floating-point numbers. */
#define mc_max(x,y)					max((mc_int_t)(x), (mc_int_t)(y))
/** @brief Clip integer value to the specified range. */
#define mc_clip(x, low, high)		mc_min(mc_max((x), (low)), (high))
/** @brief Clip floating-point value to the specified range. */
#define mc_fclip(x, low, high)		mc_fmin(mc_fmax((x), (low)), (high))
/** @brief Compute the cube root of x */
#define mc_cbrt(x)					cbrt((mc_fp_t)(x))
/** @brief Type cast a floating-point value to integer value. */
#define mc_int(x)					convert_int(x)
/** @brief Round a floating-point value towards minus infinity. */
#define mc_round(x)					roundf(x)
/** @brief Round a floating-point value towards the closest integer. */
#define mc_floor(x)					floor(x)
/** @brief Checks if a floating-point number is finite. */
#define mc_isfinite(x)				isinf(x)

/**
 * @} // end @addtogroup mc_types_constants_math
 */
/*########### End basic data types, constants and math functions  ############*/


/*############### Start 64-bit atomic increment declarations #################*/
/**
* @addtogroup mc_atomic Atomic operations
* @{
*/

#if MC_USE_64_BIT_PACKET_COUNTER || defined(__DOXYGEN__)
/**
 * @brief Software implemantation of 64-bit atomic increment operation using
 *        OpenCL 32-bit atomic operations.
 * @param[in] ptr	Pointer to a 64-bit unsigned integer to be incremented.
 * @return			Returns old value at the pointer location.
 */
inline uint64_t atomic_inc_uint64(__global uint64_t *ptr);

/**
 * @} // end @addtogroup mc_atomic
 */
#endif
/*################ End 64-bit atomic increment declarations ##################*/


/*#################### Start vector/matrix declarations ######################*/
/**
* @addtogroup mc_vector_and_matrix_types Vector and matrix data types
* @{
*/

/**
*@brief Data type used to describe a 3D point/vector coordinates in floating-point precision.
*@{
*/
struct MC_STRUCT_ATTRIBUTES mc_point2f_t{
	mc_fp_t x;	/**< @brief x coordinate. */
	mc_fp_t y;	/**< @brief y coordinate. */
};
/**
 * @}
 */
/** @brief 2D data point/vector coordinates data type. */
typedef struct mc_point2f_t mc_point2f_t;

/**
*@brief Data type used to describe a 3D point/vector of integer type.
*@{
*/
struct MC_STRUCT_ATTRIBUTES mc_point2_t{
	mc_int_t x;	/**< @brief x coordinate. */
	mc_int_t y;	/**< @brief y coordinate. */
};
/**
 * @}
 */
/** @brief 2D integer data point/vector coordinates data type. */
typedef struct mc_point2_t mc_point2_t;

/**
*@brief Data type used to describe a 3D point/vector coordinates in floating-point precision.
*@{
*/
struct MC_STRUCT_ATTRIBUTES mc_point3f_t{
	mc_fp_t x;	/**< @brief x coordinate. */
	mc_fp_t y;	/**< @brief y coordinate. */
	mc_fp_t z;	/**< @brief z coordinate. */
};
/**
 * @}
 */
/** @brief 3D data point/vector coordinates data type. */
typedef struct mc_point3f_t mc_point3f_t;

/**
*@brief Data type used to describe a 3D point/vector of integer type.
*@{
*/
struct MC_STRUCT_ATTRIBUTES mc_point3_t{
	mc_int_t x;	/**< @brief x coordinate. */
	mc_int_t y;	/**< @brief y coordinate. */
	mc_int_t z;	/**< @brief z coordinate. */
};
/**
 * @}
 */
/** @brief 3D integer data point/vector coordinates data type. */
typedef struct mc_point3_t mc_point3_t;

/**
*@brief Data type used to describe a 2D rectangle coordinates in floating-point precision.
*@{
*/
struct MC_STRUCT_ATTRIBUTES mc_rectf_t{
	mc_point2f_t top_left;	/**< @brief Coordinates of the top-left edge. */
	mc_fp_t width;	/**< @brief Rectangle width. */
	mc_fp_t height;	/**< @brief Rectangle height. */
};
/**
 * @}
 */
/** @brief Floating-point rectangle data type. */
typedef struct mc_rectf_t mc_rectf_t;

/**
*@brief Data type used to describe a 2D circle in floating-point precision.
*@{
*/
struct MC_STRUCT_ATTRIBUTES mc_circf_t{
	mc_point2f_t center;	/**< @brief Coordinates of the circle center. */
	mc_fp_t r;	/**< @brief Circle radius. */
};
/**
 * @}
 */
/** @brief Floating-point circle data type. */
typedef struct mc_circf_t mc_circf_t;


/**
*@brief Data type used to describe a 2D slot in floating-point precision.
*@{
*/
struct mc_slotf_t{
	mc_point2f_t center; /** < @brief Slot center. */
	mc_fp_t width;		/** < @brief Slot width. */
	mc_fp_t height; 	/** < @brief Slot height. */
};
/**
 * @}
 */
/** @brief Floating-point slot data type. */
typedef struct mc_slotf_t mc_slotf_t;

/**
*@brief Data type used to describe a 3D transformation matrix.
*@{
*/
struct mc_matrix3f_t{
	mc_fp_t a_11; /**< @brief first element of the first row */
	mc_fp_t a_12; /**< @brief second element of the first row */
	mc_fp_t a_13; /**< @brief third element of the first row */
	mc_fp_t a_21; /**< @brief first element of the second row */
	mc_fp_t a_22; /**< @brief second element of the second row */
	mc_fp_t a_23; /**< @brief third element of the second row */
	mc_fp_t a_31; /**< @brief first element of the third row */
	mc_fp_t a_32; /**< @brief second element of the third row */
	mc_fp_t a_33; /**< @brief third element of the third row */
};
/**
 * @}
 */
/** @brief 3D transformation matrix data type. */
typedef struct mc_matrix3f_t mc_matrix3f_t;

/**
 * Transfors a 3D point by the given 3D transformation matrix.
 * 
 * param[in] pT Pointer to a transformation matrix (mc_matrix3f_t).
 * param[in] pt Pointer to a point (mc_point3f_t) that will be transformed.
 * param[out] pres Pointer to a point (mc_point3f_t) that will hold the transformation result.
 */
#define transform_point3f(pT, pt, pres) \
	(pres)->x = (pT)->a_11*(pt)->x + (pT)->a_12*(pt)->y + (pT)->a_13*(pt)->z; \
	(pres)->y = (pT)->a_21*(pt)->x + (pT)->a_22*(pt)->y + (pT)->a_23*(pt)->z; \
	(pres)->z = (pT)->a_31*(pt)->x + (pT)->a_32*(pt)->y + (pT)->a_33*(pt)->z;

/**
 * @brief Computes the dot product of two input vectors.
 * 
 * param[in] pt1 Pointer to the first input vector (mc_point2f_t).
 * param[in] pt2 Pointer to the second input vector (mc_point2f_t).
 * 
 * returns Dot product of the two vectors.
 */
#define dot2f(pt1, pt2) \
	((pt1)->x*(pt2)->x + (pt1)->y*(pt2)->y)

/**
 * @brief Computes the dot product of two input vectors.
 * 
 * param[in] pt1 Pointer to the first input vector (mc_point3f_t).
 * param[in] pt2 Pointer to the second input vector (mc_point3f_t).
 * 
 * returns Dot product of the two vectors.
 */
#define dot3f(pt1, pt2) \
	((pt1)->x*(pt2)->x + (pt1)->y*(pt2)->y + (pt1)->z*(pt2)->z)

/**
 * @brief Computes the cross product of two input vectors.
 * 
 * param[in] pt1 Pointer to the first input vector (mc_point3f_t).
 * param[in] pt2 Pointer to the second input vector (mc_point3f_t).
 * param[out] pres Pointer to the output vector (mc_point3f_t)
 *                 filled with the cross product pt1 x pt2.
 * 
 * returns Dot product of the two vectors.
 */
#define cross3f(pt1, pt2, pres) \
	(pres)->x = (pt1->y*(pt2).z - (pt1)->z*(pt2)->y; \
	(pres)->y = (pt1)->z*(pt2).x - (pt1)->x*(pt2)->z; \
	(pres)->z = (pt1)->x*(pt2).y - (pt1)->y*(pt2)->x;

/**
 * @brief Normalizes vector length to unity. 
 * @param[in, out] pv Pointer to a vector normalized on return.
 */		
void point3f_normalize(mc_point3f_t *pv);

/**
 * Squared radius of a 2D point.
 * @param[in] pt Pointer to the point.
 * @return Returns squared polar radius of the point.
 */
#define point2f_r2(pt) mc_sqrt((pt)->x*(pt)->x + (pt)->y*(pt)->y)

/**
 * Radius of a 2D point.
 * @param[in] pt Pointer to the point.
 * @return Returns polar radius of the point.
 */
#define point2f_r(pt) mc_sqrt(point2f_r2(pt))

/**
 * Squared radius of a 3D point.
 * @param[in] pt Pointer to the point.
 * @return Returns squared polar radius of the point.
 */
#define point3f_r2(pt) mc_sqrt((pt)->x*(pt)->x + (pt)->y*(pt)->y)

/**
 * Radius of a 3D point.
 * @param[in] pt Pointer to the point.
 * @return Returns polar radius of the point.
 */
#define point3f_r(pt) mc_sqrt(point3f_r2(pt))

/**
 * @brief Calculates square of the Euclidean distance between two points.
 * @param[in] pT1 Pointer to the first point.
 * @param[in] pT2 Pointer to the second point.
 * @return Returns the square of Euclidean distance between points T1 and T2.
 */
inline mc_fp_t point3f_distance_squared(
		const mc_point3f_t * pT1, const mc_point3f_t * pT2);

/**
 * @brief Calculates the Euclidean distance between two points.
 * @param[in] pT1 Pointer to the first point.
 * @param[in] pT2 Pointer to the second point.
 * @return Returns the square of Euclidean distance between points T1 and T2.
 */
#define point3f_distance(pT1, pT2) \
		mc_sqrt(point3f_distance_squared((pT1), (pT2)))

/**
 * @brief Calculates square of the Euclidean distance between two 2D points.
 * @param[in] x1 Coordinate x of the first point.
 * @param[in] y1 Coordinate y of the first point.
 * @param[in] x2 Coordinate x the second point.
 * @param[in] x2 Coordinate y the second point.
 * @return Returns the square of Euclidean distance between (x1, y1) and (x2, y2).
 */
inline mc_fp_t xyf_distance_squared(const mc_fp_t x1, mc_fp_t y1,
		const mc_fp_t x2, const mc_fp_t y2);

/**
 * @brief Checks if a rectangle contains the given point.
 * @param[in] prectf Pointer to a rectangle (mc_rectf_t).
 * @param[in] ppoint2f Pointer to a 2D point (mc_point2f_t).
 * @return Returns nonzero if rectangle contains the given point.
 */
#define rectf_contains_pointf(prectf, ppoint2f) \
		rectf_contains_xy(prectf, (ppoint2f)->x, ppoint2f->y)

/**
 * @brief Checks if a rectangle contains the given point.
 * @param[in] prectf Pointer to a rectangle (mc_rectf_t).
 * @param[in] x Coordinate x of the point (mc_fp_t).
 * @param[in] y Coordinate y of the point (mc_fp_t).
 * @return Returns nonzero if rectangle contains the given point.
 */
#define rectf_contains_xy(prectf, x, y) \
		rectf_contains_ex( \
			(prectf)->top_left.x, (prectf)->top_left.y, \
			(prectf)->top_left.width, (prectf)->top_left.height, \
			x, y \
		)

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
		mc_fp_t x, mc_fp_t y);

/**
 * @brief Checks if a rectangle contains the given point.
 * @param[in] prectf Pointer to a rectangle (mc_rectf_t).
 * @param[in] ppoint2f Pointer to a 2D point (mc_point2f_t).
 * @return Returns nonzero if rectangle contains the given point.
 */
#define circf_contains_pointf(pcircf, ppoint2f) \
	circlf_contains_ex( \
		(pcircf)->center.x, (pcircf)->center.y, (pcircf)->r, \
		(ppointf2)->x, (ppoint2f)->y \
	)

/**
 * @brief Checks if a rectangle contains the given point.
 * @param[in] prectf Pointer to a rectangle (mc_rectf_t).
 * @param[in] x Coordinate x of the point (mc_fp_t).
 * @param[in] y Coordinate y of the point (mc_fp_t).
 * @return Returns nonzero if rectangle contains the given point.
 */
#define circf_contains_xy(pcircf, x, y) \
	circlf_contains_ex( \
		(pcircf)->center.x, (pcircf)->center.y, (pcircf)->r, \
		x, y \
	)

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
	mc_fp_t center_x, mc_fp_t center_y, mc_fp_t r, mc_fp_t x, mc_fp_t y);

/**
 * @brief Checks if slot contains a point.
 * @param[in] pslot Pointer to a slot object.
 * @param[in] pt Pointer to a 2D the point (mc_point2f_t).
 * @returns Returns nonzero if the point is within the slot.
 */
#define slotf_contains_pointf(pslot, pt) \
	slot_contains_xy(pslot, (pt)->x, (pt)->y)

/**
 * @brief Checks if slot contains a point.
 * @param[in] pslot Pointer to a slot object.
 * @param[in] x Coordinate x of the point.
 * @param[in] y Coordinate y of the point.
 * @returns Returns nonzero if the point is within the slot.
 */
#define slotf_contains_xy(pslot, x, y) \
	slotf_contains_ex( \
		(pslot)->center.x, (pslot)->center.y, \
		(pslot)->width, (pslot)->height, x, y \
	)

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
	mc_fp_t width, mc_fp_t height, mc_fp_t x, mc_fp_t y);

/**
 * @} // end @addtogroup mc_vector_and_matrix_types
 */
/*##################### End vector/matrix declarations #######################*/


/*################### Start boundary physics declarations ####################*/
/**
* @addtogroup mc_boundary_crossing Boundary crossing
* @{
*/

/**
 * Cosine of the critical angle of incidence beyond which the incident
 * beam is reflected at the boundary n1 => n2.
 * @param[in] n1	Refractive index of the material on the incident side of the 
 * 					boundary.
 * @param[in] n2	Refractive index of the material across the boundary.
 * @return			Critical angle cosine.
 */
inline mc_fp_t cos_critical(mc_fp_t n1, mc_fp_t n2);

/**
 * @brief Computes reflectance for the given boundary conditions.
 * @param[in] n1 Refractive index of the incident layer.
 * @param[in] n2 Refractive index of the layer across the boundary.
 * @param[in] cosIncidence Incidence angle (with respect to the z axis) cosine.
 * @param[in] cosCritical Critical angle cosine for the interface n1 => n2.
 * @return Returns the calculated reflectance probability from [0.0, 1.0].
 */
inline mc_fp_t reflectance(
		mc_fp_t n1, mc_fp_t n2, mc_fp_t cosIncidence, mc_fp_t cosCritical);

/**
 * @brief Computes reflectance for the given boundary conditions by
 * using data from the secondary side.
 * @param[in] n1 Refractive index of the incident layer.
 * @param[in] n2 Refractive index of the layer across the boundary.
 * @param[in] cosIncidence2 Incidence angle (with respect to z axis) cosine.
 * @return Returns the calculated reflectance probability from [0.0, 1.0].
 */
inline mc_fp_t reflectance_cos2(mc_fp_t n1, mc_fp_t n2, mc_fp_t cosIncidence2);


/**
 * Compute cosine of the critical incidence angle beyond which the incident
 * beam is reflected at the boundary of materials with refractive indices
 * n1 and n2.
 * @param[in] n1 Refractive index of the incident medium.
 * @param[in] n2 Refractive index of the material across the boundary.
 * @return Cosine of the incident angle beyond which the incident beam is
 *         reflected at the given boundary.
 */
inline mc_fp_t cos_critical(mc_fp_t n1, mc_fp_t n2);

/**
 * Computes propagation direction of the reflected beam for the given
 * incident propagation direction and boundary normal.
 * @param[in] p Incident propagation direction.
 * @param[in] n Boundary normal.
 * @param[out] r Reflected beam propagation direction filled in on return.
 *
 * @note Reflected beam propagation direction is computed as p - 2*n*(p*n)
 */
inline mc_point3f_t *reflect(
		mc_point3f_t const *p, mc_point3f_t const *n, mc_point3f_t *r);

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
		mc_fp_t n1, mc_fp_t n2, mc_fp_t cos1, mc_point3f_t *r);

/**
 * Computes propagation direction of the refracted beam for the given
 * incident propagation direction, boundary normal and refractive indices
 * of the materials on each side of the boundary.
 * @param[in] p Incident propagation direction.
 * @param[in] n Boundary normal (pointing outwards or inwards).
 * @param[in] n1 Refractive index of the material on the incident side of the boundary.
 * @param[in] n2 Refractive index of the material across the boundary.
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
inline mc_point3f_t *refract(
		mc_point3f_t const *p, mc_point3f_t const *n,
		mc_fp_t n1, mc_fp_t n2, mc_point3f_t *r);

/**
 * @} // end @addtogroup mc_boundary_crossing
 */
/*#################### End boundary physics declarations #####################*/

#endif /* #define __MC_BASE_H */


#ifndef __MCML_H
#define __MCML_H

/* #include "mcbase.template.h" */

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
/** @brief Keep the configuration data of the advanced geometry at top surface of the sample in constant memory. */
#define __mc_geometry_mem		__constant
/** @brief Keep the configuration data of the advanced geometry at bottom surface of the sample in constant memory. */
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


/*###### Start surface reflectance/transmittance detector declarations #######*/
/**
* @addtogroup mc_surface_reflectance_transmittance Photon packet accumulator
* @{
*/

/**
* @brief Data type describing the surface reflectance / transmittance detector.
* @{
*/
struct McDetector;
/**
 * @}
 */
/** @brief Surface reflectance / transmittance detector. */
typedef struct McDetector McDetector;

/* Surface reflectance / transmittance detector structure definition goes here - DO NOT EDIT! */
/* START_DETECTOR_DECLARATION_BLOCK */
{{ detector.declaration or '' }}
/* END_DETECTOR_DECLARATION_BLOCK */

/**
 * @} // end @addtogroup mc_surface_reflectance_transmittance
 */ 
/*####### End surface reflectance/transmittance detector declarations ########*/


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
	
	#if MC_USE_TOP_SURFACE_LAYOTUT || defined(__DOXYGEN__)
		__mc_geometry_mem McComplexSurfaceTop const *top_geometry;		/**< Advanced geometry at the sample top. */
	#endif
	#if MC_USE_BOTTOM_SURFACE_LAYOUT || defined(__DOXYGEN__)
		__mc_geometry_mem McComplexSurfaceBottom const *bottom_geometry;	/**< Advanced geometry at the sample bottom. */
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

	#if MC_USE_DETECTOR || defined(__DOXYGEN__)
		__mc_detector_mem const McDetector *detector;	/**< @brief Reflectance/transmittance detector configuration data. */
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
 * @brief Sets the z coordiante of the photon packet position to the specified value.
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
#define mcsim_reverse_direction_z(psim, pdir) \
	{(psim)->state.direction.z = -(paim)->state.direction.z;}

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
 * @brief Evaluates to the current photon packet weight converted to integer.
 * @param[in] psim Pointer to a simulator instance.  
 */
#define mcsim_int_weight(psim) \
	((psim)->state.weight*MC_INT_ACCUMULATOR_K + FP_0p5)

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
* @brief Evaluates to a pointer to the top layer. 
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_top_layer(psim) (&(psim)->layers[0])

/** 
* @brief Evaluates to a pointer to the bottom layer. 
* @param[in] psim Pointer to a simulator instance.
*/
#define mcsim_bottom_layer(psim) (&(psim)->layers[mcsim_layer_count(psim) - 1])

/**
 * @brief Evaluates to a pointer to the photon packet source object (::McSource).
 * @param[in] psim Pointer to a simulator instance.  
 */
#define mcsim_source(psim) ((psim)->source)

/**
 * @brief Evaluates to the advanced source geometry at the top of the layer stack.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_top_geometry(psim)	((psim)->top_geometry)

/**
 * @brief Evaluates to the advanced source geometry at the bottom of the layer stack.
 * @param[in] psim Pointer to a simulator instance.
 */
#define mcsim_bottom_geometry(psim)	((psim)->bottom_geometry)

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

#if MC_USE_FP_LUT || defined(__DOXYGEN__)
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
 * @brief Evaluates to a pointer to the surface reflectance / transmittance
 *        detector configuration structure.
 * @param[in] psim Pointer to a simulator instance.  
 */
#define mcsim_detector(psim) ((psim)->detector)

/**
 * @brief Deposit the photon packet to the accumulator.
 * 
 * @param mcsim Pointer to a simulator instance.
 *
 * @note The source code of this function is implemented in related python modules.
 */
inline void mcsim_detector_deposit(McSim *mcsim);

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
 * @brief Deposit the photon packet to the accumulator.
 * 
 * @param mcsim Pointer to a simulator instance.
 *
 * @note The source code of this function is implemented in related python modules.
 */
inline void mcsim_detector_deposit(McSim *mcsim);

/**
* @brief Evaluates to a pointer to the fluence configuration structure.
* @param[in] psim Pointer to a simulator instance. 
*/
#define mcsim_fluence(psim) ((psim)->fluence)

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
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] event_count Number of already traced events for this photon packet.
 * @return int Return nonzero to increment the event count.
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
* @brief Evaluates to a pointer to the requested user parameter.
* @param[in] psim Pointer to a simulator instance. 
* @param[in] index User parameter index.
* @note The value of index is checked for valid range.
*/
#define mcSimUserParameter(psim, index) ((psim)->userParameters[index])

/**
* @brief Evaluates to a pointer to the first element of the user-defined buffer.
* @param[in] psim Pointer to a simulator instance. 
*/
#define mcSimUserBuffer(psim) ((psim)->userBuffer)

#define INDENT "  "

#if MC_ENABLE_DEBUG

#if MC_USE_64_BIT_INTEGER
	#define FMT_INT		"l"
	#define FMT_UINT	"lu"
#else
	#define FMT_INT		"d"
	#define FMT_UINT	"u"
#endif

#if MC_USE_64_BIT_PACKET_COUNTER
	#define FMT_CNT	"lu"
#else
	#define FMT_CNT	"u"
#endif

#define dbg_print_status(psim, label) \
	printf("Status: " label "\n"); \
	printf("  Packet index/id: %" FMT_CNT "\n", \
		mcsim_packet_index(psim)); \
	printf("  Position: %.6f, %.6f, %.6f\n", mcsim_position_x(psim), \
			mcsim_position_y(psim), mcsim_position_z(psim)); \
	printf("  Direction: %.6f, %.6f, %.6f\n", mcsim_direction_x(psim), \
			mcsim_direction_y(psim), mcsim_direction_z(psim)); \
	printf("  Weight: %.6f\n", mcsim_weight(psim));

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
		dbg_print_pf(&(player)->pf);

#define dbg_print(a) printf(a)
#define dbg_print_float(label, value) \
	printf(label " %.6" FMT_FP "\n", (mc_fp_t)value);
#define dbg_print_int(label, value) \
	printf(label " %" FMT_INT "\n", (mc_int_t)value);
#define dbg_print_cnt(label, value) \
	printf(label " %" FMT_CNT "\n", (mc_cnt_t)value);
#else

#define dbg_print_status(psim, label) ;
#define dbg_print(a) ;
#define dbg_print_float(label, value) ;
#define dbg_print_int(label, value) ;
#define dbg_print_cnt(label, value) ;

#endif

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
/**@brief Photon packet has escaped the sample layers. */
#define MC_ESCAPED		2
/**
 * @} // end @addtogroup mc_boundary_interaction_outcome
 */

/**
 * @brief Handles layer boundary interactions (refraction/reflection).
 * @param[in] psim Pointer to a simulator instance.
 * @param[in] Next layer nexLayerIndex.
 * @return Returns MC_REFLECTED if the photon packet is reflected from the boundary,
 *         MC_REFRECTED if the photon packet is refracted across the layer boundary,
 *         or MC_ESCAPED if the photon packet escapes the medium through the
 *         topmost or bottommost layers.
 */
inline mc_int_t mcsim_boundary(McSim *psim, mc_int_t nexLayerIndex);

/**
 * @} // end @addtogroup mc_boundary_crossing
 */
/*################ End layer boundary handler declarations ###################*/


/*############# Start advanced sample surface geometry declaration ###########*/
/**
* @addtogroup mc_advanced_surface_geometry Advanced sample surface geometry
* @{
*/

/** @brief This value should be returned by ::mcsim_top_geometry_handler and
 *         ::mcsim_bottom_geometry_handler, if the boundary requires
 *         standard handling with the provided parameters. */
#define MC_SURFACE_GEOMETRY_CONTINUE -1

/** @brief Data type describing advanced geometry at the sample top surface. */
struct McComplexSurfaceTop;
typedef struct McComplexSurfaceTop McComplexSurfaceTop;
/* Definition of the structure of top sample surface geometry goes here - DO NOT EDIT! */
/* START_TOP_GEOMETRY_DECLARATION_BLOCK */
{{ top_geometry.declaration or '' }}
/* END_TOP_GEOMETRY_DECLARATION_BLOCK */

/** @brief Data type describing advanced geometry at the sample bottom surface. */
struct McComplexSurfaceBottom;
typedef struct McComplexSurfaceBottom McComplexSurfaceBottom;
/* Definition of the structure of top bottom surface geometry goes here - DO NOT EDIT! */
/* START_BOTTOM_GEOMETRY_DECLARATION_BLOCK */
{{ bottom_geometry.declaration or '' }}
/* END_BOTTOM_GEOMETRY_DECLARATION_BLOCK */

/**
 * @brief Handles advanced geometry at the top (z = 0) sample surface.
 * 
 * @param[in] psim Pointer to a simulator instance.
 * @param[in, out] n2 Initialized with the refractive index of the surrounding
 *                    medium. Update with a custom value of the refractive
 *                    index on return.
 * @param[out]        Update with the cosine of the critical incidence angle
 *                    for the boundary transition.
 * @return Must return ::MC_SURFACE_GEOMETRY_CONTINUE if the surface
 *         boundary should be processed using the returned values of n2 and cc.
 *         Must return MC_REFLACTED or MC_REFRACTED if the boundary was
 *         fully processed by the call (including updating of the photon
 *         packet propagation direction and weight, and current layer index
 *         if required), the value will be returned by mcsim_boundary.
 */
inline mc_int_t mcsim_top_geometry_handler(McSim *psim, mc_fp_t *n2, mc_fp_t *cc);

/**
 * @brief Handles advanced geometry at the bottom sample surface (z = sample_thickness).
 * 
 * @param[in] psim Pointer to a simulator instance.
 * @param[in, out] n2 Initialized with the refractive index of the surrounding
 *                    medium. Update with a custom value of the refractive
 *                    index on return.
 * @param[out]        Update with the cosine of the critical incidence angle
 *                    for the boundary transition.
 * @return Must return ::MC_SURFACE_GEOMETRY_CONTINUE if the surface
 *         boundary should be processed using the returned values of n2 and cc.
 *         Must return MC_REFLACTED or MC_REFRACTED if the boundary was
 *         fully processed by the call (including updating of the photon
 *         packet propagation direction and weight, and current layer index
 *         if required), the value will be returned by mcsim_boundary.
 */
inline mc_int_t mcsim_bottom_geometry_handler(McSim *psim, mc_fp_t *n2, mc_fp_t *cc);

/**
 * @} // end @addtogroup mc_advanced_surface_geometry
 */
/*############## End advanced sample surface geometry declaration ############*/


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
 * @note Scattering phase functions (their source code) are implemented in python modules.
 */
inline mc_fp_t mcsim_sample_pf(McSim *psim, mc_fp_t *azimuth);

/**
 * @} // end @addtogroup mc_scattering
 */
/*################## End scattering handling declarations ####################*/

/*################## Start linear lookup table declarations ##################*/
/**
* @addtogroup mc_linear_lookup_table Linear lookup table
* @{
*/

/**
 * @brief Linear lookup table configuration data.
 * @details Use the mcsim_sample_linear_lut macro function to sample
 * 			values from the lookup table by linear interpolation.
 * @{
*/
struct MC_STRUCT_ATTRIBUTES mc_fp_lut_t{
	/** @brief Location of the first element in the lookup table. */
	mc_fp_t first;
	/** @brief Inverse of the distance between the first and last element of the lookup table. */
	mc_fp_t inv_span;
	/**< @brief The number of elements in the lookup table. */
	mc_int_t n;
	/**< @brief Offset of the first element in the lookup table. */
	mc_int_t offset;
};
/** @} */
typedef struct mc_fp_lut_t mc_fp_lut_t;

/**
 * @brief Computes value at the given location of a lookup table using linear
 *        interpolation.
 * @note Function returns 0 for any location outside of the lookup table domain.
 * 
 * @param[in, out] mcsim	Pointer to a simulator satate.
 * @param[in] x				Location at which to compute the value of the lookup table.
 * @param[in] plut 			Pointer to a linear lut instance.
 * 
 * @return[in] Value at the given location obtained by linear interpolation.
 */
#define mcsim_sample_linear_lut(mcsim, x, plut) \
	mcsim_sample_linear_lut_(mcsim, ((x) - (plut)->first)*(plut)->inv_span, \
							(plut)->n, (plut)->offset)

/**
 * @brief Computes the value at the given relative location in the lokup table.
 * @note Function returns 0 for any location outside of the lookup table domain.
 * 
 * @param[in, out] mcsim	Pointer to a simulator satate.
 * @param[in] x_rel			Relative location ([0.0, 1.0]).
 * @param[in] n				Size of the lookup table.
 * @param[in] offset		Offset of the first lookup table element.
 * 
 * @return[in] Value at the given location obtained by linear interpolation.
 */
inline mc_fp_t mcsim_sample_linear_lut_(
		McSim *mcsim, mc_fp_t x_rel, mc_int_t n, mc_int_t offset);

/**
 * @} // end @addtogroup mc_linear_lookup_table
 */
/*################## End linear lookup table declarations ####################*/

#endif /* #define __MCML_H */


#ifndef __MCBASE_H
#define __MCBASE_H

/* #include "mcbase.template.h" */

/*############## Start 64-bit atomic increment implementation ################*/
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
		high = atomic_inc(ui32_ptr + 1) + 1;
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
 * @param[in] cosCritical2 Critical angle cosine for the interface n1 => n2.
 * @return Returns the calculated reflectance probability from [0.0, 1.0].
 */
inline mc_fp_t reflectance_cos2(mc_fp_t n1, mc_fp_t n2, mc_fp_t cos2){
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
	if (sin1 <= FP_1){
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
 * Computes propagation direction of the refracted beam for the given
 * incident propagation direction, boundary normal and refractive indices
 * of the materials on each side of the boundary.
 * @param[in] p Incident propagation direction.
 * @param[in] n Boundary normal (pointing outwards or inwards).
 * @param[in] n1 Refractive index of the material on the incident side of the boundary.
 * @param[in] n2 Refractive index of the material across the boundary.
 * @param[in] cos1 Incidence angle cosine.
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

/*#################### End boundary physics implementation ###################*/

#endif /* #define __MCBASE_H */


/* #include "mcml.template.h" */

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
		psim->state.rngStateX = (psim->state.rngStateX & (uint64_t)0xFFFFFFFFUL)*
			(psim->state.rngStateA) + (psim->state.rngStateX >> 32);

		/* Generate a random number (0,1). */
		return mc_fdiv(
			((double)(1 | (psim->state.rngStateX & (uint64_t)0xFFFFFFFFFFFFFUL))), 
			(double)0x10000000000000UL
		);
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
	psim->state.rngStateX = (psim->state.rngStateX & (uint64_t)0xFFFFFFFFUL)*
		(psim->state.rngStateA) + (psim->state.rngStateX >> 32);

	/* Generate a random number (0,1). */
	return mc_fdiv(
		((float)(1 | (unsigned int)(psim->state.rngStateX))), 
		(float)0x100000000UL
	);
};
/*############### End random number generator implementation #################*/

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

	/* get the default critical cosine */
	/* cos_critical = (dir->z < 0) ? 
		mc_layer_cc_bottom(layer) : mc_layer_cc_top(layer); */

	#if MC_USE_TOP_SURFACE_LAYOTUT
	/* advance geometry at the top sample surface */
	if(nextLayerIndex == 0)
	{
		mc_int_t top_res = mcsim_top_geometry_handler(psim, &n2, &cos_critical);
		if(top_res != MC_SURFACE_GEOMETRY_CONTINUE){
			/* boundary handling done/completed by geometry */
			return top_res;
		}
	}else
	#endif
	#if MC_USE_BOTTOM_SURFACE_LAYOUT
	if(nextLayerIndex == mcsim_layer_count(psim) - 1)
	{
		mc_int_t bottom_res = mcsim_bottom_geometry_handler(psim, &n2, &cos_critical);
		if (bottom_res != MC_SURFACE_GEOMETRY_CONTINUE){
			/* boundary handling done/completed by geometry */
			return bottom_res;
		}
	}
	else
	#endif
	{
		/* advanced geometry is off */
		/* get the critical cosine for the boundary */
		if (dir->z < FP_0)
			cos_critical = mc_layer_cc_top(currentLayer);
		else
			cos_critical = mc_layer_cc_bottom(currentLayer);

		n2 = mc_layer_n(nextLayer);
	};

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
/*################ Start photon packet trace implementation ##################*/


/*################ Start scattering handling implementations ##################*/
/**
* @brief Called when the photon packet needs to be scattered.
* @param[in, out] psim Simulator satate.
* @details Function computes and updates the packet propagation direction
*          and associated states photon packet states.
*/
inline void mcsim_scatter(McSim *psim){
	mc_fp_t fi, sinFi, cosFi;
	mc_fp_t cosTheta, sinTheta;
	mc_fp_t px, k;
	mc_fp_t sinTheta_cosFi, sinTheta_sinFi;

	/* get the photon packet propagation direction - manipulate directly */
	mc_point3f_t *dir = mcsim_direction(psim);

	/* sample the scattering phase functions */
	cosTheta = mcsim_sample_pf(psim, &fi);

	sinTheta = mc_sqrt(FP_1 - cosTheta*cosTheta);

	mc_sincos(fi, &sinFi, &cosFi);

	sinTheta_cosFi = sinTheta*cosFi;
	sinTheta_sinFi = sinTheta*sinFi;

	px = dir->x;

	if(mc_fabs(dir->z) >= FP_COS_0){
		dir->x = sinTheta_cosFi;
		dir->y = sinTheta_sinFi;
		dir->z = mc_fcopysign(cosTheta, dir->z*cosTheta); 
	}else{
		k = mc_sqrt(FP_1 - dir->z*dir->z);

		dir->x = mc_fdiv(sinTheta_cosFi*px*dir->z - sinTheta_sinFi*dir->y, k) +
			px*cosTheta;
		dir->y = mc_fdiv(sinTheta_cosFi*dir->y*dir->z + sinTheta_sinFi*px, k) +
			dir->y*cosTheta;
		dir->z = (-sinTheta_cosFi)*k + dir->z*cosTheta;
	};

	/* Single precision can lose unity vector length. */
	#if !defined(MC_DOUBLE_PRECISION)
		point3f_normalize(dir);
	#endif
};
/*################# End scattering handling implementation ###################*/


/*################# Start linear lookup table implementation #################*/
#if MC_USE_FP_LUT || defined(__DOXYGEN__)

#define print_fp_lut(prefix, plut) \
	printf(prefix " (first=%d, n=%d, offset=%d, inv_span=%e)", \
		(plut)->first, (plut)->inv_span, (plut)->n, (plut)->offset);

/**
 * @brief Computes value at the given location of a floating-point lookup table
 *        using linear interpolation.
 * @note Function returns 0 for any location outside of the lookup table domain.
 * 
 * @param[in, out] mcsim	Pointer to a simulator satate.
 * @param[in] x				Location at which to compute the value of the lookup table.
 * @param[in] plut 			Pointer to a linear lookup table instance.
 * 
 * @return[in] Value at the given location obtained by linear interpolation.
 */
#define mcsim_sample_linear_flut(mcsim, x, plut) \
	mcsim_sample_linear_flut_(mcsim, ((x) - (plut)->first)*(plut)->inv_span, \
							(plut)->n, (plut)->offset)

/**
 * @brief Computes value at the given relative location of a floating-point
 *        lookup table using linear interpolation.
 * @note Function returns 0 for any location outside of the lookup table domain.
 * 
 * @param[in, out] mcsim	Pointer to a simulator satate.
 * @param[in] xrel			Relative location (from 0.0 to 1.0) at which to
 * 							compute the value of the lookup table.
 * @param[in] plut 			Pointer to a linear lookup table instance.
 * 
 * @return[in] Value at the given location obtained by linear interpolation.
 */
#define mcsim_relsample_linear_fp_lut(mcsim, xrel, plut) \
	mcsim_sample_linear_flut_(mcsim, (xrel), (plut)->n, (plut)->offset)

/**
 * @brief Computes the value at the given relative location in the lokup table
 *        of floating-point data.
 * @note Function returns 0 for any location outside of the lookup table domain.
 * 
 * @param[in, out] mcsim	Pointer to a simulator satate.
 * @param[in] x_rel			Relative location ([0.0, 1.0]).
 * @param[in] n				Size of the lookup table.
 * @param[in] offset		Offset of the first lookup table element.
 * 
 * @return[in] Value at the given location obtained by linear interpolation.
 */
inline mc_fp_t mcsim_sample_linear_flut_(
		McSim *mcsim, mc_fp_t x_rel, mc_int_t n, mc_int_t offset){
	mc_fp_t fp_index, w;
	mc_int_t index;

	if (x_rel >= FP_0 && x_rel <= FP_1){
		fp_index = x_rel*(n - 1);
		index = mc_int(fp_index);
		w = fp_index - index;
		return mcsim_fp_lut_array(mcsim)[offset + index]*(FP_1 - w) +
			mcsim_fp_lut_array(mcsim)[offset + mc_min(index + 1, n - 1)]*w;
	}
	return FP_0;
};

#endif /* #if MC_USE_FP_LUT */
/*################# End linear lookup table implementation ###################*/


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


/*################ Start photon packet trace implementation ##################*/
/* User-defined surface reflectance/transmittance detector implementation
   goes here - DO NOT EDIT! */
/* START_DETECTOR_IMPLEMENTATION_BLOCK */
{{ detector.implementation or '' }}
/* END_DETECTOR_IMPLEMENTATION_BLOCK */
/*################# End photon packet trace implementation ###################*/


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


/*############ Start top sample surface geometry implementation ##############*/
/* User-defined top sample surface geometry implementation
   goes here - DO NOT EDIT! */
/* START_TOP_GEOMETRY_IMPLEMENTATION_BLOCK */
{{ top_geometry.implementation or '' }}
/* END_TOP_GEOMETRY_IMPLEMENTATION_BLOCK */
/*############# End top sample surface geometry implementation ###############*/


/*########### Start bottom sample surface geometry implementation ############*/
/* User-defined bottom sample surface geometry implementation
   goes here - DO NOT EDIT! */
/* START_BOTTOM_GEOMETRY_IMPLEMENTATION_BLOCK */
{{ bottom_geometry.implementation or '' }}
/* END_BOTTOM_GEOMETRY_IMPLEMENTATION_BLOCK */
/*############ End bottom sample surface geometry implementation #############*/


/* __BOTTOM_GEOMETRY_IMPLEMENTATION_BLOCK__ */


/*#################### Start Monte Carlo OpenCL kernel #######################*/
__kernel void McKernel(
	mc_cnt_t n_packets,
	__global mc_cnt_t *n_packets_done,

	__global uint32_t *num_kernels,

	mc_fp_t mc_rmax,

	__global uint64_t *rng_state_x,
	__global uint32_t const *rng_state_a,	// must nut be constant - insufficient memory on GPU

	uint32_t num_layers,
	__constant McLayer const *layers,

	__mc_source_mem McSource const *source,
	__mc_geometry_mem McComplexSurfaceTop const *top_geometry,
	__mc_geometry_mem McComplexSurfaceBottom const *bottom_geometry,

	__mc_trace_mem const McTrace *trace,

	__mc_fluence_mem McFluence const *fluence,

	__mc_detector_mem McDetector const *detector,

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
		
		#if MC_USE_TOP_SURFACE_LAYOTUT
			,top_geometry			/**< Advanced geometry at the sample top. */
		#endif
		#if MC_USE_BOTTOM_SURFACE_LAYOUT
			,bottom_geometry		/**< Advanced geometry at the sample bottom. */
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

		#if MC_USE_DETECTOR
			,detector				/* McDetector : Reflectance/transmittance configuration structure. */
		#endif

		,integer_buffer		/* __global mc_int_t *fp_buffer: Common integer buffer. */
		,float_buffer		/* __global mc_fp_t *fp_buffer: Common floating-point buffer. */
		,accumulator_buffer	/* __global mc_accu_t *accumulator_buffer: Common accumulator buffer. */
	};

	/* make unused variables void to surpresss compiler warnings */
	#if !MC_USE_TOP_SURFACE_LAYOTUT
		(void)top_geometry;
	#endif
	#if !MC_USE_BOTTOM_SURFACE_LAYOUT
		(void)bottom_geometry;
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

	#if !MC_USE_DETECTOR
		(void)detector;
	#endif

	/* check if photon packets have to be launched by this thread */
	if((sim.state.photon_index = pkt_cnt_atomic_inc(n_packets_done)) < n_packets){

		atomic_inc(num_kernels);

		#if MC_ENABLE_DEBUG
		/* print debug information ony for the first photon packet */
		if (sim.state.photon_index == 0){

			printf("\nSimulation:\n"
					INDENT "n_packets: %" FMT_CNT "\n"
					INDENT "n_packets_done: %" FMT_CNT "\n"
					INDENT "num_layers: %u\n",
				n_packets, *n_packets_done, num_layers);

			for (int i=0; i < sim.num_layers; ++i){
				dbg_print_int("Layer:", i);
				dbg_print_layer(&sim.layers[i], INDENT);
			};

			dbg_print_source(source);

			#if MC_USE_TRACE
				dbg_print_trace(trace);
			#endif
			#if MC_USE_FLUENCE
				dbg_print_fluence(fluence);
			#endif
			#if MC_USE_DETECTOR
				dbg_print_detector(detector);
			#endif
			printf("Test:\n\t(float)0x100000000UL: %.9f\n"
					"\n\tInitial rng states X: %lu A: %u"
					"\n\tTest mcsim_random() %.9f"
					"\n\tNew rng states X: %lu A: %u\n",
					(float)0x100000000,
					sim.state.rngStateX, sim.state.rngStateA,
					mcsim_random(&sim),
					sim.state.rngStateX, sim.state.rngStateA);
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
		while(!done){

			/* generate a new step */
			step = -mc_log(mcsim_random(&sim))*
				mc_layer_inv_mut(mcsim_current_layer(&sim));
			step = mc_fmin(step, FP_1);

			/* initialize the next layer index with the current layer index */ 
			nexLayerIndex = mcsim_current_layer_index(&sim);

			/* check if the photon packet hits the top layer boundary and if so,
				adjust the step size - go only to the layer boundary */
			if (mcsim_position_z(&sim) + step*mcsim_direction_z(&sim) <
					mc_layer_top(mcsim_current_layer(&sim))){
				--nexLayerIndex;
				if (mc_fabs(mcsim_direction_z(&sim)) != FP_0){
					step = mc_fdiv( mc_layer_top(mcsim_current_layer(&sim)) - 
						mcsim_position_z(&sim), mcsim_direction_z(&sim) );
				};
			};
			/* check if the photon packet hits the bottom layer boundary and if so,
				adjust the step size - go only to the layer boundary */
			if (mcsim_position_z(&sim) + 
					step*mcsim_direction_z(&sim) >=
					mc_layer_bottom(mcsim_current_layer(&sim))){
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

			/* process boundary hit or handle absorption */
			if (nexLayerIndex != mcsim_current_layer_index(&sim)){
				/* deal with the boundary hit ... returns MC_REFRACTED or MC_REFLECTED */
				mcsim_boundary(&sim, nexLayerIndex);
				/* check if photon escaped the sample */
				if ( (mcsim_current_layer_index(&sim) <= 0) || 
					(mcsim_current_layer_index(&sim) >= mcsim_layer_count(&sim) - 1) ){
						#if MC_USE_DETECTOR
							mcsim_detector_deposit(&sim);
						#endif
						
						done = true; /* photon escaped the sample */
					}

				#if MC_USE_TRACE == MC_USE_TRACE_ALL
					mcsim_trace_this_event(&sim);
				#endif
			}else{
				/* Do absorption only when no layer boundary has been hit.*/
				deposit = mcsim_weight(&sim)*
					mc_layer_mua_inv_mut(mcsim_current_layer(&sim));
				mcsim_adjust_weight(&sim, deposit);
				#if MC_USE_FLUENCE
					/* update the fluence data in fluence mode */
					mcsim_fluence_deposit(&sim, deposit);
				#endif

				/* Scatter the photon packet, */
				mcsim_scatter(&sim);
				/* Perform survival lottery if required */
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
				#if MC_USE_TRACE == MC_USE_TRACE_ALL
					mcsim_trace_this_event(&sim);
				#endif
			}

			/* check if photon escaped the predefined simulation domain */
			if (point3f_distance_squared(mcsim_position(&sim), &src_pos) >
					mc_rmax*mc_rmax){
				done = true;
			};

			if(done){
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

				if((sim.state.photon_index = pkt_cnt_atomic_inc(n_packets_done)) < n_packets){
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
