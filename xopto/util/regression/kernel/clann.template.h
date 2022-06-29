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

#ifndef __MC_OPTIONS_H
#define __MC_OPTIONS_H

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
	/** @brief Enables stdout debug information on CPU devices. */
	#define MC_ENABLE_DEBUG						FALSE
#endif
 
/**
 * @addtogroup mc_option Simulation options
 * @{
 */

#if !defined(MC_METHOD) || defined(__DOXYGEN__)
	/** @brief Set the Stepping method. Allowed values are
	 *         0 (albedo weight - default),
	 *         1 (albedo rejection) or
	 *         2 (microscopic beer-lambert). */
	#define MC_METHOD							0
#endif

#if !defined(MC_USE_DOUBLE_PRECISION) || defined(__DOXYGEN__)
	/** @brief Enable double floating-point precision arithmetics. 
	 *  @note Double precision arithmetic is significantly slower than
	 *        the default single precision arithmetics */
	#define MC_USE_DOUBLE_PRECISION				FALSE
#endif

#if !defined(MC_USE_64_BIT_SIZE_T) || defined(__DOXYGEN__)
	/** @brief Enable 64-bit defaut size type. Note that this setting is 
     *         independent of the OpenCL size_t type. */
	#define MC_USE_64_BIT_SIZE_T				FALSE
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
	#define MC_USE_64_BIT_INTEGER				FALSE
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

#if !defined(MC_FLOAT_LUT_ARRAY_MEMORY) || defined(__DOXYGEN__)
	/** @brief Integer lookup table memory space: must be __global or __constant. */
	#define MC_FLOAT_LUT_ARRAY_MEMORY			__global
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

#if !defined(MC_USE_FLUENCE_CACHE) || defined(__DOXYGEN__)
	/** @brief Enables fluence accumulator cache. */
	#define MC_USE_FLUENCE_CACHE				FALSE
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

/** @brief Albedo Weight Monte Carlo simulation method. */
#define ALBEDO_WEIGHT				0
/** @brief Albedo Rejection Monte Carlo simulation method. */
#define ALBEDO_REJECTION			1
/** @brief Microscopic Ber-Lambert Monte Carlo simulation method. */
#define MICROSCOPIC_BEER_LAMBERT	2

/**
 * @} // end @addtogroup mc_simulator_options
 */
/*################### End simulator configuration options ####################*/


/*############## Start suppot for OpenCL compiler directives  ################*/
/**
 * @addtogroup mc_opencl_pragma OpenCL compiler directives
 * @{
 */

#if __OPENCL_C_VERSION__ >= 200
	/** @brief Standard OpenCL pragma for unrolling loops valid since 2.0. */
	#define pragma_unroll_hint(n)		__attribute__((opencl_unroll_hint(n)))
#else
	#if defined(__NV_CL_C_VERSION)
		/** @brief Loop unrolling for Nvidia devices */
		#if defined(cl_nv_pragma_unroll)
			/* Enable Nvidia OpenCL extension for unrolling loops. */
			#pragma OPENCL EXTENSION cl_nv_pragma_unroll : enable
			#define pragma_message(x)			_Pragma(#x) 
			#define pragma_unroll_hint(n) 		pragma_message("unroll" n)
		#else
			/** @brief Loop unrolling extension not supported! */
			#define pragma_unroll_hint(n)	;
		#endif
	#elif defined(INTELFPGA_CL) || defined(cl_amd_device_attribute_query)
		/** @brief Loop unrolling for Intel and AMD OpenCL devices. */
		#define pragma_message(x)			_Pragma(#x) 
		#define pragma_unroll_hint(n) 		pragma_message("unroll" n)
	#else
		/** @brief Unknown OpenCL device ... no loop unrolling pragma */
		#define pragma_unroll_hint(n)	;
	#endif
#endif

/**
 * }@ // end mc_opencl_pragma
 */
/*############### End suppot for OpenCL compiler directives  #################*/

#endif /* #define __MC_OPTIONS_H */


#ifndef __MC_TYPES_H
#define __MC_TYPES_H


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
	/** @brief Using a 64-bit unsigned integer as a packet counter when
	 * ::MC_USE_64_BIT_PACKET_COUNTER is TRUE. */
	typedef ulong mc_cnt_t;

	/** @brief Maximum value that can be represented by the counter type. */
	#define MC_CNT_MAX ((ulong)0xFFFFFFFFFFFFFFFFul)
#else
	/** @brief Using a 32-bit unsigned integer as a packet counter when
	 * ::MC_USE_64_BIT_PACKET_COUNTER is NOT TRUE. */
	typedef uint mc_cnt_t;

	/** @brief Maximum value that can be represented by the counter type. */
	#define MC_CNT_MAX ((uint)0xFFFFFFFF)
#endif

#if MC_USE_64_BIT_SIZE_T || defined(__DOXYGEN__)
    /* Using 64-bit default size type. */
    typedef ulong mc_size_t;
	typedef ulong2 mc_size2_t;
	typedef ulong3 mc_size3_t;
	typedef ulong4 mc_size4_t;
	typedef ulong8 mc_size8_t;
	typedef ulong16 mc_size16_t;
#else
    /* Using 32-bit default size type. */
    typedef uint mc_size_t;
	typedef uint2 mc_size2_t;
	typedef uint3 mc_size3_t;
	typedef uint4 mc_size4_t;
	typedef uint8 mc_size8_t;
	typedef uint16 mc_size16_t;
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
	/** @brief Single precision inverse EPS. */
	#define FP_INV_EPS		FP_LITERAL(4503599627370496.0)
	/** @brief Maximum integer (4503599627370495) that can be represented by a double precision floating-point number. */
	#define FP_MAX_INT		((uint64_t)0xFFFFFFFFFFFFFul)
#else
	/** @brief Single precision EPS. */
	#define FP_EPS			FP_LITERAL(1.1920929e-07)
	/** @brief Single precision inverse EPS. */
	#define FP_INV_EPS		FP_LITERAL(8388608.0)
	/** @brief Maximum integer (8388607) that can be represented by a single precision floating-point number. */
	#define FP_MAX_INT		(0x7FFFFF)
#endif

/** @brief Minimum radius for logscale radial accumulators. */
#define FP_RMIN			FP_LITERAL(1e-12)
/** @brief Minimum radius for logscale optical path length accumulators. */
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
/** @brief A floating-point constant: two point zero. */
#define FP_4			FP_LITERAL(4.0)
/** @brief A floating-point constant: half pi. */
#define FP_HALF_PI		FP_LITERAL(1.5707963267948966)
/** @brief A floating-point constant: pi. */
#define FP_PI			FP_LITERAL(3.141592653589793)
/** @brief A floating-point constant: 2 times pi. */
#define FP_2PI			FP_LITERAL(6.283185307179586)
/** @brief A floating-point constant: cos(30deg). */
#define FP_COS_30		FP_LITERAL(0.8660254037844386)
/** @brief Cosine of 90 deg. */
#define FP_COS_90		FP_0
/** @brief Cosine of 90 deg. */
#define FP_COS_0		(FP_1 - FP_COS_90)
/** @brief floating-point infinity constant. */
#define FP_INF			INFINITY
/** @brief Conversion from radians to degrees. */
#define FP_RAD2DEG		FP_LITERAL(57.2957795130823229)
/** @brief Conversion from degrees to radians. */
#define FP_DEG2RAD		FP_LITERAL(0.017453292519943295)
/** @brief Speed of light in vacuum (m/s). */
#define FP_C			FP_LITERAL(299792458.0)
/** @brief Inverse of the speed of light in vacuum (m/s). */
#define FP_INV_C		FP_LITERAL(3.3356409519815204e-09)

#if MC_USE_NATIVE_MATH || defined(__DOXYGEN__)
	/** @brief Native tangent function. */
	#define mc_tan(x)					native_tan(x)
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
#elif MC_USE_HALF_MATH
	/** @brief Native half precision tangent function. */
	#define mc_tan(x)					half_tan(x)
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
	/** @brief Standard tangent function. */
	#define mc_tan(x)					tan(x)
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

/** @brief Arcus sin(x). */
#define mc_asin(x)					asin((mc_fp_t)(x))
/** @brief Arcus cos(x). */
#define mc_acos(x)					acos((mc_fp_t)(x))
/** @brief Arcus tan(x). */
#define mc_atan(x)					atan((mc_fp_t)(x))
/** @brief Hyperbolic tangent function. */
#define mc_tanh(x)					tanh(x)
/** @brief Arcus tangens(y/x). */
#define mc_atan2(y, x)				atan2((mc_fp_t)(y), (mc_fp_t)(x))
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
#define mc_clip(x, low, high)		clamp((mc_int_t)x, (mc_int_t)low, (mc_int_t)high)
/** @brief Clip floating-point value to the specified range. */
#define mc_fclip(x, low, high)		clamp((mc_fp_t)x, (mc_fp_t)low, (mc_fp_t)high)
/** @brief Compute the square of x */
#define mc_fsquare(x)				((mc_fp_t)(x)*(mc_fp_t)(x))
/** @brief Compute the cube root of x */
#define mc_cbrt(x)					cbrt((mc_fp_t)(x))
/** @brief Type cast a floating-point value to integer value. */
#define mc_int(x)					convert_int(x)
/** @brief Type cast a floating-point value to unsigned integer value. */
#define mc_uint(x)					convert_uint(x)
/** @brief Round a floating-point value towards minus infinity. */
#define mc_round(x)					roundf(x)
/** @brief Round a floating-point value towards the closest integer. */
#define mc_floor(x)					floor(x)
/** @brief Checks if a floating-point number is finite. */
#define mc_isfinite(x)				(!isinf(x))

/**
 * @} // end @addtogroup mc_types_constants_math
 */
/*########### End basic data types, constants and math functions  ############*/

#endif /* __MC_TYPES_H */


/*################## Start simulator configuration options ###################*/
/*################### End simulator configuration options ####################*/


/*########## Start basic data types, constants and math functions  ###########*/
/*########### End basic data types, constants and math functions  ############*/