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
#if !defined(MC_USE_BALLISTIC_KERNEL) || defined(__DOXYGEN__)
	/** @brief Switches to a ballistic implementation of the Monte Carlo.
	 *  #note Ballistic stepping is about 2-3 fold faster but produces 
	 *	noisy results and should be used only as a reference */
	#define MC_USE_BALLISTIC_KERNEL				FALSE
#endif

#if !defined(MC_USE_DOUBLE_PRECISION) || defined(__DOXYGEN__)
	/** @brief Enable double floating-point precision arithmetics. 
	 *  @note Double precision arithmetic is significatly slower than
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

#if !defined(MC_FLOAT_LUT_ARRAY_MEMORY) || defined(__DOXYGEN__)
	/** @brief Integer lookup table memory space: must be __global or __constant. */
	#define MC_FLOAT_LUT_ARRAY_MEMORY				__global
#endif

#if !defined(MC_USE_PACKED_STRUCTS) || defined(__DOXYGEN__)
	/** @brief Define to TRUE to force packed structurMC_USE_64_BIT_INTsctaes. */
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
/*############## Start suppot for OpenCL compiler directives  ################*/


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

	#define accumulator_deposit(paddress, weight) \
		accu_64_deposit_32((__global mc_accu_t *)(paddress), (uint)(weight))
#else
	typedef uint mc_accu_t;
	#define MC_ACCU_MAX ((uint)0xFFFFFFFF)

	/** @brief Deposit to 32-bit weight to 32-bit accumulators. */
	#define accumulator_deposit(paddress, weight) \
		atomic_add((__gloabal mc_accu_t *)(paddress), (mc_accu_t)(weight))
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

#if MC_USE_64_BIT_SIZE_T || defined(__DOXYGEN__)
    /* Using 64-bit default size type. */
    typedef ulong mc_size_t;
#else
    /* Using 32-bit default size type. */
    typedef uint mc_size_t;
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

#if MC_USE_NATIVE_MATH || defined(__DOXYGEN__)
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


/*################ Start 64-bit atomic accumulator deposit ###################*/
/**
* @addtogroup mc_atomic_accumulator Atomic accumulator deposit
* @{
*/

/**
 * @brief Atomic deposit a 32-bit integer weight to 64-bit accumulator
 * @param[in] address Accumulator addres.
 * @param[in] weight Weight to deposit / add to the accumulator.
 */
inline void accu_64_deposit_32(volatile __global ulong *address, uint weight);

/**
 * @} // end @addtogroup mc_atomic_accumulator
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


/*#################### Start debug support declarations ######################*/
/**
* @addtogroup mc_debug_support Debug support
* @{
*/

/** @brief Default indent string. */
#define INDENT "  "

#if MC_USE_64_BIT_INTEGER
	/** @brief Printf format specifier for default signed integer type (mc_int_t). */ 
	#define FMT_INT		"l"
	/** @brief Printf format specifier for default unsigned integer type (mc_uint_t). */
	#define FMT_UINT	"lu"
#else
	/** @brief Printf format specifier for default signed integer type (mc_uint_t). */
	#define FMT_INT		"d"
	/** @brief Printf format specifier for default unsigned integer type (mc_uint_t). */
	#define FMT_UINT	"u"
#endif

#if MC_USE_64_BIT_SIZE_T
	/** @brief Printf format specifier for default size type (mc_size_t). */ 
	#define FMT_SIZE_T	"lu"
#else
	/** @brief Printf format specifier for default size type (mc_size_t). */ 
	#define FMT_SIZE_T	"u"
#endif

#if MC_USE_64_BIT_PACKET_COUNTER
	/** @brief Printf format specifier for default counter type (mc_cnt_t). */ 
	#define FMT_CNT	"lu"
#else
	/** @brief Printf format specifier for default counter type (mc_cnt_t). */ 
	#define FMT_CNT	"u"
#endif

#if MC_ENABLE_DEBUG || defined(__DOXYGEN__)

/** 
 * @brief Print a formated message.
 * param[in] msg Message to print."
 */
#define dbg_printf	printf

/** 
 * @brief Print a message.
 * param[in] msg Message to print."
 */
#define dbg_print(msg) printf(msg "\n")

/** 
 * @brief Print a floating-point (mc_fp_t) value with label.
 * param[in] label Label that will prefix the floating point value."
 * param[in] value The floating-point value to print.
 */
#define dbg_print_float(label, value) \
	printf(label " %.6f\n", (mc_fp_t)(value))

/** 
 * @brief Print a size type (mc_size_t) value with label.
 * param[in] label Label that will prefix the size type value."
 * param[in] value The size value to print.
 */
#define dbg_print_size_t(label, value) \
	printf(label " %" FMT_SIZE_T "\n", (mc_size_t)(value))

/** 
 * @brief Print an integer (mc_int_t) value with label.
 * param[in] label Label that will prefix the size type value."
 * param[in] value The integer value to print.
 */
#define dbg_print_int(label, value) \
	printf(label " %" FMT_INT "\n", (mc_int_t)(value))

/** 
 * @brief Print an unsigned integer (mc_uint_t) value with label.
 * param[in] label Label that will prefix the size type value."
 * param[in] value The integer value to print.
 */
#define dbg_print_uint(label, value) \
	printf(label " %" FMT_UINT "\n", (mc_uint_t)(value))

/** 
 * @brief Print a counter type (mc_cnt_t) value with label.
 * param[in] label Label that will prefix the size type value."
 * param[in] value The counter value to print.
 */
#define dbg_print_cnt(label, value) \
	printf(label " %" FMT_CNT "\n", (mc_cnt_t)(value))

/**
 * @brief Print a mc_point2f_t value.
 * param[in] label Label that will prefix the size type value."
 * param[in] pvalue Pointer to a mc_point2f_t value / structure.
 */
#define dbg_print_point2f(label, pvalue) \
	printf(label "(%.6f, %.6f)\n", (mc_fp_t)(pvalue)->x, (mc_fp_t)(pvalue)->y);

/**
 * @brief Print a mc_point3f_t value.
 * param[in] label Label that will prefix the size type value."
 * param[in] pvalue Pointer to a mc_point3f_t value / structure.
 */
#define dbg_print_point3f(label, pvalue) \
	printf(label " (%.6f, %.6f, %.6f)\n", (mc_fp_t)(pvalue)->x, \
		(mc_fp_t)(pvalue)->y, (mc_fp_t)(pvalue)->z)

/** 
 * @brief Print a mc_point2_t value.
 * param[in] label Label that will prefix the size type value."
 * param[in] pvalue Pointer to a mc_point2_t value / structure.
 */
#define dbg_print_point2(label, pvalue) \
	printf(label " (%" FMT_INT ", %" FMT_INT ")\n", \
		(mc_int_t)(pvalue)->x, (mc_int_t)(pvalue)->y)

/** 
 * @brief Print a mc_point3_t value.
 * param[in] label Label that will prefix the size type value."
 * param[in] pvalue Pointer to a mc_point3_t value / structure.
 */
#define dbg_print_point3(label, pvalue) \
	printf(label " (%" FMT_INT ", %" FMT_INT ", %" FMT_INT ")\n", \
		(mc_int_t)(pvalue)->x, (mc_int_t)(pvalue)->y, (mc_int_t)(pvalue)->z)

#define dbg_print_matrix2f(label, pm) \
	printf(label " [[%.6f, %.6f], [%.6f, %.6f]]\n", \
		(pm)->a_11, (pm)->a_12, \
		(pm)->a_21, (pm)->a_22)

#define dbg_print_matrix3f(label, pm) \
	printf(label " [[%.6f, %.6f, %.6f], [%.6f, %.6f, %.6f], [%.6f, %.6f, %.6f]]\n", \
		(pm)->a_11, (pm)->a_12, (pm)->a_13, \
		(pm)->a_21, (pm)->a_22, (pm)->a_23, \
		(pm)->a_31, (pm)->a_32, (pm)->a_33)

#else

#define dbg_printf (void)
#define dbg_print(a) ;
#define dbg_print_float(label, value) ;
#define dbg_print_size_t(label, value) ;
#define dbg_print_int(label, value) ;
#define dbg_print_uint(label, value) ;
#define dbg_print_cnt(label, value) ;
#define dbg_print_point3f(label, pvalue) ;
#define dbg_print_point3(label, pvalue) ;
#define dbg_print_point2f(label, pvalue) ;
#define dbg_print_point2(label, pvalue) ;
#define dbg_print_matrix2f(label, pm) ;
#define dbg_print_matrix3f(label, pm) ;

#endif

/**
 * @} // end @addtogroup mc_simulator_core
 */
/*##################### End debug support declarations #######################*/


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
*@brief Data type used to describe a 3D point/vector of size type.
*@{
*/
struct MC_STRUCT_ATTRIBUTES mc_point3s_t{
	mc_size_t x;	/**< @brief x coordinate. */
	mc_size_t y;	/**< @brief y coordinate. */
	mc_size_t z;	/**< @brief z coordinate. */
};
/**
 * @}
 */
/** @brief 3D point of size_t data type. */
typedef struct mc_point3s_t mc_point3s_t;

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
*@brief Data type used to describe a 2D transformation matrix.
*@{
*/
struct mc_matrix2f_t{
	mc_fp_t a_11; /**< @brief first element of the first row */
	mc_fp_t a_12; /**< @brief second element of the first row */
	mc_fp_t a_21; /**< @brief first element of the second row */
	mc_fp_t a_22; /**< @brief second element of the second row */
};
/**
 * @}
 */
/** @brief 2D transformation matrix data type. */
typedef struct mc_matrix2f_t mc_matrix2f_t;

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
 * @brief Reverse the direction of a 3D vector.
 * 
 * param[in] pt1 Pointer to input vector (mc_point3f_t).
 * param[in] pt2 Pointer to output (reversed input) vector (mc_point3f_t).
 * 
 * returns Dot product of the two vectors.
 */
#define reverse3f(v1, v2) \
	{(v2)->x = -(v1)->x; (v2)->y = -(v1)->y; (v2)->z = -(v1)->z;}

/**
 * @brief Reverse the direction of a 2D vector.
 * 
 * param[in] pt1 Pointer to input vector (mc_point3f_t).
 * param[in] pt2 Pointer to output (reversed input) vector (mc_point3f_t).
 * 
 * returns Dot product of the two vectors.
 */
#define reverse2f(v1, v2) \
	{(v2)->x = -(v1)->x; (v2)->y = -(v1)->y;}

/**
 * Transfors the x coordinate of a 3D point by the given 3D transformation
 * matrix.
 * 
 * param[in] pT Pointer to a transformation matrix (mc_matrix3f_t).
 * param[in] pt Pointer to a point (mc_point3f_t) that will be transformed.
 */
#define transform_point3f_x(pT, pt) \
	((pT)->a_11*(pt)->x + (pT)->a_12*(pt)->y + (pT)->a_13*(pt)->z)

/**
 * Transfors the y coordinate of a 3D point by the given 3D transformation
 * matrix.
 * 
 * param[in] pT Pointer to a transformation matrix (mc_matrix3f_t).
 * param[in] pt Pointer to a point (mc_point3f_t) that will be transformed.
 */
#define transform_point3f_y(pT, pt) \
	((pT)->a_21*(pt)->x + (pT)->a_22*(pt)->y + (pT)->a_23*(pt)->z)

/**
 * Transfors the z coordinate of a 3D point by the given 3D transformation
 * matrix.
 * 
 * param[in] pT Pointer to a transformation matrix (mc_matrix3f_t).
 * param[in] pt Pointer to a point (mc_point3f_t) that will be transformed.
 */
#define transform_point3f_z(pT, pt) \
	((pT)->a_31*(pt)->x + (pT)->a_32*(pt)->y + (pT)->a_33*(pt)->z)

/**
 * Transfors a 3D point by the given 3D transformation matrix.
 * 
 * param[in] pT Pointer to a transformation matrix (mc_matrix3f_t).
 * param[in] pt Pointer to a point (mc_point3f_t) that will be transformed.
 * param[out] pres Pointer to a point (mc_point3f_t) that will hold the transformation result.
 */
#define transform_point3f(pT, pt, pres) \
	(pres)->x = transform_point3f_x(pT, pt); \
	(pres)->y = transform_point3f_y(pT, pt); \
	(pres)->z = transform_point3f_z(pT, pt);

/**
 * Transfors the x coordinate of a 2D point by the given 3D transformation
 * matrix.
 * 
 * param[in] pT Pointer to a transformation matrix (mc_matrix2f_t).
 * param[in] pt Pointer to a point (mc_point2f_t) that will be transformed.
 */
#define transform_point2f_x(pT, pt) \
	((pT)->a_11*(pt)->x + (pT)->a_12*(pt)->y)

/**
 * Transfors the y coordinate of a 3D point by the given 3D transformation
 * matrix.
 * 
 * param[in] pT Pointer to a transformation matrix (mc_matrix2f_t).
 * param[in] pt Pointer to a point (mc_point2f_t) that will be transformed.
 */
#define transform_point2f_y(pT, pt) \
	((pT)->a_21*(pt)->x + (pT)->a_22*(pt)->y)

/**
 * Transfors a 2D point by the given 3D transformation matrix.
 * 
 * param[in] pT Pointer to a transformation matrix (mc_matrix2f_t).
 * param[in] pt Pointer to a point (mc_point2f_t) that will be transformed.
 * param[out] pres Pointer to a point (mc_point2f_t) that will hold the transformation result.
 */
#define transform_point2f(pT, pt, pres) \
	(pres)->x = transform_point2f_x(pT, pt); \
	(pres)->y = transform_point2f_y(pT, pt); \

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
inline mc_fp_t point2f_distance_squared(
		const mc_point2f_t * pT1, const mc_point2f_t * pT2);

/**
 * @brief Calculates the Euclidean distance between two points.
 * @param[in] pT1 Pointer to the first point.
 * @param[in] pT2 Pointer to the second point.
 * @return Returns the square of Euclidean distance between points T1 and T2.
 */
#define point2f_distance(pT1, pT2) \
		mc_sqrt(point2f_distance_squared((pT1), (pT2)))

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
	mc_size_t n;
	/**< @brief Offset of the first element in the lookup table. */
	mc_size_t offset;
};
/** @} */
typedef struct mc_fp_lut_t mc_fp_lut_t;

#if !defined(fp_linear_lut_rel_sample) || defined(__DOXYGEN__)
	/**
	 * @brief Sample a floating-point lookup table with a relative index
	 * from 0 to 1. If the index is outside of the valid range, the output
	 * value (pvalue) is not changed / assigned.
	 * @param[in] buffer  The common lookup table buffer.
	 * @param[in] plut    Pointer to a mc_fp_lut_t lookup table instance.
	 * @param in where    A floating-point sampling location from 0.0 to 1.0.
	 * @param[out] pvalue Pointer to a floating-point value for the sample.
	 */
	#define fp_linear_lut_rel_sample(buffer, plut, where, pvalue) \
		{ \
			mc_fp_t _fp_index = (where)*((plut)->n - 1); \
			mc_size_t _index1 = (mc_size_t)(_fp_index + FP_0p5); \
			if (_index1 >= 0 && _index1 < (plut)->n){ \
				mc_fp_t _w2 = _fp_index - mc_floor(_fp_index); \
				mc_size_t _index2 = mc_clip(_index1 + 1, 0, (plut)->n - 1); \
				*pvalue = (buffer)[(plut)->offset + _index1]*(FP_1 - _w2) + \
					(buffer)[(plut)->offset + _index2]*_w2; \
			}; \
		};
#endif

#if !defined(fp_linear_lut_sample) || defined(__DOXYGEN__)
	/**
	 * @brief Sample a floating-point lookup table with a floating-point index
	 * from 0 to n - 1, where n is the size of lookup table.
	 * If the index is outside of the valid range, the output
	 * value (pvalue) is not changed / assigned.
	 * @param[in] buffer  The common lookup table buffer.
	 * @param[in] plut    Pointer to a mc_fp_lut_t lookup table instance.
	 * @param in where    A floating-point sampling location from 0.0 to 1.0.
	 * @param[out] pvalue Pointer to a floating-point value for the sample.
	 */
	#define fp_linear_lut_sample(buffer, plut, where, pvalue) \
		{ \
			mc_fp_t _fp_index = (where); \
			mc_size_t _index1 = (mc_size_t)(_fp_index + FP_0p5); \
			if (_index1 >= 0 && _index1 < (plut)->n){ \
				mc_fp_t _w2 = _fp_index - mc_floor(_fp_index); \
				mc_size_t _index2 = mc_clip(_index1 + 1, 0, (plut)->n - 1); \
				*pvalue = (buffer)[(plut)->offset + _index1]*(FP_1 - _w2) + \
					(buffer)[(plut)->offset + _index2]*_w2; \
			}; \
		};
#endif

#if MC_ENABLE_DEBUG || defined(__DOXYGEN__)
/**
 * @brief Debug print of linear floating point lookup table parameters.
 * @param[in] prefix Cons prefix string.
 * @param plut Pointer to a linear lookup table instance.
 */
#define dbg_print_fp_lut(prefix, plut) \
	printf(prefix " (first=%.6f, inv_span=%.6f, n=%" FMT_SIZE_T ", offset=%" FMT_SIZE_T ")", \
		(plut)->first, (plut)->inv_span, (plut)->n, (plut)->offset);
#else
#define dbg_print_fp_lut(prefix, plut) ;
#endif
/*################### End linear lookup table declarations ###################*/


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
		mc_fp_t n1, mc_fp_t n2, mc_point3f_t *r);


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
inline int refract_safe(
		mc_point3f_t const *p, mc_point3f_t const *n,
		mc_fp_t n1, mc_fp_t n2, mc_point3f_t *r);

/**
 * @} // end @addtogroup mc_boundary_crossing
 */
/*#################### End boundary physics declarations #####################*/



/*############### Start random number generator declarations #################*/
/**
* @addtogroup mc_random Random number generator
* @{
*/

/**
 * @brief Generates a single precision random number from [0.0, 1.0] and
 *        update the generator state.
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
inline mc_fp_t fp_random_single(unsigned long *x, unsigned a);

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
	inline mc_fp_t fp_random_double(unsigned long *x, unsigned a);

	/** @brief Using double precision random generator when ::MC_USE_DOUBLE_PRECISION is defined. **/
	#define fp_random	fp_random_double
#else
	/** @brief Using single precision random generator when ::MC_USE_DOUBLE_PRECISION is NOT TRUE. **/
	#define fp_random	fp_random_single
#endif

/**
 * @} // end @addtogroup mc_random
 */
/*################ End random number generator declarations ##################*/



#endif /* #define __MC_BASE_H */
