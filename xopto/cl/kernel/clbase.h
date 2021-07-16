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

#ifndef CLBASE_DEFINES_
#define CLBASE_DEFINES_

/*########## Start basic data types, constants and math functions  ###########*/
/**
 * @addtogroup cl_types_constants_math Math, Types and Constants
 * @{
 */

/** @brief 8-bit signed integer type. */
typedef char int8_t;
typedef char2 int8v2_t;
typedef char3 int8v3_t;
typedef char4 int8v4_t;
typedef char8 int8v8_t;
typedef char16 int8v16_t;
/** @brief 8-bit unsigned integer type. */
typedef unsigned char uint8_t;
typedef uchar2 uint8v2_t;
typedef uchar3 uint8v3_t;
typedef uchar4 uint8v4_t;
typedef uchar8 uint8v8_t;
typedef uchar16 uint8v16_t;
/** @brief 16-bit signed integer type. */
typedef short int16_t;
typedef short2 int16v2_t;
typedef short3 int16v3_t;
typedef short4 int16v4_t;
typedef short8 int16v8_t;
typedef short16 int16v16_t;
/** @brief 16-bit unsigned integer type. */
typedef ushort uint16_t;
typedef ushort2 uint16v2_t;
typedef ushort3 uint16v3_t;
typedef ushort4 uint16v4_t;
typedef ushort8 uint16v8_t;
typedef ushort16 uint16v16_t;
/** @brief 32-bit signed integer type. */
typedef int int32_t;
typedef int2 int32v2_t;
typedef int3 int32v3_t;
typedef int4 int32v4_t;
typedef int8 int32v8_t;
typedef int16 int32v16_t;
/** @brief 32-bit unsigned signed integer type. */
typedef uint uint32_t;
typedef uint2 uint32v2_t;
typedef uint3 uint32v3_t;
typedef uint4 uint32v4_t;
typedef uint8 uint32v8_t;
typedef uint16 uint32v16_t;
/** @brief 64-bit signed integer type. */
typedef long int64_t;
typedef long2 int64v2_t;
typedef long3 int64v3_t;
typedef long4 int64v4_t;
typedef long8 int64v8_t;
typedef long16 int64v16_t;
/** @brief 64-bit unsigned integer type. */
typedef ulong uint64_t;
typedef ulong2 uint64v2_t;
typedef ulong3 uint64v3_t;
typedef ulong4 uint64v4_t;
typedef ulong8 uint64v8_t;
typedef ulong16 uint64v16_t;

#if USE_DOUBLE_PRECISION || defined(__DOXYGEN__)
	/** @brief  Define a double precision floating-point literal when ::USE_DOUBLE_PRECISION is TRUE. */
	#define FP_LITERAL(value) value
	/** @brief Using double precision floating-point. */
	typedef double fp_t;
	typedef double2 fpv2_t;
	typedef double3 fpv3_t;
	typedef double4 fpv4_t;
	typedef double8 fpv8_t;
	typedef double16 fpv16_t;
#else
	/** @brief  Define a single precision floating-point literal when ::USE_DOUBLE_PRECISION is NOT TRUE. */
	#define FP_LITERAL(value) value##f
	/** @brief Using single precision floating-point. */
	typedef float fp_t;
	typedef float2 fpv2_t;
	typedef float3 fpv3_t;
	typedef float4 fpv4_t;
	typedef float8 fpv8_t;
	typedef float16 fpv16_t;
#endif

/** @brief A short way to define a floating-point literal using the default precision. */
#define FPL(value) FP_LITERAL(value)

#if USE_64_BIT_INTEGER || defined(__DOXYGEN__)
	/** @brief Using 64-bit signed integers when ::USE_64_BIT_INTEGER is TRUE. */
	typedef int64_t int_t;
	/** @brief Using 64-bit signed integer vector[2] when ::USE_64_BIT_INTEGER is TRUE. */
	typedef int64v2_t intv2_t;
	/** @brief Using 64-bit signed integer vector[3] when ::USE_64_BIT_INTEGER is TRUE. */
	typedef int64v3_t intv3_t;
	/** @brief Using 64-bit signed integer vector[4] when ::USE_64_BIT_INTEGER is TRUE. */
	typedef int64v4_t intv4_t;
	/** @brief Using 64-bit signed integer vector[8] when ::USE_64_BIT_INTEGER is TRUE. */
	typedef int64v8_t intv8_t;
	/** @brief Using 64-bit signed integer vector[16] when ::USE_64_BIT_INTEGER is TRUE. */
	typedef int64v16_t intv16_t;
	
	/** @brief Using 64-bit unsigned integers when ::USE_64_BIT_INTEGER is TRUE. */
	typedef uint64_t uint_t;
	/** @brief Using 64-bit unsigned integer vector[2] when ::USE_64_BIT_INTEGER is TRUE. */
	typedef uint64v2_t uintv2_t;
	/** @brief Using 64-bit unsigned integer vector[3] when ::USE_64_BIT_INTEGER is TRUE. */
	typedef uint64v3_t uintv3_t;
	/** @brief Using 64-bit unsigned integer vector[4] when ::USE_64_BIT_INTEGER is TRUE. */
	typedef uint64v4_t uintv4_t;
	/** @brief Using 64-bit unsigned integer vector[8] when ::USE_64_BIT_INTEGER is TRUE. */
	typedef uint64v8_t uintv8_t;
	/** @brief Using 64-bit unsigned integer vector[16] when ::USE_64_BIT_INTEGER is TRUE. */
	typedef uint64v16_t uintv16_t;
#else
	/** @brief Using 32-bit signed integers when ::USE_64_BIT_INTEGER is NOT TRUE. */
	typedef int32_t int_t;
	/** @brief Using 32-bit signed integer vector[2] when ::USE_64_BIT_INTEGER is NOT TRUE. */
	typedef int32v2_t intv2_t;
	/** @brief Using 32-bit signed integer vector[3] when ::USE_64_BIT_INTEGER is NOT TRUE. */
	typedef int32v3_t intv3_t;
	/** @brief Using 32-bit signed integer vector[4] when ::USE_64_BIT_INTEGER is NOT TRUE. */
	typedef int32v4_t intv4_t;
	/** @brief Using 32-bit signed integer vector[8] when ::USE_64_BIT_INTEGER is NOT TRUE. */
	typedef int32v8_t intv8_t;
	/** @brief Using 32-bit signed integer vector[16] when ::USE_64_BIT_INTEGER is NOT TRUE. */
	typedef int32v16_t intv16_t;

	/** @brief Using 32-bit unsigned integers when ::USE_64_BIT_INTEGER is NOT TRUE. */
	typedef uint32_t uint_t;
	/** @brief Using 32-bit unsigned integer vector[2] when ::USE_64_BIT_INTEGER is NOT TRUE. */
	typedef uint32v2_t uintv2_t;
	/** @brief Using 32-bit unsigned integer vector[3] when ::USE_64_BIT_INTEGER is NOT TRUE. */
	typedef uint32v3_t uintv3_t;
	/** @brief Using 32-bit unsigned integer vector[4] when ::USE_64_BIT_INTEGER is NOT TRUE. */
	typedef uint32v4_t uintv4_t;
	/** @brief Using 32-bit unsigned integer vector[8] when ::USE_64_BIT_INTEGER is NOT TRUE. */
	typedef uint32v8_t uintv8_t;
	/** @brief Using 32-bit unsigned integer vector[16] when ::USE_64_BIT_INTEGER is NOT TRUE. */
	typedef uint32v16_t uintv16_t;
#endif

/** @brief Default floating-point type. */
typedef fp_t cl_fp_t;
/** @brief Default floating-point vector[2] type. */
typedef fpv2_t cl_fpv2_t;
/** @brief Default floating-point vector[3] type. */
typedef fpv3_t cl_fpv3_t;
/** @brief Default floating-point vector[4] type. */
typedef fpv4_t cl_fpv4_t;
/** @brief Default floating-point vector[8] type. */
typedef fpv8_t cl_fpv8_t;
/** @brief Default floating-point vector[16] type. */
typedef fpv16_t cl_fpv16_t;

/** @brief Default signed integer type. */
typedef int_t cl_int_t;
/** @brief Default signed integer vector[2] type. */
typedef intv2_t cl_intv2_t;
/** @brief Default signed integer vector[3] type. */
typedef intv3_t cl_intv3_t;
/** @brief Default signed integer vector[4] type. */
typedef intv4_t cl_intv4_t;
/** @brief Default signed integer vector[8] type. */
typedef intv8_t cl_intv8_t;
/** @brief Default signed integer vector[16] type. */
typedef intv16_t cl_intv16_t;

/** @brief Default unsigned integer type. */
typedef uint_t cl_uint_t;
/** @brief Default unsigned integer vector[2] type. */
typedef uintv2_t cl_uintv2_t;
/** @brief Default unsigned integer vector[3] type. */
typedef uintv3_t cl_uintv3_t;
/** @brief Default unsigned integer vector[4] type. */
typedef uintv4_t cl_uintv4_t;
/** @brief Default unsigned integer vector[8] type. */
typedef uintv8_t cl_uintv8_t;
/** @brief Default unsigned integer vector[16] type. */
typedef uintv16_t cl_uintv16_t;

#if USE_DOUBLE_PRECISION || defined(__DOXYGEN__)
	/** @brief Double precision EPS. */
	#define FP_EPS			FP_LITERAL(2.220446049250313e-16)
	/** @brief Maximum integer (4503599627370495) that can be represented by a double precision floating-point number. */
	#define FP_MAXINT		((uint64_t)0xFFFFFFFFFFFFFul)
#else
	/** @brief Single precision EPS. */
	#define FP_EPS			FP_LITERAL(1.1920929e-07)
	/** @brief Maximum integer (8388607) that can be represented by a single precision floating-point number. */
	#define FP_MAX_INT		(0x7FFFFF)
#endif

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
/** @brief floating-point infinity constant. */
#define FP_INF			INFINITY

#if USE_NATIVE_MATH || defined(__DOXYGEN__)
	/** @brief Native sine function. */
	#define cl_sin(x)					native_sin(x)
	/** @brief Native cosine function. */
	#define cl_cos(x)					native_cos(x)
	/** @brief Native division function. */
	#define cl_fdiv(a,b)				native_divide((cl_fp_t)(a), (cl_fp_t)(b))
	/** @brief Native reciprocal. */
	#define cl_reciprocal(x)			native_recip((cl_fp_t)(x))
	/** @brief Native natural logarithm function. */
	#define cl_log(x)					native_log(x)
	/** @brief Native square root function. */
	#define cl_sqrt(x)					native_sqrt((cl_fp_t)(x))
	/** @brief Native reciprocal square root function. */
	#define cl_rsqrt(x)					native_rsqrt((cl_fp_t)(x))
	/** @brief Native power function, where x >= 0. */
	#define cl_pow(x, y)				native_powr(x, y)
	/** @brief Standard exponential function, where x >= 0. */
	#define cl_exp(x)					native_exp(x)
	/** @brief Native Simultaneous sine and cosine computation. */
	#define cl_sincos(fi, ps, pc)		(*(ps)=sincos(fi, (pc)))
	/* Native sin and cos functions seem to be highly inacurate on many platforms. */
	/*#define cl_sincos(fi, ps, pc)		{*(ps) = native_sin(fi); *(pc) = native_cos(fi);} */
#elif USE_HALF_MATH
	/** @brief Native half precision sine function. */
	#define cl_sin(x)					half_sin(x)
	/** @brief Native half precision cosine function. */
	#define cl_cos(x)					half_cos(x)
	/** @brief Native half precision division function. */
	#define cl_fdiv(a,b)				half_divide((cl_fp_t)(a), (cl_fp_t)(b))
	/** @brief Native half precision reciprocal. */
	#define cl_reciprocal(x)			half_recip((cl_fp_t)(x))
	/** @brief Native half precision natural logarithm function. */
	#define cl_log(x)					half_log(x)
	/** @brief Native half precision square root function. */
	#define cl_sqrt(x)					half_sqrt((cl_fp_t)(x))
	/** @brief Native half precision reciprocal square root function. */
	#define cl_rsqrt(x)					half_rsqrt((cl_fp_t)(x))
	/** @brief Native half precision power function, where x >= 0. */
	#define cl_pow(x, y)				half_powr(x, y)
	/** @brief Half precision exponential function, where x >= 0. */
	#define cl_exp(x)					half_exp(x)
	/** @brief Half precision simultaneous sine and cosine computation. */
	#define cl_sincos(fi, ps, pc)		(*(ps)=sincos(fi, (pc)))
	/* Half sin and cos functions are far too inacurate. */
	/* #define cl_sincos(fi, ps, pc)	{*(ps) = half_sin(fi); *(pc) = half_cos(fi);} */
#else
	/** @brief Standard sine function. */
	#define cl_sin(x)					sin(x)
	/** @brief Standard cosine function. */
	#define cl_cos(x)					cos(x)
	/** @brief Standard single precision division function. */
	#define cl_fdiv(a,b)				((cl_fp_t)(a)/(cl_fp_t)(b))
	/** @brief Native single precision reciprocal. */
	#define cl_reciprocal(x)			cl_fdiv((cl_fp_t)(1), (cl_fp_t)(x))
	/** @brief Standard single precision reciprocal square root function. */
	#define cl_rsqrt(x)					rsqrt((cl_fp_t)(x))
	/** @brief Standard single precision square root function. */
	#define cl_sqrt(x)					sqrt((cl_fp_t)(x))
	/** @brief Standard natural logarithm function. */
	#define cl_log(x)					log(x)
	/** @brief Standard power function, where x >= 0. */
	#define cl_pow(x, y)				powr(x, y)
	/** @brief Standard exponential function, where x >= 0. */
	#define cl_exp(x)					exp(x)
	/** @brief Simultaneous sine and cosine computation. */
	#define cl_sincos(fi, ps, pc)		(*(ps)=sincos(fi, (pc)))
#endif

/** @brief Copy the sign of a floating-point number. */
#define cl_fcopysign(to, from)		(copysign(to, from))
/** @brief Evaluate to integer sign of a floating-point number. */
#define cl_fsign(x)					(((x) >= FP_0) ? 1 : -1)
/** @brief Evaluate to integer sign of an integer number. */
#define cl_sign(x)					(((x) >= 0) ? 1 : -1)
/** @brief Absolute value of a floating-point number. */
#define cl_fabs(x)					fabs((cl_fp_t)(x))
/** @brief Minimum of two floating-point numbers. */
#define cl_fmin(x,y)				fmin((cl_fp_t)(x), (cl_fp_t)(y))
/** @brief Maximum of two floating-point numbers. */
#define cl_fmax(x,y)				fmax((cl_fp_t)(x), (cl_fp_t)(y))
/** @brief Minimum of two integer numbers. */
#define cl_min(x,y)					min((cl_int_t)(x), (cl_int_t)(y))
/** @brief Maximum of two floating-point numbers. */
#define cl_max(x,y)					max((cl_int_t)(x), (cl_int_t)(y))
/** @brief Clip integer value to the specified range. */
#define cl_clip(x, low, high)		cl_min(cl_max((x), (low)), (high))
/** @brief Clip floating-point value to the specified range. */
#define cl_fclip(x, low, high)		cl_fmin(cl_fmax((x), (low)), (high))
/** @brief Compute the cube root of x */
#define cl_cbrt(x)					cbrt((cl_fp_t)(x))
/** @brief Type cast a floating-point value to integer value. */
#define cl_int(x)					convert_int(x)
/** @brief Round a floating-point value towards minus infinity. */
#define cl_round(x)					roundf(x)
/** @brief Round a floating-point value towards the closest integer. */
#define cl_floor(x)					floor(x)
/** @brief Checks if a floating-point number is finite. */
#define cl_isfinite(x)				isinf(x)

/**
 * @} // end @addtogroup cl_types_constants_math
 */
/*########### End basic data types, constants and math functions  ############*/


/*###################### Start random number generator #######################*/
/**
 * @addtogroup cl_rng OpneCL random number generators
 * @{
 */

/**
 * @brief Generates a random number from open interval (0.0, 1.0) given a state
 *	of the generator in terms of two unsigned integers. The two states/seeds 
 *	should be initialized independently for each thread. See Python module
 *	clrng for generating pairs of unsigned integer seeds.
 * @param[in, out] rngStateX Pointer to 64-bit unsigned number/state
 *					 updated on return.
 * @param[in, out] rngStateA Pointer to 32-bit unsigned number/state.
 * @return A single random number from (0.0, 1.0).
 * 
 * @details George Marsaglia's Random Number Generator. A single precision 
 *	floating point number has a 23-bit mantisa, hence only integers from 
 *	[0, 2^23 - 1] can be represented without loss of information.
 */
 inline float mlrngoof(uint64_t *rngStateX, uint32_t rngStateA){
	*rngStateX = (*rngStateX & (uint64_t)0xFFFFFFFFUL)*
		(rngStateA) + (*rngStateX >> 32);
	
	/* Generate a random number from open interval (0.0, 1.0). */
	return cl_fdiv( ((float)(1 | ((uint32_t)(*rngStateX) & 0x7FFFFF) )), 
		(float)0x800000);
};

/**
 * @brief Generates a random number from open/closed interval (0.0, 1.0] given 
 *	a state of the generator in terms of two unsigned integers. The two 
 *	states/seeds  should be initialized independently for each thread. 
 *	See Python module clrng for generating pairs of unsigned integer seeds.
 * @param[in, out] rngStateX Pointer to 64-bit unsigned number/state
 *					 updated on return.
 * @param[in, out] rngStateA Pointer to 32-bit unsigned number/state.
 * @return A single random number from (0.0, 1.0).
 * 
 * @details George Marsaglia's Random Number Generator. A single precision 
 *	floating point number has a 23-bit mantisa, hence only integers from 
 *	[0, 2^23 - 1] can be represented without loss of information.
 */
 inline float mlrngocf(uint64_t *rngStateX, uint32_t rngStateA){
	*rngStateX = (*rngStateX & (uint64_t)0xFFFFFFFFUL)*
		(rngStateA) + (*rngStateX >> 32);
	
	/* Generate a random number from open interval (0.0, 1.0). */
	return cl_fdiv( ((float)(1 | ((uint32_t)(*rngStateX) & 0x7FFFFF) )), 
		(float)0x7FFFFF);
};

/**
 * @brief Generates a random number from closed/open interval [0.0, 1.0) given 
 *	a state of the generator in terms of two unsigned integers. The two 
 *	states/seeds  should be initialized independently for each thread. 
 *	See Python module clrng for generating pairs of unsigned integer seeds.
 * @param[in, out] rngStateX Pointer to 64-bit unsigned number/state
 *					 updated on return.
 * @param[in, out] rngStateA Pointer to 32-bit unsigned number/state.
 * @return A single random number from (0.0, 1.0).
 * 
 * @details George Marsaglia's Random Number Generator. A single precision 
 *	floating point number has a 23-bit mantisa, hence only integers from 
 *	[0, 2^23 - 1] can be represented without loss of information.
 */
 inline float mlrngcof(uint64_t *rngStateX, uint32_t rngStateA){
	*rngStateX = (*rngStateX & (uint64_t)0xFFFFFFFFUL)*
		(rngStateA) + (*rngStateX >> 32);
	
	/* Generate a random number from open interval (0.0, 1.0). */
	return cl_fdiv( ((float)(((uint32_t)(*rngStateX) & 0x7FFFFF) )), 
		(float)0x800000);
};

/**
 * @brief Generates a random number from closed/closed interval [0.0, 1.0] given 
 *	a state of the generator in terms of two unsigned integers. The two 
 *	states/seeds  should be initialized independently for each thread. 
 *	See Python module clrng for generating pairs of unsigned integer seeds.
 * @param[in, out] rngStateX Pointer to 64-bit unsigned number/state
 *					 updated on return.
 * @param[in, out] rngStateA Pointer to 32-bit unsigned number/state.
 * @return A single random number from (0.0, 1.0).
 * 
 * @details George Marsaglia's Random Number Generator. A single precision 
 *	floating point number has a 23-bit mantisa, hence only integers from 
 *	[0, 2^23 - 1] can be represented without loss of information.
 */
 inline float mlrngccf(uint64_t *rngStateX, uint32_t rngStateA){
	*rngStateX = (*rngStateX & (uint64_t)0xFFFFFFFFUL)*
		(rngStateA) + (*rngStateX >> 32);
	
	/* Generate a random number from open interval (0.0, 1.0). */
	return cl_fdiv( ((float)(((uint32_t)(*rngStateX) & 0x7FFFFF) )), 
		(float)0x7FFFFF);
};

#if defined(ML_DOUBLE_PRECISION)
	/**
	 * @brief Generates a random number from open interval (0.0, 1.0) given a state 
	 *	of the generator in terms of two unsigned integers. The two states/seeds 
	 *	should be initialized independently for each thread. See Python module
	 *	clrng for generating pairs of unsigned integer seeds.
	 * @param[in, out] rngStateX Pointer to 64-bit unsigned number/state
	 *					updated on return.
	 * @param[in, out] rngStateA Pointer to 32-bit unsigned number/state.
	 * @return A double precision random number from (0.0, 1.0).
	 * 
	 * @details George Marsaglia's Random Number Generator. A double precision 
	 *	floating point number has a 52-bit mantisa, hence only integers from 
	 *	[0, 2^52 - 1] can be represented without loss of information.
	 */
	 inline double mlrngood(uint64_t *rngStateX, uint32_t rngStateA){
		*rngStateX = (*rngStateX & (uint64_t)0xFFFFFFFFUL)*
			(rngStateA) + (*rngStateX >> 32);
		
		/* Generate a random number from open interval (0.0, 1.0). */
		return cl_fdiv( ((double)(1 | (*rngStateX & (uint64_t)0xFFFFFFFFFFFFFUL))), 
			(double)0x10000000000000UL);
	};

	/**
	 * @brief Generates a random number from open/closed interval (0.0, 1.0] given
	 *	 a state of the generator in terms of two unsigned integers. The two  
	 *	states/seeds should be initialized independently for each thread.
	 *	See Python module clrng for generating pairs of unsigned integer seeds.
	 * @param[in, out] rngStateX Pointer to 64-bit unsigned number/state
	 *					updated on return.
	 * @param[in, out] rngStateA Pointer to 32-bit unsigned number/state.
	 * @return A double precision random number from (0.0, 1.0).
	 * 
	 * @details George Marsaglia's Random Number Generator. A double precision 
	 *	floating point number has a 52-bit mantisa, hence only integers from 
	 *	[0, 2^52 - 1] can be represented without loss of information.
	 */
	 inline double mlrngocd(uint64_t *rngStateX, uint32_t rngStateA){
		*rngStateX = (*rngStateX & (uint64_t)0xFFFFFFFFUL)*
			(rngStateA) + (*rngStateX >> 32);
		
		/* Generate a random number from open interval (0.0, 1.0). */
		return cl_fdiv( ((double)(1 | (*rngStateX & (uint64_t)0xFFFFFFFFFFFFFUL))), 
			(double)0xFFFFFFFFFFFFFUL);
	};

	/**
	 * @brief Generates a random number from closed/open interval [0.0, 1.0) given
	 *	 a state of the generator in terms of two unsigned integers. The two  
	 *	states/seeds should be initialized independently for each thread.
	 *	See Python module clrng for generating pairs of unsigned integer seeds.
	 * @param[in, out] rngStateX Pointer to 64-bit unsigned number/state
	 *					updated on return.
	 * @param[in, out] rngStateA Pointer to 32-bit unsigned number/state.
	 * @return A double precision random number from (0.0, 1.0).
	 * 
	 * @details George Marsaglia's Random Number Generator. A double precision 
	 *	floating point number has a 52-bit mantisa, hence only integers from 
	 *	[0, 2^52 - 1] can be represented without loss of information.
	 */
	 inline double mlrngcod(uint64_t *rngStateX, uint32_t rngStateA){
		*rngStateX = (*rngStateX & (uint64_t)0xFFFFFFFFUL)*
			(rngStateA) + (*rngStateX >> 32);
		
		/* Generate a random number from open interval (0.0, 1.0). */
		return cl_fdiv( ((double)((*rngStateX & (uint64_t)0xFFFFFFFFFFFFFUL))), 
			(double)0x10000000000000UL);
	};

	/**
	 * @brief Generates a random number from closed/closed interval [0.0, 1.0] given
	 *	 a state of the generator in terms of two unsigned integers. The two  
	 *	states/seeds should be initialized independently for each thread.
	 *	See Python module clrng for generating pairs of unsigned integer seeds.
	 * @param[in, out] rngStateX Pointer to 64-bit unsigned number/state
	 *					updated on return.
	 * @param[in, out] rngStateA Pointer to 32-bit unsigned number/state.
	 * @return A double precision random number from (0.0, 1.0).
	 * 
	 * @details George Marsaglia's Random Number Generator. A double precision 
	 *	floating point number has a 52-bit mantisa, hence only integers from 
	 *	[0, 2^52 - 1] can be represented without loss of information.
	 */
	 inline double mlrngccd(uint64_t *rngStateX, uint32_t rngStateA){
		*rngStateX = (*rngStateX & (uint64_t)0xFFFFFFFFUL)*
			(rngStateA) + (*rngStateX >> 32);
		
		/* Generate a random number from open interval (0.0, 1.0). */
		return cl_fdiv( ((double)((*rngStateX & (uint64_t)0xFFFFFFFFFFFFFUL))), 
			(double)0xFFFFFFFFFFFFFUL);
	};
#endif /* ML_DOUBLE_PRECISION */

/**
 * @} // @addtogroup cl_rng OpneCL random number generators
 */
/*####################### End random number generator ########################*/

#endif /* CLBASE_DEFINES_ */
