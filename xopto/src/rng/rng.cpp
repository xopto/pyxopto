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

#if defined(_MSC_VER)
    // #include "stdafx.h"
#endif
#include <stdint.h>
#include <math.h>

#if defined(EXPORT_FIND_PRIMES)
	#include "gmp.h"
#endif

// see: http://stackoverflow.com/questions/2164827/explicitly-exporting-shared-library-functions-in-linux
// details: http://gcc.gnu.org/wiki/Visibility

#if defined(_MSC_VER)
    //  Microsoft 
    // #define EXPORT __declspec(dllexport)
    // #define IMPORT __declspec(dllimport)
    #ifdef RNG_EXPORTS
		#define RNG_API __declspec(dllexport) 
	#else
		#define RNG_API __declspec(dllimport) 
	#endif
#elif defined(__GNUC__)
    //  GCC
    // #define EXPORT __attribute__((visibility("default")))
    // #define IMPORT
    #ifdef RNG_EXPORTS
		#define RNG_API __attribute__((visibility("default")))
	#else
		#define RNG_API
	#endif
#else
    //  do nothing and hope for the best?
    #define EXPORT
    #define IMPORT
    #pragma warning Unknown dynamic link import/export semantics.
#endif

#ifdef __cplusplus
extern "C" {
#endif

RNG_API int init_RNG(uint64_t *x, uint32_t *a, uint32_t *fora,
	const uint32_t n_rng, uint64_t xinit){

	uint32_t begin = fora[0];

	// Here we set up a loop, using the first multiplier in the file to generate x's and c's
	// There are some restictions to these two numbers:
	// 0<=c<a and 0<=x<b, where a is the multiplier and b is the base (2^32)
	// also [x,c]=[0,0] and [b-1,a-1] are not allowed.

	//Make sure xinit is a valid seed (using the above mentioned restrictions)
	if ( (xinit == 0ULL) |
			(((uint32_t)(xinit >> 32)) >= (begin - 1)) |
			(((uint32_t)xinit) >= 0xffffffffUL)){
		//xinit (probably) not a valid seed! (we have excluded a few unlikely exceptions)
		return 1;
	}

	for (uint32_t i = 0; i < n_rng; i++){
		a[i] = fora[i + 1];
		x[i] = 0;
		while ((x[i] == 0) |
				(((uint32_t)(x[i] >> 32)) >= (fora[i + 1] - 1)) |
				(((uint32_t)x[i]) >= 0xffffffffUL)){
		
			//generate a random number
			xinit = (xinit & 0xffffffffULL)*(begin) + (xinit >> 32);

			//calculate c and store in the upper 32 bits of x[i]
			x[i] = (uint32_t)floor( (((double)((uint32_t)xinit))/(double)0x100000000ULL)
									*fora[i + 1]); //Make sure 0<=c<a
			x[i] = x[i] << 32;

			//generate a random number and store in the lower 32 bits of x[i] (as the initial x of the generator)
			xinit = (xinit & 0xffffffffULL)*(begin) + (xinit >> 32);//x will be 0<=x<b, where b is the base 2^32
			x[i] += (uint32_t)xinit;
		}
	}
	return 0;
};

RNG_API int rand_cc(float *data, uint32_t n, uint64_t* x, uint32_t a)
{
	// Generate random numbers from [0,1].
	for (uint32_t i = 0; i < n; ++i){
		*x = (*x & 0xffffffffULL)*(a) + (*x >> 32);
		data[i] = ((float)((unsigned int)(*x) & 0x7FFFFF))/(float)0x7FFFFF;
	}
	return 0;
};

RNG_API int rand_co(float *data, uint32_t n, uint64_t* x, uint32_t a)
{
	// Generate random numbers from [0,1).
	for (uint32_t i = 0; i < n; ++i){
		*x = (*x & 0xffffffffULL)*(a) + (*x >> 32);
		data[i] = ((float)((unsigned int)(*x) & 0x7FFFFF))/(float)0x800000;
	}
	return 0;
};

RNG_API int rand_oc(float *data, uint32_t n, uint64_t* x, uint32_t a){

	// Generate random numbers from (0,1].
	for (uint32_t i = 0; i < n; ++i){
		*x = (*x & 0xffffffffULL)*(a) + (*x >> 32);
		data[i] = ((float)(1UL | ((unsigned int)(*x) & 0x7FFFFF))) / (float)0x7FFFFF;
	}
	return 0;
};

RNG_API int rand_oo(float *data, uint32_t n, uint64_t* x, uint32_t a)
{
	// Generate random numbers from [0,1).
	for (uint32_t i = 0; i < n; ++i){
		*x = (*x & 0xffffffffULL)*(a) + (*x >> 32);
		data[i] = ((float)(1UL | ((unsigned int)(*x) & 0x7FFFFF))) / (float)0x800000;
	}
	return 0;
};

#if defined(EXPORT_FIND_PRIMES)
	static int isprime(mpz_t n){
		int i;
		int kind;
		int test_list[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41};
		for(i = 0; i < 13; i++){
			kind = mpz_probab_prime_p(n, test_list[i]);
			if(kind == 0)
				return 0;
		}
		return 1;
	};

	RNG_API int find_primes(
			unsigned int n, uint32_t *vec_a, uint64_t *vec_n1, uint64_t *vec_n2){
		mpz_t n1, n2, a;
		unsigned int i=0;

		mpz_init(n1);
		mpz_init(n2);
		mpz_init(a);
		mpz_set_str(a, "4294967118", 0);

		printf("Generating %d primes\n", n);

		while(i < n){
			mpz_mul_2exp(n2, a, 32);
			mpz_sub_ui(n2, n2, 1lu);
			if(isprime(n2)){
				mpz_sub_ui(n1, n2, 1lu);
				mpz_fdiv_q_2exp(n1, n1, 1lu);
				if(isprime(n1)){
					//We have found our safeprime, calculate a and print to file
					mpz_export (&vec_a[i], NULL, 1, sizeof(uint32_t), 0, 0, a);
					mpz_export (&vec_n2[i], NULL, 1, sizeof(uint64_t), 0, 0, n2);
					mpz_export (&vec_n1[i], NULL, 1, sizeof(uint64_t), 0, 0, n1);
					printf("\rProcessing %d/%d", ++i, n);
				}
			}
			mpz_sub_ui(a, a, 1lu);
		}

		printf("\nDone\n");

		return i;
	}
#endif


#if defined(_MSC_VER)
	BOOL APIENTRY DllMain(
			HMODULE hModule,
			DWORD  ul_reason_for_call,
			LPVOID lpReserved)
	{
		switch (ul_reason_for_call)
		{
		case DLL_PROCESS_ATTACH:
		case DLL_THREAD_ATTACH:
		case DLL_THREAD_DETACH:
		case DLL_PROCESS_DETACH:
			break;
		}
		return TRUE;
	}

#endif

#ifdef __cplusplus
};
#endif
