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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "gmp.h"

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

int main(int argc, char *argv[]){
	FILE* file = NULL;
	mpz_t n1, n2, a;
	int i=0;

	mpz_init(n1);
	mpz_init(n2);
	mpz_init(a);
	mpz_set_str(a, "4294967118", 0);

	int n = 15000;

	if (argc >= 2 && argv[1]){
		n = atoi( argv[1] );
	};

	if (argc >= 3 && argv[2]){
		file = fopen(argv[2], "w");
	}else{
		file = fopen("mc_primes.txt", "w");
	}

	if (!file){
		printf("Failed to open output file!");
		exit(1);
	}

    printf("Generating %d primes\n", n);

	while(i < n){
		mpz_mul_2exp(n2, a, 32);
		mpz_sub_ui(n2, n2, 1lu);
		if(isprime(n2)){
			mpz_sub_ui(n1, n2, 1lu);
			mpz_fdiv_q_2exp(n1, n1, 1lu);
			if(isprime(n1)){
				//We have found our safeprime, calculate a and print to file
				mpz_out_str (file, 10, a);
				fprintf(file, " ");
				mpz_out_str (file, 10, n2);
				fprintf(file, " ");
				mpz_out_str (file, 10, n1);
				fprintf(file, "\n");
				printf("\rProcessing %d/%d", ++i, n);
			 }
		}
		mpz_sub_ui(a, a, 1lu);
	}

    printf("\nDone\n");

	fclose(file);
	exit (0);
};
