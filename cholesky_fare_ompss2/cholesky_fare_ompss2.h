/*
 * Copyright (C) 2021  Jimmy Aguilar Mena
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CHOLESKY_OMP_MPI_H
#define CHOLESKY_OMP_MPI_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>

#include "benchmarks_ompss.h"


//! This is to initialize blocks[i][j]
static inline void fill_block(const size_t ts, double block[ts][ts],
                              const size_t i, const size_t j, const size_t dim
) {
	myassert(block);

	const size_t seed = (i > j) ? i * ts + j : j * ts + i;
	struct drand48_data status;       	// using re-entrant version for rng
	srand48_r(seed, &status);
	double rnd1, rnd2;

	for (size_t k = 0; k < ts; ++k) {
		if (i == j) {
			drand48_r(&status, &rnd1);
			drand48_r(&status, &rnd2);
			block[k][k] = rnd1 * rnd2 + dim;

			for (size_t l = k + 1; l < ts; ++l) {
				drand48_r(&status, &rnd1);
				drand48_r(&status, &rnd2);

				const double val = rnd1 * rnd2;
				block[k][l] = val;
				block[l][k] = val;
			}

		} else {
			for (size_t l = 0; l < ts; ++l) {
				drand48_r(&status, &rnd1);
				drand48_r(&status, &rnd2);

				if (i > j) {
					block[k][l] = rnd1 * rnd2;
				} else {
					block[l][k] = rnd1 * rnd2;
				}
			}
		}
	}
}

static inline bool compare_blocks(
	const size_t ts,
	double B1[ts][ts],
	double B2[ts][ts]
) {
	for (size_t k = 0; k < ts; ++k) {
		for (size_t l = 0; l < ts; ++l) {
			if (B1[k][l] != B2[k][l]) {
				printf("Check failed for B[%zu][%zu] = %lf expected: %lf\n",
				       k, l,
				       B1[k][l], B2[k][l]);
				return false;
			}
		}
	}
	return true;
}

static inline void print_block(const size_t ts, double A[ts][ts])
{

	for (size_t i = 0; i < ts; ++i) {
		for (size_t j = 0; j < ts; ++j) {
			printf("%5.2f ", (float)A[i][j]);
		}
		printf("\n");
	}
	fflush(stdout);
}

static inline void oss_potrf(int ts, double A[ts][ts], int k, int y, int x, int prvanim)
{
	inst_blas_kernel(prvanim, BLAS_POTRF, k, y, x);
	assert(A != NULL);

	static int INFO;
	static const char L = 'L';

	dpotrf_(&L, &ts, (double *)A, &ts, &INFO);
	inst_blas_kernel(prvanim, BLAS_NONE, k,y,x);
}

static inline void oss_trsm(int ts, double A[ts][ts], double B[ts][ts], int k, int y, int x, int prvanim)
{
	inst_blas_kernel(prvanim, BLAS_TRSM, k,y,x);
	assert(A != NULL);
	assert(B != NULL);


    char LO = 'L', TR = 'T', NU = 'N', RI = 'R';
    double DONE = 1.0;
    dtrsm_(&RI, &LO, &TR, &NU, &ts, &ts, &DONE,
	       (double *)A, &ts,
	       (double *)B, &ts);

	inst_blas_kernel(prvanim, BLAS_NONE, k,y,x);
}

static inline void oss_gemm(
	int ts,
	double A[ts][ts],
	double B[ts][ts],
	double C[ts][ts],
	int k, int y, int x, int prvanim
) {
	inst_blas_kernel(prvanim, BLAS_GEMM, k,y,x);
	assert(A != NULL);
	assert(B != NULL);
	assert(C != NULL);

    const char TR = 'T', NT = 'N';
    double DONE = 1.0, DMONE = -1.0;
    dgemm_(&NT, &TR, &ts, &ts, &ts, &DMONE,
	       (double *)A, &ts,
	       (double *)B, &ts, &DONE,
	       (double *)C, &ts);

	inst_blas_kernel(prvanim, BLAS_NONE, k,y,x);
}

static inline void oss_syrk(int ts, double A[ts][ts], double B[ts][ts], int k, int y, int x, int prvanim)
{
	inst_blas_kernel(prvanim, BLAS_SYRK, k,y,x);
	assert(A != NULL);
	assert(B != NULL);

    static char LO = 'L', NT = 'N';
    static double DONE = 1.0, DMONE = -1.0;
    dsyrk_(&LO, &NT, &ts, &ts, &DMONE,
	       (double *)A, &ts, &DONE,
	       (double *)B, &ts);
	inst_blas_kernel(prvanim, BLAS_NONE, k,y,x);
}

#ifdef __cplusplus
}
#endif


#endif /* CHOLESKY_OMP_MPI_H */
