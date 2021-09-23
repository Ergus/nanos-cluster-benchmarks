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

#include <mkl.h>
#include <mpi.h>
#include <omp.h>
#include <extrae.h>

#include "benchmarks_ompss.h"

#define PRVANIM_NONE 0
#define PRVANIM_POTRF 1
#define PRVANIM_TRSM 2
#define PRVANIM_GEMM 3
#define PRVANIM_SYRK 4
#define PRVANIM_EVENT 9200042

#if __WITH_EXTRAE

#define BLAS_EVENT 9910003

#define BLAS_EVENT_VALUES						\
	EVENT(BLAS_NONE)							\
	EVENT(BLAS_POTRF)							\
	EVENT(BLAS_TRSM)							\
	EVENT(BLAS_GEMM)							\
	EVENT(BLAS_SYRK)

enum blas_values_t {
#define EVENT(evt) evt,
	BLAS_EVENT_VALUES
#undef EVENT
	BLAS_NEVENTS
};

void register_blas_events()
{
	extrae_type_t event = BLAS_EVENT;

	unsigned nvalues = BLAS_NEVENTS;

	static extrae_value_t blas_values[BLAS_NEVENTS] = {
#define EVENT(evt) (extrae_value_t) evt,
		BLAS_EVENT_VALUES
#undef EVENT
	};

	static char *blas_names[BLAS_NEVENTS] = {
#define EVENT(evt) #evt,
		BLAS_EVENT_VALUES
#undef EVENT
	};

	inst_define_event_type(&event, "blas_event", &nvalues, blas_values, blas_names);
}

void inst_prvanim(int kernel, int k, int y, int x)
{
	nanos6_instrument_event(PRVANIM_EVENT, kernel);
	nanos6_instrument_event(PRVANIM_EVENT, k + 1000000);
	nanos6_instrument_event(PRVANIM_EVENT, y + 2000000);
	nanos6_instrument_event(PRVANIM_EVENT, x + 3000000);
}

void register_prvanim_events()
{
	// extrae_type_t event = PRVANIM_KERNEL;
	// unsigned int nvalues = 0;
	// inst_define_event_type(&event, "prvanim kernel", &nvalues, NULL, NULL);
	// event = PRVANIM_K;
	// inst_define_event_type(&event, "prvanim k", &nvalues, NULL, NULL);
	// event = PRVANIM_Y;
	// inst_define_event_type(&event, "prvanim y", &nvalues, NULL, NULL);
	// event = PRVANIM_X;
	// inst_define_event_type(&event, "prvanim x", &nvalues, NULL, NULL);
}


#else // __WITH_EXTRAE

#define BLAS_EVENT 0
#define register_blas_events()
#define register_prvanim_events()
#define inst_prvanim(kernel, k, y, x)

#endif // __WITH_EXTRAE

void dgemm_(const char *transa, const char *transb, int *l, int *n, int *m, double *alpha,
	const void *a, int *lda, void *b, int *ldb, double *beta, void *c, int *ldc);

void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha,
	double *a, int *lda, double *b, int *ldb);

void dsyrk_(char *uplo, char *trans, int *n, int *k, double *alpha, double *a, int *lda,
	double *beta, double *c, int *ldc);


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
	if (prvanim) {
		inst_prvanim(PRVANIM_POTRF, k,y,x);
	}
	inst_event(BLAS_EVENT, BLAS_POTRF);
	assert(A != NULL);

	static int INFO;
	static const char L = 'L';

	dpotrf_(&L, &ts, (double *)A, &ts, &INFO);
	inst_event(BLAS_EVENT, BLAS_NONE);
	if (prvanim) {
		inst_prvanim(PRVANIM_NONE, k,y,x);
	}
}

static inline void oss_trsm(int ts, double A[ts][ts], double B[ts][ts], int k, int y, int x, int prvanim)
{
	if (prvanim) {
		inst_prvanim(PRVANIM_TRSM, k,y,x);
	}
	inst_event(BLAS_EVENT, BLAS_TRSM);
	assert(A != NULL);
	assert(B != NULL);


    char LO = 'L', TR = 'T', NU = 'N', RI = 'R';
    double DONE = 1.0;
    dtrsm_(&RI, &LO, &TR, &NU, &ts, &ts, &DONE,
	       (double *)A, &ts,
	       (double *)B, &ts);

	inst_event(BLAS_EVENT, BLAS_NONE);
	if (prvanim) {
		inst_prvanim(PRVANIM_NONE, k,y,x);
	}
}

static inline void oss_gemm(
	int ts,
	double A[ts][ts],
	double B[ts][ts],
	double C[ts][ts],
	int k, int y, int x, int prvanim
) {
	if (prvanim) {
		inst_prvanim(PRVANIM_GEMM, k,y,x);
	}
	inst_event(BLAS_EVENT, BLAS_GEMM);
	assert(A != NULL);
	assert(B != NULL);
	assert(C != NULL);

    const char TR = 'T', NT = 'N';
    double DONE = 1.0, DMONE = -1.0;
    dgemm_(&NT, &TR, &ts, &ts, &ts, &DMONE,
	       (double *)A, &ts,
	       (double *)B, &ts, &DONE,
	       (double *)C, &ts);

	inst_event(BLAS_EVENT, BLAS_NONE);
	if (prvanim) {
		inst_prvanim(PRVANIM_NONE, k,y,x);
	}
}

static inline void oss_syrk(int ts, double A[ts][ts], double B[ts][ts], int k, int y, int x, int prvanim)
{
	if (prvanim) {
		inst_prvanim(PRVANIM_SYRK, k,y,x);
	}
	inst_event(BLAS_EVENT, BLAS_SYRK);
	assert(A != NULL);
	assert(B != NULL);

    static char LO = 'L', NT = 'N';
    static double DONE = 1.0, DMONE = -1.0;
    dsyrk_(&LO, &NT, &ts, &ts, &DMONE,
	       (double *)A, &ts, &DONE,
	       (double *)B, &ts);
	inst_event(BLAS_EVENT, BLAS_NONE);
	if (prvanim) {
		inst_prvanim(PRVANIM_NONE, k,y,x);
	}
}

static inline void wait(MPI_Request *comm_req)
{
	int comm_comp = 0;

	do {
		MPI_Test(comm_req, &comm_comp, MPI_STATUS_IGNORE);
	} while (!comm_comp);
}


#ifdef __cplusplus
}
#endif


#endif /* CHOLESKY_OMP_MPI_H */
