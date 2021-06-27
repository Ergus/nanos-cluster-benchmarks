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

#include "benchmarks_ompss.h"

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

#else // __WITH_EXTRAE

#define BLAS_EVENT 0
#define register_blas_events()

#endif // __WITH_EXTRAE

void dgemm_(const char *transa, const char *transb, int *l, int *n, int *m, double *alpha,
	const void *a, int *lda, void *b, int *ldb, double *beta, void *c, int *ldc);

void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha,
	double *a, int *lda, double *b, int *ldb);

void dsyrk_(char *uplo, char *trans, int *n, int *k, double *alpha, double *a, int *lda,
	double *beta, double *c, int *ldc);


static inline void oss_potrf(int ts, double A[ts][ts])
{
	inst_event(BLAS_EVENT, BLAS_POTRF);
	assert(A != NULL);

	static int INFO;
	static const char L = 'L';

	dpotrf_(&L, &ts, (double *)A, &ts, &INFO);
	inst_event(BLAS_EVENT, BLAS_NONE);
}

static inline void oss_trsm(int ts, double A[ts][ts], double B[ts][ts])
{
	inst_event(BLAS_EVENT, BLAS_TRSM);
	assert(A != NULL);
	assert(B != NULL);

    char LO = 'L', TR = 'T', NU = 'N', RI = 'R';
    double DONE = 1.0;
    dtrsm_(&RI, &LO, &TR, &NU, &ts, &ts, &DONE,
	       (double *)A, &ts,
	       (double *)B, &ts);

	inst_event(BLAS_EVENT, BLAS_NONE);
}

static inline void oss_gemm(
	int ts,
	double A[ts][ts],
	double B[ts][ts],
	double C[ts][ts]
) {
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
}

static inline void oss_syrk(int ts, double A[ts][ts], double B[ts][ts])
{
	inst_event(BLAS_EVENT, BLAS_SYRK);
	assert(A != NULL);
	assert(B != NULL);

    static char LO = 'L', NT = 'N';
    static double DONE = 1.0, DMONE = -1.0;
    dsyrk_(&LO, &NT, &ts, &ts, &DMONE,
	       (double *)A, &ts, &DONE,
	       (double *)B, &ts);
	inst_event(BLAS_EVENT, BLAS_NONE);
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
