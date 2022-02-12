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

#ifndef BENCHMARKS_OMPSS_H
#define BENCHMARKS_OMPSS_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <libgen.h>  // basename
#include <assert.h>
#include <syscall.h>
#include <errno.h>

#include <numaif.h>

#include <nanos6.h>
#include <nanos6/debug.h>

#include "cmacros/macros.h"
#include "ArgParserC/argparser.h"

	void __print_task(const double * const mat,
	                  const size_t rows, const size_t cols,
	                  const char prefix[64], const char name[64]
	) {
		#pragma oss task in(mat[0; rows * cols]) label("matrix_print")
		{
			__print(mat, rows, cols, prefix, name);
		}
	}

#define printmatrix_task(mat, rows, cols, prefix)	\
	__print_task(mat, rows, cols, prefix, #mat)


	int get_numa_from_address(void *ptr)
	{
		int status;
		int numa_node = -1;
		const int ret = get_mempolicy(&numa_node, NULL, 0, ptr, MPOL_F_NODE | MPOL_F_ADDR);
		/* const int ret = move_pages(0 , 1, &ptr, NULL, &numa_node, 0); */

		if (ret != 0) {
			int errnum = errno;
			perror("Numa error: ");
			return -1;
		}

		return numa_node;
	}


	// Declare some blas routines.
	#include <mkl.h>
	#include <limits.h>

	void dcopy_(const int *n, const double *dx, const int *incx, double *dy, const int *incy);

	void dgemv_ (const char *trans, const int *m, const int *n,
	             const double *alpha, const double *A, const int *lda,
	             const double *x, const int *incx,
	             const double *beta, double *y, const int *incy);

	void dgemm_(const char *transa, const char *transb, int *l, int *n, int *m,
	            double *alpha, const void *a, int *lda, void *b, int *ldb,
	            double *beta, void *c, int *ldc);

	void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n,
	            double *alpha, double *a, int *lda, double *b, int *ldb);

	void dsyrk_(char *uplo, char *trans, int *n, int *k,
	            double *alpha, double *a, int *lda,
	            double *beta, double *c, int *ldc);


#if __WITH_EXTRAE // #####################

#include <extrae.h>
#include "extrae_user_events.h"

	typedef extrae_type_t inst_type_t;
	typedef extrae_value_t inst_value_t;

#define inst_define_event_type(type,name,nvalues,values,descriptions)	\
	Extrae_define_event_type(type,name,nvalues,values,descriptions)
#define inst_event(evt, val)					\
	Extrae_event(evt, val)

#define BLAS_EVENT 9910003
#define PRVANIM_EVENT 9200042

#define BLAS_EVENT_VALUES						\
	EVENT(BLAS_NONE)							\
	EVENT(BLAS_POTRF)							\
	EVENT(BLAS_TRSM)							\
	EVENT(BLAS_GEMM)							\
	EVENT(BLAS_GEMV)							\
	EVENT(BLAS_COPY)							\
	EVENT(BLAS_SYRK)

	enum blas_values_t {
		#define EVENT(evt) evt,
		BLAS_EVENT_VALUES
		#undef EVENT
		BLAS_NEVENTS
	};

	void register_blas_events()
	{
		inst_type_t event = BLAS_EVENT;

		unsigned nvalues = BLAS_NEVENTS;

		static extrae_value_t blas_values[BLAS_NEVENTS] = {
			#define EVENT(evt) (inst_value_t) evt,
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

	void inst_blas_kernel(bool emmit, int kernel, int k, int y, int x)
	{
		inst_event(BLAS_EVENT, kernel);

		#ifdef __WITH_PRVANIM
		if (!emmit) {
			return;
		}
		nanos6_instrument_event(PRVANIM_EVENT, kernel);
		nanos6_instrument_event(PRVANIM_EVENT, k + 1000000);
		nanos6_instrument_event(PRVANIM_EVENT, y + 2000000);
		nanos6_instrument_event(PRVANIM_EVENT, x + 3000000);
		#endif // __WITH_PRVANIM
	}

#else // __WITH_EXTRAE // #####################

	typedef size_t inst_type_t;
	typedef size_t inst_value_t;

#define inst_define_event_type(type,name,nvalues,values,descriptions)
#define inst_event(evt, val)

#define BLAS_EVENT 0
#define register_blas_events()
#define inst_blas_kernel(emmit, kernel, k, y, x)

#endif // __WITH_EXTRAE // #####################

#ifdef __cplusplus
}
#endif

#endif // BENCHMARKS_OMPSS_H
