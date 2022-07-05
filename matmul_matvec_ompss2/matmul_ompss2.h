/*
 * Copyright (C) 2019  Jimmy Aguilar Mena
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

#ifndef MATVEC_H
#define MATVEC_H

#ifdef __cplusplus
extern "C" {
#endif

#include "benchmarks_ompss.h"

#if ISMATVEC
	void matmul_base(const double *A, const double *B, double * const C,
	                 size_t ts, size_t dim, size_t colsBC
	) {
		#if BLAS == 0
		inst_event(USER_EVENT, USER_MATVEC);
		for (size_t i = 0; i < ts; ++i) {
			C[i] = 0.0;

			for (size_t j = 0; j < dim; ++j) {
				C[i] += A[i * dim + j] * B[j];
			}
		}
		inst_event(USER_EVENT, USER_NONE);

		#elif BLAS == 1

		inst_blas_kernel(false, BLAS_GEMV, 0, 0, 0);

		myassert(dim < (size_t) INT_MAX);

		const char TR = 'T';
		const int M = (int) dim;
		const int N = (int) ts;
		const double alpha = 1.0;
		const double beta = 0.0;
		const int incx = 1;

		dgemv_(&TR, &M, &N, &alpha, A, &M, B, &incx, &beta, C, &incx);
		inst_blas_kernel(false, BLAS_NONE, 0, 0, 0);
		#else // BLAS
		#error No valid BLAS value
		#endif // BLAS
	}
#else
	void matmul_base(const double *A, const double *B, double * const C,
	                 size_t ts, size_t dim, size_t colsBC
	) {
		#if BLAS == 0
		inst_event(USER_EVENT, USER_MATMUL);
		for (size_t i = 0; i < ts; ++i) {
			for (size_t k = 0; k < colsBC; ++k)
				C[i * colsBC + k] = 0.0;

			for (size_t j = 0; j < dim; ++j) {
				const double temp = A[i * dim + j];

				for (size_t k = 0; k < colsBC; ++k) {
					C[i * colsBC + k] += (temp * B[j * colsBC + k]);
				}
			}
		}
		inst_event(USER_EVENT, USER_NONE);

		#elif BLAS == 1

		inst_blas_kernel(false, BLAS_GEMM, 0, 0, 0);

		const char TA = 'N';
		const char TB = 'N';
		const int M = (int) dim;
		const int N = (int) ts;
		const int K = (int) colsBC;
		const double ALPHA = 1.0;
		const int LDA = M;
		const int LDB = K;
		const double BETA = 0.0;
		const int LDC = M;

		dgemm_(&TA, &TB, &M, &N, &K, &ALPHA,
		       B, &LDB,
		       A, &LDA, &BETA,
		       C, &LDC);

		inst_blas_kernel(false, BLAS_NONE, 0, 0, 0);

		#else // BLAS
		#error No valid BLAS value
		#endif // BLAS
	}
#endif


	void free_matrix(double *mat, size_t size)
	{
		nanos6_dfree(mat, size * sizeof(double));
	}

	bool validate(const double *A, const double *B, double *C,
	              size_t dim, size_t colsBC
	) {
		bool success = false;

		#pragma oss task in(A[0; dim * dim])			\
			in(B[0; dim * colsBC])						\
			in(C[0; dim * colsBC])						\
			inout(success) label("validate")
		{
			success = true;

			double *expected = (double *) malloc(dim * colsBC * sizeof(double));
			matmul_base(A, B, expected, dim, dim, colsBC);

			for (size_t i = 0; i < dim; ++i) {
				for (size_t j = 0; j < colsBC; ++j) {

					if (abs(C[i * colsBC + j] - expected[i * colsBC + j]) > 1e-12) {
						success = false;
						break;
					}
				}
			}

			free(expected);
		}
		#pragma oss taskwait
		return success;
	}


	double *alloc_init(const size_t rows, size_t cols, size_t ts, bool init)
	{
		const size_t numNodes = nanos6_get_num_cluster_nodes();
		myassert(rows >= ts);              // at least 1 portion per task
		myassert(rows / ts >= numNodes);   // at least 1 task / node.
		modcheck(rows, ts);

		const size_t size = cols * rows;

		double *ret =
			(double *) nanos6_dmalloc(size * sizeof(double),
			                          nanos6_equpart_distribution, 0, NULL);
		myassert(ret != NULL);

		if (init) { // Initialize to random??

			const size_t rowsPerNode = rows / numNodes;

			for (size_t i = 0; i < rows; i += rowsPerNode) { // loop nodes

#if WITHNODE == 1
				const int nodeid = i / rowsPerNode;
#else
				const int nodeid = nanos6_cluster_no_hint;
#endif
				#pragma oss task weakout(ret[i * cols; rowsPerNode * cols]) \
					node(nodeid) label("initalize_weak")
				{
					for (size_t j = i; j < i + rowsPerNode; j += ts) { // loop tasks

						#pragma oss task out(ret[j * cols; ts * cols])		\
							node(nanos6_cluster_no_offload) label("initalize_slice")
						{
							struct drand48_data drand_buf;
							srand48_r(j, &drand_buf);
							double x;

							const size_t elems = ts * cols;

							for (size_t k = 0; k < elems; ++k) {
								drand48_r(&drand_buf, &x);
								ret[j * cols + k] = x;
							}
						}
					}
				}
			}

		}
		return ret;
	}


#ifdef __cplusplus
}
#endif

#endif
