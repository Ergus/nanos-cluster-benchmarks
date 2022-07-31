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

#ifndef MATMUL_OMPSS2_H
#define MATMUL_OMPSS2_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <libgen.h>

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

#else // ISMATVEC => ISMATMUL

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
#endif // ISMATVEC


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

				const int nodeid
					= (WITHNODE ? i / rowsPerNode : nanos6_cluster_no_hint);

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

	// checkpoint and restart functions.
	double *alloc_restart(const size_t rows, size_t cols, int process, int id)
	{
		myassert(process != 0);
		const size_t size = cols * rows;

		double *ret = (double *) nanos6_dmalloc(
			size * sizeof(double), nanos6_equpart_distribution, 0, NULL
		);
		myassert(ret != NULL);

		int rc = nanos6_serialize(ret, size * sizeof(double), process, id, false);
		myassert(rc == 0);

		return ret;
	}

	void checkpoint_matrix(
		double *mat, const size_t rows, size_t cols, int process, int id
	) {
		myassert(process != 0);
		myassert(id != 0);
		const size_t size = cols * rows;

		int rc = nanos6_serialize(mat, size * sizeof(double), process, id, true);
		myassert(rc == 0);
	}

//=================== Conditionals after here ==================================

#if TASKTYPE == 0 // strong flat

	void matvec_tasks_ompss2(const double *A, const double *B, double *C,
	                         size_t ts, size_t dim, size_t colsBC, size_t it
	) {
		if (it == 0) {
			printf("# %s strong flat tasks ompss2\n",
			       (ISMATVEC ? "matvec" : "matmul"));
		}

		myassert(ts <= dim);
		modcheck(dim, ts);

		const size_t numNodes = nanos6_get_num_cluster_nodes();
		const size_t rowsPerNode = dim / numNodes;

		for (size_t i = 0; i < dim; i += ts) {
			const int nodeid = (WITHNODE ? i / rowsPerNode : nanos6_cluster_no_hint);

			#pragma oss task in(A[i * dim; ts * dim])	\
				in(B[0; dim * colsBC])					\
				out(C[i * colsBC; ts * colsBC])			\
				node(nodeid) label("strongmatvec")
			matmul_base(&A[i * dim], B, &C[i * colsBC], ts, dim, colsBC);
		}
	}

#elif TASKTYPE == 1 // strong nested

	void matvec_tasks_ompss2(const double *A, const double *B, double *C,
	                         size_t ts, size_t dim, size_t colsBC, size_t it
	) {
		if (it == 0) {
			printf("# %s strong nested task ompss2 withnode=%d\n",
			       (ISMATVEC ? "matvec" : "matmul"), WITHNODE);
		}

		myassert(ts <= dim);
		modcheck(dim, ts);

		const size_t numNodes = nanos6_get_num_cluster_nodes();

		const size_t rowsPerNode = dim / numNodes;
		myassert(ts <= rowsPerNode);
		modcheck(rowsPerNode, ts);

		for (size_t i = 0; i < dim; i += rowsPerNode) {

			const int nodeid = (WITHNODE ? i / rowsPerNode : nanos6_cluster_no_hint);

			#pragma oss task in(A[i * dim; rowsPerNode * dim])			\
				in(B[0; dim * colsBC])									\
				out(C[i * colsBC; rowsPerNode * colsBC])				\
				node(nodeid) label("weakmatvec")
			{
				for (size_t j = i; j < i + rowsPerNode; j += ts) {
					#pragma oss task in(A[j * dim; ts * dim])			\
						in(B[0; dim * colsBC])							\
						out(C[j * colsBC; ts * colsBC])					\
						node(nanos6_cluster_no_offload) label("strongmatvec")
					matmul_base(&A[j * dim], B, &C[j * colsBC], ts, dim, colsBC);
				}
			}
		}
	}

#elif TASKTYPE == 2 // taskfor

	void matvec_tasks_ompss2(const double *A, const double *B, double *C,
	                         size_t nchunks, size_t dim, size_t colsBC, size_t it
	) {
		if (it == 0) {
			printf("# %s weak taskfor ompss2 withnode=%d\n",
			       (ISMATVEC ? "matvec" : "matmul"), WITHNODE);
		}

		const size_t numNodes = nanos6_get_num_cluster_nodes();
		myassert(nchunks <= dim);
		modcheck(dim, nchunks);

		const size_t ts = dim / numNodes / nchunks;
		myassert(ts * numNodes * nchunks == dim);

		const size_t rowsPerNode = dim / numNodes;
		myassert(ts <= rowsPerNode);
		modcheck(rowsPerNode, ts);

		for (size_t i = 0; i < dim; i += rowsPerNode) {

			const int nodeid = (WITHNODE ? i / rowsPerNode : nanos6_cluster_no_hint);

			#pragma oss task weakin(A[i * dim; rowsPerNode * dim])		\
				weakin(B[0; dim * colsBC])								\
				weakout(C[i * colsBC; rowsPerNode * colsBC])			\
				node(nodeid)											\
				firstprivate(ts) label("weakmatvec_task")
			{
				#pragma oss task for in(A[i * dim; rowsPerNode * dim])	\
					in(B[0; dim * colsBC])								\
					out(C[i * colsBC; rowsPerNode * colsBC])			\
					firstprivate(ts) chunksize(ts)						\
					node(nanos6_cluster_no_offload)						\
					label("taskfor_matvec")
				for (size_t j = i; j < i + rowsPerNode; j += ts) {
					matmul_base(&A[j * dim], B, &C[j * colsBC], ts, dim, colsBC);
				}
			}
		}
	}

#elif TASKTYPE == 3 // weak tasks

	void matvec_tasks_ompss2(const double *A, const double *B, double *C,
	                         size_t ts, size_t dim, size_t colsBC, size_t it
	) {
		if (it == 0) {
			printf("# %s weak task FETCHTASK=%d withnode=%d\n",
			       (ISMATVEC ? "matvec" : "matmul"), FETCHTASK, WITHNODE);
		}

		myassert(ts <= dim);
		modcheck(dim, ts);

		const size_t numNodes = nanos6_get_num_cluster_nodes();

		const size_t rowsPerNode = dim / numNodes;
		myassert(ts <= rowsPerNode);
		modcheck(rowsPerNode, ts);

		for (size_t i = 0; i < dim; i += rowsPerNode) {

			const int nodeid = (WITHNODE ? i / rowsPerNode : nanos6_cluster_no_hint);

			#pragma oss task weakin(A[i * dim; rowsPerNode * dim])		\
				weakin(B[0; dim * colsBC])								\
				weakout(C[i * colsBC; rowsPerNode * colsBC])			\
				node(nodeid) label("weakmatvec")
			{
				#if FETCHTASK == 1
				#pragma oss task in(A[i * dim; rowsPerNode * dim])		\
					in(B[0; dim * colsBC])								\
					out(C[i * colsBC; rowsPerNode * colsBC])			\
					node(nanos6_cluster_no_offload) label("fetchtask")
				{
				}
				#endif

				for (size_t j = i; j < i + rowsPerNode; j += ts) {
					#pragma oss task in(A[j * dim; ts * dim])			\
						in(B[0; dim * colsBC])							\
						out(C[j * colsBC; ts * colsBC])					\
						node(nanos6_cluster_no_offload) label("strongmatvec")
					{
						matmul_base(&A[j * dim], B, &C[j * colsBC], ts, dim, colsBC);
					}
				}
			}
		}
	}

#endif // TASKTYPE

#ifdef __cplusplus
}
#endif

#endif // MATMUL_OMPSS2_H
