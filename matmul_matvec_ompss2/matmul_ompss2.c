/*
 * Copyright (C) 2019  Jimmy Aguilar Mena
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <libgen.h>

#include "ArgParserC/argparser.h"

#include "matmul_ompss2.h"

#if ISMATVEC
void matmul_base(const double *A, const double *B, double * const C,
                 size_t ts, size_t dim, size_t colsBC
) {
	inst_blas_kernel(false, BLAS_DGEMV, 0, 0, 0);

	myassert(dim < (size_t) INT_MAX);

	const char TR = 'T';
	const int M = (int) dim;
	const int N = (int) ts;
	const double alpha = 1.0;
	const double beta = 0.0;
	const int incx = 1;

	dgemv_(&TR, &M, &N, &alpha, A, &M, B, &incx, &beta, C, &incx);

	// for (size_t i = 0; i < ts; ++i) {
	// 	C[i] = 0.0;

	// 	for (size_t j = 0; j < dim; ++j) {
	// 		C[i] += A[i * dim + j] * B[j];
	// 	}
	// }

	inst_blas_kernel(false, BLAS_NONE, 0, 0, 0);
}
#else
void matmul_base(const double *A, const double *B, double * const C,
                 size_t ts, size_t dim, size_t colsBC
) {
	inst_blas_kernel(false, BLAS_DGEMM, 0, 0, 0);

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

	// for (size_t i = 0; i < ts; ++i) {
	// 	for (size_t k = 0; k < colsBC; ++k)
	// 		C[i * colsBC + k] = 0.0;

	// 	for (size_t j = 0; j < dim; ++j) {
	// 		const double temp = A[i * dim + j];

	// 		for (size_t k = 0; k < colsBC; ++k) {
	// 			C[i * colsBC + k] += (temp * B[j * colsBC + k]);
	// 		}
	// 	}
	// }

	inst_blas_kernel(false, BLAS_NONE, 0, 0, 0);
}
#endif

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

		#pragma oss task in(A[i * dim; rowsPerNode * dim])				\
			in(B[0; dim * colsBC])										\
			out(C[i * colsBC; rowsPerNode * colsBC])					\
			node(nodeid) label("weakmatvec")
		{
			for (size_t j = i; j < i + rowsPerNode; j += ts) {
				#pragma oss task in(A[j * dim; ts * dim])				\
					in(B[0; dim * colsBC])								\
					out(C[j * colsBC; ts * colsBC])						\
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

		#pragma oss task weakin(A[i * dim; rowsPerNode * dim])			\
			weakin(B[0; dim * colsBC])									\
			weakout(C[i * colsBC; rowsPerNode * colsBC])				\
			node(nodeid)												\
			firstprivate(ts) label("weakmatvec_task")
		{
			#pragma oss task for in(A[i * dim; rowsPerNode * dim])		\
				in(B[0; dim * colsBC])									\
				out(C[i * colsBC; rowsPerNode * colsBC])				\
				firstprivate(ts) chunksize(ts)							\
				node(nanos6_cluster_no_offload)							\
				label("taskfor_matvec")
			for (size_t j = i; j < i + rowsPerNode; j += ts) {
				matmul_base(&A[j * dim], B, &C[j * colsBC], ts, dim, colsBC);
			}
		}
	}
}

#elif TASKTYPE == 3

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

		#pragma oss task weakin(A[i * dim; rowsPerNode * dim])			\
			weakin(B[0; dim * colsBC])									\
			weakout(C[i * colsBC; rowsPerNode * colsBC])				\
			node(nodeid) label("weakmatvec")
		{
#if FETCHTASK == 1
			#pragma oss task in(A[i * dim; rowsPerNode * dim])			\
				in(B[0; dim * colsBC])									\
				out(C[i * colsBC; rowsPerNode * colsBC])				\
				node(nanos6_cluster_no_offload) label("fetchtask")
			{
			}
#endif
			for (size_t j = i; j < i + rowsPerNode; j += ts) {
				#pragma oss task in(A[j * dim; ts * dim])				\
					in(B[0; dim * colsBC])								\
					out(C[j * colsBC; ts * colsBC])						\
					node(nanos6_cluster_no_offload) label("strongmatvec")
				{
					matmul_base(&A[j * dim], B, &C[j * colsBC], ts, dim, colsBC);
				}
			}
		}
	}
}

#endif // TASKTYPE

int main(int argc, char* argv[])
{
	init_args(argc, argv);

	const char *PREFIX = basename(argv[0]);
	const size_t ROWS = create_cl_size_t ("Rows");
	const size_t TS = create_cl_size_t ("Tasksize");
	const size_t ITS = create_optional_cl_size_t ("Iterations", 1);
	const int PRINT = create_optional_cl_int ("Print", 0);

	printf("# Initializing data\n");
	timer ttimer = create_timer("Total_time");

	const size_t colsBC = (ISMATVEC == 1 ? 1 : ROWS);

	double *A = alloc_init(ROWS, ROWS, TS, true);    // this initialized by blocks ts x rows
	double *B = alloc_init(ROWS, colsBC, TS, true);  // this splits the array in ts
	double *C = alloc_init(ROWS, colsBC, TS, false);
	#pragma oss taskwait

	printf("# Starting algorithm\n");
	timer atimer = create_timer("Algorithm_time");

	for (int i = 0; i < ITS; ++i) {
		matvec_tasks_ompss2(A, B, C, TS, ROWS, colsBC, i);
	}
	#pragma oss taskwait

	stop_timer(&atimer);

	printf("# Finished algorithm...\n");
	stop_timer(&ttimer);

	if (PRINT) {
		printmatrix_task(A, ROWS, ROWS, PREFIX);
		printmatrix_task(B, ROWS, colsBC, PREFIX);
		printmatrix_task(C, ROWS, colsBC, PREFIX);

		if (PRINT > 1) {
			const bool valid = validate(A, B, C, ROWS, colsBC);
			printf("# Verification: %s\n", (valid ? "Success" : "Failed"));
		}

		#pragma oss taskwait
	}

	free_matrix(A, ROWS * ROWS);
	free_matrix(B, ROWS * colsBC);
	free_matrix(C, ROWS * colsBC);

	create_reportable_int("worldsize", nanos6_get_num_cluster_nodes());
	create_reportable_int("cpu_count", nanos6_get_num_cpus());
	create_reportable_int("namespace_enabled", nanos6_get_namespace_is_enabled());
	create_reportable_string("nanos6_version", nanos6_get_runtime_version());

	report_args();
	free_args();

	return 0;
}
