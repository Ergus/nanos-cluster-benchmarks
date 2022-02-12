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

#include "jacobi_ompss2.h"

void init_AB_taskfor(double *A, double *B, const size_t dim, size_t ts)
{
	const size_t numNodes = nanos6_get_num_cluster_nodes();
	myassert(dim >= ts);
	myassert(dim / ts >= numNodes);
	modcheck(dim, ts);
	modcheck(dim / numNodes, ts);

	const size_t rowsPerNode = dim / numNodes;

	for (size_t i = 0; i < dim; i += rowsPerNode) { // loop nodes

		const int nodeid = i / rowsPerNode;

		#pragma oss task weakout(A[i * dim; rowsPerNode * dim])			\
			weakout(B[i; rowsPerNode])									\
			firstprivate(ts) node(nodeid) wait label("init_AB_weak")
		{
			#pragma oss task for out(A[i * dim; rowsPerNode * dim])		\
				out(B[i; rowsPerNode])									\
				firstprivate(ts) chunksize(ts)							\
				node(nanos6_cluster_no_offload)							\
				label("init_AB_taskfor")
			for (size_t j = i; j < i + rowsPerNode; ++j) { // loop tasks

				struct drand48_data drand_buf;
				srand48_r(j, &drand_buf);
				double x;

				double cum = 0.0, sum = 0.0;
				for (size_t l = 0; l < dim; ++l) {
					drand48_r(&drand_buf, &x);
					A[j * dim + l] = x;
					cum += fabs(x);
					sum += x;
				}
				// Diagonal element condition.
				const double valjj = A[j * dim + j];
				if (signbit(valjj)) {
					A[j * dim + j] = valjj - cum;
					B[j] = sum - cum;
				} else {
					A[j * dim + j] = valjj + cum;
					B[j] = sum + cum;
				}
			}
		}
	}
}

void init_x_taskfor(double *x, const size_t dim, size_t ts, double val)
{
	const size_t numNodes = nanos6_get_num_cluster_nodes();
	myassert(dim >= ts);
	modcheck(dim, ts);

	const size_t rowsPerNode = dim / numNodes;

	for (size_t i = 0; i < dim; i += rowsPerNode) { // loop nodes

		int nodeid = i / rowsPerNode;

		#pragma oss task weakout(x[i; rowsPerNode])			\
			node(nodeid)									\
			firstprivate(ts) wait label("init_vec_weak")
		{
			#pragma oss task for out(x[i; rowsPerNode])					\
				firstprivate(ts) chunksize(ts)							\
				node(nanos6_cluster_no_offload)							\
				label("init_vec_taskfor")
			for (size_t j = i; j < i + rowsPerNode; ++j) { // loop tasks
				x[j] = val;
			}
		}
	}
}

void jacobi_modify_taskfor(double *A, double *b, size_t dim, size_t ts)
{
	const size_t numNodes = nanos6_get_num_cluster_nodes();
	myassert(dim >= ts);
	modcheck(dim, ts);

	const size_t rowsPerNode = dim / numNodes;

	for (size_t i = 0; i < dim; i += rowsPerNode) { // loop nodes

		int nodeid = i / rowsPerNode;

		#pragma oss task weakinout(A[i * dim; rowsPerNode * dim])	\
			weakinout(b[i; rowsPerNode])							\
			node(nodeid)											\
			firstprivate(ts) wait label("jacobi_modify_weak")
		{
			#pragma oss task for inout(A[i * dim; rowsPerNode * dim])	\
				inout(b[i; rowsPerNode])								\
				firstprivate(ts) chunksize(ts)							\
				node(nanos6_cluster_no_offload)							\
				label("jacobi_modify_taskfor")
			for (size_t j = i; j < i + rowsPerNode; ++j) { // loop tasks
				const double iAjj = 1 / fabs(A[j * dim + j]);

				for (size_t l = 0; l < dim; ++l) {
					A[j * dim + l] = (j == l) ? 0 : -1. * A[j * dim + l] * iAjj;
				}
				b[j] *= iAjj;
			}
		}
	}
}


void jacobi_taskfor_ompss2(const double *A, const double *B,
                           double *xin, double *xout,
                           size_t nchunks, size_t dim, size_t it
) {
	if (it == 0){
		printf("# jacobi taskfor");
	}

	const size_t numNodes = nanos6_get_num_cluster_nodes();
	myassert(nchunks <= dim);
	modcheck(dim, nchunks);

	const size_t ts = dim / numNodes / nchunks;
	myassert(ts * numNodes * nchunks == dim);

	const size_t rowsPerNode = dim / numNodes;
	myassert(ts <= rowsPerNode);
	modcheck(rowsPerNode, ts);

	if (numNodes > 0) {
		#pragma oss task inout(xin[0; dim]) node(0) label("trick2")
		{
		}
	}

	for (int i = 0; i < dim; i += rowsPerNode) {

		const int nodeid = i / rowsPerNode;

		#pragma oss task weakin(A[i * dim; rowsPerNode * dim])		\
			weakin(xin[0; dim])											\
			weakin(B[i; rowsPerNode])									\
			weakout(xout[i; rowsPerNode])								\
			node(nodeid)												\
			firstprivate(ts) wait label("weakjacobi_task")
		{
			#pragma oss task for in(A[i * dim; rowsPerNode * dim])	\
				in(xin[0; dim])											\
				in(B[i; rowsPerNode])									\
				out(xout[i; rowsPerNode])								\
				firstprivate(ts) chunksize(ts)							\
				node(nanos6_cluster_no_offload)							\
				label("jacobi_taskfor")
			for (size_t j = i; j < i + rowsPerNode; ++j) {
				jacobi(&A[j * dim], &B[j], xin, &xout[j], 1, dim);
			}
		}
	}
}

int main(int argc, char* argv[])
{
	init_args(argc, argv);

	const char *PREFIX = basename(argv[0]);
	const size_t ROWS = create_cl_size_t ("Rows");
	const size_t TS = create_cl_size_t ("Tasksize");
	const size_t ITS = create_optional_cl_size_t ("Iterations", 1);
	const int PRINT = create_optional_cl_int ("Print", 0);
	myassert(ITS > 0);

	printf("# Initializing data\n");
	timer ttimer = create_timer("Total_time");

	double *A = nanos6_dmalloc(ROWS * ROWS * sizeof(double),
	                           nanos6_equpart_distribution, 0, NULL);

	double *B = nanos6_dmalloc(ROWS * sizeof(double),
	                           nanos6_equpart_distribution, 0, NULL);

	double *x1 = nanos6_dmalloc(ROWS * sizeof(double),
	                            nanos6_equpart_distribution, 0, NULL);

	double *x2 = nanos6_dmalloc(ROWS * sizeof(double),
	                            nanos6_equpart_distribution, 0, NULL);

	init_AB_taskfor(A, B, ROWS, TS);
	jacobi_modify_taskfor(A, B, ROWS, TS);

	init_x_taskfor(x1, ROWS, TS, 0.0);
	init_x_taskfor(x2, ROWS, TS, 0.0);

	#pragma oss taskwait

	printf("# Starting algorithm\n");
	timer atimer = create_timer("Algorithm_time");

	double *xin = NULL;
	double *xout = NULL;

	for (int i = 0; i < ITS; ++i) {
		xin = (i % 2 == 0) ? x1 : x2;
		xout = (i % 2 == 0) ? x2 : x1;

		jacobi_taskfor_ompss2(A, B, xin, xout, TS, ROWS, i);
	}
	#pragma oss taskwait

	stop_timer(&atimer);

	printf("# Finished algorithm...\n");
	stop_timer(&ttimer);

	if (PRINT) {
		printmatrix_task(A, ROWS, ROWS, PREFIX);
		printmatrix_task(B, ROWS, 1, PREFIX);
		printmatrix_task(xout, ROWS, 1, PREFIX);
		#pragma oss taskwait
	}

	nanos6_dfree(A, ROWS * ROWS * sizeof(double));
	nanos6_dfree(B, ROWS * sizeof(double));
	nanos6_dfree(x1, ROWS * sizeof(double));
	nanos6_dfree(x2, ROWS * sizeof(double));

	create_reportable_int("worldsize", nanos6_get_num_cluster_nodes());
	create_reportable_int("cpu_count", nanos6_get_num_cpus());
	create_reportable_int("namespace_enabled", nanos6_get_namespace_is_enabled());
	create_reportable_string("nanos6_version", nanos6_get_runtime_version());

	report_args();
	free_args();

	return 0;
}
