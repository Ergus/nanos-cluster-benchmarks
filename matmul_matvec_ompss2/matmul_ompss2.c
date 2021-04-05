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

#if ISSTRONG == 1

void matvec_tasks(const double *A, const double *B, double *C,
                  size_t ts, size_t dim, size_t colsBC, size_t it
) {
	if (it == 0)
		printf("# matvec_strong_flat \n");
	myassert(ts <= dim);
	modcheck(dim, ts);

	const size_t numNodes = nanos6_get_num_cluster_nodes();
	const size_t rowsPerNode = dim / numNodes;

	for (size_t i = 0; i < dim; i += ts) {
		int nodeid = i / rowsPerNode;

		#pragma oss task in(A[i * dim; ts * dim]) \
			in(B[0; colsBC * dim]) \
			out(C[i * colsBC; ts * colsBC]) \
			node(nodeid) label("strongmatvec")
		matmul_base(&A[i * dim], B, &C[i], ts, dim, colsBC);
	}
}

#elif ISSTRONG == 2

void matvec_tasks(const double *A, const double *B, double *C,
                  size_t ts, size_t dim, size_t colsBC, size_t it
) {
	if (it == 0)
		printf("# matvec_strong_interleave\n");
	myassert(ts <= dim);
	modcheck(dim, ts);

	const size_t numNodes = nanos6_get_num_cluster_nodes();

	for (size_t i = 0; i < dim; i += ts) {
		int nodeid = i % numNodes;

		#pragma oss task in(A[i * dim; ts * dim]) \
			in(B[0; colsBC * dim]) \
			out(C[i * colsBC; ts * colsBC]) \
			node(nodeid) label("strongmatvec")
		matmul_base(&A[i * dim], B, &C[i], ts, dim, colsBC);
	}
}


#elif ISSTRONG == 3

void matvec_tasks(const double *A, const double *B, double *C,
                  size_t ts, size_t dim, size_t colsBC, size_t it
) {
	if (it == 0)
		printf("# matvec_strong_nested\n");
	myassert(ts <= dim);
	modcheck(dim, ts);

	const size_t numNodes = nanos6_get_num_cluster_nodes();

	const size_t rowsPerNode = dim / numNodes;
	myassert(ts <= rowsPerNode);
	modcheck(rowsPerNode, ts);

	for (size_t i = 0; i < dim; i += rowsPerNode) {

		int nodeid = i / rowsPerNode;

		#pragma oss task in(A[i * dim; rowsPerNode * dim])				\
			in(B[0; dim * colsBC])										\
			out(C[i; rowsPerNode * colsBC]) node(nodeid) label("weakmatvec")
		{
			for (size_t j = i; j < i + rowsPerNode; j += ts) {
				#pragma oss task in(A[j * dim; ts * dim])				\
					in(B[0; dim * colsBC])								\
					out(C[j; ts * colsBC])								\
					node(nanos6_cluster_no_offload) label("strongmatvec")
				matmul_base(&A[j * dim], B, &C[j], ts, dim, colsBC);
			}
		}
	}
}

#else // ISSTRONG == 0

#if FETCHTASK == 0
#define THECOND 0
#elif FETCHTASK == 1
#define THECOND 1
#elif FETCHTASK == 2
#define THECOND it < 1
#else  // FETCHTASK
#error FETCHTASK value not valid.
#endif // FETCHTASK

void matvec_tasks(const double *A, const double *B, double *C,
                  size_t ts, size_t dim, size_t colsBC, size_t it
) {
	if (it == 0)
		printf("# matvec_weak_node FETCHTASK %d\n", FETCHTASK);
	myassert(ts <= dim);
	modcheck(dim, ts);

	const size_t numNodes = nanos6_get_num_cluster_nodes();

	const size_t rowsPerNode = dim / numNodes;
	myassert(ts <= rowsPerNode);
	modcheck(rowsPerNode, ts);

	for (size_t i = 0; i < dim; i += rowsPerNode) {

#if WITH_NODE == 1
		int nodeid = i / rowsPerNode;
#else
		int nodeid = nanos6_cluster_no_hint;
#endif

		#pragma oss task weakin(A[i * dim; rowsPerNode * dim])			\
			weakin(B[0; dim * colsBC])									\
			weakout(C[i; rowsPerNode * colsBC])							\
			node(nodeid) label("weakmatvec")
		{
			if (THECOND) {
				#pragma oss task in(A[i * dim; rowsPerNode * dim])		\
					in(B[0; dim * colsBC])								\
					out(C[i; rowsPerNode * colsBC])						\
					node(nanos6_cluster_no_offload) label("fetchtask")
				{
				}
			}

			for (size_t j = i; j < i + rowsPerNode; j += ts) {
				#pragma oss task in(A[j * dim; ts * dim])		\
					in(B[0; dim * colsBC])						\
					out(C[j; ts * colsBC])						\
					node(nanos6_cluster_no_offload) label("strongmatvec")
				matmul_base(&A[j * dim], B, &C[j], ts, dim, colsBC);
			}
		}
	}
}
#endif // ISSTRONG

int main(int argc, char* argv[])
{
	init_args(argc, argv);

	const char *PREFIX = basename(argv[0]);
	const int ROWS = create_cl_int ("Rows");
	const int TS = create_cl_int ("Tasksize");
	const int ITS = create_optional_cl_int ("Iterations", 1);
	const int PRINT = create_optional_cl_int ("Print", 0);

	printf("# Initializing data\n");
	timer ttimer = create_timer("Total_time");

	const size_t colsBC = (ISMATVEC == 1 ? 1 : ROWS);

	double *A = alloc_init(ROWS, ROWS, TS, true);      // this initialized by blocks ts x rows
	double *B = alloc_init(ROWS, colsBC, TS, true);    // this splits the array in ts
	double *C = alloc_init(ROWS, colsBC, TS, false);
	#pragma oss taskwait

	printf("# Starting algorithm\n");
	timer atimer = create_timer("Algorithm_time");

	for (int i = 0; i < ITS; ++i) {
		matvec_tasks(A, B, C, TS, ROWS, colsBC, i);
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

	report_args();
	free_args();

	return 0;
}
