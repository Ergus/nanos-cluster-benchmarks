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

#include "ArgParserC/argparser.h"

#include "matvec.h"

#ifdef ISSTRONG

void matvec_tasks(const double *A, const double *B, double *C,
                  size_t ts, size_t dim, size_t colsBC, size_t it)
{
	myassert(ts <= dim);
	modcheck(dim, ts);

	for (size_t i = 0; i < dim; i += ts) {
		#pragma oss task in(A[i * dim; dim * ts]) \
			in(B[0; dim * colsBC]) \
			out(C[i * colsBC; ts * colsBC]) label("weakmatvec")
		matmul_base(&A[i * dim], C, &C[i], ts, dim, colsBC);
	}
}

#else // ISSTRONG

#if FETCHTASK == 0
#define THECOND 0
#elif FETCHTASK == 1
#define THECOND 1
#elif FETCHTASK == 2
#define THECOND it < 1
#else
#error FETCHTASK value not valid.
#endif

void matvec_tasks(const double *A, const double *B, double *C,
                  size_t ts, size_t dim, size_t colsBC, size_t it)
{
	const size_t numNodes = nanos6_get_num_cluster_nodes();
	myassert(ts <= dim);
	modcheck(dim, ts);

	const size_t rowsPerNode = dim / numNodes;
	myassert(ts <= rowsPerNode);
	modcheck(rowsPerNode, ts);

	for (size_t i = 0; i < dim; i += rowsPerNode) {

		#pragma oss task weakin(A[i * dim; rowsPerNode * dim])		\
			weakin(B[0; dim * colsBC])								\
			weakout(C[i; rowsPerNode * colsBC]) label("weakmatvec")
		{
			if (THECOND) {
				#pragma oss task in(A[i * dim; rowsPerNode * dim])		\
					in(B[0; dim * colsBC])								\
					out(C[i; rowsPerNode * colsBC]) label("fetchtask")
				{
				}
			}

			for (size_t j = i; j < i + rowsPerNode; j += ts) {
				#pragma oss task in(A[j * dim; ts * dim])		\
					in(B[0; dim * colsBC])						\
					out(C[j; ts * colsBC]) label("strongmatvec")
				matmul_base(&A[j * dim], B, &C[j], ts, dim, colsBC);
			}
		}
	}
}

#endif // ISSTRONG

#if ISMATVEC
#define PREFIX "matvec"
#else
#define PREFIX "matmul"
#endif

int main(int argc, char* argv[])
{
	init_args(argc, argv);

	const int ROWS = create_cl_int ("Rows");
	const int TS = create_cl_int ("Tasksize");
	const int ITS = create_optional_cl_int ("Iterations", 1);
	const int PRINT = create_optional_cl_int ("Print", 0);

	printf("# Initializing data\n");
	timer ttimer = create_timer("Total time");

	const size_t colsBC = (ISMATVEC ? 1 : ROWS);

	double *A = alloc_init(ROWS, ROWS, TS);      // this initialized by blocks ts x rows
	double *B = alloc_init(ROWS, colsBC, TS);    // this splits the array in ts
	double *C = alloc_init(ROWS, colsBC, ROWS);  // This one initializes all the arrays
	#pragma oss taskwait

	printf("# Starting algorithm\n");
	timer atimer = create_timer("Algorithm time");

	for (int i = 0; i < ITS; ++i)
		matvec_tasks(A, B, C, TS, ROWS, colsBC, i);
	#pragma oss taskwait

	stop_timer(&atimer);

	printf("# Finished algorithm...\n");
	stop_timer(&ttimer);

	if (PRINT) {
		printmatrix_task(A, ROWS, ROWS, PREFIX);
		printmatrix_task(B, ROWS, colsBC, PREFIX);
		printmatrix_task(C, ROWS, colsBC, PREFIX);

		if (PRINT > 1) {
			const bool valid = validate(A, B, C, ROWS, ROWS);
			printf("# Verification: %s\n", (valid ? "Success" : "Failed"));
		}
	}

	free_matrix(A, ROWS * ROWS);
	free_matrix(B, ROWS * colsBC);
	free_matrix(C, ROWS * colsBC);

	create_reportable_int("worldsize", nanos6_get_num_cluster_nodes());
	create_reportable_int("cpu_count", count_sched_cpus());

	const double performance = ITS * ROWS * ROWS * 2000.0 / getNS_timer(&atimer);
	create_reportable_double("performance", performance);

	report_args();
	free_args();

	return 0;
}
