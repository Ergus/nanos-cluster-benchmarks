/*
 * Copyright (C) 2022  Jimmy Aguilar Mena
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

// Call this program like: srun ... ./executable DIM TS 1 2 4 8 16 32 64

// Where all the numbers after the TS are the dimensions to resize before. next
// iteration. Only matmul/size will be executed every time.

int main(int argc, char* argv[])
{
	init_args(argc, argv);

	const char *PREFIX = basename(argv[0]);
	const size_t ROWS = create_cl_size_t ("Rows");
	const size_t TS = create_cl_size_t ("Tasksize");

	const size_t initialSize = nanos6_get_num_cluster_nodes();

	inst_register_events();  // Register the events in the instrumentation

	timer ttimer = create_timer("Total_time");

	const size_t colsBC = (ISMATVEC == 1 ? 1 : ROWS);

	printf("# Initializing data\n");
	double *A = alloc_init(ROWS, ROWS, TS, true);    // this initialized by blocks ts x rows
	double *B = alloc_init(ROWS, colsBC, TS, true);  // this splits the array in ts
	double *C = alloc_init(ROWS, colsBC, TS, false);
	#pragma oss taskwait

	printf("# Starting algorithm\n");

	// This is from the ArgParser API to access the remaining extra arguments
	// The pointer is to the original argv array, so don't change it.
	char **rest = NULL;
	int nrest = get_rest_args(&rest);

	assert(rest != NULL || nrest == 0);

	for (int n = 0;  n < nrest; ++n) {
		const int oldsize = nanos6_get_num_cluster_nodes();
		const int newsize = strtol(rest[n], NULL, 0);
		assert(newsize >= initialSize);

		const int delta = newsize - oldsize;

		if (delta == 0) {
			continue;
		}

		printf("# Resize start: %d %d %d\n", oldsize, newsize, delta);

		struct timespec startRes, endRes, end1, end2;

		getTime(&startRes);
		nanos6_cluster_resize(delta);
		getTime(&endRes);

		matvec_tasks_ompss2(A, B, C, TS, ROWS, colsBC, 0);
		#pragma oss taskwait noflush
		getTime(&end1);

		matvec_tasks_ompss2(A, B, C, TS, ROWS, colsBC, 0);
		#pragma oss taskwait noflush
		getTime(&end2);

		const struct timespec deltaRes = diffTime(&startRes, &endRes);
		const struct timespec delta1 = diffTime(&endRes, &end1);
		const struct timespec delta2 = diffTime(&end1, &end2);

		printf("## Resize:%d->%d:%d %lg %lg %lg\n",
		       oldsize, newsize, delta,
		       getNS(&deltaRes), getNS(&delta1), getNS(&delta2));
	}

	printf("# Finished algorithm and resizes...\n");
	stop_timer(&ttimer);

	timer ftimer = create_timer("Free_time");
	free_matrix(A, ROWS * ROWS);
	free_matrix(B, ROWS * colsBC);
	free_matrix(C, ROWS * colsBC);
	stop_timer(&ftimer);

	create_reportable_int("cpu_count", nanos6_get_num_cpus());
	create_reportable_int("namespace_enabled", nanos6_get_namespace_is_enabled());
	create_reportable_string("nanos6_version", nanos6_get_runtime_version());

	report_args();
	free_args();

	return 0;
}
