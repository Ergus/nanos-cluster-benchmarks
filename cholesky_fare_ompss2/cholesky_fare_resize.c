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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.	 If not, see <http://www.gnu.org/licenses/>.
 */

#include "cholesky_fare.h"


int main(int argc, char *argv[])
{
	init_args(argc, argv);

	const char *PREFIX = basename(argv[0]);
	const size_t ROWS = create_cl_size_t("Rows");
	const size_t TS = create_cl_size_t("Tasksize");
	modcheck(ROWS, TS);

	const size_t initialSize = nanos6_get_num_cluster_nodes();

	inst_register_events();  // Register the events in the instrumentation

	printf("# Initializing data\n");

	timer ttimer = create_timer("Total_time");

	// Allocate memory
	const size_t nt = ROWS / TS;

	int (*block_rank)[nt] = nanos6_lmalloc(nt * nt * sizeof(int));

	double (*A)[nt][TS][TS] = nanos6_dmalloc(ROWS * ROWS * sizeof(double),
	                                         nanos6_equpart_distribution, 0, NULL);
	assert(A != NULL);

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

		struct timespec startRes, endRes, endInit, endCho;

		getTime(&startRes);
		nanos6_cluster_resize(delta);
		getTime(&endRes);

		printf("# Resize done: %d %d %d\n", oldsize, newsize, delta);

		get_block_rank_task(nt, block_rank);

		#pragma oss task in(block_rank[0:nt][0:nt])					\
			out(A[0;nt][0;nt][0;TS][0;TS])							\
			node(nanos6_cluster_no_offload) label("init_wrapper")
		cholesky_init_task(nt, TS, A, NULL, block_rank);

		#pragma oss taskwait
		getTime(&endInit);

		// ===========================================
		// printf("# Starting algorithm\n");

		// cholesky_ompss2(nt, TS, A, block_rank, /* prvanim */ true);
		// #pragma oss taskwait
		// getTime(&endCho);

		// ===========================================

		const struct timespec deltaRes = diffTime(&startRes, &endRes);
		const struct timespec deltaInit = diffTime(&endRes, &endInit);
		// const struct timespec deltaCho = diffTime(&endInit, &endCho);

		printf("## Resize:%d->%d:%d %lg %lg\n",
		       oldsize, newsize, delta, getNS(&deltaRes), getNS(&deltaInit));

	}
	printf("# Finished algorithm\n");
	stop_timer(&ttimer);

	#pragma oss taskwait

	create_reportable_int("Iterations", 1);
	create_reportable_int("worldsize", nanos6_get_num_cluster_nodes());
	create_reportable_int("cpu_count", nanos6_get_num_cpus());
	create_reportable_int("namespace_enabled", nanos6_get_namespace_is_enabled());
	create_reportable_string("nanos6_version", nanos6_get_runtime_version());

	report_args();

	assert(A);
	nanos6_dfree(A, ROWS * ROWS * sizeof(double));

	assert(block_rank);
	nanos6_lfree(block_rank, nt * nt * sizeof(int));

	free_args();

	return 0;
}
