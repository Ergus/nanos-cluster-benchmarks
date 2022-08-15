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
	const int CHECK = create_optional_cl_int("Check", 0);
	modcheck(ROWS, TS);

	inst_register_events();  // Register the events in the instrumentation

	printf("# Initializing data\n");;

	timer ttimer = create_timer("Total_time");

	// Allocate memory
	const size_t nt = ROWS / TS;

	int (*block_rank)[nt] = nanos6_lmalloc(nt * nt * sizeof(int));
	get_block_rank_task(nt, block_rank);

	double (*A)[nt][TS][TS] = nanos6_dmalloc(ROWS * ROWS * sizeof(double),
	                                         nanos6_equpart_distribution, 0, NULL);
	assert(A != NULL);

	double (*Ans)[nt][TS][TS] = NULL;
	if (CHECK) {
		Ans = nanos6_dmalloc(ROWS * ROWS * sizeof(double),
		                     nanos6_equpart_distribution, 0, NULL);
		assert(Ans != NULL);
	}

	#pragma oss taskwait

	// ===========================================
	// WARMUP: The first iteration is slow, likely because of the
	// overhead of thread creation, which seems to impact MPI. Extra
	// threads are needed for the weak tasks, which stay alive for
	// autowait. Note: we want them to stay alive also so that there
	// is more opportunity to connect later tasks in the namespace.
	cholesky_init_task(nt, TS, A, Ans, block_rank);
	#pragma oss taskwait

	// ===========================================
	printf("# Starting warmup\n");

	timer atimer_warmup = create_timer("Warmup_time");
	cholesky_ompss2(nt, TS, A, block_rank, /* prvanim */ false);
	#pragma oss taskwait

	stop_timer(&atimer_warmup);
	printf("# Finished warmup\n");

	// ===========================================
	// ACTUAL CALCULATION
	cholesky_init_task(nt, TS, A, Ans, block_rank);
	#pragma oss taskwait

	// ===========================================
	printf("# Starting algorithm\n");

	timer atimer = create_timer("Algorithm_time");
	cholesky_ompss2(nt, TS, A, block_rank, /* prvanim */ true);
	#pragma oss taskwait

	stop_timer(&atimer);
	// ===========================================

	printf("# Finished algorithm\n");

	if (CHECK) {
		timer stimer = create_timer("Single_time");
		cholesky_single(nt, TS, Ans);
		#pragma oss taskwait
		stop_timer(&stimer);
	}

	stop_timer(&ttimer);

	// Verification
	if (Ans != NULL) {
		#pragma oss task in(A[0;nt][0;nt][0;TS][0;TS]) \
			in(Ans[0;nt][0;nt][0;TS][0;TS])					\
			node(nanos6_cluster_no_offload) label("check")
		{
			bool match = true;
			for (size_t i = 0; i < nt && match; i++) {
				for (size_t j = 0; j < nt && match; j++) {
					match = compare_blocks(TS, A[i][j], Ans[i][j]);

					if (!match) {
						fprintf(stderr, "# Check failed in block A[%zu][%zu]\n",
						        i, j);
					}
				}
			}
		}
	}

	#pragma oss taskwait

	create_reportable_int("Iterations", 1);
	create_reportable_int("worldsize", nanos6_get_num_cluster_nodes());
	create_reportable_int("cpu_count", nanos6_get_num_cpus());
	create_reportable_int("namespace_enabled", nanos6_get_namespace_is_enabled());
	create_reportable_string("nanos6_version", nanos6_get_runtime_version());

	report_args();

	// Release memory
	if (Ans != NULL) {
		nanos6_dfree(Ans, ROWS * ROWS * sizeof(double));
	}

	assert(A);
	nanos6_dfree(A, ROWS * ROWS * sizeof(double));

	assert(block_rank);
	nanos6_lfree(block_rank, nt * nt * sizeof(int));

	free_args();

	return 0;
}
