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

#include "benchmarks_ompss.h"
#include "matmul_ompss2.h"

int main(int argc, char* argv[])
{
	init_args(argc, argv);

	const char *PREFIX = basename(argv[0]);
	const size_t ROWS = create_cl_size_t ("Rows");
	const size_t TS = create_cl_size_t ("Tasksize");
	const size_t ITS = create_cl_size_t ("Iterations");
	const int PRINT = create_optional_cl_int ("Print", 0);
	const size_t CHECKPOINT = create_optional_cl_size_t("Checkpoint", 0);
	const size_t RESTART = create_optional_cl_size_t("Restart", 0);
	const size_t TIME_S = create_optional_cl_size_t("Epoch_secs", 0);
	const size_t TIME_NS = create_optional_cl_size_t("Epoch_nsecs", 0);

	// Measure the time to start the process and initialize the runtime
	// This requires to pass $(date "+%s %N") for the last two arguments.
	if (TIME_S != 0 && TIME_NS != 0) {
		struct timespec t1 = {.tv_sec = TIME_S, .tv_nsec = TIME_NS};

		struct timespec t2;
		clock_gettime(CLOCK_REALTIME, &t2);

		struct timespec delta = diffTime(&t1, &t2);
		create_reportable_double("srun_time", getNS(&delta));
	}

	inst_register_events();  // Register the events in the instrumentation

	timer ttimer = create_timer("Total_time");

	const size_t colsBC = (ISMATVEC == 1 ? 1 : ROWS);

	double *A = NULL, *B  = NULL, *C = NULL;

	if (RESTART != 0) {
		printf("# Recovering checkpoint: %zu\n", RESTART);
		timer rtimer = create_timer("Restart_time");
		A = alloc_restart(ROWS, ROWS, RESTART, 1);    // this initialized by blocks ts x rows
		B = alloc_restart(ROWS, colsBC, RESTART, 2);  // this splits the array in ts
		C = alloc_restart(ROWS, colsBC, RESTART, 3);
		#pragma oss taskwait
		stop_timer(&rtimer);
	} else {
		printf("# Initializing data\n");
		A = alloc_init(ROWS, ROWS, TS, true);    // this initialized by blocks ts x rows
		B = alloc_init(ROWS, colsBC, TS, true);  // this splits the array in ts
		C = alloc_init(ROWS, colsBC, TS, false);
		#pragma oss taskwait
	}

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

	if (CHECKPOINT != 0) {
		timer ctimer = create_timer("Checkpoint_time");
		checkpoint_matrix(A, ROWS, ROWS, CHECKPOINT, 1);
		checkpoint_matrix(B, ROWS, colsBC, CHECKPOINT, 2);
		checkpoint_matrix(C, ROWS, colsBC, CHECKPOINT, 3);
		#pragma oss taskwait
		stop_timer(&ctimer);
	}

	timer ftimer = create_timer("Free_time");
	free_matrix(A, ROWS * ROWS);
	free_matrix(B, ROWS * colsBC);
	free_matrix(C, ROWS * colsBC);
	stop_timer(&ftimer);

	create_reportable_int("worldsize", nanos6_get_num_cluster_nodes());
	create_reportable_int("cpu_count", nanos6_get_num_cpus());
	create_reportable_int("namespace_enabled", nanos6_get_namespace_is_enabled());
	create_reportable_string("nanos6_version", nanos6_get_runtime_version());

	report_args();
	free_args();

	return 0;
}
