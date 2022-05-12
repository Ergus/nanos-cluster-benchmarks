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

#include "cholesky_fare_ompss2.h"

void print_matrix_task(const size_t nt, const size_t ts,
                       double A[nt][nt][ts][ts], int nodeid
) {
	#pragma oss task in(A[0:nt][0:nt][0;ts][0;ts]) \
		node(nodeid) label("print_matrix")
	{
		for (size_t i = 0; i < nt; i++) {
			for (size_t k = 0; k < ts; ++k) {
				for (size_t j = 0; j < nt; ++j) {
					for (size_t l = 0; l < ts; ++l) {
						printf("%5.2f ", (float)A[i][j][k][l]);
					}
					printf(" ");
				}
				printf("\n");
			}
			printf(" \n");
		}
		printf("--------\n");
		fflush(stdout);
	}
}

void get_block_rank_task(const size_t nt, int block_rank[nt][nt])
{
	const size_t np = nanos6_get_num_cluster_nodes();

	#pragma oss task out(block_rank[0:nt][0:nt])				\
		node(nanos6_cluster_no_offload) label("get_block_rank")
	{
		int row = np, col = np;

		if (np != 1) {
			while (1) {
				row = row / 2;
				if (row * col == np) {
					break;
				}
				col = col / 2;
				if (row * col == np) {
					break;
				}
			}
		}
		dbprintf("# row = %d, col = %d\n", row, col);

		size_t tmp_rank = 0, offset = 0;
		for (size_t i = 0; i < nt; i++) {
			for (size_t j = 0; j < nt; j++) {
				block_rank[i][j] = tmp_rank + offset;
				++tmp_rank;
				if (tmp_rank >= col)
					tmp_rank = 0;
			}
			tmp_rank = 0;
			offset = (offset + col >= np) ? 0 : offset + col;
		}
	}
}

void cholesky_init_task(const size_t nt, const size_t ts,
                        double A[nt][nt][ts][ts], double Ans[nt][nt][ts][ts],
                        int block_rank[nt][nt]
) {
	const size_t dim = nt * ts;

	for (size_t i = 0; i < nt; ++i) {
		for (size_t j = 0; j < nt; ++j) {
			int nodeij = block_rank[i][j];

			if (Ans != NULL) {
				#pragma oss task out(Ans[i][j][0;ts][0;ts])			\
					node(nanos6_cluster_no_offload) label("init_Ans")
				fill_block(ts, Ans[i][j], i, j, dim);

				#pragma oss task weakin(Ans[i][j][0;ts][0;ts])		\
					weakout(A[i][j][0;ts][0;ts])					\
					node(nodeij) label("weak_copy_Ans")
				{
					#pragma oss task in(Ans[i][j][0;ts][0;ts])		\
						out(A[i][j][0;ts][0;ts])					\
						node(nanos6_cluster_no_offload) label("copy_Ans")
					memcpy(A[i][j], Ans[i][j], ts * ts * sizeof(double));
				}

			} else {
				#pragma oss task out(A[i][j][0;ts][0;ts])	\
					node(nodeij) label("init_A")
				fill_block(ts, A[i][j], i, j, dim);
			} // pragma
		} // for j
	} // for i
}

void cholesky_single(const size_t nt, const size_t ts,
                     double A[nt][nt][ts][ts]
) {
	for (size_t k = 0; k < nt; ++k) {
		#pragma oss task inout(A[k][k][0;ts][0;ts])					\
			node(nanos6_cluster_no_offload) label("single_potrf")
		oss_potrf(ts, A[k][k], k,k,k,0);

		for (size_t i = k + 1; i < nt; ++i) {
			#pragma oss task in(A[k][k][0;ts][0;ts])			\
				inout(A[k][i][0;ts][0;ts])								\
				node(nanos6_cluster_no_offload) label("single_trsm")
			oss_trsm(ts, A[k][k], A[k][i], k,k,i,0);
		}

		for (size_t i = k + 1; i < nt; ++i) {
			for (size_t j = k + 1; j < i; ++j) {
				#pragma oss task in(A[k][i][0;ts][0;ts])			\
					in(A[k][j][0;ts][0;ts])							\
					inout(A[j][i][0;ts][0;ts])						\
					node(nanos6_cluster_no_offload) label("single_gemm")
				oss_gemm(ts, A[k][i], A[k][j], A[j][i], k,j,i,0);
			}

			#pragma oss task in(A[k][i][0;ts][0;ts])					\
				inout(A[i][i][0;ts][0;ts])								\
				node(nanos6_cluster_no_offload) label("single_syrk")
			oss_syrk(ts, A[k][i], A[i][i], k,i,i,0);
		} // for i
	} // for k
}

#if TYPEID == 0

void cholesky_ompss2(const size_t nt, const size_t ts,
                     double A[nt][nt][ts][ts],
                     int block_rank[nt][nt],
					 int prvanim
) {
	printf("# cholesky weak\n");
	const size_t np = nanos6_get_num_cluster_nodes();

	for (size_t k = 0; k < nt; ++k) {

		int nodekk = block_rank[k][k];

		#pragma oss task weakinout(A[k][k][0;ts][0;ts])	\
			node(nodekk) label("weak_potrf") priority(nt-1-k)
		{
			#pragma oss task inout(A[k][k][0;ts][0;ts]) \
				node(nanos6_cluster_no_offload) label("potrf") priority(nt-1-k)
			{
				(void)nt;
				oss_potrf(ts, A[k][k], k,k,k, prvanim);
			}
		}

		// Order by block_rank
		for(size_t v = 0; v < np; v++) {
			size_t p = np-1-v;

			for (size_t i = k + 1; i < nt; ++i) {
				int nodeki = block_rank[k][i];
				if (nodeki == p) {

					#pragma oss task weakin(A[k][k][0;ts][0;ts])	\
						weakinout(A[k][i][0;ts][0;ts])				\
						node(nodeki) label("weak_trsm") priority(nt-1-k)
					{
						#pragma oss task in(A[k][k][0;ts][0;ts])	\
							inout(A[k][i][0;ts][0;ts])				\
							node(nanos6_cluster_no_offload) label("trsm") priority(nt-1-k)
						{
							(void)nt;
							oss_trsm(ts, A[k][k], A[k][i], k,k,i, prvanim);
						}
					}
				}
			}
		}

		// Order by block_rank
		for(size_t v = 0; v < np; v++) {
			size_t p = np-1-v;

			for (size_t i = k + 1; i < nt; ++i) {
				for (size_t j = k + 1; j < i; ++j) {
					int nodeji = block_rank[j][i];
					if (nodeji == p) {

						#pragma oss task weakin(A[k][i][0;ts][0;ts])		\
							weakin(A[k][j][0;ts][0;ts])						\
							weakinout(A[j][i][0;ts][0;ts])					\
							node(nodeji) label("weak_gemm") priority(nt-1-j)
						{
							#pragma oss task in(A[k][i][0;ts][0;ts])		\
								in(A[k][j][0;ts][0;ts])						\
								inout(A[j][i][0;ts][0;ts])					\
								node(nanos6_cluster_no_offload) label("gemm") priority(nt-1-j)
							{
								(void)nt;
								oss_gemm(ts, A[k][i], A[k][j], A[j][i], k,j,i, prvanim);
							}
						}
					}
				} // for j

				int nodeii = block_rank[i][i];
				if (nodeii == p) {

					#pragma oss task weakin(A[k][i][0;ts][0;ts])	\
						weakinout(A[i][i][0;ts][0;ts])				\
						node(nodeii) label("weak_syrk") priority(nt-1-i)
					{
						int p = (i == k+1) ? 1 : 0;
						#pragma oss task in(A[k][i][0;ts][0;ts])	\
							inout(A[i][i][0;ts][0;ts])				\
							node(nanos6_cluster_no_offload) label("syrk") priority(nt-1-i)
						{
							(void)nt;
							oss_syrk(ts, A[k][i], A[i][i], k,i,i, prvanim);
						}
					}
				}
			} // for i
		} // for k
	} // for v
}

#elif TYPEID == 1

void cholesky_ompss2(const size_t nt, const size_t ts,
                     double A[nt][nt][ts][ts],
                     int block_rank[nt][nt],
					 int prvanim
) {
	printf("# cholesky strong\n");

	for (size_t k = 0; k < nt; ++k) {

		int nodekk = block_rank[k][k];

		#pragma oss task inout(A[k][k][0;ts][0;ts]) \
			node(nodekk) label("potrf")
		oss_potrf(ts, A[k][k],  k,k,k, prvanim);

		for (size_t i = k + 1; i < nt; ++i) {
			int nodeki = block_rank[k][i];

			#pragma oss task in(A[k][k][0;ts][0;ts])	\
				inout(A[k][i][0;ts][0;ts])				\
				node(nodeki) label("trsm")
			oss_trsm(ts, A[k][k], A[k][i], k,k,i, prvanim);
		}

		for (size_t i = k + 1; i < nt; ++i) {
			for (size_t j = k + 1; j < i; ++j) {
				int nodeji = block_rank[j][i];

				#pragma oss task in(A[k][i][0;ts][0;ts])			\
					in(A[k][j][0;ts][0;ts])							\
					inout(A[j][i][0;ts][0;ts])						\
					node(nodeji) label("gemm")
				oss_gemm(ts, A[k][i], A[k][j], A[j][i], k,j,i, prvanim);
			}

			int nodeii = block_rank[i][i];

			#pragma oss task in(A[k][i][0;ts][0;ts])		\
				inout(A[i][i][0;ts][0;ts])					\
				node(nodeii) label("syrk")
			oss_syrk(ts, A[k][i], A[i][i], k, i,i, prvanim);
		} // for i
	} // for k
}

#else  // TYPE
#error Cholesky type value not valid.
#endif // TYPE

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
