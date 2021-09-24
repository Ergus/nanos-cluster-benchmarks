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

#include <unistd.h>
#include "cholesky_fare_ompss2_taskfor.h"

struct matrix_info {
	size_t nt;                     // number of blocks
	size_t ts;                     // block_size
	size_t np;                     // Number of processes
	size_t prows, pcols;           // processes cols and rows
	size_t blocks_per_node;        // blocks per node.
};

// node block i, j
size_t get_block_node(const struct matrix_info *info, size_t i, size_t j)
{
	return (i % info->prows) * info->pcols + j % info->pcols;
}

// local index for block i, j in node
size_t get_block_local_index(const struct matrix_info *info, size_t i, size_t j)
{
	return (i / info->prows) * (info->nt / info->pcols) + (j / info->pcols);
}

// global index for block i, j
size_t get_block_global_index(const struct matrix_info *info, size_t i, size_t j)
{
	const size_t node = get_block_node(info, i, j);
	const size_t blocks_per_node = info->nt * info->nt / info->np;
	const size_t local_index = get_block_local_index(info, i, j);

	return node * blocks_per_node + local_index;
}

// Round up x to be congruent to y modulo mod
// e.g. round_up(5, 3, 10) = 13, because the smallest number
// >= 5 that is congruent to 3 modulo 10 is 13
int round_up(int x, int y, int mod)
{
	// Assume 0 <= y < mod
	int delta = (x + mod - y) % mod;
	return delta ? x + (mod-delta) : x;
}

/* void print_matrix_task( */
/* 	size_t nt, size_t ts, */
/* 	double A[nt * nt][ts][ts], */
/* 	int nodeid */
/* 	const struct matrix_info info, */
/* ) { */
/* 	size_t nt = info.nt; */
/* 	size_t ts = info.ts; */

/* 	#pragma oss task in(A[0:nt][0:nt][0;ts][0;ts]) in(info[0]) \ */
/* 		node(nodeid) label("print_matrix") */
/* 	{ */
/* 		for (size_t i = 0; i < nt; i++) { */
/* 			for (size_t k = 0; k < ts; ++k) { */
/* 				for (size_t j = 0; j < nt; ++j) { */
/* 					const size_t idx = get_block_global_index(info, i, j); */
/* 					double lA[ts][ts] = A[idx]; */

/* 					for (size_t l = 0; l < ts; ++l) { */
/* 						printf("%5.2f ", (float)lA[k][l]); */
/* 					} */
/* 					printf(" "); */
/* 				} */
/* 				printf("\n"); */
/* 			} */
/* 			printf(" \n"); */
/* 		} */
/* 		printf("--------\n"); */
/* 		fflush(stdout); */
/* 	} */
/* } */

void get_block_info(int ROWS, int TS, struct matrix_info *info)
{
	info->nt = ROWS / TS;
	info->ts = TS;
	info->np = nanos6_get_num_cluster_nodes();
	modcheck(info->nt, info->np);

	info->blocks_per_node = info->nt * info->nt / info->np;

	size_t lrow = info->np;
	size_t lcol = info->np;

	// get num rows and columns.
	if (info->np != 1) {
		while (1) {
			lrow = lrow / 2;
			if (lrow * lcol == info->np) {
				break;
			}
			lcol = lcol / 2;
			if (lrow * lcol == info->np) {
				break;
			}
		}
	}

	info->prows = lrow;
	info->pcols = lcol;

	dbprintf("# row = %d, col = %d\n", lrow, lcol);
}


void cholesky_init_task(const struct matrix_info *pinfo,
                        double A[pinfo->nt * pinfo->nt][pinfo->ts][pinfo->ts],
                        double Ans[pinfo->nt * pinfo->nt][pinfo->ts][pinfo->ts]
) {
	struct matrix_info info = *pinfo;
	const size_t dim = pinfo->nt * pinfo->ts;
	const size_t nt = pinfo->nt;
	const size_t ts = pinfo->ts;
	const size_t np = pinfo->np;

	for (size_t i = 0; i < nt; ++i) {
		for (size_t j = 0; j < nt; ++j) {
			int nodeij = get_block_node(pinfo, i, j);
			double (*Aij)[ts] = A[ get_block_global_index(&info, i, j) ];

			if (Ans != NULL) {
				double (*Ansij)[ts] = Ans[ get_block_global_index(&info, i, j) ];

				#pragma oss task out(Ansij[0;ts][0;ts])			\
					node(nanos6_cluster_no_offload) label("init_Ans")
				fill_block(ts, Ansij, i, j, dim);

				#pragma oss task weakin(Ansij[0;ts][0;ts])		\
					weakout(Aij[0;ts][0;ts])					\
					node(nodeij) label("weak_copy_Ans")
				{
					#pragma oss task in(Ansij[0;ts][0;ts])		\
						out(Aij[0;ts][0;ts])					\
						node(nanos6_cluster_no_offload) label("copy_Ans")
					memcpy(Aij, Ansij, ts * ts * sizeof(double));
				}

			} else {
				#pragma oss task out(Aij[0;ts][0;ts])	\
					node(nodeij) label("init_A")
				fill_block(ts, Aij, i, j, dim);
			} // pragma
		} // for j
	} // for i
}

void cholesky_single(
	const struct matrix_info *pinfo,
	double A[pinfo->nt * pinfo->nt][pinfo->ts][pinfo->ts]
) {
	const size_t dim = pinfo->nt * pinfo->ts;

	const size_t nt = pinfo->nt;
	const size_t ts = pinfo->ts;
	const size_t np = pinfo->np;

	for (size_t k = 0; k < nt; ++k) {

		double (*Akk)[ts] = A[ get_block_global_index(pinfo, k, k) ];

		#pragma oss task inout(Akk[0;ts][0;ts])						\
			node(nanos6_cluster_no_offload) label("single_potrf")
		oss_potrf(ts, Akk, k,k,k, 0);

		for (size_t i = k + 1; i < nt; ++i) {
			double (*Aki)[ts] = A[get_block_global_index(pinfo, k, i)];

			#pragma oss task in(Akk[0;ts][0;ts])						\
				inout(Aki[0;ts][0;ts])									\
				node(nanos6_cluster_no_offload) label("single_trsm")
			oss_trsm(ts, Akk, Aki, k,k,i,0);

		}

		for (size_t i = k + 1; i < nt; ++i) {
			double (*Aki)[ts] = A[get_block_global_index(pinfo, k, i)];
			double (*Aii)[ts] = A[get_block_global_index(pinfo, i, i)];

			#pragma oss task in(Aki[0;ts][0;ts])					\
				inout(Aii[0;ts][0;ts])								\
				node(nanos6_cluster_no_offload) label("single_syrk")
			oss_syrk(ts, Aki, Aii, k,i,i,0);

			for (size_t j = k + 1; j < i; ++j) {

				double (*Akj)[ts] = A[get_block_global_index(pinfo, k, j)];
				double (*Aji)[ts] = A[get_block_global_index(pinfo, j, i)];

				#pragma oss task in(Aki[0;ts][0;ts])					\
					in(Akj[0;ts][0;ts])									\
					inout(Aji[0;ts][0;ts])								\
					node(nanos6_cluster_no_offload) label("single_gemm")
				oss_gemm(ts, Aki, Akj, Aji, k,j,i,0);
			}
		} // for i
	} // for k
}

inline int min(int x, int y)
{
	return (x < y) ? x : y;
}

inline int max(int x, int y)
{
	return (x < y) ? y : x;
}

void cholesky_ompss2(
	const struct matrix_info *pinfo,
	double A[pinfo->nt * pinfo->nt][pinfo->ts][pinfo->ts]
) {
	int prvanim = 1;
	const size_t dim = pinfo->nt * pinfo->ts;
	const size_t nt = pinfo->nt;
	const size_t ts = pinfo->ts;
	const size_t np = pinfo->np;

	printf("# cholesky weak\n");

	// Try a copy to trick mercurium.
	const struct matrix_info info = *pinfo;

	const size_t pcols = info.pcols;
	const size_t npcols = info.nt / info.pcols;
	const size_t blocks_per_node = info.blocks_per_node;
	int delay[64];

	// Size of each taskfor in terms of the whole matrix (the number
	// of elements in an actual taskfor is chop / pcols)
	int chop_size = nt / 2;
	int num_chops = nt / chop_size;

	assert((chop_size % pcols) == 0); // chop should be a multiple of pcols

	for (size_t k = 0; k < nt; ++k) {

		// POTRF
		const int nodekk = get_block_node(&info, k, k);
		int idx_kk = get_block_global_index(&info, k, k);
		double (*Akk)[ts] = A[ get_block_global_index(&info, k, k) ];

		#pragma oss task weakinout(Akk[0;ts][0;ts])	\
			node(nodekk) label("weak_potrf") priority(nt)
		{
			#pragma oss task inout(Akk[0;ts][0;ts]) \
				node(nanos6_cluster_no_offload) label("potrf") priority(nt)
				{
				 (void)nt;
				 oss_potrf(ts, Akk, k,k,k, prvanim);
				// printf("%ld. POTRF on node %d      inout %d\n", k, nodekk, idx_kk);
				}
		}

		// Row k trsm
		// Calculate global start index and size for each node in the k-th row
		int first_i_n[8], start_p[8], len_p[8], ofs_p[8];
		for (size_t x = 0; x < info.pcols; x++) {
			first_i_n[x] = round_up(k+1, x, pcols); 
			start_p[x] = get_block_global_index(&info, k, first_i_n[x]);
			ofs_p[x] = get_block_global_index(&info, k, x);
			len_p[x] = (nt - first_i_n[x] + pcols-1) / pcols;
		}

		// Chop up the new row to get a start and end global index for each chop
		int start_global[num_chops+1][pcols];
		for(int x = 0; x < pcols; x++) {
			for(int chop = 0; chop < num_chops; chop++) {
				int first_i_chop = x + chop*chop_size;
				start_global[chop][x] = get_block_global_index(&info, k, first_i_chop);
			}
			start_global[num_chops][x] = get_block_global_index(&info, k+1, x);
		}

// Chop containing a given column
#define CHOP(col)            ((col)/chop_size)
// Start of a given chop: column number on current node
#define START_CHOP_COL(chop, row,pcol) max((chop)*chop_size + pcol, round_up(row+1, pcol, pcols))
// End of a given chop: column number on current node
#define END_CHOP_COL(chop,pcol) (((chop)+1)*chop_size + pcol)
// Start of a given chop: global index on current node
#define START_CHOP_IDX(chop,row,pcol) (get_block_global_index(&info, row, pcol) + (START_CHOP_COL(chop,row,pcol))/pcols)
// Length of a given chop: elements to process
#define CHOP_LEN(chop,row,pcol) ((END_CHOP_COL(chop,pcol)-START_CHOP_COL(chop,row,pcol))/pcols)


		const size_t nodek0 = get_block_node(&info, k, 0);        // node of first block

		for (size_t x = 0; x < pcols; ++x) {
			const size_t node = nodek0 + x;
			if (len_p[x] > 0) {
				#pragma oss task weakin(Akk[0;ts][0;ts])			\
					weakinout(A[start_p[x];len_p[x]][0;ts][0;ts])			\
					node(node) label("weak_trsm") priority(2+nt-1-k)
				{
					// Chop all columns
					int first_i = first_i_n[x];
					for(int chop = CHOP(first_i); chop < num_chops; chop++) {
						int start_col = START_CHOP_COL(chop,k,x);
						int start_idx = START_CHOP_IDX(chop,k,x);
						int len = CHOP_LEN(chop,k,x);
						#pragma oss task for in(Akk[0;ts][0;ts]) inout(A[start_idx;len][0;ts][0;ts]) \
									node(nanos6_cluster_no_offload) label("trsm") priority(nt-1-k)
						for (int w=0; w<len; w++) {
							int i = start_col + w * pcols;
							(void)nt; (void)start_idx; (void)len_p;
							double (*Aki)[ts] = A[start_idx+w];
							oss_trsm(ts, Akk, Aki, k,k,i, prvanim);
						}
					}
				}
			}
		}

		// info of first block in row k + 1
		const size_t nodek10 = get_block_node(&info, k + 1, 0);
		const size_t lidxk10 = get_block_local_index(&info, k + 1, 0);

		// Node rows
		for (size_t y = 0; y < info.prows; ++y) {
			// Node columns
			for (size_t x = 0; x < info.pcols; ++x) {
				int node = y * info.pcols + x;

				int first_j = round_up(k+1, y, info.prows);

				const size_t first
					= node * blocks_per_node                  // start of this
					+ lidxk10                                 // offset 
					+ ((node < nodek10) ? npcols : 0);

				const int count = (node + 1) * blocks_per_node - first;

				if (first >= info.np * blocks_per_node || count <= 0) {
					continue;
				}

				#pragma oss task weakin( {A[start_p[i];len_p[i]][0;ts][0;ts], i=0;pcols} ) \
					weakinout(A[first; count][0;ts][0;ts])					\
					node(node) label("weak_syrk") priority(2+nt-1-k) weakinout(delay[node])
				{
					 nanos6_set_early_release(nanos6_no_wait);
					for (size_t j = first_j; j < nt; ++j) {
						int nodejj = get_block_node(&info, j, j);

						if (nodejj == node) {
							double (*Ajj)[ts] = A[ get_block_global_index(&info, j, j) ];
							double (*Akj)[ts] = A[ get_block_global_index(&info, k, j) ];

							// Identify chop containing the Akj element
							int chop = CHOP(j);
							int start_idx = START_CHOP_IDX(chop,k,x);
							int len = CHOP_LEN(chop,k,x);

							// Take whole row on relevant node to avoid splitting dependency
							#pragma oss task in(A[start_idx;len][0;ts][0;ts]) \
								inout(Ajj[0;ts][0;ts])				\
								node(nanos6_cluster_no_offload) label("syrk") priority(nt-1-j)
							{
								(void)nt; (void)start_p; (void)len_p;
								oss_syrk(ts, Akj, Ajj, k,j,j, prvanim);
							}
						}
					}
					
					int hack_val = round_up(num_chops,0,2) + 1; // and round up to next odd number (assuming 2 sockets)
					int hack_count = 0;

					for (size_t j = first_j; j < nt+1; j += info.prows) {
						int first_i = round_up(j+1, x, info.pcols);
						if (first_i >= nt) {
							// No work for this node
							continue;
						}

						// All GEMM operations in this j iteration reference
						// A[k,j]; get the pointer to it and the whole chop containing it
						double (*Akj)[ts] = A[ get_block_global_index(&info, k, j) ];
						int chop_kj = CHOP(j);
						int colkj = j % pcols;
						int start_chop_kj = START_CHOP_IDX(chop_kj,k,colkj);
						int len_chop_kj = CHOP_LEN(chop_kj,k,colkj);

						for(int chop = CHOP(first_i); chop < num_chops; chop++) {

							// Inout elements to update on row j
							int start_idx = START_CHOP_IDX(chop,j,x);
							int start_i = START_CHOP_COL(chop,j,x);
							int len = CHOP_LEN(chop,j,x);

							// Whole chop for ki
							int start_chop_ki = START_CHOP_IDX(chop,k,x);
							int len_chop_ki = CHOP_LEN(chop,k,x);

							// First element accessed by ki
							int start_idx_ki = get_block_global_index(&info, k, start_i);

							// Hack scheduling: we introduce a delay after the first two rows
							// which will force the scheduler to schedule the TRSMs for the next
							// k earlier. Otherwise all the GEMM taskfors occupy the cores so
							// the TRSMs get scheduled at the end of the whole calculation.
							if (hack_count < hack_val) {
								// First two have delay as input
								#pragma oss task for in(A[start_chop_ki;len_chop_ki][0;ts][0;ts]) \
									in(A[start_chop_kj;len_chop_kj][0;ts][0;ts]) \
									inout(A[start_idx;len][0;ts][0;ts])					\
									node(nanos6_cluster_no_offload) label("gemm1") priority(nt-1-j) in(delay[node])
									for (int w=0; w<len; w++) {
										int i = start_i + w*pcols;
										(void)nt; (void)delay;
										assert(start_idx_ki +w >= start_chop_ki && start_idx_ki+w < start_chop_ki+len_chop_ki);
										assert(get_block_global_index(&info,k,j) >= start_chop_kj);
										assert(get_block_global_index(&info,k,j) < start_chop_kj + len_chop_kj);
										double (*Aki)[ts] = A[start_idx_ki + w];
										double (*Aji)[ts] = A[start_idx + w];
										oss_gemm(ts, Aki, Akj, Aji, k,j,i,prvanim);
									}
							} else if (hack_count == hack_val) {
								// Next in parallel with delay
								#pragma oss task for in(A[start_chop_ki;len_chop_ki][0;ts][0;ts]) \
									in(A[start_chop_kj;len_chop_kj][0;ts][0;ts]) \
									inout(A[start_idx;len][0;ts][0;ts])					\
									node(nanos6_cluster_no_offload) label("gemm2") priority(nt-1-j)
									for (int w=0; w<len; w++) {
										int i = start_i + w*pcols;
										(void)nt; (void)delay;
										double (*Aki)[ts] = A[start_idx_ki + w];
										double (*Aji)[ts] = A[start_idx + w];
										oss_gemm(ts, Aki, Akj, Aji, k,j,i,prvanim);
									}
									#pragma oss task inout(delay[node])
									{
										(void)delay;
										usleep(1);
									}
								} else {
								// Last after the delay
								#pragma oss task for in(A[start_chop_ki;len_chop_ki][0;ts][0;ts]) \
									in(A[start_chop_kj;len_chop_kj][0;ts][0;ts]) \
									inout(A[start_idx;len][0;ts][0;ts])					\
									node(nanos6_cluster_no_offload) label("gemm3") priority(nt-1-j) in(delay[node])
									for (int w=0; w<len; w++) {
										int i = start_i + w*pcols;
										(void)nt; (void)delay;
										double (*Aki)[ts] = A[start_idx_ki + w];
										double (*Aji)[ts] = A[start_idx + w];
										oss_gemm(ts, Aki, Akj, Aji, k,j,i,prvanim);
									}
								}
								hack_count ++;
							}
						}
					}
				}
			}
	} // for k
}

int main(int argc, char *argv[])
{
	init_args(argc, argv);

	const char *PREFIX = basename(argv[0]);
	const int ROWS = create_cl_int("Rows");
	const int TS = create_cl_int("Tasksize");
	const int CHECK = create_optional_cl_int("Check", 0);
	modcheck(ROWS, TS);

	printf("# Initializing data\n");;

	struct matrix_info info;
	get_block_info(ROWS, TS, &info);

	timer ttimer = create_timer("Total_time");

	double (*A)[TS][TS] = nanos6_dmalloc(ROWS * ROWS * sizeof(double),
	                                     nanos6_equpart_distribution, 0, NULL);

	myassert(A != NULL);

	double (*Ans)[TS][TS] = NULL;
	if (CHECK) {
		Ans = nanos6_dmalloc(ROWS * ROWS * sizeof(double),
		                     nanos6_equpart_distribution, 0, NULL);
		myassert(Ans != NULL);
	}

	#pragma oss taskwait

	// ===========================================
	// WARMUP: The first iteration is slow, likely because of the
	// overhead of thread creation, which seems to impact MPI. Extra
	// threads are needed for the weak tasks, which stay alive for
	// autowait. Note: we want them to stay alive also so that there
	// is more opportunity to connect later tasks in the namespace.
	cholesky_init_task(&info, A, Ans);
	#pragma oss taskwait

	// ===========================================
	printf("# Starting warmup\n");

	timer atimer_warmup = create_timer("Warmup_time");
	cholesky_ompss2(&info, A);
	#pragma oss taskwait

	stop_timer(&atimer_warmup);
	// ===========================================

	printf("# Finished warmup\n");

#if 0
	{
		int i,j;
		int NT = ROWS / TS;
		for(i=0;i<NT;i++) {
			for(j=0; j<NT; j++) {
				int idx = get_block_global_index(&info, i,j);
				double *Aij = A[idx];
				printf("%p ", Aij);
			}
			printf("\n");
		}
	}
#endif

	// ===========================================
	// ACTUAL CALCULATION
	cholesky_init_task(&info, A, Ans);


#if 0 // not needed as taskwaits unfragment even if different writeid
	// Distribute neatly by node
	for (int node=0; node < nanos6_get_num_cluster_nodes(); node++) {
		int bn = info.blocks_per_node;
		double (*An)[TS] = A[node*bn];
		// printf("node %d %p\n", node, An);
		#pragma oss task inout(An[0;bn]) node(node)
		{
		}
	}
#endif

	#pragma oss taskwait

	// ===========================================
	printf("# Starting algorithm\n");
	for (int node=0; node < nanos6_get_num_cluster_nodes(); node++) {
		#pragma oss task node(node)
		nanos6_stats_clear();
	}
	timer atimer = create_timer("Algorithm_time");
	cholesky_ompss2(&info, A);
	#pragma oss taskwait

	stop_timer(&atimer);
	// ===========================================

	printf("# Finished algorithm\n");

	if (CHECK) {
		timer stimer = create_timer("Single_time");
		cholesky_single(&info, Ans);
		#pragma oss taskwait
		stop_timer(&stimer);
	}

	stop_timer(&ttimer);

	// Verification
	if (Ans != NULL) {
		size_t nt = info.nt;
		size_t ts = info.ts;

		#pragma oss task in(A[0;nt*nt][0;TS][0;TS])			\
			in(Ans[0;nt*nt][0;TS][0;TS])					\
			node(nanos6_cluster_no_offload) label("check")
		{
			bool match = true;
			for (size_t i = 0; i < nt && match; i++) {
				for (size_t j = 0; j < nt && match; j++) {
					double (*Aij)[ts] = A[ get_block_global_index(&info, i, j) ];
					double (*Ansij)[ts] = Ans[ get_block_global_index(&info, i, j) ];

					match = compare_blocks(TS, Aij, Ansij);

					if (!match) {
						fprintf(stderr, "# Check failed in block A[%zu][%zu]\n",
						        i, j);
						assert(0);
					}
				}
			}
		}
	}

	#pragma oss taskwait

	create_reportable_int("worldsize", nanos6_get_num_cluster_nodes());
	create_reportable_int("cpu_count", nanos6_get_num_cpus());
	create_reportable_int("namespace_enabled", nanos6_get_namespace_is_enabled());

	report_args();

	// Release memory
	if (Ans != NULL) {
		nanos6_dfree(Ans, ROWS * ROWS * sizeof(double));
	}

	myassert(A);
	nanos6_dfree(A, ROWS * ROWS * sizeof(double));

	free_args();

	return 0;
}
