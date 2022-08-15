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


void cholesky_memory_init_task(const struct matrix_info *pinfo,
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

void cholesky_memory_single(
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

void cholesky_memory_ompss2(
	const struct matrix_info *pinfo,
	double A[pinfo->nt * pinfo->nt][pinfo->ts][pinfo->ts]
) {
	printf("# cholesky memory\n");

	const size_t dim = pinfo->nt * pinfo->ts;
	const size_t nt = pinfo->nt;
	const size_t ts = pinfo->ts;
	const size_t np = pinfo->np;

	// Try a copy to trick mercurium.
	const struct matrix_info info = *pinfo;

	const size_t pcols = info.pcols;
	const size_t npcols = info.nt / info.pcols;
	const size_t blocks_per_node = info.blocks_per_node;

	for (size_t k = 0; k < nt; ++k) {

		const int nodekk = get_block_node(&info, k, k);

		double (*Akk)[ts] = A[ get_block_global_index(&info, k, k) ];

		#pragma oss task weakinout(Akk[0;ts][0;ts])	\
			node(nodekk) label("weak_potrf")
		{
			#pragma oss task inout(Akk[0;ts][0;ts]) \
				node(nanos6_cluster_no_offload) label("potrf")
			oss_potrf(ts, Akk, k,k,k, 0);
		}

		const size_t nodek0 = get_block_node(&info, k, 0);  // node of first block

		// Row k trsm
		for (size_t it = 0; it < pcols; ++it) {

			const size_t node = nodek0 + it;

			const size_t start = get_block_global_index(&info, k, it);

			#pragma oss task weakin(Akk[0;ts][0;ts])			   \
				weakinout(A[start;npcols][0;ts][0;ts])			   \
				node(node) label("weak_trsm")
			{
				for (size_t i = k + 1; i < nt; ++i) {
					int nodeki = get_block_node(&info, k, i);

					if (nodeki != node)
						continue;

					double (*Aki)[ts] = A[ get_block_global_index(&info, k, i) ];

					#pragma oss task in(Akk[0;ts][0;ts])			\
						inout(Aki[0;ts][0;ts])						\
						node(nanos6_cluster_no_offload) label("trsm")
					oss_trsm(ts, Akk, Aki, k,k,i, 0);
				}
			}
		}

		const size_t idxk0 = get_block_global_index(&info, k, 0);     // idx of first block

		// info of first block in row k + 1
		const size_t nodek10 = get_block_node(&info, k + 1, 0);
		const size_t lidxk10 = get_block_local_index(&info, k + 1, 0);

		for (size_t node = 0; node < info.np; ++node) {

			const size_t first
				= node * blocks_per_node                  // start of this
				+ lidxk10                                 // offset
				+ ((node < nodek10) ? npcols : 0);

			const int count = (node + 1) * blocks_per_node - first;

			if (first >= info.np * blocks_per_node || count <= 0) {
				continue;
			}


			#pragma oss task weakin( {A[idxk0 + i*blocks_per_node; npcols][0;ts][0;ts], i=0;pcols} ) \
				weakinout(A[first; count][0;ts][0;ts])					\
				node(node) label("weak_syrk")
			{
				for (size_t i = k + 1; i < nt; ++i) {
					int nodeii = get_block_node(&info, i, i);

					if (nodeii == node) {
						double (*Aii)[ts] = A[ get_block_global_index(&info, i, i) ];
						double (*Aki)[ts] = A[ get_block_global_index(&info, k, i) ];

						#pragma oss task in(Aki[0;ts][0;ts])	\
							inout(Aii[0;ts][0;ts])				\
							node(nanos6_cluster_no_offload) label("syrk")
						oss_syrk(ts, Aki, Aii, k,i,i, 0);
					}


					for (size_t j = k + 1; j < i; ++j) {
						int nodeji = get_block_node(&info, j, i);

						if (nodeji == node) {

							double (*Aki)[ts] = A[ get_block_global_index(&info, k, i) ];
							double (*Akj)[ts] = A[ get_block_global_index(&info, k, j) ];
							double (*Aji)[ts] = A[ get_block_global_index(&info, j, i) ];

							#pragma oss task in(Aki[0;ts][0;ts])		\
								in(Akj[0;ts][0;ts])						\
								inout(Aji[0;ts][0;ts])					\
								node(nanos6_cluster_no_offload) label("gemm")
							oss_gemm(ts, Aki, Akj, Aji, k, j, i, 0);
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
	const size_t ROWS = create_cl_size_t("Rows");
	const size_t TS = create_cl_size_t("Tasksize");
	const int CHECK = create_optional_cl_int("Check", 0);
	modcheck(ROWS, TS);

	inst_register_events();  // Register the events in the instrumentation

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
	cholesky_memory_init_task(&info, A, Ans);
	#pragma oss taskwait

	// ===========================================
	printf("# Starting warmup\n");

	timer atimer_warmup = create_timer("Warmup_time");
	cholesky_memory_ompss2(&info, A);
	#pragma oss taskwait

	stop_timer(&atimer_warmup);
	printf("# Finished warmup\n");

	// ===========================================
	// ACTUAL CALCULATION
	cholesky_memory_init_task(&info, A, Ans);
	#pragma oss taskwait

	// ===========================================
	printf("# Starting algorithm\n");

	timer atimer = create_timer("Algorithm_time");
	cholesky_memory_ompss2(&info, A);
	#pragma oss taskwait

	stop_timer(&atimer);
	// ===========================================

	printf("# Finished algorithm\n");

	if (CHECK) {
		timer stimer = create_timer("Single_time");
		cholesky_memory_single(&info, Ans);
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

	myassert(A);
	nanos6_dfree(A, ROWS * ROWS * sizeof(double));

	free_args();

	return 0;
}
