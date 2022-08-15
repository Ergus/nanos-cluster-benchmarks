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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CHOLESKY_OMP_MPI_H
#define CHOLESKY_OMP_MPI_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>

#include "benchmarks_ompss.h"

void oss_potrf(int ts, double A[ts][ts], int k, int y, int x, int prvanim)
{
	inst_blas_kernel(prvanim, BLAS_POTRF, k, y, x);
	assert(A != NULL);

	static int INFO;
	static const char L = 'L';

	dpotrf_(&L, &ts, (double *)A, &ts, &INFO);
	inst_blas_kernel(prvanim, BLAS_NONE, k,y,x);
}

void oss_trsm(int ts, double A[ts][ts], double B[ts][ts], int k, int y, int x, int prvanim)
{
	inst_blas_kernel(prvanim, BLAS_TRSM, k,y,x);
	assert(A != NULL);
	assert(B != NULL);


    char LO = 'L', TR = 'T', NU = 'N', RI = 'R';
    double DONE = 1.0;
    dtrsm_(&RI, &LO, &TR, &NU, &ts, &ts, &DONE,
	       (double *)A, &ts,
	       (double *)B, &ts);

	inst_blas_kernel(prvanim, BLAS_NONE, k,y,x);
}

void oss_gemm(
	int ts,
	double A[ts][ts],
	double B[ts][ts],
	double C[ts][ts],
	int k, int y, int x, int prvanim
) {
	inst_blas_kernel(prvanim, BLAS_GEMM, k,y,x);
	assert(A != NULL);
	assert(B != NULL);
	assert(C != NULL);

    const char TR = 'T', NT = 'N';
    double DONE = 1.0, DMONE = -1.0;
    dgemm_(&NT, &TR, &ts, &ts, &ts, &DMONE,
	       (double *)A, &ts,
	       (double *)B, &ts, &DONE,
	       (double *)C, &ts);

	inst_blas_kernel(prvanim, BLAS_NONE, k,y,x);
}

void oss_syrk(int ts, double A[ts][ts], double B[ts][ts], int k, int y, int x, int prvanim)
{
	inst_blas_kernel(prvanim, BLAS_SYRK, k,y,x);
	assert(A != NULL);
	assert(B != NULL);

    static char LO = 'L', NT = 'N';
    static double DONE = 1.0, DMONE = -1.0;
    dsyrk_(&LO, &NT, &ts, &ts, &DMONE,
	       (double *)A, &ts, &DONE,
	       (double *)B, &ts);
	inst_blas_kernel(prvanim, BLAS_NONE, k,y,x);
}



void print_matrix_task(const size_t nt, const size_t ts,
                       double A[nt][nt][ts][ts], int nodeid
) {
	#pragma oss task in(A[0;nt][0;nt][0;ts][0;ts])	\
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

	#pragma oss task out(block_rank[0;nt][0;nt])				\
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

//! This is to initialize blocks[i][j]
void fill_block(const size_t ts, double block[ts][ts],
                const size_t i, const size_t j, const size_t dim
) {
	myassert(block);

	const size_t seed = (i > j) ? i * ts + j : j * ts + i;
	struct drand48_data status;       	// using re-entrant version for rng
	srand48_r(seed, &status);
	double rnd1, rnd2;

	for (size_t k = 0; k < ts; ++k) {
		if (i == j) {
			drand48_r(&status, &rnd1);
			drand48_r(&status, &rnd2);
			block[k][k] = rnd1 * rnd2 + dim;

			for (size_t l = k + 1; l < ts; ++l) {
				drand48_r(&status, &rnd1);
				drand48_r(&status, &rnd2);

				const double val = rnd1 * rnd2;
				block[k][l] = val;
				block[l][k] = val;
			}

		} else {
			for (size_t l = 0; l < ts; ++l) {
				drand48_r(&status, &rnd1);
				drand48_r(&status, &rnd2);

				if (i > j) {
					block[k][l] = rnd1 * rnd2;
				} else {
					block[l][k] = rnd1 * rnd2;
				}
			}
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

bool compare_blocks(const size_t ts, double B1[ts][ts], double B2[ts][ts])
{
	for (size_t k = 0; k < ts; ++k) {
		for (size_t l = 0; l < ts; ++l) {
			if (B1[k][l] != B2[k][l]) {
				printf("Check failed for B[%zu][%zu] = %lf expected: %lf\n",
				       k, l,
				       B1[k][l], B2[k][l]);
				return false;
			}
		}
	}
	return true;
}

void print_block(const size_t ts, double A[ts][ts])
{

	for (size_t i = 0; i < ts; ++i) {
		for (size_t j = 0; j < ts; ++j) {
			printf("%5.2f ", (float)A[i][j]);
		}
		printf("\n");
	}
	fflush(stdout);
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


#ifdef __cplusplus
}
#endif


#endif /* CHOLESKY_OMP_MPI_H */
