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


#include "cholesky_utils.h"

#include "ArgParserC/argparser.h"

//! Initialize the blocked matrix

void cholesky(const size_t nblocks,
              const size_t bsize,
              double A[nblocks][nblocks][bsize][bsize]
) {
	const size_t numNodes = nanos6_get_num_cluster_nodes();
	myassert(nblocks >= numNodes);
	modcheck(nblocks, numNodes);

	const size_t blocks_per_node = nblocks / numNodes;

	for (size_t i = 0; i < nblocks; ++i) {
		int nodei = i / blocks_per_node;

		#pragma oss task weakinout(A[i][i][0;bsize][0;bsize])	\
			node(nodei) label("weak_potrf")
		{
			oss_potrf(bsize, A[i][i]);    // Diagonal Block Factorization
		}

		#pragma oss task weakin(A[i][i][0;bsize][0;bsize])		\
			weakinout(A[i][i+1:nblocks-1][0;bsize][0;bsize])	\
			node(nodei) label("weak_trsm")
		{
			for (size_t j = i + 1; j < nblocks; ++j)      // Triangular Systems
				oss_trsm(bsize, A[i][i], A[i][j]);
		}

		for (size_t j = i + 1; j < nblocks; ++j) {    // Update Trailing Matrix
			int nodej = j / blocks_per_node;

			#pragma oss task weakin(A[i][j][0;bsize][0;bsize])		\
				weakout(A[j][j][0;bsize][0;bsize])					\
				node(nodej) label("weak_syrk")
			{
				oss_syrk(bsize, A[i][j], A[j][j]);
			}

			#pragma oss task weakin(A[i][i+1:j][0;bsize][0;bsize])		\
				weakinout(A[i+1:j-1][i+1:nblocks-1][0;bsize][0;bsize])	\
				node(nodei) label("weak_gemm")
			{
				for (size_t k = i + 1; k < j; ++k) {
					int nodek = k / blocks_per_node;
					#pragma oss task weakin(A[i][j][0;bsize][0;bsize])	\
						weakin(A[i][k][0;bsize][0;bsize])				\
						weakout(A[k][j][0;bsize][0;bsize])				\
						node(nodek) label("weak_gemm")
					oss_gemm(bsize, A[i][j], A[i][k], A[k][j]);
				}
			}
		}
	}
}

int main(int argc, char **argv)
{
	init_args(argc, argv);

	const int DIM = create_cl_int("Dimension");
	const int BSIZE = create_cl_int("Block_size");
	const int CHECK = create_optional_cl_int("Check", 0);

	myassert(DIM >= BSIZE);
	modcheck(DIM, BSIZE);

	//========= End Command Line ===============

	const size_t nblocks = DIM / BSIZE;          // Number of blocks
	const size_t dim2 = DIM * DIM * sizeof(double);
	printf("# Initializing data\n");
	timer ttimer = create_timer("Total_time");

	//======= Allocate matrices ===============
	double (*matrix)[nblocks][BSIZE][BSIZE] = NULL;
	matrix = nanos6_dmalloc(dim2, nanos6_equpart_distribution, 0, NULL);
	assert(matrix);

	init_matrix(nblocks, BSIZE, matrix);

	double (*original)[DIM] = NULL;
	double (*factorized)[DIM] = NULL;

	//========= Init matrix ===================
	#pragma oss taskwait

	if (CHECK) {
		original = nanos6_dmalloc(dim2, nanos6_equpart_distribution, 0, NULL);
		assert(original);

		factorized = nanos6_dmalloc(dim2, nanos6_equpart_distribution, 0, NULL);
		assert(factorized);

		task_blocked2flat(nblocks, BSIZE, DIM, matrix, original);

		write_matrix_block("orig_block.txt", nblocks, BSIZE, matrix);
		write_matrix_flat("orig_flat.txt", DIM, original);
	}

	// ===========================================
	printf("# Starting algorithm in process\n");
	timer atimer = create_timer("Algorithm_time");

	cholesky(nblocks, BSIZE, matrix);
	#pragma oss taskwait

	stop_timer(&atimer);

	//======== Check if set =====================
	if (CHECK) {
		// Allocate new matrix to transform it from tiled to flat
		// Transform from tile (matrix) to flat (factorized_matrix)
		task_blocked2flat(nblocks, BSIZE, DIM, matrix, factorized);

		printf("# Writting factorized\n");
		write_matrix_block("fact_block.txt", nblocks, BSIZE, matrix);
		write_matrix_flat("fact_flat.txt", DIM, factorized);

		printf("# Checking the correctness of the factorization...\n");
		const double EPS = BLAS_dfpinfo(blas_eps);
		check_factorization(DIM, original, factorized, DIM, EPS);

		#pragma oss taskwait
		nanos6_dfree(factorized, dim2);
		nanos6_dfree(original, dim2);
	}

	stop_timer(&ttimer);
	nanos6_dfree(matrix, dim2);

	report_args();
	free_args();

	return 0;
}

