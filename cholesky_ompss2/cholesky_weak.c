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

void cholesky(const size_t nblocks, const size_t bsize,
              double A[nblocks][nblocks][bsize][bsize]
) {
	for (size_t i = 0; i < nblocks; ++i) {
		oss_potrf(bsize, A[i][i]);    // Diagonal Block Factorization

		#pragma oss task weakin(A[i][i][0;bsize][0;bsize])		\
			weakinout(A[i][i+1:nblocks-1][0;bsize][0;bsize])
		{
			for (size_t j = i + 1; j < nblocks; ++j)      // Triangular Systems
				oss_trsm(bsize, A[i][i], A[i][j]);
		}

		#pragma oss task weakin(A[i][i+1:nblocks-1][0;bsize][0;bsize])	\
			weakinout(A[i+1:nblocks-1][i+1:nblocks-1][0;bsize][0;bsize])
		{
			for (size_t j = i + 1; j < nblocks; ++j) {    // Update Trailing Matrix
				oss_syrk(bsize, A[i][j], A[j][j]);
				#pragma oss task weakin(A[i][i+1:j][0;bsize][0;bsize])	\
					weakinout(A[i+1:j-1][i+1:nblocks-1][0;bsize][0;bsize])
				{
					for (size_t k = i + 1; k < j; ++k)
						oss_gemm(bsize, A[i][j], A[i][k], A[k][j]);
				}
			}
		}
	}
}

int main(int argc, char **argv)
{
	init_args(argc, argv);

	const int dim = create_cl_int ("Dimension");
	const int bsize = create_cl_int ("Block_size");
	const int check = create_optional_cl_int ("Check", 0);

	myassert(dim >= bsize);
	modcheck(dim, bsize);

	//========= End Command Line ===============

	const size_t nblocks = dim / bsize;          // Number of blocks
	const size_t dim2 = dim * dim * sizeof(double);
	printf("# Initializing data\n");
	timer ttimer = create_timer("Total time");

	//======= Allocate matrices ===============
	double (*matrix)[nblocks][bsize][bsize] = NULL;
	matrix = nanos6_dmalloc(dim2, nanos6_equpart_distribution, 0, NULL);
	assert(matrix);

	init_matrix(nblocks, bsize, matrix);

	double (*original)[dim] = NULL;
	double (*factorized)[dim] = NULL;

	//========= Init matrix ===================
	#pragma oss taskwait

	if (check) {
		original = nanos6_dmalloc(dim2, nanos6_equpart_distribution, 0, NULL);
		assert(original);

		factorized = nanos6_dmalloc(dim2, nanos6_equpart_distribution, 0, NULL);
		assert(factorized);

		task_blocked2flat(nblocks, bsize, dim, matrix, original);

		write_matrix_block("orig_block.txt", nblocks, bsize, matrix);
		write_matrix_flat("orig_flat.txt", dim, original);
	}

	// Algorithm
	printf("# Executing the factorization...\n");
	timer atimer = create_timer("Algorithm time");

	cholesky(nblocks, bsize, matrix);
	#pragma oss taskwait

	stop_timer(&atimer);

	const double performance = dim * dim * dim * 3000.0 / getNS_timer(&atimer);

	create_reportable_double("performance", performance);

	//======== Check if set =====================
	if (check) {
		// Allocate new matrix to transform it from tiled to flat
		// Transform from tile (matrix) to flat (factorized_matrix)
		task_blocked2flat(nblocks, bsize, dim, matrix, factorized);

		printf("# Writting factorized\n");
		write_matrix_block("fact_block.txt", nblocks, bsize, matrix);
		write_matrix_flat("fact_flat.txt", dim, factorized);

		printf("# Checking the correctness of the factorization...\n");
		const double EPS = BLAS_dfpinfo(blas_eps);
		check_factorization(dim, original, factorized, dim, EPS);

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

