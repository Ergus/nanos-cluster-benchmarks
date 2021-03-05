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


//############## Blas wrappers #####################

void oss_potrf(const size_t bsize, double A[bsize][bsize])
{
	#pragma oss task inout(A[0;bsize][0;bsize])			\
		node(nanos6_cluster_no_offload) label("potrf")
	{
		int error = LAPACKE_dpotrf(
			LAPACK_COL_MAJOR,
			'L',
			bsize,
			&A[0][0],
			bsize
		);

		assert(error == 0);
	}
}

void oss_trsm(const size_t bsize, double A[bsize][bsize], double B[bsize][bsize])
{
	#pragma oss task in(A[0;bsize][0;bsize])			\
		inout(B[0;bsize][0;bsize])						\
		node(nanos6_cluster_no_offload) label("trsm")
	cblas_dtrsm(
		CblasColMajor,
		CblasRight,
		CblasLower,
		CblasTrans,
		CblasNonUnit,
		bsize, bsize, 1.0,
		&A[0][0], bsize,
		&B[0][0], bsize
	);
}

void oss_syrk(const size_t bsize, double A[bsize][bsize], double B[bsize][bsize])
{
	#pragma oss task in(A[0;bsize][0;bsize])			\
		inout(B[0;bsize][0;bsize])						\
		node(nanos6_cluster_no_offload) label("syrk")
	cblas_dsyrk(
		CblasColMajor,
		CblasLower,
		CblasNoTrans,
		bsize, bsize, -1.0,
		&A[0][0], bsize,
		1.0,
		&B[0][0], bsize
	);
}

void oss_gemm(const size_t bsize,
              double A[bsize][bsize],
              double B[bsize][bsize],
              double C[bsize][bsize]
) {
	#pragma oss task in(A[0;bsize][0;bsize])			\
		in(B[0;bsize][0;bsize])							\
		inout(C[0;bsize][0;bsize])						\
		node(nanos6_cluster_no_offload) label("gemm")
	cblas_dgemm(
		CblasColMajor,
		CblasNoTrans,
		CblasTrans,
		bsize, bsize, bsize, -1.0,
		&A[0][0], bsize,
		&B[0][0], bsize,
		1.0,
		&C[0][0], bsize
	);
}

//######## Check functions ########

void check_factorization(const size_t n,
                         double A1[n][n], double A2[n][n],
                         const size_t lda, const double eps
) {
	#pragma oss task in(A1[0;n][0;n]) in(A2[0;n][0;n]) label("check_factorization")
	{
		printf("# ===================================\n");
		printf("# Checking the Cholesky Factorization\n");

		const size_t len = n * n * sizeof(double);
		double *Residual = (double *) malloc(len);
		double *L1       = (double *) malloc(len);
		double *L2       = (double *) malloc(len);
		memset(L1, 0, len);
		memset(L2, 0, len);

		const double alpha = 1.0;
		LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', n, n, &A1[0][0], lda, Residual, n);
		LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'L', n, n, &A2[0][0], lda, L1, n);
		LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'L', n, n, &A2[0][0], lda, L2, n);

		cblas_dtrmm(
			CblasColMajor,
			CblasRight,
			CblasLower,
			CblasTrans,
			CblasNonUnit,
			n, n, alpha, L1, n, L2, n
		);

		// Compute the Residual A - L'*L
		// daxpy: y := a*x + y  (with alpha = -1, x:= L'*L
		cblas_daxpy( n*n, -1.0, L2, 1, Residual, 1 );

		const double Rnorm = LAPACKE_dlange(LAPACK_COL_MAJOR, 'I', n, n, Residual, n);
		const double Anorm = LAPACKE_dlange(LAPACK_COL_MAJOR, 'I', n, n, &A1[0][0], n);

		const int ok = isnan(Rnorm/(Anorm*n*eps)) || (Rnorm/(Anorm*n*eps) > 60.0);

		printf("# ||L'L-A||_oo/(||A||_oo.n.eps) = %e \n", Rnorm / (Anorm * n * eps));
		printf("# Factorization is %s... ! \n", (ok ? "SUSPICIOUS" : "CORRECT"));
		printf("# ===================================\n");

		free(Residual);
		free(L1);
		free(L2);
	}
}

void write_matrix_flat(const char filename[64], size_t ld, double matrix[ld][ld])
{
	#pragma oss task in(matrix[0;ld][0;ld]) label("write")
	{
		printf("# Writing matrix in file %s\n", filename);
		FILE *fd = fopen(filename, "w");
		assert(fd);

		for (size_t i = 0; i < ld; ++i) {
			for(size_t j = 0; j < ld; ++j){
				fprintf(fd,"%lf ", matrix[i][j]);
			}
			fprintf(fd,"\n");
		}
		fclose(fd);
	}
}

void write_matrix_block(const char filename[64],
                        const size_t nblocks, const size_t bsize,
                        double matrix[nblocks][nblocks][bsize][bsize]
) {
	#pragma oss task in(matrix[0;nblocks][0;nblocks][0;bsize][0;bsize]) label("write_blocked")
	{
		printf("# Writing blocked matrix in file %s\n", filename);
		FILE *fd = fopen(filename, "w");
		assert(fd);

		for (size_t i = 0; i < nblocks; ++i) {          // block row
			for (size_t k = 0; k < bsize; ++k) {        // line row within block
				for (size_t j = 0; j < nblocks; ++j) {  // block column
					for (size_t l = 0; l < bsize; ++l) {
						fprintf(fd,"%lf ", matrix[i][j][k][l]);
					}
				}
				fprintf(fd,"\n");
			}
		}
		fclose(fd);
	}
}

//! This is to initialize non-diagonal blocks (assuming 1 block/task)
void task_fill_block_ij(const size_t nblocks, const size_t bsize,
                        double matrix[nblocks][nblocks][bsize][bsize],
                        const size_t i, const size_t j
) {
	assert(i != j);

	#pragma oss task out(matrix[i][j][0;bsize][0;bsize]) label("init_block_ij")
	{
		const size_t seed = (i > j) ? i * nblocks + j : j * nblocks + i;
		struct drand48_data status;       	// using re-entrant version for rng
		srand48_r(seed, &status);
		double rnd1, rnd2;

		for (size_t k = 0; k < bsize; ++k) {
			for (size_t l = 0; l < bsize; ++l) {
				drand48_r(&status, &rnd1);
				drand48_r(&status, &rnd2);

				if (i > j) {
					matrix[i][j][k][l] = rnd1 * rnd2;
				} else {
					matrix[i][j][l][k] = rnd1 * rnd2;
				}
			}
		}
	}
}

//! This is to initialize diagonal blocks (assuming 1 block/task)
void task_fill_block_ii(const size_t nblocks, const size_t bsize,
                        double matrix[nblocks][nblocks][bsize][bsize],
                        const size_t i
) {
	#pragma oss task out(matrix[i][i][0;bsize][0;bsize]) label("init_block_ii")
	{
		const size_t seed = i * nblocks + i, ld = nblocks * bsize;
		double rnd1, rnd2;
		struct drand48_data status;     	// using re-entrant version for rng

		srand48_r(seed, &status);

		for (size_t k = 0; k < bsize; ++k) {
			drand48_r(&status, &rnd1);
			drand48_r(&status, &rnd2);
			matrix[i][i][k][k] = rnd1 * rnd2 + ld;

			for (size_t l = k+1; l < bsize; ++l) {
				drand48_r(&status, &rnd1);
				drand48_r(&status, &rnd2);

				const double val = rnd1 * rnd2;
				matrix[i][i][k][l] = val;
				matrix[i][i][l][k] = val;
			}
		}
	}
}

//! Convert flat matrix to blocked
void task_flat2blocked(size_t nblocks, size_t bsize,
                       double flat[nblocks*bsize][nblocks*bsize],
                       double blocked[nblocks][nblocks][bsize][bsize]
) {
	#pragma oss task weakin(flat[0;nblocks*bsize][0;nblocks*bsize]) \
		weakout(blocked[0;nblocks][0;nblocks][0;bsize][0;bsize]) label("flat2blocked")
	{
		const size_t ld = nblocks * bsize;

		for (size_t i = 0; i < nblocks; ++i) {
			#pragma oss task weakin(flat[i*bsize; bsize][0;ld]) \
				weakout(blocked[i][0; nblocks][0; bsize][0; bsize])
			{
				for (size_t j = 0; j < nblocks; ++j) {
					#pragma oss task in(flat[i*bsize; bsize][j*bsize; bsize])	  \
						out(blocked[i][j][0; bsize][0; bsize])
					{
						for (size_t ii = 0; ii < bsize; ++ii)
							for (size_t jj = 0; jj < bsize; ++jj)
								blocked[i][j][ii][jj] = flat[i * bsize + ii][j * bsize + jj];
					}
				}
			}
		}
	}
}

//! Convert blocked matrix to flat
void task_blocked2flat(size_t nblocks, size_t bsize, size_t ld,
                       double blocked[nblocks][nblocks][bsize][bsize],
                       double flat[nblocks*bsize][nblocks*bsize]
) {
	#pragma oss task weakin(blocked[0; nblocks][0; nblocks][0; bsize][0; bsize]) \
		weakout(flat[0; ld][0; ld]) label("blocked2flat")
	{
		for (size_t i = 0; i < nblocks; ++i) {
			#pragma oss task weakin(blocked[i][0;nblocks][0;bsize][0;bsize]) \
				weakout(flat[i*bsize; bsize][0;ld]) label("blocked2flat_slice")
			{
				for (size_t j = 0; j < nblocks; ++j) {
					#pragma oss task in(blocked[i][j][0;bsize][0;bsize]) \
						out(flat[i*bsize; bsize][j*bsize; bsize]) label("blocked2flat_block")
					{
						for (size_t ii = 0; ii < bsize; ++ii)
							for (size_t jj = 0; jj < bsize; ++jj)
								flat[i * bsize + ii][j * bsize + jj] = blocked[i][j][ii][jj];
					}
				}
			}
		}
	}
}

//! Copy blocked matrix (respecting the blocked possible distribution)
void task_copy_matrix_blocked(size_t nblocks, size_t bsize,
                              double mbin[nblocks][nblocks][bsize][bsize],
                              double mbout[nblocks][nblocks][bsize][bsize]
) {
	#pragma oss task weakin(mbin[0;nblocks][0;nblocks][0;bsize][0;bsize]) \
		weakout(mbout[0;nblocks][0;nblocks][0;bsize][0;bsize]) label("copy_blocked")
	{
		const size_t bsize2 = bsize * bsize;

		for (size_t i = 0; i < nblocks; ++i) {
			#pragma oss task weakin(mbin[i][0;nblocks][0;bsize][0;bsize]) \
				weakout(mbout[i][0;nblocks][0;bsize][0;bsize]) label("copy_blocked_slice")
			{
				for (size_t j = 0; j < nblocks; ++j) {
					#pragma oss task in(mbin[i][j][0;bsize][0;bsize])	  \
						out(mbout[i][j][0;bsize][0;bsize]) label("copy_blocked_block")
					{
						memcpy(&mbout[i][j][0][0], &mbin[i][j][0][0], bsize2 * sizeof(double));
					}
				}
			}
		}
	}
}

void init_matrix(const size_t nblocks, const size_t bsize,
                 double matrix[nblocks][nblocks][bsize][bsize])
{
	for (size_t i = 0; i < nblocks; ++i) {
		task_fill_block_ii(nblocks, bsize, matrix, i);

		for (size_t j = i + 1; j < nblocks; ++j) {
			task_fill_block_ij(nblocks, bsize, matrix, i, j);
			task_fill_block_ij(nblocks, bsize, matrix, j, i);
		}
	}
}
