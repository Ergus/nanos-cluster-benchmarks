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


#ifndef CHOLESKY_UTILS_H
#define CHOLESKY_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>

#include <unistd.h>
#include <fcntl.h>
#include <assert.h>
#include <math.h>

#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "cmacros/macros.h"
#include "mkl.h"

#include <nanos6.h>
#include <nanos6/debug.h>

//############## Blas wrappers #####################

void oss_potrf(const size_t bsize, double A[bsize][bsize], size_t i);

void oss_trsm(const size_t bsize,
              double A[bsize][bsize],
              double B[bsize][bsize]);

void oss_syrk(const size_t bsize, double A[bsize][bsize], double B[bsize][bsize]);

void oss_gemm(const size_t bsize, double A[bsize][bsize],
              double B[bsize][bsize], double C[bsize][bsize]);

//######## Check functions ########

void check_factorization(const size_t n,
                         double A1[n][n], double A2[n][n],
						 const size_t lda, const double eps);

void write_matrix_flat(const char filename[64], size_t ld, double matrix[ld][ld]);

void write_matrix_block(const char filename[64], const size_t nblocks, const size_t bsize,
                          double matrix[nblocks][nblocks][bsize][bsize]);

// ####### Fill Blocks ############
void task_fill_block_ij(const size_t nblocks, const size_t bsize,
                        double matrix[nblocks][nblocks][bsize][bsize],
                        const size_t i, const size_t j);

void task_fill_block_ii(const size_t nblocks, const size_t bsize,
                        double matrix[nblocks][nblocks][bsize][bsize],
                        const size_t i);

// ########### Reshape #############
void task_flat2blocked(size_t nblocks, size_t bsize,
                  double flat[nblocks*bsize][nblocks*bsize],
                  double blocked[nblocks][nblocks][bsize][bsize]);

void task_blocked2flat(size_t nblocks, size_t bsize, size_t ld,
                       double blocked[nblocks][nblocks][bsize][bsize],
                       double flat[nblocks*bsize][nblocks*bsize]);

void task_copy_matrix_blocked(size_t nblocks, size_t bsize,
                         double in[nblocks][nblocks][bsize][bsize],
                         double out[nblocks][nblocks][bsize][bsize]);

void init_matrix(const size_t nblocks, const size_t bsize,
                 double matrix[nblocks][nblocks][bsize][bsize]);

//############## Blass wrappers #####################

enum blas_cmach_type {
	blas_base      = 151,
	blas_t         = 152,
	blas_rnd       = 153,
	blas_ieee      = 154,
	blas_emin      = 155,
	blas_emax      = 156,
	blas_eps       = 157,
	blas_prec      = 158,
	blas_underflow = 159,
	blas_overflow  = 160,
	blas_sfmin     = 161
};

static inline double BLAS_dpow_di(double x, int n)
{
	double rv = 1.0;

	if (n < 0) {
		n = -n;
		x = 1.0 / x;
	}

	for (; n; n >>= 1, x *= x) {
		if (n & 1) rv *= x;
	}

	return rv;
}

static inline double BLAS_dfpinfo(enum blas_cmach_type cmach)
{
	const double b = 2.0;
	const int    t = 53, m = -1021; // l = 1024,
	const double eps = BLAS_dpow_di(b, -t);
	const double r   = BLAS_dpow_di(b, m-1);

	//double o = ((1.0 - eps) * BLAS_dpow_di(b, l-1)) * b;

	switch (cmach) {
	case blas_eps:   return eps;
	case blas_sfmin: return r;
	default:
		fprintf(stderr, "%s %d %d %d\n", __func__, -1, cmach, 0);
		abort();
		break;
	}

	return 0.0;
}

#endif // CHOLESKY_UTILS_H
