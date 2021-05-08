/*
 * Copyright (C) 2019  Jimmy Aguilar Mena
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

#ifndef MATVEC_H
#define MATVEC_H

#ifdef __cplusplus
extern "C" {
#endif

#include "benchmarks_ompss.h"

	void free_matrix(double *mat, size_t size)
	{
		nanos6_dfree(mat, size * sizeof(double));
	}


	void matmul_base(const double *A, const double *B, double * const C,
	                 size_t lrowsA, size_t dim, size_t colsBC)
	{
		for (size_t i = 0; i < lrowsA; ++i) {
			for (size_t k = 0; k < colsBC; ++k)
				C[i * colsBC + k] = 0.0;

			for (size_t j = 0; j < dim; ++j) {
				const double temp = A[i * dim + j];

				for (size_t k = 0; k < colsBC; ++k) {
					C[i * colsBC + k] += (temp * B[j * colsBC + k]);
				}
			}
		}
	}

	bool validate(const double *A, const double *B, double *C,
	              size_t dim, size_t colsBC)
	{
		bool success = false;

		#pragma oss task in(A[0; dim * dim])			\
			in(B[0; dim * colsBC])						\
			in(C[0; dim * colsBC])						\
			inout(success) label("validate")
		{
			success = true;

			double *expected = (double *) malloc(dim * colsBC * sizeof(double));
			matmul_base(A, B, expected, dim, dim, colsBC);

			for (size_t i = 0; i < dim; ++i) {
				for (size_t j = 0; j < colsBC; ++j) {

					if (abs(C[i * colsBC + j] - expected[i * colsBC + j]) > 1e-12) {
						success = false;
						break;
					}
				}
			}

			free(expected);
		}
		#pragma oss taskwait
		return success;
	}

#ifdef __cplusplus
}
#endif

#endif
