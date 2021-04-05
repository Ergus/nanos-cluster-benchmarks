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

#include <nanos6/debug.h>

#include "cmacros/macros.h"

	void init_AB(double *A, double *B, const size_t dim, size_t ts)
	{
		const size_t numNodes = nanos6_get_num_cluster_nodes();
		myassert(dim >= ts);
		myassert(dim / ts >= numNodes);
		modcheck(dim, ts);
		modcheck(dim / numNodes, ts);

		const size_t rowsPerNode = dim / numNodes;

		for (size_t i = 0; i < dim; i += rowsPerNode) { // loop nodes

			const int nodeid = i / rowsPerNode;

			#pragma oss task weakout(A[i * dim; rowsPerNode * dim])	\
				weakout(B[i; rowsPerNode])							\
				node(nodeid) wait label("init_AB_weak")
			{
				for (size_t j = i; j < i + rowsPerNode; j += ts) { // loop tasks

					#pragma oss task out(A[j * dim; ts * dim])			\
						out(B[j; ts])									\
						node(nanos6_cluster_no_offload) label("init_AB_strong")
					{
						struct drand48_data drand_buf;
						srand48_r(j, &drand_buf);
						double x;

						for (size_t k = j; k < j + ts; ++k) {
							double cum = 0.0, sum = 0.0;
							for (size_t l = 0; l < dim; ++l) {
								drand48_r(&drand_buf, &x);
								A[k * dim + l] = x;
								cum += fabs(x);
								sum += x;
							}
							// Diagonal element condition.
							const double valkk = A[k * dim + k];
							if (signbit(valkk)) {
								A[k * dim + k] = valkk - cum;
								B[k] = sum - cum;
							} else {
								A[k * dim + k] = valkk + cum;
								B[k] = sum + cum;
							}
						}
					}
				}
			}
		}
	}

	void init_x(double *x, const size_t dim, size_t ts, double val)
	{
		const size_t numNodes = nanos6_get_num_cluster_nodes();
		myassert(dim >= ts);
		modcheck(dim, ts);

		const size_t rowsPerNode = dim / numNodes;

		for (size_t i = 0; i < dim; i += rowsPerNode) { // loop nodes

			int nodeid = i / rowsPerNode;

			#pragma oss task weakout(x[i; rowsPerNode])		\
				node(nodeid) wait label("init_vec_weak")
			{
				for (size_t j = i; j < i + rowsPerNode; j += ts) { // loop tasks

					#pragma oss task out(x[j; ts])						\
						node(nanos6_cluster_no_offload) label("init_vec_strong")
					{
						for (size_t k = j; k < j + ts; ++k) {
							x[k] = val;
						}
					}
				}
			}
		}
	}

	void jacobi_modify(double *A, double *b, size_t dim, size_t ts)
	{
		const size_t numNodes = nanos6_get_num_cluster_nodes();
		myassert(dim >= ts);
		modcheck(dim, ts);

		const size_t rowsPerNode = dim / numNodes;

		for (size_t i = 0; i < dim; i += rowsPerNode) { // loop nodes

			int nodeid = i / rowsPerNode;

			#pragma oss task weakinout(A[i * dim; rowsPerNode * dim])	\
				weakinout(b[i; rowsPerNode])							\
				node(nodeid) wait label("jacobi_modify_weak")
			{
				for (size_t j = i; j < i + rowsPerNode; j += ts) { // loop tasks

					#pragma oss task inout(A[j * dim; ts * dim])		\
						inout(b[j; ts])									\
						node(nanos6_cluster_no_offload) label("jacobi_modify_strong")
					{
						for (size_t k = j; k < j + ts; ++k) {
							const double Akk = fabs(A[k * dim + k]);

							for (size_t l = 0; l < dim; ++l) {
								if (l == k) {
									A[k * dim + l] = 0;
								} else {
									A[k * dim + l] = - (A[k * dim + l] /  Akk);
								}
							}
							b[k] /= Akk;
						}
					}
				}
			}
		}
	}

	void free_matrix(double *mat, size_t size)
	{
		nanos6_dfree(mat, size * sizeof(double));
	}


	void matmul_base(const double *A, const double *B, double * const C,
	                 size_t lrowsA, size_t dim)
	{
		for (size_t i = 0; i < lrowsA; ++i) {
			C[i] = 0.0;

			for (size_t j = 0; j < dim; ++j) {
				const double temp = A[i * dim + j];

				C[i] += (temp * B[j]);
			}
		}
	}

	void __print_task(const double * const mat,
	                  const size_t rows, const size_t cols,
	                  const char prefix[64], const char name[64])
	{
		#pragma oss task in(mat[0; rows * cols]) label("matrix_print")
		{
			__print(mat, rows, cols, prefix, name);
		}
	}

#define printmatrix_task(mat, rows, cols, prefix) \
	__print_task(mat, rows, cols, prefix, #mat)

#ifdef __cplusplus
}
#endif

#endif
