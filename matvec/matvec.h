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

#ifdef NANOS6
#include <nanos6/debug.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#define min(x,y) ((x) < (y)) ? (x) : (y)
#define max(x,y) ((x) > (y)) ? (x) : (y)

#define myassert(cond) {										\
		if (!(cond)) {											\
			fprintf(stderr, "%s%s:%u Assertion `%s' failed.\n", \
			        __func__, __FILE__, __LINE__, #cond);		\
			abort();											\
		}														\
	}

double *alloc_init(const size_t rows, size_t cols, size_t TS)
{
	myassert(rows >= TS);
	myassert(rows % TS == 0);

	const size_t size = cols * rows;

	double *ret = (double *) nanos6_dmalloc(size * sizeof(double),
	                                        nanos6_equpart_distribution, 0, nullptr);
	myassert (ret);

	const size_t piece = cols * TS;

	for (size_t i = 0; i < rows; i += TS) {
		double *first = &ret[i * cols];

		#pragma oss task out(first[0; piece]) label(initalize_vector)
		{
			struct drand48_data drand_buf;
			srand48_r(i, &drand_buf);
			double x;

			for (size_t j = 0; j < piece; ++j) {
				drand48_r(&drand_buf, &x);
				first[j] = x;
			}
		}
	}

	return ret;
}

void free_matrix(double *mat, size_t size)
{
	nanos6_dfree(mat, size * sizeof(double));
}

// Multiply 
void matvec_base(const double *A, const double *x, double *b,
                 size_t rows, size_t cols)
{
	for (size_t i = 0; i < rows; ++i) {
		b[i] = 0.0;
		for (size_t j = 0; j < cols; ++j)
			b[i] += A[i * cols + j] * x[j];
	}
}

void matvec_print2d(const double *A, size_t rows, size_t cols,
                    const char filename[128])
{
	#pragma oss task in(A[0; rows * cols]) label(matrix_print)
	{
		FILE *fp = fopen(filename, "w+");
		myassert(fp);

		for (size_t i = 0; i < rows; ++i)
			for (size_t j = 0; j < cols; ++j)
				fprintf(fp, "%.3lf%s",
				        A[i * cols + j],
				        (j == cols - 1) ? "\n" : " ");

		fclose(fp);
	}
}

void matvec_print1d(const double *vec, size_t size, const char filename[128])
{
	#pragma oss task in(vec[0; size]) label(print_vector)
	{
		FILE *fp = fopen(filename, "w+");
		myassert(fp);

		for (size_t i = 0; i < size; ++i)
			fprintf(fp, "%.3lf%s", vec[i],
			        (i < size - 1) ? " " : "\n");
		fclose(fp);
	}
}

bool validate(const double *A, const double *b, double *x, size_t rows, size_t cols)
{
	bool success = true;

	#pragma oss task in(A[0; rows * cols]) in(b[0; cols]) in(x[0; rows])\
		inout(success) label(validate)
	{
		success = true;

		double *expected = (double *) malloc(rows * sizeof(double));
		matvec_base(A, b, expected, rows, cols);

		for (size_t i = 0; i < rows; ++i) {
			if (abs(x[i] - expected[i]) > 1e-12) {
				success = false;
				break;
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
