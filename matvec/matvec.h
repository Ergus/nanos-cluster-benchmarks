#ifndef MATVEC_H
#define MATVEC_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef NANOS6
#include <nanos6/debug.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define min(x,y) ((x) < (y)) ? (x) : (y)
#define max(x,y) ((x) > (y)) ? (x) : (y)

#define myassert(cond) {						\
		if (!cond) {						\
			fprintf(stderr, "%s%s:%u Assertion `%s' failed.\n", \
			        __func__, __FILE__, __LINE__, #cond);	\
			abort();					\
		}							\
	}

double *alloc_init(const size_t rows, size_t cols, size_t TS)
{
	myassert(rows > TS == 0);
	myassert(rows % TS == 0);

	const size_t size = cols * TS;

	double *ret = (double *) nanos6_dmalloc(size * sizeof(double),
	                                        nanos6_equpart_distribution, 0, nullptr);
	myassert (ret);

	for (size_t i = 0; i < rows; i += TS) {
		double *first = &ret[i * cols];

		#pragma oss task out(first[0; size]) label(initalize_vector)
		{
			struct drand48_data drand_buf;
			srand48_r(i, &drand_buf);
			double x;

			for (size_t j = 0; j < size; ++j) {
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

bool validate(const double *A, const double *b, double *x,
               size_t rows, size_t cols)
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
