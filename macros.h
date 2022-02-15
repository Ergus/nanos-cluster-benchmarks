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

#ifndef benchmarks_h
#define benchmarks_h

#ifdef __cplusplus
extern "C" {
#endif

// This before ANY include
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

// C headers
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <sched.h>

#ifdef _OPENMP
#include <omp.h>
#else // _OPENMP

#ifndef _OMPSS_2
#warning "Compiling without OpenMP or OmpSs"
#endif // _OMPSS_2

#define omp_get_thread_num() 0
#define omp_set_num_threads(var) {}
#define omp_set_dynamic(NUM)
#define omp_get_max_threads() 1

#define omp_set_schedule(...)
#define omp_sched_static
#endif // _OPENMP

#define imin(x, y) (((x) < (y)) ? (x) : (y))
#define imax(x, y) (((x) > (y)) ? (x) : (y))

#define frand()(4.*(double)rand()/(RAND_MAX)-2.) //uniform rng in [-2,2]

#define printme() {										   \
		fprintf(stderr,"Func: %s in %s:%d (process %s)\n", \
		        __PRETTY_FUNCTION__, __FILE__, __LINE__,   \
		        getenv("OMPI_COMM_WORLD_RANK"));		   \
	}

#define modcheck(a, b){	  \
		const int themod=(a) % (b); \
		if(themod) { \
			fprintf(stderr,"Error: %s %% %s = %d\n", #a, #b, themod); \
			printme(); \
			exit(EXIT_FAILURE); \
		} \
	}

#define myassert(cond) {										\
		if (!(cond)) {											\
			fprintf(stderr, "%s%s:%u Assertion `%s' failed.\n", \
			        __func__, __FILE__, __LINE__, #cond);		\
			exit(EXIT_FAILURE);									\
		}														\
	}

#ifndef NDEBUG
#define dbprintf(...) fprintf(stderr, __VA_ARGS__)
#define dbprintf_if(cond, ...) if (cond) fprintf(stderr, __VA_ARGS__)
#else
#define dbprintf(...)
#define dbprintf_if(...)
#endif

	inline FILE *get_file(const char prefix[64], const char name[64],
	               const char *restrict mode)
	{
		FILE *fp = stdout;
		if (prefix != NULL) {
			char filename[256];
			sprintf(filename,"%s_%s.mat", prefix, name);
			fp = fopen(filename, mode);
			myassert(fp);
		}
		return fp;
	}

	inline void print_matrix_header(FILE *fp, const char name[64],
	                                const size_t rows, const size_t cols
	) {
		fprintf(fp, "# name: %s\n", name);
		fprintf(fp, "# type: matrix\n");
		fprintf(fp, "# rows: %lu\n", rows);
		fprintf(fp, "# columns: %lu\n", cols);
	}

	inline void print_matrix_data(FILE *fp, const double * const mat,
	                              const size_t rows, const size_t cols
	) {
		for (size_t i = 0; i < rows; ++i) {
			for(size_t j = 0; j < cols; ++j) {
				fprintf(fp, "%3.8lf ", mat[i * cols + j]);
			}
			fprintf(fp,"\n");
		}

	}

	static inline void __print(const double * const mat,
	                           const size_t rows, const size_t cols,
	                           const char prefix[64], const char name[64]
	) {
		FILE *fp = get_file(prefix, name, "w+");
		print_matrix_header(fp, name, rows, cols);
		print_matrix_data(fp, mat, rows, cols);

		if (fp != stdout) {
			fclose(fp);
		}
	}

#define printmatrix(mat, rows, cols, prefix)		\
	__print(mat, rows, cols, prefix, #mat)

	static inline void block_init(double * const __restrict__ array,
	                              const size_t rows, const size_t cols, int seed
	) {
		const size_t fullsize = rows * cols;

		struct drand48_data drand_buf;
		srand48_r(seed == 0 ? 1 : seed, &drand_buf);
		double x;

		for(size_t i = 0; i < fullsize; ++i) {
			drand48_r(&drand_buf, &x);
			array[i] = x;
		}
	}

	static inline int count_sched_cpus()
	{
		cpu_set_t mask;
		if (sched_getaffinity(0, sizeof(mask), &mask) == 0)
			return CPU_COUNT(&mask);
		return -1;
	}

	static inline int copy_file(const char in[], const char out[])
	{
		FILE *fin = fopen(in, "r"),
			*fout = (out == NULL) ? stdout : fopen(out, "w");

		if (fin == NULL || fout == NULL)
		{
			fprintf(stderr, "Error coping file in=%p out=%p\n",
			        fin, fout);
			return -1;
		}

		char ch;
		while( ( ch = fgetc(fin) ) != EOF ) {
			fputc(ch, fout);
		}

		fclose(fin);

		if (fout != NULL && fout != stdout) {
			fclose(fout);
		}

		return 0;
	}

#ifdef __cplusplus
}
#endif

#endif // benchmarks_h
