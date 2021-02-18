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

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_set_num_threads(var) {}
#define omp_set_dynamic(NUM)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <sched.h>

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
#define dprintf(...) fprintf(stderr,__VA_ARGS__)
#else
#define dprintf(...) {}
#endif


	static inline void __print(const double * const mat,
	                           const size_t rows, const size_t cols,
	                           const char prefix[64], const char name[64])
	{
		char filename[256];
		sprintf(filename,"%s_%s.mat", prefix, name);
		FILE *fp = fopen(filename, "w+");
		myassert(fp);

		fprintf(fp, "# name: %s\n", name);
		fprintf(fp, "# type: matrix\n");
		fprintf(fp, "# rows: %lu\n", rows);
		fprintf(fp, "# columns: %lu\n", cols);

		for (size_t i = 0; i < rows; ++i) {
			for(size_t j = 0; j < cols; ++j) {
				fprintf(fp, "%3.8lf ", mat[i * cols + j]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
	}

#define printmatrix(mat, rows, cols, prefix) \
	__print(mat, rows, cols, prefix, #mat)


	static inline void matrix_init(double * const __restrict__ array,
	                               const size_t rows, const size_t cols, int seed)
	{
		const size_t fullsize = rows * cols;

#ifdef _OPENMP
		#pragma omp parallel
#endif
		{
			size_t i;

			struct drand48_data drand_buf;
			srand48_r(seed + omp_get_thread_num(), &drand_buf);
			double x;

#ifdef _OPENMP
			#pragma omp for
#endif
			for(i = 0; i < fullsize; ++i) {
				drand48_r(&drand_buf, &x);
				array[i] = x;
			}
		}
	}

	static inline int count_sched_cpus()
	{
		cpu_set_t mask;
		if (sched_getaffinity(0, sizeof(mask), &mask) == 0)
			return CPU_COUNT(&mask);
		return -1;
	}

#ifdef __cplusplus
}
#endif

#endif // benchmarks_h
