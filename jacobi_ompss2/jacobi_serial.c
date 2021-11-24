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

#include "jacobi_ompss2.h"

extern void jacobi_base(
	const double * __restrict__ A,
	double Bi,
	const double * __restrict__ xin,
	double * __restrict__ xouti, size_t dim
);


void init_AB_serial(double *A, double *B, const size_t dim)
{
	for (size_t i = 0; i < dim; ++i) { // loop nodes

		struct drand48_data drand_buf;
		srand48_r(i, &drand_buf);
		double x;

		double cum = 0.0, sum = 0.0;
		for (size_t l = 0; l < dim; ++l) {
			drand48_r(&drand_buf, &x);
			A[i * dim + l] = x;
			cum += fabs(x);
			sum += x;
		}
		// Diagonal element condition.
		const double valii = A[i * dim + i];
		if (signbit(valii)) {
			A[i * dim + i] = valii - cum;
			B[i] = sum - cum;
		} else {
			A[i * dim + i] = valii + cum;
			B[i] = sum + cum;
		}
	}
}


void init_x_serial(double *x, const size_t dim, double val)
{
	for (size_t i = 0; i < dim; ++i) { // loop nodes
		x[i] = val;
	}
}


void jacobi_modify_serial(double *A, double *b, size_t dim)
{
	for (size_t i = 0; i < dim; ++i) { // loop nodes

		const double Aii = fabs(A[i * dim + i]);

		for (size_t l = 0; l < dim; ++l) {
			A[i * dim + l] = (l != i) ? (- A[i * dim + l] /  Aii) : 0;
		}
		b[i] /= Aii;
	}
}


void jacobi_tasks_serial(const double *A, const double *B, double *xin, double *xout,
                         size_t dim, size_t it
) {
	if (it == 0)
		printf("# jacobi_tasks_serial\n");

	for (size_t i = 0; i < dim; ++i) {
		inst_event(9910002, dim);

		jacobi_base(&A[i * dim], B[i], xin, &xout[i], dim);

		inst_event(9910002, 0);
	}
}

int main(int argc, char* argv[])
{
	init_args(argc, argv);

	const char *PREFIX = basename(argv[0]);
	const size_t ROWS = create_cl_size_t ("Rows");

	printf("# Initializing data\n");
	timer ttimer = create_timer("Total_time");

	double *A = malloc(ROWS * ROWS * sizeof(double));
	double *B = malloc(ROWS * sizeof(double));
	double *x1 = malloc(ROWS * sizeof(double));
	double *x2 = malloc(ROWS * sizeof(double));

	init_AB_serial(A, B, ROWS);
	jacobi_modify_serial(A, B, ROWS);
	init_x_serial(x1, ROWS, 0.0);
	init_x_serial(x2, ROWS, 0.0);

	printf("# Starting algorithm\n");
	timer atimer = create_timer("Algorithm_time");

	double *xin = NULL;
	double *xout = NULL;

	for (int i = 0; i < 2; ++i) {
		xin = (i % 2 == 0) ? x1 : x2;
		xout = (i % 2 == 0) ? x2 : x1;

		jacobi_tasks_serial(A, B, xin, xout, ROWS, i);
	}

	stop_timer(&atimer);

	printf("# Finished algorithm...\n");
	stop_timer(&ttimer);

	/* printmatrix(A, ROWS, ROWS, PREFIX); */
	/* printmatrix(B, ROWS, 1, PREFIX); */
	/* printmatrix(xout, ROWS, 1, PREFIX); */

	free(A);
	free(B);
	free(x1);
	free(x2);

	create_reportable_int("worldsize", nanos6_get_num_cluster_nodes());
	create_reportable_int("cpu_count", nanos6_get_num_cpus());
	create_reportable_int("namespace_enabled", nanos6_get_namespace_is_enabled());

	report_args();
	free_args();

	return 0;
}
