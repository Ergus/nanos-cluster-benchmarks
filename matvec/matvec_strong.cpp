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

#include <iostream>
#include <cstdlib>

#include "argparser.h"

#include "matvec.h"

void matvec_tasks(const double *A, const double *b, double *x,
            size_t rows, size_t cols, size_t TS)
{
	assert(TS <= rows);
	myassert (rows % TS == 0);

	for (size_t i = 0; i < rows; i += TS) {
		#pragma oss task in(A[i * cols; cols * TS]) in(b[0; cols]) out(x[i; TS])
		matvec_base(&A[i * cols], b, &x[i], TS, cols);
	}
	#pragma oss taskwait
}

int main(int argc, char* argv[])
{
	init_args (argc, argv);

	const int ROWS = create_cl_int ("Rows");
	const int COLS = create_cl_int ("Columns");
	const int TS = create_cl_int ("Task size");
	const int its = create_optional_cl_int ("iterations", 1);
	const int print = create_optional_cl_int ("print", 1);

	std::cout << "Initializing data" << std::endl;
	timer *ttimer = create_timer("Total time");

	double *A = alloc_init(ROWS, COLS, TS);   // This initialized by blocks TS x cols
	double *b = alloc_init(ROWS, 1, TS);      // this splits the array in TS
	double *x = alloc_init(COLS, 1, COLS);    // This one initializes all the arrays
	#pragma oss taskwait

	std::cout << "Starting algorithm" << std::endl;
	timer *atimer = create_timer("Algorithm time");

	matvec_tasks(A, x, b, ROWS, COLS, TS);

	free_timer(atimer);

	std::cout << "Finished algorithm..." << std::endl;

	if (print) {
		matvec_print2d(A, ROWS, COLS, "matrix_A.mat");
		matvec_print1d(b, ROWS, "vector_b.mat");
		matvec_print1d(x, COLS, "vector_x.mat");

		const bool valid = validate (A, b, x, ROWS, COLS);

		std::cout << "Verification: " << (valid ? "Success" : "Failed") << std::endl;
		std::cout << "Done printing results..." << std:: endl;
	}
	free_timer(ttimer);

	free_matrix(A, ROWS * COLS);
	free_matrix(x, COLS);
	free_matrix(b, ROWS);

	const double performance = its * ROWS * COLS * 2000.0 ; // (double) atimer;

	create_reportable_double ("performance", performance);
	report_args ();
	free_args ();

	return 0;
}
