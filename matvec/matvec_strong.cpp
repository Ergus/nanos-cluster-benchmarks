/*
 * Copyright (C) 2019  Jimmy Aguilar Mena
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>

#include "ArgParserC/argparser.h"

#include "matvec.h"

void matvec_tasks_strong(const double *A, const double *b, double *x,
                         int rows, int ts, size_t it)
{
	assert(ts <= rows);
	myassert (rows % ts == 0);

	for (size_t i = 0; i < rows; i += ts) {
		#pragma oss task in(A[i * rows; rows * ts]) in(b[0; rows]) out(x[i; ts])
		matvec_base(&A[i * rows], b, &x[i], ts, rows);
	}
}

int main(int argc, char* argv[])
{
	init_args(argc, argv);

	const int ROWS = create_cl_int ("Rows");
	const int ts = create_cl_int ("Tasksize");
	const int its = create_optional_cl_int ("iterations", 1);
	const int print = create_optional_cl_int ("print", 0);

	std::cout << "Initializing data" << std::endl;
	timer *ttimer = create_timer("Total time");

	double *A = alloc_init(ROWS, ROWS, ts); // this initialized by blocks ts x rows
	double *b = alloc_init(ROWS, 1, ts);    // this splits the array in ts
	double *x = alloc_init(ROWS, 1, ROWS);  // This one initializes all the arrays
	#pragma oss taskwait

	std::cout << "Starting algorithm" << std::endl;
	timer *atimer = create_timer("Algorithm time");

	for (int i = 0; i < its; ++i)
		matvec_tasks_strong(A, x, b, ROWS, ts, i);
	#pragma oss taskwait

	stop_timer(atimer);

	std::cout << "Finished algorithm..." << std::endl;

	if (print) {
		matvec_print2d(A, ROWS, ROWS, "matrix_A.mat");
		matvec_print1d(b, ROWS, "vector_b.mat");
		matvec_print1d(x, ROWS, "vector_x.mat");

		const bool valid = validate(A, x, b, ROWS, ROWS);

		std::cout << "Verification: "
		          << (valid ? "Success" : "Failed")
		          << std::endl;
	}
	stop_timer(ttimer);

	free_matrix(A, ROWS * ROWS);
	free_matrix(x, ROWS);
	free_matrix(b, ROWS);

	const double performance = its * ROWS * ROWS * 2000.0 / getNS_timer(atimer);

	create_reportable_double("performance", performance);
	report_args();
	free_timer(atimer);
	free_timer(ttimer);
	free_args();

	return 0;
}
