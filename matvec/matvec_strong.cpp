#include "CommandLineParameter.hpp"
#include "Report.hpp"
#include "Reportable.hpp"
#include "Timer.hpp"

#include <iostream>

#include <cstdlib>

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

int main(int argc, char* argv[]) {

	CommandLine::initialize(argc, argv, "matvec", "cluster");
	CommandLineParameter<long> ROWS("ROWS", "Number of rows", "rows of matrix");
	CommandLineParameter<long> COLS("COLS", "Number of columns", "columns of matrix");
	CommandLineParameter<long> TS("TS", "Task size", "task size");
	OptionalCommandLineParameter<long> its("its", 1, "iterations", "inner repetitions");
	OptionalCommandLineParameter<bool> print("print", 0, "print matrices", "print matrices");
	CommandLine::validate();

	std::cout << "Initializing data" << std::endl;
	Timer ttime("ttime", "Total execution time");

	double *A = alloc_init(ROWS, COLS, TS);   // This initialized by blocks TS x cols
	double *b = alloc_init(ROWS, 1, TS);      // this splits the array in TS
	double *x = alloc_init(COLS, 1, COLS);    // This one initializes all the arrays
	#pragma oss taskwait

	std::cout << "Starting algorithm" << std::endl;

	Timer atimer("atime", "Algorithm execution time");

	matvec_tasks(A, x, b, ROWS, COLS, TS);

	ttime.stop();

	std::cout << "Finished algorithm..." << std::endl;

	if (print) {
		matvec_print2d(A, ROWS, COLS, "matrix_A.mat");
		matvec_print1d(b, ROWS, "vector_b.mat");
		matvec_print1d(x, COLS, "vector_x.mat");

		const bool valid = validate (A, b, x, ROWS, COLS);

		std::cout << "Verification: " << (valid ? "Success" : "Failed") << std::endl;
		std::cout << "Done printing results..." << std:: endl;
	}

	free_matrix(A, ROWS * COLS);
	free_matrix(x, COLS);
	free_matrix(b, ROWS);

	const double performance = its * ROWS * COLS * 2000.0 / (double) atimer;

	ReportEntry<double> passes("RESULT", "performance",
	                           "Megaflops per second", performance, "MFlops/s");

	Report::emit();

	return 0;
}
