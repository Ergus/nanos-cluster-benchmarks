
#ifndef benchmarks_h
#define benchmarks_h

// C headers
#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_set_num_threads(var) {}
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

#define printme() {	  \
		fprintf(stderr,"Func: %s in %s:%d (process %s)\n", \
		        __PRETTY_FUNCTION__, __FILE__, __LINE__, \
		        getenv("OMPI_COMM_WORLD_RANK")); \
	}

#define modcheck(a, b){	  \
		const int themod=(a) % (b); \
		if(themod) { \
			fprintf(stderr,"Error: %s %% %s = %d\n", #a, #b, themod); \
			printme(); \
			MPI_Abort(MPI_COMM_WORLD, -1); \
		} \
	}

#define myassert(cond) {										\
		if (!(cond)) {											\
			fprintf(stderr, "%s%s:%u Assertion `%s' failed.\n", \
			        __func__, __FILE__, __LINE__, #cond);		\
			abort();											\
		}														\
	}

#ifndef NDEBUG
#define dprintf(...) fprintf(stderr,__VA_ARGS__)
#else
#define dprintf(...) {}
#endif

#ifdef __cplusplus

#include <string>
#include <sstream>

template<typename... Args>
inline bool string_in(std::string var, std::string first, Args... args)
{
	return (var==first) || string_in(var,args...);
}

template<>
inline bool string_in(std::string var, std::string first)
{
	return (var==first);
}

template<class T>
inline std::string toString(const T& t, bool *ok = NULL)
{
	std::ostringstream stream;
	stream << t;
	if(ok != NULL)
		*ok = (stream.fail() == false);
	return stream.str();
}

#endif // __cplusplus

#endif // benchmarks_h
