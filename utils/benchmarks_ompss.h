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

#ifndef BENCHMARKS_OMPSS_H
#define BENCHMARKS_OMPSS_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <libgen.h>  // basename
#include <assert.h>
#include <syscall.h>
#include <errno.h>

#include <numaif.h>

#include <nanos6.h>
#include <nanos6/debug.h>

#include "cmacros/macros.h"
#include "ArgParserC/argparser.h"

// Define extrae if macro is set.
#if __WITH_EXTRAE
#include "extrae_user_events.h"

	typedef extrae_type_t inst_type_t;
	typedef extrae_value_t inst_value_t;
#define inst_define_event_type(type,name,nvalues,values,descriptions)	\
	Extrae_define_event_type(type,name,nvalues,values,descriptions)
#define inst_event(evt, val) Extrae_event(evt, val)
#else // __WITH_EXTRAE
	typedef size_t inst_type_t;
	typedef size_t inst_value_t;
#define inst_define_event_type(type,name,nvalues,values,descriptions)
#define inst_event(evt, val)
#endif // __WITH_EXTRAE


	void __print_task(const double * const mat,
	                  const size_t rows, const size_t cols,
	                  const char prefix[64], const char name[64]
	) {
		#pragma oss task in(mat[0; rows * cols]) label("matrix_print")
		{
			__print(mat, rows, cols, prefix, name);
		}
	}

#define printmatrix_task(mat, rows, cols, prefix)	\
	__print_task(mat, rows, cols, prefix, #mat)


	int get_numa_from_address(void *ptr)
	{
		int status;
		int numa_node = -1;
		const int ret = get_mempolicy(&numa_node, NULL, 0, ptr, MPOL_F_NODE | MPOL_F_ADDR);
		/* const int ret = move_pages(0 , 1, &ptr, NULL, &numa_node, 0); */

		if (ret != 0) {
			int errnum = errno;
			perror("Numa error: ");
			return -1;
		}

		return numa_node;
	}


#ifdef __cplusplus
}
#endif

#endif // BENCHMARKS_OMPSS_H
