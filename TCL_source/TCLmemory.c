/*
 *
 *
 *        _/       _/   _/       _/    _/_/_/_/_/   _/_/_/          _/
 *       _/       _/   _/_/     _/    _/           _/    _/       _/  _/
 *      _/       _/   _/ _/    _/    _/           _/      _/    _/    _/
 *     _/       _/   _/  _/   _/    _/_/_/_/     _/      _/   _/      _/ 
 *    _/       _/   _/   _/  _/    _/           _/      _/   _/_/_/_/_/  
 *   _/       _/   _/    _/ _/    _/           _/      _/   _/      _/          
 *   _/     _/    _/     _/_/    _/           _/     _/    _/      _/          
 *    _/_/_/     _/       _/    _/_/_/_/_/   _/_/_/_/     _/      _/   
 *
 *
 *   UNEDA - The Universal Engine for Decision Analysis
 *
 *   Website: https://people.dsv.su.se/~mad/UNEDA
 *   GitHub:  https://github.com/uneda-cda/UNEDA
 *
 *   Licensed under CC BY 4.0: https://creativecommons.org/licenses/by/4.0/.
 *   Provided "as is", without warranty of any kind, express or implied.
 *   Reuse and modifications are encouraged, with proper attribution.
 *
 *
 *
 *                     UNEDA Tree Core Layer (TCL)
 *                     --------------------------- 
 *
 *    +----- o o o -----------------------------------------------+
 *    |    o       o             Prof. Mats Danielson             |
 *    |   o  STHLM  o            DECIDE Research Group            |
 *    |   o         o   Dept. of Computer and Systems Sciences    |
 *    |   o   UNI   o            Stockholm University             |
 *    |    o       o     PO Box 1203, SE-164 25 Kista, SWEDEN     |
 *    +----- o o o -----------------------------------------------+
 *
 *                Copyright (c) 2012-2025 Mats Danielson
 *                     Email: mats.danielson@su.se
 *
 */

/*
 *   File: TCLmemory.c
 *
 *   Purpose: memory allocation functions
 *
 *
 *   Functions exported outside TCL
 *   ------------------------------
 *   mem_alloc
 *   mem_free
 *   mem_map
 *   mem_exit
 *
 *   Functions outside module, inside TCL
 *   ------------------------------------
 *   NONE
 *
 *   Functions internal to module
 *   ----------------------------
 *   NONE
 *
 */

#include "TCLinternal.h"
#include <malloc.h>

#define MAX_CHUNKS 10000

struct mem_chunk {
	void *ptr;
	size_t size;
	char *type;
	char *source;
	};

static struct mem_chunk mem_tab[MAX_CHUNKS];
static int n_chunks = 0;


 /*********************************************************
  *
  *  Memory allocation functions
  *
  *********************************************************/

void *mem_alloc(size_t size, char *type, char *source) {
	void *mem_ptr;

	/* Check input parameters */
	if (n_chunks >= MAX_CHUNKS)
		return NULL;
	/* Get memory chunk */
	mem_ptr = malloc(size);
	if (!mem_ptr)
		return NULL;
	/* Fill entries */
	mem_tab[n_chunks].ptr = mem_ptr;
	mem_tab[n_chunks].size = size;
	mem_tab[n_chunks].type = type;
	mem_tab[n_chunks++].source = source;
	return mem_ptr;
	}


rcode mem_free(void* mem_ptr) {
	int i;

	/* Check input parameters */
	for (i=0; i<n_chunks; i++)
		if (mem_ptr == mem_tab[i].ptr)
			break;
	if (i == n_chunks)
		return TCL_INPUT_ERROR;
	/* Move last entry up to fill hole */
	mem_tab[i].ptr = mem_tab[--n_chunks].ptr;
	mem_tab[i].size = mem_tab[n_chunks].size;
	mem_tab[i].type = mem_tab[n_chunks].type;
	mem_tab[i].source = mem_tab[n_chunks].source;
	free(mem_ptr);
	return TCL_OK;
	}


void mem_map() {
	int i;

	/* Print all allocated memory chunks */
	for (i=0; i<n_chunks; i++)
		fprintf(stderr,"0x%p  %s  size: %ld  source: %s\n",mem_tab[i].ptr,
					 mem_tab[i].type,(long)mem_tab[i].size,mem_tab[i].source);
	}


rcode mem_exit() {

	/* Check no memory chunks remaining */
	if (n_chunks) {
		mem_map();
		return TCL_MEMORY_LEAK;
		}
	return TCL_OK;
	}
