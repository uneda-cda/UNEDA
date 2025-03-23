/*
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
 *   Copyright (c) 2012-2025  Prof. Mats Danielson, Stockholm University
 *
 *   Website: https://people.dsv.su.se/~mad/UNEDA
 *   GitHub:  https://github.com/uneda-cda/UNEDA
 *
 *   Licensed under CC BY 4.0: https://creativecommons.org/licenses/by/4.0/.
 *   Provided "as is", without warranty of any kind, express or implied.
 *   Reuse and modifications are encouraged, with proper attribution.
 *
 */

/*
 *   File: uct.h
 *
 *   Purpose: the header file for UCT (the UNEDA Core Tester).
 *
 */

/* BATCH MODE SWITCH: 0=INTERACTIVE  1=BATCH  2=PERFORMANCE */
/* Removed from public release on June 6, 2025 */
#define BATCH 0

#if BATCH == 0 // Interactive
#define NO_BATCH_MODE
#endif

#if BATCH == 1 // Verification sequence
#define BATCH_MODE
#define VERIFICATION
#endif

#if BATCH == 2 // Performance timing
#define BATCH_MODE
#define PROFILE
#endif

/* PATENT SWITCH: 0=NO PATENT  1=PATENT */
/* Will be inactivated on June 6, 2025 */
#define PATENT 1

#if PATENT == 0
/* Don't use patented code */
#define NO_DDT_PATENT
#endif

#if PATENT == 1
/* Have license, use patented code */
#define DDT_PATENT
#define TCL_PATENT
#endif

/* Copied from UNEDA-DTL */
#define DTL_MAIN 7
#define DTL_FUNC 21
#define DTL_TECH 1

/* Make directives */
#include "make.h"
#ifdef THINKC
#include <console.h>
#endif
#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>
#ifdef UNIX
#include <unistd.h>
#else
#include <io.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#ifdef PERFORMANCE
#include <Perf.h>
#endif
/* UNEDA header files */
#include "DTLparameters.h"
#include "TCL.h"
/* UCT header files */
#include "global.h"
#include "prototype.h"

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

/* Input modes */
#define ALT_MODE 1
#define V_MODE 2
#define PW_MODE 3

#define EPS 1E-5

#ifndef max
#define max(a,b) (a >= b ? a : b)
#define min(a,b) (a <= b ? a : b)
#endif

/* Console parameters */
#define N_CONSOLE_ROWS 72
	
struct user_frame {
	char frame_name[FNSIZE+2];
	char alt_name[MAX_ALTS+1][32];
	struct d_frame *df;
	int multilevel;
	int n_crit;
	int loaded;
	int dirty;
	double v_lo,v_up;
	int weak;
	};


/* from ufile.c */
int read_ufile(char *fn, char *folder, struct user_frame *uf);
int write_ufile(char *fn, char *folder, struct user_frame *uf);
int backup_ufile(char *fn, char *folder);

/* from random.c */
void init_random();
double xrandom();
double drandom(double lobo, double upbo);
int irandom(int lobo, int upbo);

