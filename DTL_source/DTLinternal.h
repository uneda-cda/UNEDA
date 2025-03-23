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
 *                   UNEDA Decision Tree Layer (DTL)
 *                   -------------------------------
 *
 *    +----- o o o ------------------------------------------------+
 *    |    o       o              Prof. Mats Danielson             |
 *    |   o  STHLM  o             DECIDE Research Group            |
 *    |   o         o    Dept. of Computer and Systems Sciences    |
 *    |   o   UNI   o             Stockholm University             |
 *    |    o       o      PO Box 1203, SE-164 25 Kista, SWEDEN     |
 *    +----- o o o ------------------------------------------------+
 *
 *                Copyright (c) 2012-2025 Mats Danielson
 *                     Email: mats.danielson@su.se
 *
 */

/*
 *   File: DTLinternal.h
 *
 *   Purpose: Internal header file for the DTL library.
 *            Must be included by all DTL modules.
 *
 *
 *   DTL layers
 *   ----------
 *   Layer 0: Add-in functions above DTL proper (CAR, SML)
 *   Layer 1: API functions that handle the smx semaphore
 *   Layer 2: Support functions below the API (invisible)
 *
 *   If a function belongs to another layer than the obvious
 *   one, this is indicated by a comment in the beginning.
 *
 */


/*  DTL package C language note
 *  ---------------------------
 *
 *  ANSI C, ISO C, and Standard C are standards for the C
 *  programming language published by ANSI, ISO, and IEC.
 *  The names refer to the original version of the standard
 *  known as C89. The first standard for C was published by
 *  ANSI. That document was subsequently adopted by ISO/IEC
 *  and later revisions published by ISO/IEC have been adopted
 *  by ANSI. The standard was completed in 1989 and ratified
 *  as ANSI X3.159-1989 Programming Language C. This version
 *  of the language is mostly referred to as C89 or ANSI C in
 *  order to distinguish it from C90 which was ratified by
 *  ISO/IEC as ISO/IEC 9899:1990 with only formatting changes.
 *  Thus, the terms C89 and C90 refer to the same language
 *  and the terminology of C89 is being used in this package
 *  since it was conceived in 1994 and continuously evolved.
 *
 *  One reason to stay with C89 is that most C compilers are
 *  actually C++ compilers having C as a proper subset. But
 *  this subset is quite often only C89 with parts of C99,
 *  with one reason being that C and C++ later diverged.
 *  Thus, portability is ensured by sticking to using C89.
 *
 *  In reality, there are a few exceptions to the C89 rule.
 *  They have mostly to do with function declarations.
 *
 *  The following additional rules ensure optimal portability:
 *  1. One-line comments from C99 using // are allowed
 *  2. K&R-style functions from C89 are disallowed
 *  3. Anything invalidated in C99 is disallowed, such as e.g.
 *     o Implicit int function declarations
 *     o Function declarations without parameter specifications
 *     o Array type in struct without size specification
 *
 *  This should ensure that the DTL/TCL libraries will compile
 *  on any well-known platform using C compilers up to C23 and
 *  C++ compilers up to C++23, thus being rather future-safe.
 */


 /*****************************************************
  *
  *  Standard C library includes
  *
  *****************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <limits.h>
#include <signal.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <errno.h>
#include <setjmp.h>


 /*****************************************************
  *
  *  Exported calls for testing and debugging
  *
  *****************************************************/

#include "DTLdebug.h"


 /*****************************************************
  *
  *  MS VC++ compiler directive
  *
  *****************************************************/

#ifdef _MSC_VER
/* This is deliberate C89 code since that always compiles under C++.
 * Stop nagging about sprintf etc. (which are correct in C89/C99) */
#pragma warning(disable:4996)
#endif


 /********************************************************
  *
  *  Configuration parameters
  *
  ********************************************************/

#define ZOR_MML // only for ZOR debugger

#ifdef ZOR_MML
/* Enable output of DTL trace and debug info to console window */
#define CONSOLE
#endif

/* Enable internal test of error paths */
#define noINTERNAL_TEST

/* trunc() is in C99 but not in MS VC++ 6.0 */
#define _no_trunc


 /*****************************************************
  *
  *  Declarations local to DTL
  *
  *****************************************************/

#include "TCL.h"

#define DTL_MAIN 7
#define DTL_FUNC 21
#define DTL_TECH 1

#define DTLF_SIZE  224 // file name size
#define FILE_EXISTS 17 // posix error code

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#ifndef max
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#endif

/* Frame type checks */
#define PS (uf->frame_type == PS_FRAME)
#define DM (uf->frame_type == DM_FRAME)
#define PM (uf->frame_type == PM_FRAME)

/* Constants for B-normal calculations */
#define PI 3.1415926535897932384626433832795028841
#define DELTAPI 1.13799131882385 // 2.0*((4.0-PI)/2.0)^(2.0/3.0)


 /*****************************************************
  *
  *  Local data structures
  *
  *****************************************************/

struct user_frame {
	/* Permanent - preserved on file save */
	int frame_type;
	int frame_nbr;
	char frame_name[FNSIZE+1];
	int n_alts;
	int n_crit;
	int n_sh; // not yet saved due to compatibility with old file format
	/* Dynamic - not preserved on file save */
	int load_crit;
	struct d_frame *df;
	struct d_frame *df_list[MAX_CRIT+1];
	/* Dynamic - not preserved on file save since
	 * belonging to an add-in package (CAR or AS) */
	int WP_autogen[MAX_CRIT+1]; // CAR
	int V_n_rels[MAX_CRIT+1];   // CAR
	double av_min[MAX_CRIT+1];  // AS
	double av_max[MAX_CRIT+1];  // AS
	};

struct bn_rec {
	/* B-normal distribution parameters */
	int valid;
	double location;
	double scale2;
	double gr_min;
	double gr_max;
	double alpha;
	};


 /*****************************************************
  *
  *  Global data within DTL
  *
  *****************************************************/

extern int frame_loaded;
extern int dtl_error_count;
extern int dtl_trace_count;

extern struct user_frame *uf;
extern struct user_frame *uf_list[];

extern int cst_on;
extern int cst_ext;
extern FILE* cst;
extern char strg[],msg[];
extern char dtl_buffer[];


 /*****************************************************
  *
  *  Internal DTL calls - not exported
  *
  *****************************************************/

// DTLinternal.c
int dtl_is_init();
rcode dtl_kernel_error();
rcode dtl_error(rcode drc);
rcode call(rcode rc, char *proc);
int f_NaN(double test_num);
int f_NaN2(double test_num);
int imod_(int x, int y);
double fmod_(double x, double y);
int trace_ifperr(int type, int rv);
void shut_down(int tag);
void handle_illegal_ptr(void *ptr, int tag);
void _certify_ptr(void *ptr, int tag);
int get_days();
rcode cst_open();
void cst_close();
void cst_log(char *msg);
void cst_trace(char *msg);
rcode load_PV_stmt(int crit, struct user_stmt_rec *ustmt, struct stmt_rec *stmt, char type);
rcode load_W_stmt(struct user_w_stmt_rec *ustmt, struct stmt_rec *stmt);
struct user_frame *new_uf();
struct user_frame *get_uf(int index);
int dispose_uf(int index);
rcode load_df0(int crit);
rcode check_df0(int crit);
rcode load_df00(int crit);
rcode load_df1(int crit);
rcode check_df1(int crit);
#ifndef _trunc
double trunc(double d);
#endif

// DTLbnormal.c
double sgn(double x);
double n_cdf(double x);
double owens_t(double x, double alpha);
double bn_cdf(double val, double mean, double var, double alpha);
double bn_inv_cdf(double cdf, double mean, double var, double alpha);
double b_delta(double skew);

// DTLframe.c
rcode dtl_dispose_frame(int ufnbr);

// DTLmisc.c
int dtl_node2crit(int node);
int dtl_crit2node(int crit);
rcode dtl_is_shadow_crit(int crit);

// DTLeval.c
void sort_b(int order[], double maxmin[], int start, int stop, bool max);
void eval_cache_invalidate();
rcode evaluate_frame(int crit, int method, int Ai, int Aj, e_matrix e_result);
rcode evaluate_frameset(int crit, int method, int Ai, int Aj, e_matrix e_result);
rcode dtl_ev_to_cdf(int crit, double ev_level, double *mass);
rcode dtl_cdf_to_ev(int crit, double belief_level, double *lobo, double *upbo);

// DTLwbase.c
rcode dtl_set_W_check(h_vector lobox, h_vector mbox, h_vector upbox);
rcode dtl_set_W_mbox_auto(h_vector lobox, h_vector upbox);
int dtl_W_node_parents(int node1, int node2);
int dtl_W_nbr_of_siblings(int node);
int dtl_pure_W_tree();
int dtl_real_W_crit(int node);
int dtl_nbr_W_midpoints();

// DTLpbase.c
rcode dtl_set_P_check(int crit, h_matrix lobox, h_matrix mbox, h_matrix upbox);
rcode dtl_set_P_mbox_auto(int crit, h_matrix lobox, h_matrix upbox);
int dtl_P_node_parents(int crit, int alt, int node1, int node2);
int dtl_P_nbr_of_siblings(int crit, int alt, int node);
int dtl_nbr_P_midpoints(int crit);

// DTLvbase.c
rcode dtl_set_V_check(int crit, h_matrix lobox, h_matrix mbox, h_matrix upbox);
rcode dtl_set_V_mbox_rels(int crit, int V_n_rels, h_matrix lobox, h_matrix upbox);
int dtl_nbr_V_midpoints(int crit);
double dtl_mid_to_modal(double lobo, double mid, double upbo, int oor);
#define dtl_real_V_node DTI_real_V_node
rcode dtl_check_V_modality(int crit, int Ai, int Aj);

// DTLautoscale.c
rcode dtl_copy_AV_crit_scale(int cr_from, int cr_to);

// DTLdominance.c
rcode dtl_get_dominance(int crit, int Ai, int Aj, double *cd_value, int *d_order);

// TCL.h
extern t_matrix t2f,t2r,t2i,r2t,i2t;
extern tn_row f2r,f2i,r2f,i2f,i2end;

// SML.h
int sml_is_open();


 /*****************************************************
  *
  *  Abort handling
  *
  *****************************************************/

extern int dtl_abort_request;
extern int dtl_abort_cst;

#define dtl_abort_init() \
	dtl_abort_cst = cst_on; \
	dtl_abort_request = FALSE

#define dtl_abort_check() \
	if (dtl_abort_request) { \
		cst_on = dtl_abort_cst; \
		if (cst_on) { \
			sprintf(msg," *** DTL::%s user abort ***\n",dtl_func); \
			cst_log(msg); \
			} \
		_smx_end(); \
		return DTL_USER_ABORT; \
		}


 /*****************************************************
  *
  *  Single thread mechanism - prohibits reentrancy
  *
  *****************************************************/

/* Single thread mutex semaphore: 99.999% perfect. There are small
 * time-critical holes connected to calling _smx_end() twice in an
 * error procedure and having a thread scheduled in between allowing
 * entry of another thread thus being prematurely marked non-busy
 * by the second _smx_end() intended for the first caller.
 *
 * smx calling conventions
 * -----------------------
 * DTL_ calls are API calls and must end the smx semaphore
 *
 * dtl_ calls are API support calls (often static) and could end the smx semaphore
 *
 * DTI_ calls are test and debug API calls and must end the smx semaphore
 *
 * Functions with other names should not end the smx semaphore themselves but rather
 * rely on the callers on top of them in the call stack ending them in time due.
 *
 * There are a few exceptions to this rule, most notably evaluate_frame and
 * evaluate_frameset in DTLeval.c. */

extern int smx_busy;
extern char *dtl_func;

#define _smx_begin(fn) \
{	if (smx_busy) \
		return DTL_BUSY; \
	dtl_func = fn; \
	_init_assert(); \
	smx_busy = TRUE; }

#define _smx_continue _smx_begin

#define _smx_end() \
{	dtl_func = "NULL"; \
	smx_busy = FALSE; }

#define _smx_name(fn) \
{	dtl_func = fn; }

#define _smx_noname() \
{	dtl_func = "NULL"; }


 /***********************************************************
  *
  *  Assert = in-line conditional code guard. Log error and
  *  return an error code to the caller. Compare this to the
  *  built-in C assert which exits the program - too harsh.
  *
  ***********************************************************/

extern jmp_buf assert_envir;
extern int a_tag;

/* Must be a macro inline, not a function because there
 * is nothing to return to on the stack after a longjmp */

#define _init_assert() \
	if (a_tag = setjmp(assert_envir)) { \
		sprintf(msg," DTL::%.10s runtime assert failed (%03d)\n",strcmp(dtl_func,"NULL")?dtl_func:"SYS",a_tag); \
		if (cst_on) \
			cst_log(msg); \
		else \
			cst_trace(msg); \
		return dtl_error(DTL_ASSERT_FAILED); \
		} \

/* Important: _dtl_assert() must be preceeded by _smx_begin() that
 * in turn contains _init_assert(), else nowhere to jump -> SIGSEGV */

#define _dtl_assert(cond,tag) \
	if (!(cond)) \
		longjmp(assert_envir,tag)


 /************************************************************
  *
  *  Guarded math operators
  *  ----------------------
  *
  *  These operators are for testing out new algorithms,
  *  not for permanent use in released library versions.
  *
  *  Three possible versions of div, rem, mod, etc:
  *  1. inline macro -> write trace on error and continue
  *  2. inline macro with longjmp -> will abort execution
  *  3. function call -> yields inefficient stack handling
  *  (only version 1 is implemented in DTL for now)
  *
  *  Two possible versions of iabs:
  *  1. Saturated integer abs: interpret INT_MIN as -INT_MAX
  *     (2-comp hardware does bitwise invert followed by +1)
  *  2. NaN integer abs: INT_MIN as NaN/overflow indicator
  *  (both versions are implemented in DTL)
  *
  *  Three possible versions of fabs:
  *  1. Detect and guard from NaN (INF & IND) exceptions
  *  2. Standard C library fabs
  *  3. f_abs defined to circumvent a compiler bug
  *  (all three versions exist in DTL)
  *
  *  NOTE: In C89, integer division is partially undefined.
  *        Negative numbers can be either rounded toward 0
  *        or toward negative infinity. As a consequence, %
  *        is equally partially undefined. Example: -5/2=-2
  *        and -5%2=-1 are allowed as are -5/2=-3 and -5%2=1.
  *        In C99 & C11, they are always -5/2=-2 and -5%2=-1.
  *
  *        Already MS VC++ 6.0 (which is pre-C99) has the C99
  *        versions of div and rem, even though other parts
  *        of C99 are still not implemented in VS 2019 (while
  *        the complete C99 is implemented in gcc 4.5 and on).
  *
  *  NOTE: In C89-C11, fmod is not mod at all but rather rem
  *        and this is not expected to change in C17 or C23.
  *        (C17 is a tidy-up of ISO C and C23 is yet unknown).
  *
  *  NOTE: These macro functions were much used during the
  *        development of the library but removed once the.
  *        correctness of the functions were certified.
  *        Some instances still remain in use to this day.
  *
  ************************************************************/

#define  idiv(x,y) ((y)==0?trace_ifperr(0,(x)==0?1:0):(x)/(y))           // x/0 takes precedence
#define _idiv(x,y) ((y)==0?(x)==0?1:trace_ifperr(0,0):(x)/(y))           // x/x takes precedence
#define  irem(x,y) ((y)==0?trace_ifperr(1,0):(x)%(y))                    // x irem 0 is 0 (0<=r<y & y->0 => r->0)
#define  imod(x,y) ((y)==0?trace_ifperr(2,0):imod_(x,y))                 // x imod 0 is 0 (0<=r<y & y->0 => r->0)
#define _imod(x,y) ((y)==0?trace_ifperr(2,x):imod_(x,y))                 // x imod 0 is x (x-r=n*y => x=r)
#define  emod(x,y) ((y)==0?trace_ifperr(2,0):imod_(x,abs(y)))            // Euclidean imod = true clock function
#define  iabs(x)   ((x)==INT_MIN?trace_ifperr(3,INT_MAX):((x)<0?-(x):(x))) // iabs(-iabs(x)) != x for INT_MIN
#define _iabs(x)   ((x)==INT_MIN?trace_ifperr(3,0):((x)<0?-(x):(x)))     // iabs controlled NaN handling
#define _isgn(x)   ((x)<0?-1:(x?1:0))                                    // isgn(x) is 0 for x=0
#define  fdiv(x,y) ((y)==0.0?(double)trace_ifperr(0,(x)==0.0?1:0):(x)/(y))     // x/0 takes precedence
#define _fdiv(x,y) ((y)==0.0?((x)==0.0?1.0:(double)trace_ifperr(0,0)):(x)/(y)) // x/x takes precedence
#define _frem(x,y) ((y)==0.0?(double)trace_ifperr(1,0):fmod(x,y))        // C fmod is not fmod but frem (!!)
#define _fmod(x,y) ((y)==0.0?(double)trace_ifperr(2,0):fmod_(x,y))       // the classic fmod is here instead
#define  gmod(x,y) ((y)==0.0?(double)trace_ifperr(2,0):fmod_(x,fabs(y))) // Euclidean fmod = true clockwork function
#define _gmod(x,y) ((y)==0.0?(x)+trace_ifperr(2,0):fmod_(x,fabs(y)))     // Euclidean fmod = true clockwork function
#define _fabs(x)   (f_NaN(x)?(double)trace_ifperr(3,0):fabs(x))          // detect NaN (INF & IND)
#define _fsgn(x)   (f_NaN(x)?(double)trace_ifperr(4,0):((x)<0.0?-1.0:(x?1.0:0.0)))     // _fsgn(x) yields 0.0 for x=0.0
#define _sqrt(x)   (_fsgn(x)<0.0?((x)>-DTL_EPS?0.0:(double)trace_ifperr(5,0)):sqrt(x)) // _sqrt(x) forgives rounding errors
