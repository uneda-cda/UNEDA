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
 *   File: DTLmisc.c
 *
 *   Purpose: DTL management + TCL management calls
 *
 *
 *   Functions exported outside DTL
 *   ------------------------------
 *   DTL_abort
 *   DTL_init
 *   DTL_exit
 *   DTL_get_release
 *   DTL_get_release_long
 *   DTL_get_capacity
 *   DTL_get_J_properties
 *   DTL_nbr_of_W_stmts
 *   DTL_nbr_of_P_stmts
 *   DTL_nbr_of_V_stmts
 *   DTL_nbr_of_weights
 *   DTL_nbr_of_crit
 *   DTL_nbr_of_alts
 *   DTL_total_cons
 *   DTL_total_nodes
 *   DTL_nbr_of_cons
 *   DTL_nbr_of_nodes
 *   DTL_error/2
 *   DTL_u_error/2
 *
 *   Functions outside module, debug use
 *   -----------------------------------
 *   DTI_get_API_type
 *   DTI_list_all_frames
 *   DTI_is_tree
 *   DTI_tree_structure
 *   DTI_node2crit
 *   DTI_nbr_of_P_mbox
 *   DTI_nbr_of_V_mbox
 *
 *   Functions outside module, inside DTL
 *   ------------------------------------
 *   dtl_node2crit
 *
 *   Functions internal to module
 *   ----------------------------
 *   user_abort
 *   invalid_ptr
 *   draw_tree/2/3
 *
 */

#include "DTL.h"
#include "DTLinternal.h"
// fp exceptions
#ifdef _MSC_VER // MS Win
#include <float.h>
#else // Unix/Mac
#include <fenv.h>
#endif

 /********************************************************
  *
  *  Internal data
  *
  ********************************************************/

static int dtl_init = FALSE;


 /*********************************************************
  *
  *  JSON property structure
  *
  *********************************************************/

#define JPROP_NBR 11

static struct J_entry {
  char* name; // identifier (label)
  int value;  // content (number)
	} J_son;

static const struct J_entry J_props[JPROP_NBR] = {
//{"lib conf",  LIB_CONF},
  {"rel main",  DTL_MAIN},
  {"rel func",  DTL_FUNC},
  {"rel tech",  DTL_TECH},
  {"max frames",MAX_FRAMES-1}, // 1 reserved for internal use
  {"max crit",  MAX_CRIT},
  {"max alts",  MAX_ALTS},
  {"max cons",  MAX_CONS},
  {"max copa",  MAX_COPA},
  {"max nodes", MAX_NODES},
  {"max nopa",  MAX_NOPA},
  {"max stmts", MAX_STMTS-1} // 1 reserved for internal use
	};

static const char *LIB_CONF =
#ifdef C_SML
"SML"
#else
"DTL"
#endif
;

 /*********************************************************
  *
  *  DTL library management
  *
  *********************************************************/

void DTLAPI DTL_abort() {

	/* Must be called by a different thread in SAME address space.
	   If you are another process, use the signal SIGINT instead. */
	dtl_abort_request = TRUE;
	if (cst_on)
		cst_log("DTL_abort(0)\n");
	}


static void user_abort(int sig) {

	/* If you are another process, use this abort call
	   via an installed interrupt handler instead */
	dtl_abort_request = TRUE;
	if (cst_on) {
		sprintf(msg,"DTL_abort(%d)\n",sig);
		cst_log(msg);
		}
	signal(sig,user_abort); // reinstall handler
	}


static void invalid_ptr(int sig) {

	/* Called by signal SIGSEGV when an illegal address
	 * is not caught by certify_ptr but slipped through.
	 *
	 * Note: signals during load, most often indicating
	 * insufficient virtual memory, are not let through
	 * to the application; they will stay within the OS. */
	handle_illegal_ptr(NULL,sig);
	/* Do not reinstall handler, the signal will just
	 * reoccur when the instruction is restarted */
	_smx_end();
	exit(EXIT_FAILURE);
	}


static void calc_bug(int sig) {
	unsigned fpe_flag;

	/* Called by signal SIGFPE when an illegal floating point
	 * operation is not caught by macros but slipped through.
	 * On some non-MS systems, also for _integer_ div_by_0.
	 *
	 * Flags are implementation-dependent, but
	 * 0 is always no_flag (bitmasks since C11).
	 * MSVC: 0x04 = overflow (_EM_OVERFLOW)
	 *       0x08 = div_by_0 (_EM_ZERODIVIDE)
	 *       0x10 = invalid  (_EM_INVALID)
	 * gcc:  0x01 = invalid  (FE_INVALID)
	 *       0x04 = div_by_0 (FE_DIVBYZERO)
	 *       0x08 = overflow (FE_OVERFLOW)
	 * IEEE754 does not assign values/bits to errors.
	 */
#ifdef _MSC_VER
	fpe_flag = _clearfp();
#else // fetestexcept not in MSVC
	fpe_flag = fetestexcept(FE_ALL_EXCEPT);
#endif
	if (!fpe_flag) // flags cleared on entry
		fpe_flag = (unsigned)(sig<<3);
	/* Log error and close the run */
	sprintf(msg," Little man lost count (0x%02X)\n",fpe_flag);
	shut_down(sig);
	_smx_end();
	exit(EXIT_FAILURE);
	}


int dtl_is_init() {

	// status of DTL package
	return dtl_init;
	}


rcode DTLAPI DTL_init() {
	rcode rc;
	int i;

	/* Begin single thread semaphore */
	_smx_begin("INIT");
	/* Set up log file */
	if ((rc = cst_open()) == DTL_INTERNAL_ERROR)
		return dtl_error(rc);
	/* Log function call if log file open */
	if (cst_on)
		cst_log("DTL_init()\n");
	/* Check if function can start */
	if (dtl_init)
		return dtl_error(DTL_STATE_ERROR);
	/* Override default interrupt handlers */
	signal(SIGINT,user_abort);
	signal(SIGSEGV,invalid_ptr);
	/* Enable suitablefp exceptions before installing signal handler */
#ifdef _MSC_VER // MS Win: select triggers
	_controlfp(~(unsigned int)(_EM_ZERODIVIDE|_EM_INVALID|_EM_OVERFLOW),(unsigned int)_MCW_EM);
#else // Unix: select triggers
	feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
	fedisableexcept(FE_UNDERFLOW|FE_INEXACT);
#endif
	signal(SIGFPE,calc_bug);
	/* Initialise frame pointers (only auto in debug mode) */
	for (i=0; i<=MAX_FRAMES; i++)
		uf_list[i] = NULL;
	/* Now ready to fly */
	dtl_error_count = 0;
	dtl_trace_count = 0;
	dtl_init = TRUE;
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTL_exit() {
	rcode rc;
	int i;

	/* Begin single thread semaphore */
	_smx_begin("EXIT");
	/* Log function call */
	if (cst_on)
		cst_log("DTL_exit()\n");
	/* Check if function can start */
	if (!dtl_init)
		return dtl_error(DTL_STATE_ERROR);
	if (frame_loaded)
		return dtl_error(DTL_FRAME_IN_USE);
	/* Release resources */
	for (i=1; i<=MAX_FRAMES; i++) {
		if (uf_list[i])
			if (dtl_dispose_frame(i)) {
				_smx_continue("EXT2");
				}
		}
	rc = mem_exit();
	/* Log function result */
	if (rc)
		if (cst_on)
			cst_log("  Memory leakage detected\n");
		else
			cst_trace(" memory leakage detected\n");
	cst_close();
	dtl_init = FALSE;
	/* End single thread semaphore */
	_smx_end();
	if (rc)
		return DTL_MEMORY_LEAK;
	else
		return dtl_trace_count;
	}


 /*********************************************************
  *
  *  Structure commands
  *
  *  (Problems with size_t in VC++ 5.0 -> use unsigned int)
  *
  *********************************************************/

 /*
  * Call semantics: Obtains the release of the DTL package.
  * The format is "DM.DF.DT", where DM=DTL main version,
  * DF=DTL functional version and DT=DTL technical version.
  * The long version includes number of days of existence.
  */

void DTLAPI DTL_get_release(char *relstrg, unsigned c_size) {
	unsigned len;

	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_release(strg[%u])\n",c_size);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(relstrg,1);
	/* Collect information */
	len = sprintf(dtl_buffer,"%d.%02d.%d",DTL_MAIN,DTL_FUNC,DTL_TECH);
	if (len >= c_size) { // buffer overrun
		relstrg[0] = '\0';
		sprintf(msg," DTL_get_release(strg[%u]) buffer too short (min %u)\n",c_size,len+1);
		cst_trace(msg);
		return;
		}
	/* Transfer to caller */
	strcpy(relstrg,dtl_buffer);
	}


void DTLAPI DTL_get_release_long(char *relstrg, unsigned c_size) {
	unsigned len;

	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_release_long(strg[%u])\n",c_size);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(relstrg,1);
	/* Collect requested information */
	len = sprintf(dtl_buffer,"%d.%02d.%d [%d]",DTL_MAIN,DTL_FUNC,DTL_TECH,get_days());
	if (len >= c_size) { // buffer overrun
		relstrg[0] = '\0';
		sprintf(msg," DTL_get_release_long(strg[%u]) buffer too short (min %u)\n",c_size,len+1);
		cst_trace(msg);
		return;
		}
	/* Transfer to caller */
	strcpy(relstrg,dtl_buffer);
	}


 /*
  * Call semantics: Obtains the capacities of the DTL package.
  */

void DTLAPI DTL_get_capacity(char *capstrg, unsigned c_size) {
	unsigned len;

	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_capacity(strg[%u])\n",c_size);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(capstrg,1);
	/* Collect requested information */
	len = sprintf(dtl_buffer,"%d %d %d %d %d %d %d %d",
			MAX_FRAMES,MAX_CRIT,MAX_ALTS,MAX_NODES,MAX_NOPA,MAX_CONS,MAX_COPA,MAX_STMTS);
	if (len >= c_size) { // buffer overrun
		capstrg[0] = '\0';
		sprintf(msg," DTL_get_capacity(strg[%u]) buffer too short (min %u)\n",c_size,len+1);
		cst_trace(msg);
		return;
		}
	/* Transfer to caller */
	strcpy(capstrg,dtl_buffer);
	}


/* DTL_get_J_properties delivers a serialised JSON object
 * JSON is the JavaScript Object Notation exchange format
 * The J_strg argument is as if created by JSON.stringify()
 * Property names are encoded as strings
 * Property values (except the first) are numbers (integers)
 * The JSON object should be decoded by calling JSON.parse() */

void DTLAPI DTL_get_J_properties(char *J_strg, unsigned c_size) {
	int i;
	unsigned len;
	char prop[40]; // property identifier + value

	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_J_properties(strg[%u])\n",c_size);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(J_strg,1);
	/* Generate stringified JSON property object */
	sprintf(dtl_buffer,"{\"lib conf\":\"%s\",",LIB_CONF);
	for (i=0; i<JPROP_NBR; i++) {
		sprintf(prop,"\"%s\":%d,",J_props[i].name,J_props[i].value);
		strcat(dtl_buffer,prop);
		}
	dtl_buffer[strlen(dtl_buffer)-1] = '}'; // "," -> "}"
	len = (unsigned)strlen(dtl_buffer);
	if (len >= c_size) { // buffer overrun
		J_strg[0] = '\0';
		sprintf(msg," DTL_get_J_properties(strg[%u]) buffer too short (min %u)\n",c_size,len+1);
		cst_trace(msg);
		return;
		}
	/* Transfer to caller */
	strcpy(J_strg,dtl_buffer);
	}


void DTLAPI DTI_get_API_type(char *typestrg, unsigned c_size) {
	unsigned len;

	/* Log function call */
	if (cst_ext) {
		sprintf(msg,"DTI_get_API_type(strg[%u])\n",c_size);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(typestrg,1);
	/* Collect requested information */
	len = sprintf(dtl_buffer,"%s",
			CALL_STACK?(CALL_STACK==1?"_cdecl(C89)":(CALL_STACK==2?"_stdcall(MS_WIN)":"_non_std")):"_default");
	if (len >= c_size) { // buffer overrun
		typestrg[0] = '\0';
		sprintf(msg," DTI_get_API_type(strg[%u]) buffer too short (min %u)\n",c_size,len+1);
		cst_trace(msg);
		return;
		}
	/* Transfer to caller */
	strcpy(typestrg,dtl_buffer);
	}


#ifdef CONSOLE

/*
 * For use only by test and debug. Print structures not made
 * available on the outside of DTL due to abstraction layers.
 */

void DTLAPI DTI_list_all_frames() {
	int i,j,hit=0;

	/* Also runs unloaded - list all frames to choose from */
	for (i=1; i<=MAX_FRAMES; i++)
		if (uf_list[i]) {
			hit++;
			printf("%s-frame %d: %s %s\n",uf_list[i]->frame_type==PM_FRAME?"PM":uf_list[i]->frame_type==PS_FRAME?"PS":"DM",
					i,uf_list[i]->frame_name,i==frame_loaded?"[-loaded-]":"");
			if (uf_list[i]->frame_type == PM_FRAME)
				for (j=0; j<=MAX_CRIT; j++)
					if (uf_list[i]->df_list[j])
						printf("    crit %d: %s\n",j,uf_list[i]->df_list[j]->name);
			}
	if (!hit)
		printf("No frames\n");
	}


static int vbar[14];

static void draw_tree2(int alt, int snode, int level, int mode) {
	struct d_frame *df = uf->df;
	int i,tnode,is_at = level;

	/* Mode 0 = nodes, mode 1 = values, mode 2 = criteria numbers */
	for (tnode=df->down[alt][snode]; tnode; tnode=df->next[alt][tnode]) {
		if (!is_at)
			printf(" ");
		for (i=is_at; i<level; i++)
			if (vbar[i])
				printf(" |   ");
			else
				printf("     ");
		if (df->down[alt][tnode]) {
			printf("%2d---",mode?0:tnode);
			vbar[level] = df->next[alt][tnode];
			if (level<=12)
				draw_tree2(alt,tnode,level+1,mode);
			else
				printf("*\n");
			}
		else
			printf("%2d\n",mode>1?TCL_get_V_index(df,alt,tnode):tnode);
		is_at = 0;
		}
	}


static void draw_tree3(int alt, int snode, int level) {
	struct d_frame *df = uf->df;
	int i,tnode,is_at = level;

	/* For trees with 3-digit nodes (100 nodes or more) in mode 0 */
	for (tnode=df->down[alt][snode]; tnode; tnode=df->next[alt][tnode]) {
		if (!is_at)
			printf(" ");
		for (i=is_at; i<level; i++)
			if (vbar[i])
				printf("  |  ");
			else
				printf("     ");
		if (df->down[alt][tnode]) {
			if (tnode>9)
				printf("%3d--",tnode);
			else
				printf(" %02d--",tnode);
			vbar[level] = df->next[alt][tnode];
			if (level<=12)
				draw_tree3(alt,tnode,level+1);
			else
				printf("*\n");
			}
		else
			if (tnode>9)
				printf("%3d\n",tnode);
			else
				printf(" %02d\n",tnode);
		is_at = 0;
		}
	}


static void draw_tree(int alt, int mode) {
	struct d_frame *df = uf->df;

	if (mode)
		draw_tree2(alt,0,0,mode);
	else
		if (df->tot_cons[alt] > 99)
			draw_tree3(alt,0,0);
		else
			draw_tree2(alt,0,0,0);
	}


void DTLAPI DTI_tree_structure(int crit) {
	int i;

	if (!frame_loaded) {
		printf("No frame loaded\n");
		return;
		}
	if (load_df00(crit)) {
		if (crit < 0)
			printf("No intermediate node %d in weight tree\n",-crit);
		else
			printf("No criterion %d in frame\n",crit);
		return;
		}
	if (!(uf->df->tree)) {
		printf("Flat structure\n");
		return;
		}
	if (crit > 0) {
		printf("Event trees\n");
		for (i=1; i<=uf->df->n_alts; i++) {
			printf("Alternative %d\n",i);
			draw_tree(i,0);
			}
		printf("Value trees\n");
		for (i=1; i<=uf->df->n_alts; i++) {
			printf("Alternative %d\n",i);
			draw_tree(i,1);
			}
		}
	else if (crit < 0) {
		printf("Criteria subtree node numbers\n");
		draw_tree2(1,-crit,0,0);
		printf("Real criteria numbers\n");
		draw_tree2(1,-crit,0,2);
		return;
		}
	else {
		printf("Criteria tree node numbers (%spure)\n",dtl_pure_W_tree()?"":"im");
		draw_tree(1,0);
		printf("Real criteria numbers\n");
		draw_tree(1,2);
		}
	}

#endif


bool DTLAPI DTI_is_tree(int crit) {

	if (!frame_loaded)
		return FALSE;
	if (load_df0(crit))
		return FALSE;
	return uf->df->tree;
	}


 /*********************************************************
  *
  *  DTL number peeks
  *
  *  They return a negative number on error
  *  Not protected by smx mechanism
  *
  *********************************************************/

int DTLAPI DTL_nbr_of_W_stmts() {

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Return metric */
	if (PM)
		return uf->df_list[0]->P_base->n_stmts;
	else
		return 0;
	}


int DTLAPI DTL_nbr_of_P_stmts(int crit) {
	int i,sum;

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	if (PM) {
		/* Check input parameters */
		if ((crit < 0) || (crit > uf->n_crit))
			return DTL_CRIT_UNKNOWN;
		/* Return metric */
		if (crit) {
			if (uf->df_list[crit])
				return uf->df_list[crit]->P_base->n_stmts;
			else
				return 0;
			}
		else { // sum over all criteria
			for (sum=0,i=1; i<=uf->n_crit; i++)
				if (uf->df_list[i])
					sum += uf->df_list[i]->P_base->n_stmts;
			return sum;
			}
		}
	else // PS
		return crit==1?uf->df->P_base->n_stmts:DTL_CRIT_UNKNOWN;
	}


int DTLAPI DTI_nbr_of_P_mbox(int crit) {
	int i,stop,sum;
	struct base *P;

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	if (PM) { // also for W
	/* Check input parameters */
		if ((crit < 0) || (crit > uf->n_crit))
			return DTL_CRIT_UNKNOWN;
		/* Return metric */
		P = uf->df_list[crit]->P_base;
		stop = uf->df_list[crit]->tot_cons[0];
		}
	else { // PS
		if (crit != 1)
			return DTL_CRIT_UNKNOWN;
		P = uf->df->P_base;
		stop = uf->df->tot_cons[0];
		}
	/* Return metric */
	for (sum=0,i=1; i<=stop; i++)
		if ((P->lo_midbox[i] > -1.0) || (P->up_midbox[i] > -1.0))
			sum++;
	return sum;
	}


int DTLAPI DTL_nbr_of_V_stmts(int crit) {
	int i,sum;

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	if (PM) {
		/* Check input parameters */
		if ((crit < 0) || (crit > uf->n_crit))
			return DTL_CRIT_UNKNOWN;
		/* Return metric */
		if (crit) {
			if (uf->df_list[crit])
				return uf->df_list[crit]->V_base->n_stmts;
			else
				return 0;
				}
		else { // sum over all criteria
			for (sum=0,i=1; i<=uf->n_crit; i++)
				if (uf->df_list[i])
					sum += uf->df_list[i]->V_base->n_stmts;
			return sum;
			}
		}
	else // PS
		return crit==1?uf->df->V_base->n_stmts:DTL_CRIT_UNKNOWN;
	}


int DTLAPI DTI_nbr_of_V_mbox(int crit) {
	int i,stop,sum;
	struct base *V;

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameters */
	if (PM) {
	/* Check input parameters */
		if ((crit < 1) || (crit > uf->n_crit))
			return DTL_CRIT_UNKNOWN;
		/* Return metric */
		V = uf->df_list[crit]->V_base;
		stop = uf->df_list[crit]->tot_cons[0];
		}
	else { // PS
		if (crit != 1)
			return DTL_CRIT_UNKNOWN;
		V = uf->df->V_base;
		stop = uf->df->tot_cons[0];
		}
	/* Return metric */
	for (sum=0,i=1; i<=stop; i++)
		if ((V->lo_midbox[i] > -1.0) || (V->up_midbox[i] > -1.0))
			sum++;
	return sum;
	}


int DTLAPI DTL_nbr_of_weights() {

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Return metric */
	if (PM)
		return uf->df_list[0]->tot_cons[1];
	else // PS
		return 0;
	}


int DTLAPI DTL_nbr_of_crit() {

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Return metric */
	return uf->n_crit;
	}


int DTLAPI DTL_nbr_of_alts() {

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Return metric */
	return uf->n_alts;
	}


int DTLAPI DTL_total_cons(int crit) {
	int i,sum;

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	if (PM) {
		/* Check input parameters */
		if ((crit < -1) || (crit > uf->n_crit))
			return DTL_CRIT_UNKNOWN;
		/* Return metric */
		if (crit > -1) {
			if (uf->df_list[crit])
				return uf->df_list[crit]->n_cons[0];
			else
				return 0;
			}
		else { // sum over all criteria
			for (sum=0,i=1; i<=uf->n_crit; i++)
				if (uf->df_list[i])
					sum += uf->df_list[i]->n_cons[0];
			return sum;
			}
		}
	else // PS
		return crit==1?uf->df->n_cons[0]:DTL_CRIT_UNKNOWN;
	}


int DTLAPI DTL_total_nodes(int crit) {
	int i,sum;

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	if (PM) {
		/* Check input parameters */
		if ((crit < -1) || (crit > uf->n_crit))
			return DTL_CRIT_UNKNOWN;
		/* Return metric */
		if (crit > -1) {
			if (uf->df_list[crit])
				return uf->df_list[crit]->tot_cons[0];
			else
				return 0;
			}
		else { // sum over all criteria
			for (sum=0,i=1; i<=uf->n_crit; i++)
				if (uf->df_list[i])
					sum += uf->df_list[i]->tot_cons[0];
			return sum;
			}
		}
	else // PS
		return crit==1?uf->df->tot_cons[0]:DTL_CRIT_UNKNOWN;
	}


int DTLAPI DTL_nbr_of_cons(int crit, int alt) {

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameters */
	if ((crit < 1) || (crit > uf->n_crit))
		return DTL_CRIT_UNKNOWN;
	if ((alt < 1) || (alt > uf->n_alts))
		return DTL_ALT_UNKNOWN;
	/* Return metric */
	if (PM) {
		if (uf->df_list[crit])
			return uf->df_list[crit]->n_cons[alt];
		else
			return 0;
		}
	else // PS
		return crit==1?uf->df->n_cons[alt]:DTL_CRIT_UNKNOWN;
	}


int DTLAPI DTL_nbr_of_nodes(int crit, int alt) {

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameters */
	if ((crit < 1) || (crit > uf->n_crit))
		return DTL_CRIT_UNKNOWN;
	if ((alt < 1) || (alt > uf->n_alts))
		return DTL_ALT_UNKNOWN;
	/* Return metric */
	if (PM) {
		if (uf->df_list[crit])
			return uf->df_list[crit]->tot_cons[alt];
		else
			return 0;
		}
	else // PS
		return crit==1?uf->df->tot_cons[alt]:DTL_CRIT_UNKNOWN;
	}


int dtl_node2crit(int node) {

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	if (!PM)
		return DTL_WRONG_FRAME_TYPE;
	if (load_df0(0))
		return DTL_SYS_CORRUPT;
	return TCL_get_V_index(uf->df,1,node);
	}


int DTLAPI DTI_node2crit(int node) {

	/* Call stack refurbishing */
	return dtl_node2crit(node);
	}


int dtl_crit2node(int crit) {

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	if (!PM)
		return DTL_WRONG_FRAME_TYPE;
	if (load_df0(0))
		return DTL_SYS_CORRUPT;
	return TCL_get_tot_index(1,crit);
	}


int DTLAPI DTI_crit2node(int node) {

	/* Call stack refurbishing */
	return dtl_crit2node(node);
	}


 /*********************************************************
  *
  *  Interpret (severity of) error codes
  *
  *  There are two types of call information returned:
  *  1. Functions returning DTL_OK on success. In case
  *     of error, a (negative) error code is returned.
  *     This is the most common type of function in DTL.
  *  2. Functions returning a (positive) integer as the
  *     result of the call (sometimes this is even the
  *     main reason for calling). In case of error, an
  *     error number (negative) is overloaded instead.
  *
  *  NOTE: In order to differentiate the two series of
  *  error codes, TCL has its error range as positives.
  *  This way, both series can be handled in parallel.
  *  The TCL codes have DTL_KERNEL_ERROR added to them
  *  when returned to the DTL caller, thus fitting into
  *  the same unified error signalling system from DTL,
  *  with TCL being above DTL_KERNEL_ERROR and DTL below.
  *  Forgiveness is key: allow the user to make mistakes.
  *
  *********************************************************/

int DTLAPI DTL_error2(rcode drc) {

	/* Check if error in return code or just information
	   0 = drc contains information, result valid
	   1 = drc contains information, result invalid
	   2 = drc contains an error, result invalid */
	switch (drc) {
		// codes interpreted as information
		case DTL_DIFFERING_RANKS: // also DTL_LAST_ROW_INCOMPLETE
		case DTL_WEAK_MASS_DISTR: // also DTL_INFINITE_MASS
		case DTL_SCALE_CHANGE:
			return 0;
		case DTL_USER_ABORT:
		case DTL_FILE_UNKNOWN:
		case DTL_KERNEL_ERROR+TCL_TOO_MANY_STMTS:
			return 1;
		// codes interpreted as error if < 0, otherwise as 'ok'
		default:
			return drc<DTL_OK?2:0;
		}
	}


int DTLAPI DTL_u_error2(rcode drc) {

	/* Check if error in return code from numeric user input
	   0 = drc contains information, result valid
	   1 = drc contains information or user mistake, result invalid
	   2 = drc contains an error, result invalid */
	switch (drc) {
		// codes interpreted as user mistake
		case DTL_INCONSISTENT:
		case DTL_KERNEL_ERROR+TCL_INCONSISTENT:
			return 1;
		// codes interpreted as standard DTL error check
		default:
			return DTL_error2(drc);
		}
	}


bool DTLAPI DTL_error(rcode drc) {

	/* Check if error in return code or just information
	   FALSE = drc contains only information
	   TRUE =  drc contains an error */
	return DTL_error2(drc)>1;
	}


bool DTLAPI DTL_u_error(rcode drc) {

	/* Check if error in return code or just information
	   FALSE = drc contains only information or user mistake
	   TRUE =  drc contains an error */
	return DTL_u_error2(drc)>1;
	}
