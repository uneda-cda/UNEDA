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
 *   File: DTLinternal.c
 *
 *   Purpose: Support routines for all procedures in DecisionTree
 *            plus external error texts and system folder functions
 *
 *
 *   Functions exported outside DTL
 *   ------------------------------
 *   DTL_get_errtxt
 *   DTL_get_errtxt_p
 *   DTL_get_errtxt_i/16
 *   DTI_set_folder/16
 *   DTI_reset_folder
 *   DTI_get_folder/16
 *
 *   Functions outside module, inside DTL
 *   ------------------------------------
 *   call
 *   dtl_kernel_error
 *   dtl_error
 *   get_days
 *   imod_
 *   fmod_
 *   trunc
 *   f_NaN
 *   trace_ifperr
 *   cst_open
 *   cst_close
 *   cst_log
 *   cst_trace
 *   load_PV_stmt
 *   load_W_stmt
 *   new_uf
 *   get_uf
 *   dispose_uf
 *   load_df0
 *   check_df0
 *   load_df00
 *   load_df1
 *   check_df1
 *
 *   Functions internal to module
 *   ----------------------------
 *   get_date
 *   check_bounds
 *   load_df
 *
 *
 *   DTL layer 2: most functions are below API level
 *
 */

#include "DTL.h"
#include "DTLinternal.h"
#include "CAR.h"


 /********************************************************
  *
  *  Configuration
  *
  ********************************************************/

/* Select method to detect NaN conditions */
#ifdef _MSC_VER
#define range_NaN // C99 can handle IND and INF ranges
#endif

/* File permission bit field (read, write, exec) */
#ifndef S_IRWXU
#define S_IRWXU 0700
#endif


 /********************************************************
  *
  *  Global data (internal to DTL, not module)
  *
  ********************************************************/

int dtl_abort_request;
int dtl_abort_cst;
char *dtl_func = "NULL";
struct user_frame *uf = NULL;
struct user_frame *uf_list[MAX_FRAMES+1];
char strg[64],msg[200]; // oversized to prevent 64-bit %d & %p buffer overruns
int frame_loaded = 0;
int dtl_error_count;
int dtl_trace_count;
int smx_busy = FALSE;
jmp_buf assert_envir;
int a_tag;
char dtl_folder[DTLF_SIZE] = "";
char dtl_buffer[1000];


 /********************************************************
  *
  *  Internal data (within module)
  *
  ********************************************************/

#define FRESCATI_MOVE 5713 // move Aug 30, 2027 according to plan

/* General load status */
static rcode latest_kernel_rc = TCL_OK;


 /*********************************************************
  *
  *  Error handling - DTL error texts are in caps
  *
  *  (should reside in DTLmisc.c, not DTLinternal.c)
  *
  *********************************************************/

#define MAX_ET 30 // max errtxt size (must be less than 254)

static char *DTL_errtxt[] = {
	"DTL OK",
	"KERNEL ERROR",
	"INPUT ERROR",
	"TREE ERROR",
	"OUTPUT ERROR",
	"FRAME EXISTS",
	"FRAME UNKNOWN",
	"FRAME IN USE",
	"FRAME NOT LOADED",
	"FRAME CORRUPT",
	"WRONG FRAME TYPE",
	"WRONG STATEMENT TYPE",
	"TOO MANY CONSEQUENCES",
	"TOO MANY CRITERIA",
	"LOG FILE ERROR",
	"INCONSISTENT",
	"DIFFERING RANKS",
	"STATEMENT ERROR",
	"SYSTEM CORRUPT",
	"TOO MANY ALTERNATIVES",
	"TOO MANY NODES IN TREE",
	"CRITERION MISSING",
	"TOO FEW ALTERNATIVES",
	"USER ABORT",
	"STATE ERROR",
	"CRITERION UNKNOWN",
	"CRITERION EXISTS",
	"ALTERNATIVE UNKNOWN",
	"ALTERNATIVE MISMATCH",
	"DTL BUSY",
	"NAME MISSING",
	"NAME TOO LONG",
	"NAME EXISTS",
	"NOT ALLOWED",
	"WRONG METHOD",
	"WRONG TOLERANCE",
	"FILE/FOLDER UNKNOWN",
	"SCALE CHANGE",
	"INTERNAL ERROR",
	"WEAK MASS DISTRIBUTION",
	"MEMORY LEAK",
	"BUFFER OVERRUN",
	"ASSERT FAILED",
	"- RETURN CODE OUT OF RANGE -"
	};


 /* Get an error text for simple user feedback.
  * Handle both levels (DTL + TCL) in one call. */

char* DTLAPI DTL_get_errtxt(rcode drc) {

	if (drc>=0)
		return DTL_errtxt[0]; // DTL OK
	else if ((drc < DTL_KERNEL_ERROR) && (drc >= MAX_DTL_ERR))
		return DTL_errtxt[-drc-99]; // DTL error text
	else if (drc == DTL_KERNEL_ERROR)
		return "UNSPECIFIED KERNEL ERROR"; // kernel error text
	else if (drc > DTL_KERNEL_ERROR)
		return TCL_get_errtxt(drc-DTL_KERNEL_ERROR); // TCL error text
	else
		return DTL_errtxt[-MAX_DTL_ERR-98]; // DTL rcode out of range
	}


static char errtxt_C2P[MAX_ET+1];

char* DTLAPI DTL_get_errtxt_p(rcode drc) {
	unsigned i,et_size;

	/* Convert C error message (null-terminated) to Pascal shortstring */
	et_size = (unsigned)strlen(DTL_get_errtxt(drc));
	if (et_size > MAX_ET) { // internal error, should not occur
		sprintf(msg," DTL_get_errtxt_p(%d) buffer too short (min %u)\n",drc,et_size+1);
		cst_trace(msg);
		return ""; // equals shortstring with size 0
		}
	/* Safe to convert to length-preceded string */
	strcpy(errtxt_C2P,DTL_get_errtxt(drc));
	for (i=et_size; i; i--)
		errtxt_C2P[i] = errtxt_C2P[i-1];
	errtxt_C2P[0] = (unsigned char)et_size;
	return errtxt_C2P; // not guarded by smx semaphore
	}


rcode DTLAPI DTL_get_errtxt_i(rcode drc, char *str, unsigned *len) {
	unsigned et_size;

	/* Move error message to in-situ string (pointer originates at caller) */
	et_size = (unsigned)strlen(DTL_get_errtxt(drc));
	/* Incoming mem 'len' is sizeof, outgoing content size 'len' is excl '\0' */
	if (et_size >= *len) { // receiving string too short
		sprintf(msg," DTL_get_errtxt_i(%d) buffer size (%u) too small (min %u)\n",drc,*len,et_size+1);
		cst_trace(msg);
		*len = 0U;
		return DTL_BUFFER_OVERRUN;
		}
	/* Safe to copy string (size fits) */
	strcpy(str,DTL_get_errtxt(drc));
	*len = et_size;
	return DTL_OK;
	}


rcode DTLAPI DTL_get_errtxt_i16(rcode drc, char *str, unsigned *len, bool LE) {
	int i;

	/* Incoming string length is in 16-bit chars but used as count of
	 * 8-bit bytes in DTL call to reserve for two-fold inflation.
	 * This enables the use of 8-bit overrun tests also for 16-bit.
	 * LE flag indicates if the architecture is little- or big-endian. */
	if (DTL_get_errtxt_i(drc,str,len) < DTL_OK)
		return DTL_BUFFER_OVERRUN;
	/* Inflate error message to 16-bit in-situ */
	for (i=(int)*len; i>=0; i--) // incl null terminator
		if (LE) { // little-endian
			str[2*i] = str[i];
			str[2*i+1] = '\0';
			}
		else  { // big-endian
			str[2*i] = '\0';
			str[2*i+1] = str[i];
			}
	/* Length returned in chars, not bytes */
	return DTL_OK;
	}


 /*********************************************************
  *
  *  Post-processing of TCL kernel calls
  *
  *********************************************************/

rcode call(rcode rc, char *proc) {

	latest_kernel_rc = rc;
	if (rc)
		sprintf(msg," %.30s: %.30s\n",proc,TCL_get_errtxt(rc));
	else
		sprintf(msg," %.30s: ok\n",proc);
	if (cst_on)
		cst_log(msg);
	else if (rc && (rc != TCL_INCONSISTENT)) {
		dtl_error_count++;
		cst_trace(msg);
		}
	return rc;
	}


 /*********************************************************
  *
  *  DTL error handling (flush cache + end smx semaphore)
  *
  *********************************************************/

rcode dtl_kernel_error() {

	eval_cache_invalidate();
	_smx_end();
	return DTL_KERNEL_ERROR+latest_kernel_rc;
	}


/*
 * dtl_error: Log errors on file that are of a system nature, not user-inflicted.
 * E.g. an abort is not an error in itself but rather a user-initiated action result.
 */

rcode dtl_error(rcode drc) {

	eval_cache_invalidate();
	if (strcmp(dtl_func,"NULL")) // first log entry for an error occurrence
		sprintf(msg," DTL error: %.30s in DTL::%.10s\n",DTL_get_errtxt(drc),dtl_func);
	else // subsequent log entry for same error (from e.g. evaluate_frame)
		sprintf(msg," DTL clean-up after %.30s\n",DTL_get_errtxt(drc));
	/* Log msg on cst log file if open, else - only if error - on trace log file.
	 * Step error counter only if in unattended mode, not in development/debug mode. */
	if (cst_on && (drc<DTL_KERNEL_ERROR) && (drc!=DTL_USER_ABORT))
		cst_log(msg);
	else if (dtl_is_init() && DTL_u_error(drc)) { // NOT test/debug = runtime error
		dtl_error_count++;
		cst_trace(msg);
		}
	_smx_end();
	return drc;
	}


static void get_date(char *strg) {
	time_t now;
	struct tm *date;

	time(&now);
	date = localtime(&now);
	strftime(strg,32,"%Y-%b-%d %H:%M",date);
	}


 /*********************************************************
  *
  *  Math error handling
  *
  *********************************************************/

/* Guarded math operators: defined as macros in DTLinternal.h.
 * Invoked functions reside here: integer mod and double mod.
 * Plus also NaN/IND/INF check and int/float error handler.

 * Note that there are at least four definitions of rem/mod:
 * truncated (%), floored (mod), Euclidean & rounded (Lisp).
 * Truncated and floored coincide if the dividend and divisor
 * signs are the same. Euclidean disregards the divisor sign.
 * While rem is machine hacking and mod is software hacking,
 * the most satisfying is Euclidean with its clock property
 * and its connections to set and group theories in algebra.
 * But staying with software conventions, it is floored here.
 * Note that mod loses its cyclical property for divisor = 0.
 *
 * Reference paper from MS on their implementation of modulus arithmetic:
 * microsoft.com/en-us/research/wp-content/uploads/2016/02/divmodnote.pdf
 * But the paper fails to properly take negative numbers into account. */

int ipow(int x, int y) {

	return (int)(pow(x,y)+1.0E-8);
	}


int imod_(int x, int y) {
	int m = x%y; // y!=0 guaranteed by macro

	// keep costly ops down
	if (y==-1)
		// for INT_MIN % -1
		return 0;
	if (m) // ((x%y)+y)%y costs more
		if (x>0)
			m += min(y,0);
		else
			m += max(y,0);
	return m;
	}


double fmod_(double x, double y) {
	double m = fmod(x,y); // y!=0.0 guaranteed by macro

	// keep costly ops down
	if (m!=0.0)
		if (x>0.0)
			m += min(y,0.0);
		else
			m += max(y,0.0);
	return m;
	}


#ifndef _trunc
double trunc(double d) {

	/* trunc is in C99 and on (not in C89) but MS VC++ doesn't have it yet */
	return d>0.0?floor(d):ceil(d);
	}
#endif


/* Checking for invalid numbers of type NaN (IND or INF). There are various
 * representations and behaviours depending on compiler and implementation.
 * For example, x != x is supposed to fail for NaN but does not on some of
 * the MS VC++ platforms -> very unreliable. Range tests with max values
 * seem to be more reliable and are thus preferred in DTL. Since MS VC++
 * does not have working built-in isnan or isinf, use only NaN detection.
 * If ranges are not working, undefine range_NaN and use bitmask instead. */

#ifdef range_NaN // use range NaN check (preferred but library dependent)

int f_NaN(double test_num) {

	/* Catch both NaN and IND/INF: If IND/INF, one of the two tests will
	 * trigger the clause depending on positive or negative infinity. On
	 * the other hand, if NaN then both tests will trigger the clause. */
	return test_num > DBL_MAX || test_num < -DBL_MAX;
	}

#else // use bitmask NaN check instead (IEEE 754 dependent)

int f_NaN(double test_num) {

	/* Explicit typecast to double to know the size (assuming IEEE 754 standard).
	 * This is truly portable since the sizes of both data types are corrected for,
	 * even though some compilers don't like it and complain about non-portability. */
	unsigned mask = *((unsigned *)(&test_num)+sizeof(test_num)/sizeof(mask)-1);

	// handles both big-endian and little-endian hardware architectures
	return (mask & 0x7FF00000) == 0x7FF00000;
	}

#endif


/* Integer and floating point (ipf) error handler. NOTE: the replacement
 * value is defined as int which influences how _fmod & _gmod are defined. */

int trace_ifperr(int type, int rv) {
	char tmsg[36]; // buf size incl. %d (32 or 64 bit = max 19 digits + sign)

	switch (type) {
		case 0:  sprintf(tmsg,"div by zero (%d)",rv);  break;
		case 1:  sprintf(tmsg,"rem by zero (%d)",rv);  break;
		case 2:  sprintf(tmsg,"mod by zero (%d)",rv);  break;
		case 3:  sprintf(tmsg,"abs on NaN (IND/INF)"); break;
		case 4:  sprintf(tmsg,"sgn on NaN (IND/INF)"); break;
		case 5:  sprintf(tmsg,"sqrt on neg num");      break;
		default: sprintf(tmsg,"unknown ifperr (%d)",type);
		}
	sprintf(msg," ifperr: %.36s in DTL::%.10s\n",tmsg,dtl_func);
	cst_trace(msg); // log to trace file
	return rv; // replacement value
	}


/* The SU DTL/TCL research platform was first launched on Jan 8, 2012.
 * This day, the first CAR algorithms (MCR+CRC ver 0.01) were released,
 * having been prepared during the latter part of 2011. This marks the
 * start of the modern computational decision analysis research at SU. */

#define BOOT_DAY 15347 // Jan 8, 2012 (relative to Jan 1, 1970 = Unix day 1)
#define SEC_PER_DAY 86400

int get_days() {
	time_t now;

	time(&now);
	return (int)(now/SEC_PER_DAY)-BOOT_DAY;
	}


 /*********************************************************
  *
  *  Call sequence trace log
  *
  *********************************************************/

#define TRACE_FNAME "trace.log"  // error call tracing
#define CST_FNAME "call_seq.log" // general call logging
#define CST_FNAME_EXT "call_seq_ext.log" // serious error logging
int cst_on=FALSE;
int cst_ext=FALSE;
FILE* cst;
FILE* dtr;


/* Set system folder (and create if necessary)
 * Style 0 -> C-style ANSI string as input (null-terminated)
 * Style 1 -> Pascal-style string as input (length-preceded)
 * NOTE: This function uses the (Object) Pascal type shortstring
 *       which is the original string type in old Pascal as well,
 *       not least back in 1994 when this package was concieved.
 *       The string will be converted in situ to native C style. */

#define C_STYLE 0

rcode DTLAPI DTI_set_folder(char *folder, int style) {
	unsigned i,len;

	/* Cannot log its own invocation since this call changes
	 * the system folder in which all the log files reside
	 * and a change would pull the rug from under their feet.
	 * Thus, the call is only permitted prior to DTL open. */

	/* Check if function can start */
	_certify_ptr(folder,1);
	if (dtl_is_init())
		return DTL_STATE_ERROR;
	if (style) { // Pascal-style shortstring -> convert to C
		len = (unsigned)folder[0]; // one-byte length field
		if (len > DTLF_SIZE-2)
			return DTL_BUFFER_OVERRUN;
		for (i=0U; i<len; i++) {
			if (folder[i+1] < 0x20) // unprintable char - no good
				return DTL_INPUT_ERROR;
			folder[i] = folder[i+1]; // good char - transfer to C
			}
		folder[len] = '\0'; // null-terminate new C string
		}
	else // already C-style
		if (strlen(folder) > DTLF_SIZE-2) // also room for '/' + '\0'
			return DTL_BUFFER_OVERRUN;
	if (folder[0]) {
		/* Away folder */
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4013) // "assuming extern returning int"
		if (_mkdir(folder))
#pragma warning(pop)
#else // Unix/Linux, Mac
		if (mkdir(folder,S_IRWXU))
#endif
			if (errno != FILE_EXISTS)
				return DTL_FILE_UNKNOWN; // piggyback folder error
		strcpy(dtl_folder,folder); // unguarded (known length)
		/* Must end in slash to be a proper folder path */
		if (dtl_folder[strlen(dtl_folder)-1] != '/')
			strcat(dtl_folder,"/");
		}
	else
		/* Home folder - circumvent path settings */
		strcpy (dtl_folder,"./");
	return (rcode)strlen(dtl_folder); // for string debugging
	}


/* Corresponds to strlen for 2-byte chars */

static unsigned strlen16(char* strg) {
	unsigned b_len=0U; // byte length

	while (max(strg[b_len],strg[b_len+1])) // both endians work
		b_len += 2;
	return b_len/2; // return chars, not bytes
	}


static char folder8[DTLF_SIZE]; // C-style string

rcode DTLAPI DTI_set_folder16(char *folder16) {
	unsigned i,len;

	/* Check if function can start */
	_certify_ptr(folder16,1);
	/* Check incoming length (in chars, not bytes) */
	len = strlen16(folder16);
	if (len > DTLF_SIZE-2)
		return DTL_BUFFER_OVERRUN;
	/* Convert 2-byte to 1-byte chars regardless of endian */
	for (i=0U; i<len; i++)
		folder8[i] = (char)max(folder16[2*i],folder16[2*i+1]);
	folder8[len] = '\0'; // null-terminate new C string
	return DTI_set_folder(folder8,C_STYLE);
	}


rcode DTLAPI DTI_reset_folder() {

	return min(DTI_set_folder("./",0),DTL_OK);
	}


rcode DTLAPI DTI_get_folder(char *folder, unsigned *c_size, int style) {
	unsigned i,fo_size;

	/* Permitted in both DTL open and closed states */
	_certify_ptr(folder,1);
	_certify_ptr(c_size,2);
	/* Log function call. NOTE: If this is called after DTI_set_folder, make
	 * sure there resides an appropriate cst_log file in the new system folder
	 * if the call sequence (cst) log is enabled at DTL_init file opening.
	 * The basic trace file trace.log is created automatically by DTL_init. */
	if (cst_ext) {
		sprintf(msg,"DTI_get_folder(strg[%u],%d)\n",*c_size,style);
		cst_log(msg);
		}
	fo_size = (unsigned)strlen(dtl_folder);
	/* c_size input is sizeof, output is content excl '\0' */
	if (fo_size >= *c_size) {
		/* Folder name does not fit into receiving string */
		sprintf(msg," DTI_get_folder(strg[%u],%d) buffer too short (min %u)\n",*c_size,style,fo_size+1);
		cst_trace(msg);
		*c_size = 0U;
		return DTL_BUFFER_OVERRUN;
		}
	if (style) {
		/* Convert C string (null-terminated) to Pascal shortstring */
		if (fo_size > 0x0FF) {
			/* Length does not fit into one byte */
			*c_size = 0U;
			return DTL_BUFFER_OVERRUN;
			}
		/* Safe to convert to length-preceded shortstring */
		for (i=fo_size; i; i--)
			folder[i] = dtl_folder[i-1];
		folder[0] = (unsigned char)fo_size;
		}
	else {
		/* Update string in-situ */
		strncpy(folder,dtl_folder,*c_size-1);
		folder[*c_size-1] = '\0'; // C-style end of string
		}
	*c_size = fo_size;
	/* Log function result (again, given an appropriate log file) */
	if (cst_ext) {
		sprintf(msg,"  System folder is \"%.40s\" (%u)\n",dtl_folder,fo_size);
		cst_log(msg);
		}
	return DTL_OK;
	}


rcode DTLAPI DTI_get_folder16(char *folder, unsigned *c_size, bool LE) {
	int i;

	/* Incoming string size is in 16-bit chars but use as count of
	 * 8-bit bytes in DTL call to reserve for inflation in copying.
	 * This enables the use of 8-bit overrun tests also for 16-bit. */
	if (DTI_get_folder(folder,c_size,C_STYLE) < DTL_OK)
		return DTL_BUFFER_OVERRUN;
	/* Inflate folder name to 16-bit in-situ */
	for (i=(int)*c_size; i>=0; i--) // incl null terminator
		if (LE) { // little-endian
			folder[2*i] = folder[i];
			folder[2*i+1] = '\0';
			}
		else  { // big-endian
			folder[2*i] = '\0';
			folder[2*i+1] = folder[i];
			}
	/* Size returned in chars, not bytes */
	return DTL_OK;
	}


/* cst_log levels:
 *
 * Level 0 (OFF):     no log data written
 * Level 1 (cst_on):  input data + execution status
 * Level 2 (cst_ext): level 1 + output data
 *
 * Trace file silently opened in the background */

rcode cst_open() {
	char otxt[64],fname[DTLF_SIZE+64];
	int ext_tmp=FALSE;

	if (cst_on) // not first time, already open
		return DTL_OK;
	/* Open trace file first */
	strcpy(fname,dtl_folder);
	strcat(fname,TRACE_FNAME);
	if ((dtr = fopen(fname,"a")) == NULL)
		return DTL_INTERNAL_ERROR;
	/* Next, open log file */
	strcpy(fname,dtl_folder);
	strcat(fname,CST_FNAME);
	if ((cst = fopen(fname,"r")) == NULL)
		return DTL_LOGFILE_ERROR;
	fscanf(cst," %255s",fname);
	fclose(cst);
	/* Set cst_log type */
	if (strcmp(fname,CST_FNAME))
		if (strcmp(fname,CST_FNAME_EXT))
			return DTL_LOGFILE_ERROR;
		else
			ext_tmp = TRUE;
	strcpy(fname,dtl_folder);
	strcat(fname,CST_FNAME);
	if ((cst = fopen(fname,"a")) == NULL)
		return DTL_LOGFILE_ERROR;
	// cst_log is accepted
	if (ext_tmp) {
		fprintf(cst,"\n\n");
		fprintf(cst,"\t\t       _/_/_/      _/_/_/_/_/   _/\n");
		fprintf(cst,"\t\t      _/    _/        _/       _/\n");
		fprintf(cst,"\t\t     _/      _/      _/       _/\n");
		fprintf(cst,"\t\t    _/      _/      _/       _/\n");
		fprintf(cst,"\t\t   _/      _/      _/       _/\n");
		fprintf(cst,"\t\t  _/      _/      _/       _/\n");
		fprintf(cst,"\t\t _/     _/       _/       _/\n");
		fprintf(cst,"\t\t_/_/_/_/        _/       _/_/_/_/_/\n\n\n");
		fprintf(cst,"\t    (+)\n");
		fprintf(cst,"     +--- o &&& o --------------------------------------------+\n");
		fprintf(cst,"     |  o         o            Prof. Mats Danielson           |\n");
		fprintf(cst,"     | o   STHLM   o           DECIDE Research Group          |\n");
		fprintf(cst,"     | o           o  Dept. of Computer and Systems Sciences  |\n");
		fprintf(cst,"     | o    UNI    o           Stockholm University           |\n");
		if (get_days() > FRESCATI_MOVE) // returning to main campus Frescati, Building F
		fprintf(cst,"     |  o         o         SE-106 91 Stockholm, SWEDEN       |\n");
		else // after being in exile for 36.5 years (1991-2027) in Kista, 11 kilometres north
		fprintf(cst,"     |  o         o    PO Box 1073, SE-164 25 Kista, SWEDEN   |\n");
		fprintf(cst,"     +---- o x o ---------------------------------------------+\n\n");
		}
	fprintf(cst,"\n------CSTINIT------\n");
	get_date(otxt);
	fprintf(cst,otxt);
	fprintf(cst,"\n-------------------\n");
	sprintf(otxt,"DTL release %d.%02d.%d\n",DTL_MAIN,DTL_FUNC,DTL_TECH);
	fprintf(cst,otxt);
	sprintf(otxt,"Existed %d days\n",get_days());
	fprintf(cst,otxt);
#ifdef _MSC_VER
	sprintf(otxt,"%sVisual C++ %d.%d\n",_MSC_VER<2600?"MS ":"",_MSC_VER/100-(_MSC_VER<1900?6:5),_MSC_VER%100);
	fprintf(cst,otxt);
#endif
	fprintf(cst,"------CSTLOG-------\n");
	fflush(cst);
	cst_on = TRUE;
	cst_ext = ext_tmp;
	return DTL_OK;
	}


void cst_close() {
	char date[64];

	/* Close cst_log if open */
	if (cst_on)	{
		fprintf(cst,"------CSTEXIT------\n");
		get_date(date);
		fprintf(cst,date);
		fprintf(cst,"\n-------------------\n");
		fflush(cst);
		fclose(cst);
		cst_on = cst_ext = FALSE;
		}
	/* Close silent trace log if run-time system
	 * does not force it to close prematurely */
	if (dtl_trace_count)
		cst_trace(" Closing DTL trace file\n");
	else
		fflush(dtr);
	fclose(dtr);
	}


void cst_log(char *msg) {

	if (cst_on)	{
		fprintf(cst,msg);
		fflush(cst);
		}
	}


void cst_trace(char *msg) {
	char date[64];

	/* The trace is written regardless of cst state */
	get_date(date);
	fprintf(dtr,date);
	fprintf(dtr,msg);
	fflush(dtr);
	dtl_trace_count++;
	}


void DTLAPI DTI_raw_trace(char *msg) {

	/* The msg is written regardless of cst state */
	fprintf(dtr,msg);
	fflush(dtr);
	}


 /*********************************************************
  *
  *  User statement handling
  *
  *********************************************************/

 /*
  * check_bounds: check bounds for input statement
  */

static rcode check_bounds(struct stmt_rec *stmt, double lolim) {

	if ((stmt->lobo < lolim) || (stmt->lobo > 1.0)) {
		if (cst_on) {
			sprintf(msg," lower bound check (%.10le)\n",stmt->lobo);
			cst_log(msg);
			}
		return DTL_INPUT_ERROR;
		}
	if ((stmt->upbo < 0.0) || (stmt->upbo > 1.0)) {
		if (cst_on) {
			sprintf(msg," upper bound check (%.10le)\n",stmt->upbo);
			cst_log(msg);
			}
		return DTL_INPUT_ERROR;
		}
	if (stmt->upbo < stmt->lobo) {
		if (cst_on) {
			sprintf(msg," bounds cross check [%.10le %.10le]\n",
					stmt->lobo,stmt->upbo);
			cst_log(msg);
			}
		return DTL_INPUT_ERROR;
		}
	return DTL_OK;
	}


 /*
  * load_stmt: move user prob or util statement to internal format
  */

rcode load_PV_stmt(int crit, struct user_stmt_rec *ustmt, struct stmt_rec *stmt, char type) {
	int i;

	if (cst_on) {
		cst_log(" stmt: ");
		for (i=1; i<=ustmt->n_terms; i++) {
			if (PS)
				sprintf(msg,"%c%c%d.%d ",ustmt->sign[i]==1?'+':'-',type,
						ustmt->alt[i],ustmt->cons[i]);
			else if (PM)
				sprintf(msg,"%c%c%d.%d.%d ",ustmt->sign[i]==1?'+':'-',type,
						crit,ustmt->alt[i],ustmt->cons[i]);
			else /* error */
				sprintf(msg,"ERROR ");
			cst_log(msg);
			}
		cst_log("= ");
		if (((ustmt->lobo < 0.00) && (ustmt->lobo > -0.01)) ||
				((ustmt->upbo < 0.00) && (ustmt->upbo > -0.01)) ||
				((ustmt->lobo > 0.00) && (ustmt->lobo < 0.01)) ||
				((ustmt->upbo > 0.00) && (ustmt->upbo < 0.01)))
			sprintf(msg,"[%.3le %.3le] (%le)\n",ustmt->lobo,
					ustmt->upbo,(ustmt->upbo-ustmt->lobo));
		else
			sprintf(msg,"[%.4lf %.4lf] (%le)\n",ustmt->lobo,
					ustmt->upbo,(ustmt->upbo-ustmt->lobo));
		cst_log(msg);
		}
	if (PS && (crit != 1))
		return DTL_CRIT_UNKNOWN;
	stmt->n_terms = ustmt->n_terms;
	for (i=1; i<=ustmt->n_terms; i++)	{
		stmt->alt[i] = ustmt->alt[i];
		stmt->cons[i] = ustmt->cons[i];
		if ((ustmt->sign[i] != -1) && (ustmt->sign[i] != 1))
			return DTL_INPUT_ERROR;
		stmt->sign[i] = ustmt->sign[i];
		}
	stmt->lobo = ustmt->lobo;
	stmt->upbo = ustmt->upbo;
	if (check_bounds(stmt,ustmt->n_terms==1?0.0:-1.0))
		return DTL_INPUT_ERROR;
	return DTL_OK;
	}


 /*
  * load_W_stmt: move user weight statement to internal format
  */

rcode load_W_stmt(struct user_w_stmt_rec *ustmt, struct stmt_rec *stmt) {
	int i;

	if (cst_on) {
		cst_log(" stmt: ");
		for (i=1; i<=ustmt->n_terms; i++) {
			sprintf(msg,"%cW%d ",ustmt->sign[i]==1?'+':'-',ustmt->crit[i]);
			cst_log(msg);
			}
		cst_log("= ");
		if (((ustmt->lobo < 0.00) && (ustmt->lobo > -0.01)) ||
				((ustmt->upbo < 0.00) && (ustmt->upbo > -0.01)) ||
				((ustmt->lobo > 0.00) && (ustmt->lobo < 0.01)) ||
				((ustmt->upbo > 0.00) && (ustmt->upbo < 0.01)))
			sprintf(msg,"[%.3le %.3le] (%le)\n",ustmt->lobo,
					ustmt->upbo,(ustmt->upbo-ustmt->lobo));
		else
			sprintf(msg,"[%.4lf %.4lf] (%le)\n",ustmt->lobo,
					ustmt->upbo,(ustmt->upbo-ustmt->lobo));
		cst_log(msg);
		}
	stmt->n_terms = ustmt->n_terms;
	for (i=1; i<=ustmt->n_terms; i++) {
		stmt->alt[i] = 1;
		stmt->cons[i] = ustmt->crit[i];
		if ((ustmt->sign[i] != -1) && (ustmt->sign[i] != 1))
			return DTL_INPUT_ERROR;
		stmt->sign[i] = ustmt->sign[i];
		}
	stmt->lobo = ustmt->lobo;
	stmt->upbo = ustmt->upbo;
	if (check_bounds(stmt,ustmt->n_terms==1?0.0:-1.0))
		return DTL_INPUT_ERROR;
	return DTL_OK;
	}


 /*********************************************************
  *
  *  User frame list management
  *
  *********************************************************/

 /*
  * new_uf: create new user frame and enter into uf_list.
  * Returns pointer (NULL on failure).
  */

struct user_frame *new_uf(int ufnr) {
	int j;

	if (uf_list[ufnr])
		return NULL;
	/* Create a new user frame */
	uf_list[ufnr] = (struct user_frame *)mem_alloc(sizeof(struct user_frame),
			 "struct user_frame","new_uf");
	/* Init all fields in new user frame */
	uf_list[ufnr]->frame_type = ANY_FRAME;
	uf_list[ufnr]->frame_nbr = 0;
	uf_list[ufnr]->frame_name[0] = '\0';
	uf_list[ufnr]->n_alts = 0;
	uf_list[ufnr]->n_crit = 0;
	uf_list[ufnr]->load_crit = -1;
	uf_list[ufnr]->df = NULL;
	for (j=0; j<=MAX_CRIT; j++) {
		uf_list[ufnr]->df_list[j] = NULL;
		uf_list[ufnr]->WP_autogen[j] = FALSE;
		uf_list[ufnr]->V_n_rels[j] = 0;
		uf_list[ufnr]->av_min[j] = 0.0;
		uf_list[ufnr]->av_max[j] = 1.0;
		}
	return uf_list[ufnr];
	}


 /*
  * get_uf: fetch user frame pointer.
  * Returns pointer or NULL on failure.
  */

struct user_frame *get_uf(int index) {

	if ((index < 1) || (index > MAX_FRAMES))
		return NULL;
	else
		return uf_list[index];
	}


 /*
  * dispose_uf: delete user frame and remove from uf_list.
  * Returns list position or 0 on failure.
  */

int dispose_uf(int index) {

	/* Check index bounds */
	if ((index < 1) || (index > MAX_FRAMES))
		return 0;
	/* Check if frame exists */
	if (uf_list[index] == NULL)
		return 0;
	mem_free((void *)uf_list[index]);
	uf_list[index] = NULL;
	return index;
	}


 /*
  * load_df: load df corresponding to crit
  *
  * Accepts:
  *
  * PS  1
  * PM  0..n_crit
  */

static rcode load_df(int crit) {

	if (uf == NULL) // guard pointer
		return DTL_FRAME_NOT_LOADED;
	if ((crit < 0) || (crit > uf->n_crit))
		return DTL_CRIT_UNKNOWN;
	if (PM) {
		/* See if crit exists */
		if (uf->df_list[crit] == NULL)
			return DTL_CRIT_UNKNOWN;
		/* See if already loaded */
		if (crit == uf->load_crit)
			return DTL_OK;
		/* Unload old if necessary */
		if (uf->load_crit >= 0) {
			if (call(TCL_detach_frame(uf->df),"TCL_detach_frame"))
				return DTL_FRAME_CORRUPT;
			uf->df = NULL;
			uf->load_crit = -1;
			}
		/* Load new crit */
		if (frame_loaded) {
			if (call(TCL_attach_frame(uf->df_list[crit]),"TCL_attach_frame"))
				return DTL_FRAME_CORRUPT;
			}
		else
			return DTL_FRAME_NOT_LOADED;
		uf->df = uf->df_list[crit];
		uf->load_crit = crit;
		}
	else if (PS && crit!=1)
		return DTL_CRIT_UNKNOWN;
	return DTL_OK;
	}


 /*
  * load_df0: load df corresponding to crit
  *
  * For evaluation and W-base input (w/o PS)
  *
  * Accepts:
  *
  * PS  1
  * PM  0..n_crit
  */

rcode load_df0(int crit) {

	return load_df(crit);
	}


rcode check_df0(int crit) {

	if (uf == NULL) // guard pointer
		return DTL_FRAME_NOT_LOADED;
	if (PS && (crit==1))
		return DTL_OK;
	if (PM) {
		if ((crit < 0) || (crit > uf->n_crit))
			return DTL_CRIT_UNKNOWN;
		if (uf->df_list[crit])
			return DTL_OK;
		}
	return DTL_CRIT_UNKNOWN;
	}


 /*
  * load_df00: load df corresponding to crit OR
  * accept negative number for partial MC eval
  */

rcode load_df00(int crit) {
	rcode rc;

	if (crit>=0)
		return load_df0(crit);
	else {
		if (rc = load_df(0))
			return rc;
		if (-crit > uf->df->tot_cons[1])
			return DTL_CRIT_UNKNOWN;
		if (TCL_get_P_index(uf->df,1,-crit))
			return DTL_CRIT_UNKNOWN;
		return DTL_OK;
		}
	}


 /*
  * load_df1: load df corresponding to crit
  *
  * For P- and V-base input
  *
  * Accepts:
  *
  * PS  1
  * PM  1..n_crit
  */

rcode load_df1(int crit) {

	if (!crit)
		return DTL_CRIT_UNKNOWN;
	return load_df(crit);
	}


rcode check_df1(int crit) {

	if (uf == NULL) // guard pointer
		return DTL_FRAME_NOT_LOADED;
	if ((crit < 1) || (crit > uf->n_crit))
		return DTL_CRIT_UNKNOWN;
	if (PM && (uf->df_list[crit] == NULL))
		return DTL_CRIT_UNKNOWN;
	return DTL_OK;
	}


rcode dtl_is_shadow_crit(int crit) {

	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	if (!PM)
		return DTL_WRONG_FRAME_TYPE;
	return load_df0(crit);
	}


/* Plan: graceful degradation in four steps. No remedy is available. */

void shut_down(int tag) {

	/* Step 1: Write log entry */
	if (cst_on) {
		cst_log(msg);
		cst_log(" --> EMERGENCY SHUTDOWN\n");
		}
	else
		cst_trace(msg);
	/* Step 2: Disengage TCL (for step 3 to work properly) */
	DTL_unload_frame();
	/* Step 3: Shut down DTL entirely */
	CAR_exit();
	DTL_exit();
	/* Step 4: Shout out on stderr (hoping there is someone listening) */
#ifndef _MSC_VER
	fprintf(stderr,"Aborted\n"); // Unix csh/bash-style report
#endif
	fprintf(stderr,"DTL core process aborted: tag=%03d (see logfile)\n",tag);
	fflush(stderr); // not all runtime environments autoflush at exit
	/* Caller is expected to terminate process after its final actions */
	}


 /**************************************************************
  *
  *  Fortify against faulty calls by pointer safeguard handling.
  *  Developement started in 1994 using MPW on a Powerbook 170.
  *  Has to be reasonably portable across Windows, Unix and Mac
  *  platforms from 32-bit Win 95, old MacOS and SVR3 to Win NT
  *  servers to modern 64-bit Win 11, Debian Linux and all else.
  *
  **************************************************************/

#define LO_MEM 0x00010000UL // lower limit for useful address

void handle_illegal_ptr(void *ptr, int tag) {

	/* Accessing the violating pointer will only result in an interrupt.
	 * This code is the guard against that interrupt because the powers
	 * of interrupt handlers are to some degree undefined. Alternatively,
	 * we could just wait for the interrupt to occur, but that is not so
	 * portable across platforms. Beware of interrupt handlers' low power.
	 * (Tags in the range 011-019 are indeed uncaught MMU interrupts.)
	 *
	 * Trigger tag ranges
	 * ------------------
	 * 001-009 DTL proper
	 * 011-019 MMU
	 * 101-109 CAR
	 * 201-209 Autoscale
	 * 301-309 SML/SXL/MML
	 * 401-409 BRS
	 * 901-909 Intro */

	/* Bug in MSVC++ 6.0: %p prints A-F as if it was %P */
#ifdef _MSC_VER // bug in MSVC++ 6.0 for _printf: doesn't print 0x if zero
	sprintf(msg," Address exception in DTL::%.10s at %s%p:%03d\n",dtl_func,ptr?"":"0x",ptr,tag);
#else // Unix/Linux/Mac gcc prints "(nil)" for a null pointer
	sprintf(msg," Address exception in DTL::%.10s at %p:%03d\n",dtl_func,ptr,tag);
#endif
	shut_down(tag);
	}


/* Developed for 32-bit MacOS machines in 1994, ugly on 64-bit -> should be redone */

void _certify_ptr(void *ptr, int tag) {

	/* Make sure PTR_WIDTH is properly defined in DTL.h */
	if ((unsigned verylong)ptr < LO_MEM) {
		/* Not a proper memory address, halt DTL */
		handle_illegal_ptr(ptr,tag);
		_smx_end();
		exit(EXIT_FAILURE);
		}
	}
