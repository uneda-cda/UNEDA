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
 *   File: DTL.h
 *
 *   Purpose: header file for the DTL API
 *
 *
 *   Functions exported outside DTL
 *   ------------------------------
 *   NONE
 *
 *   Functions outside module, inside DTL
 *   ------------------------------------
 *   NONE
 *
 */


 /********************************************************
  *
  *  Package configuration parameters
  *
  ********************************************************/

/* Set the C_SML flag for UNEDA stakeholder layer */

#define noC_SML


 /*****************************************************
  *
  *  API call stack configuration
  *  ----------------------------
  *
  *  The following are the CALL_STACK options:
  *  0 = none     (default + 64-bit)
  *  1 = cdecl    (standard C89/C99)
  *  2 = stdcall  (standard MS WIN)
  *  3 = fastcall
  *
  *  The DTL library should go unlabeled as default.
  *  It is recommended to use _cdecl for C and C++
  *  32-bit DLLs and _stdcall for other languages.
  *  For 64-bit DLLs call stack options are ignored.
  *
  *****************************************************/

#define CALL_STACK 0

#if CALL_STACK == 1
#define DTLAPI __cdecl

#elif CALL_STACK == 2
#define DTLAPI __stdcall

#elif CALL_STACK == 3
#define DTLAPI __fastcall

#else // default
#define DTLAPI
#endif


 /*****************************************************
  *
  *  Pointer and long integer size configuration
  *  -------------------------------------------
  *
  *  C89, from which DTL and TCL originates, does not
  *  permit the use of long long 64-bit integers. They
  *  were not concieved of at that time. By the same
  *  token, most C compilers are in reality C++ ones.
  *  The original C++ standard incorporated C89 as a
  *  subset but not newer C standards. Thus, to stay
  *  compatible with an as large set of compilers as
  *  possible, the C89 standard is adhered to. This,
  *  however, comes with the caveat of not being able
  *  to represent 64-bit pointers as integers. Thus,
  *  the DTL package allows for another definition of
  *  very long 64-bit integers where such are provided.
  *
  *  The following are the PTR_WIDTH options:
  *  0 = unknown architecture (default)
  *  1 = 32-bit architecture  (older)
  *  2 = 64-bit architecture  (modern)
  *
  *****************************************************/

#define PTR_WIDTH 0

#if PTR_WIDTH == 2
#define verylong long long

#elif (PTR_WIDTH == 0) && defined(_MSC_VER)
#define verylong __int64

#else // default is 32 bits
#define verylong long
#endif


 /*****************************************************
  *
  *  DTL package parameters
  *
  *****************************************************/

#include "DTLparameters.h"

#define MAX_CRIT   300
#define MAX_FRAMES 101 // one for internal use

#define MAX_RESULTSTEPS 21

#define DTL_EPS 1.0E-5


 /*****************************************************
  *
  *  General macro definitions
  *
  *****************************************************/

/* String handling with safeguard to prevent buffer overrun. Beware of
 * buffer underrun if byte vectors are processed with string functions. */
#define  string(s) s,sizeof(s)         // s is already a string
#define _string(x) (char*)&x,sizeof(x) // treat x as a byte string


 /*****************************************************
  *
  *  Exported data types
  *
  *****************************************************/

typedef int rcode;
typedef int bool; // might get problems with newer MS VC++ compilers

typedef double a_vector[MAX_ALTS+1];
typedef a_vector ar_matrix[MAX_ALTS+1];
typedef int ai_vector[MAX_ALTS+1];
typedef ai_vector ai_matrix[MAX_ALTS+1];
typedef double h_vector[MAX_NOPA+1];
typedef h_vector h_matrix[MAX_ALTS+1];
typedef int o_matrix[MAX_ALTS+1][MAX_COPA+1];
typedef double e_matrix[MAX_ERESULT+1][MAX_RESULTSTEPS];
typedef int t_row[MAX_NOPA+1];
typedef t_row t_matrix[MAX_ALTS+1];
typedef double ar_col[MAX_ALTS+1];
typedef double cr_col[MAX_CRIT+1];
typedef int ai_col[MAX_ALTS+1];
typedef int ci_col[MAX_CRIT+1];
typedef double s_vector[MAX_ALTS+1];
typedef s_vector s_matrix[MAX_ERESULT+1];


 /*****************************************************
  *
  *  Return codes
  *
  *****************************************************/

#define DTL_OK                  0
#define DTL_KERNEL_ERROR     -100
#define DTL_INPUT_ERROR      -101
#define DTL_TREE_ERROR       -102
#define DTL_OUTPUT_ERROR     -103
#define DTL_FRAME_EXISTS     -104
#define DTL_FRAME_UNKNOWN    -105
#define DTL_FRAME_IN_USE     -106
#define DTL_FRAME_NOT_LOADED -107
#define DTL_FRAME_CORRUPT    -108
#define DTL_WRONG_FRAME_TYPE -109
#define DTL_WRONG_STMT_TYPE  -110
#define DTL_CONS_OVERFLOW    -111
#define DTL_CRIT_OVERFLOW    -112
#define DTL_LOGFILE_ERROR    -113
#define DTL_INCONSISTENT     -114
#define DTL_DIFFERING_RANKS  -115
#define DTL_STMT_ERROR       -116
#define DTL_SYS_CORRUPT      -117
#define DTL_ALT_OVERFLOW     -118
#define DTL_NODE_OVERFLOW    -119
#define DTL_CRIT_MISSING     -120
#define DTL_TOO_FEW_ALTS     -121
#define DTL_USER_ABORT       -122
#define DTL_STATE_ERROR      -123
#define DTL_CRIT_UNKNOWN     -124
#define DTL_CRIT_EXISTS      -125
#define DTL_ALT_UNKNOWN      -126
#define DTL_ALT_MISMATCH     -127
#define DTL_BUSY             -128
#define DTL_NAME_MISSING     -129
#define DTL_NAME_TOO_LONG    -130
#define DTL_NAME_EXISTS      -131
#define DTL_NOT_ALLOWED      -132
#define DTL_WRONG_METHOD     -133
#define DTL_WRONG_TOLERANCE  -134
#define DTL_FILE_UNKNOWN     -135
#define DTL_SCALE_CHANGE     -136
#define DTL_INTERNAL_ERROR   -137
#define DTL_WEAK_MASS_DISTR  -138
#define DTL_MEMORY_LEAK      -139
#define DTL_BUFFER_OVERRUN   -140
#define DTL_ASSERT_FAILED    -141
#define MAX_DTL_ERR DTL_ASSERT_FAILED
// piggybacked return codes
#define DTL_INFINITE_MASS       DTL_WEAK_MASS_DISTR
#define DTL_LAST_ROW_INCOMPLETE DTL_DIFFERING_RANKS

 /*****************************************************
  *
  *  Exported data structures
  *
  *****************************************************/

#define MAX_TERMS 2

struct user_stmt_rec {
	int n_terms;
	int alt[MAX_TERMS+1];
	int cons[MAX_TERMS+1];
	int sign[MAX_TERMS+1];
	double lobo;
	double upbo;
	};

struct user_w_stmt_rec {
	int n_terms;
	int crit[MAX_TERMS+1];
	int sign[MAX_TERMS+1];
	double lobo;
	double upbo;
	};

typedef struct tt_node {
	char type;
	int next;
	int down;
	} ttnode;

typedef ttnode ta_tree[MAX_NOPA+1];
typedef ta_tree tt_tree[MAX_ALTS+1];


 /*****************************************************
  *
  *  Exported constants
  *
  *****************************************************/

/* Frame types */
#define ANY_FRAME 0
#define PS_FRAME  1
#define DM_FRAME  2
#define PM_FRAME  3

/* Evaluation rules */
#define E_DELTA    0
#define E_GAMMA    4
#define E_PSI      8
#define E_DIGAMMA 12

/* Evaluation mask */
#define M_EVAL    28

/* Autoscale types */
#define ABS_SCALE  1
#define DIFF_SCALE 2
#define DIST_SCALE 3
#define REVD_SCALE 4


 /*****************************************************
  *
  *  DTL API library calls
  *
  *****************************************************/

/*** System commands ***/
rcode DTLAPI DTL_init();
rcode DTLAPI DTL_exit();
void  DTLAPI DTL_abort();

/*** Structure commands ***/
rcode DTLAPI DTL_new_PS_flat_frame(int ufnbr, int n_alts, int n_cons[]);
rcode DTLAPI DTL_new_PS_tree_frame(int ufnbr, int n_alts, int n_nodes[], tt_tree xtree);
rcode DTLAPI DTL_new_PM_flat_frame(int ufnbr, int n_crit, int n_alts);
rcode DTLAPI DTL_new_PM_tree_frame(int ufnbr, int n_alts, int n_wtnodes, ta_tree wtree);
rcode DTLAPI DTL_new_PM_crit_tree(int crit, int n_nodes[], tt_tree xtree);
rcode DTLAPI DTL_load_PM_crit(int crit, int ufnbr);
rcode DTLAPI DTL_unload_PM_crit(int crit, int new_ufnbr);
rcode DTLAPI DTL_delete_PM_crit(int crit);
rcode DTLAPI DTL_new_DM_flat_frame(int ufnbr, int n_crit, int n_alts);
rcode DTLAPI DTL_new_DM_tree_frame(int ufnbr, int n_alts, int n_wtnodes, ta_tree wtree);
rcode DTLAPI DTL_new_SM_tree_frame(int ufnbr, int mode, int n_alts, int n_sh, int n_wtnodes, ta_tree wtree);
rcode DTLAPI DTL_dispose_frame(int ufnbr);
rcode DTLAPI DTL_load_frame(int ufnbr);
rcode DTLAPI DTL_unload_frame();
rcode DTLAPI DTL_unload_frame2();
rcode DTLAPI DTL_frame_name(char *fname, unsigned c_size, int *ftype);
rcode DTLAPI DTL_load_status(int *f_loaded);
rcode DTLAPI DTL_frame_type(int ufnr, int *type);
rcode DTLAPI DTL_PM_crit_exists(int crit, int *exists);

/*** File commands ***/
rcode DTLAPI DTL_read_frame(int ufnbr, char *fn, char *folder, int mode);
rcode DTLAPI DTL_read_ddt_frame(int ufnbr, char *fn, char *folder, int mode);
rcode DTLAPI DTL_write_frame(char *fn, char *folder);

/*** Weight commands ***/
rcode DTLAPI DTL_add_W_statement(struct user_w_stmt_rec* uwstmtp);
rcode DTLAPI DTL_change_W_statement(int stmt_number, double lobo, double upbo);
rcode DTLAPI DTL_replace_W_statement(int stmt_number, struct user_w_stmt_rec* uwstmtp);
rcode DTLAPI DTL_delete_W_statement(int stmt_number);
rcode DTLAPI DTL_add_W_mid_statement(struct user_w_stmt_rec* uwstmtp);
rcode DTLAPI DTL_delete_W_mid_statement(struct user_w_stmt_rec* uwstmtp);
rcode DTLAPI DTL_set_W_box(h_vector lobox, h_vector upbox);
rcode DTLAPI DTL_set_W_mbox(h_vector lobox, h_vector upbox);
rcode DTLAPI DTL_set_W_mbox1(h_vector mbox);
rcode DTLAPI DTL_remove_W_mbox();
rcode DTLAPI DTL_get_W_hull(int global, h_vector lobo, h_vector mid, h_vector upbo);
rcode DTLAPI DTL_reset_W_base();

/*** Probability commands ***/
rcode DTLAPI DTL_add_P_statement(int crit, struct user_stmt_rec* ustmtp);
rcode DTLAPI DTL_change_P_statement(int crit, int stmt_number, double lobo, double upbo);
rcode DTLAPI DTL_replace_P_statement(int crit, int stmt_number, struct user_stmt_rec* ustmtp);
rcode DTLAPI DTL_delete_P_statement(int crit, int stmt_number);
rcode DTLAPI DTL_add_P_mid_statement(int crit, struct user_stmt_rec* ustmtp);
rcode DTLAPI DTL_delete_P_mid_statement(int crit, struct user_stmt_rec* ustmtp);
rcode DTLAPI DTL_set_P_box(int crit, h_matrix lobox, h_matrix upbox);
rcode DTLAPI DTL_set_P_mbox(int crit, h_matrix lobox, h_matrix upbox);
rcode DTLAPI DTL_set_P_mbox1(int crit, h_matrix mbox);
rcode DTLAPI DTL_remove_P_mbox(int crit);
rcode DTLAPI DTL_get_P_hull(int crit, int global, h_matrix lobo, h_matrix mid, h_matrix upbo);
rcode DTLAPI DTL_reset_P_base(int crit);

/*** Value commands ***/
rcode DTLAPI DTL_add_V_statement(int crit, struct user_stmt_rec* ustmtp);
rcode DTLAPI DTL_change_V_statement(int crit, int stmt_number, double lobo, double upbo);
rcode DTLAPI DTL_replace_V_statement(int crit, int stmt_number, struct user_stmt_rec* ustmtp);
rcode DTLAPI DTL_delete_V_statement(int crit, int stmt_number);
rcode DTLAPI DTL_add_V_mid_statement(int crit, struct user_stmt_rec* ustmtp);
rcode DTLAPI DTL_delete_V_mid_statement(int crit, struct user_stmt_rec* ustmtp);
rcode DTLAPI DTL_set_V_box(int crit, h_matrix lobox, h_matrix upbox);
rcode DTLAPI DTL_set_V_mbox(int crit, h_matrix lobox, h_matrix upbox);
rcode DTLAPI DTL_set_V_mbox1(int crit, h_matrix mbox);
rcode DTLAPI DTL_set_V_modal(int crit, int mode, h_matrix lobox, h_matrix modalx, h_matrix upbox);
rcode DTLAPI DTL_remove_V_mbox(int crit);
rcode DTLAPI DTL_get_V_hull(int crit, h_matrix lobo, h_matrix mid, h_matrix upbo);
rcode DTLAPI DTL_get_V_modal(int crit, h_matrix modal);
rcode DTLAPI DTL_check_V_modality(int crit, int Ai, int Aj);
rcode DTLAPI DTL_get_V_modality_matrix(int crit, ai_matrix modal_mx);
rcode DTLAPI DTL_reset_V_base(int crit);

/*** Automatic scale commands ***/
rcode DTLAPI DTL_set_AV_box(int crit, bool rev, bool renorm, h_matrix lobox, h_matrix upbox);
rcode DTLAPI DTL_set_AV_mbox(int crit, h_matrix lobox, h_matrix upbox);
rcode DTLAPI DTL_set_AV_mbox1(int crit, h_matrix mbox);
rcode DTLAPI DTL_set_AV_modal(int crit, int mode, bool rev, bool renorm, h_matrix lobox, h_matrix modalx, h_matrix upbox);
rcode DTLAPI DTL_get_AV_crit_scale(int crit, double *v_min, double *v_max);
rcode DTLAPI DTL_set_AV_MC_scale(double v_min, double v_max);
rcode DTLAPI DTL_copy_AV_MC_scale(int crit);
rcode DTLAPI DTL_reset_AV_MC_scale();
rcode DTLAPI DTL_get_AV_MC_scale(double *v_min, double *v_max);
rcode DTLAPI DTL_get_AV_user_vector(int crit, int type, int size, double v_val[], double av_val[]);
rcode DTLAPI DTL_get_AV_user_value(int crit, int type, double v_val, double *av_val);
rcode DTLAPI DTL_get_AV_user_intervals(int crit, int type, int size, double v_lobo[], double v_upbo[], 
		double av_lobo[], double av_upbo[]);
rcode DTLAPI DTL_get_AV_user_interval(int crit, int type, double v_lobo, double v_upbo, 
		double *av_lobo, double *av_upbo);
rcode DTLAPI DTL_get_AV_norm_vector(int crit, int type, int size, double av_val[], double v_val[]);
rcode DTLAPI DTL_get_AV_norm_value(int crit, int type, double av_val, double *v_val);
rcode DTLAPI DTL_get_AV_norm_intervals(int crit, int type, int size, double av_lobo[], double av_upbo[], 
		double v_lobo[], double v_upbo[]);
rcode DTLAPI DTL_get_AV_norm_interval(int crit, int type, double av_lobo, double av_upbo, 
		double *v_lobo, double *v_upbo);
rcode DTLAPI DTL_check_AV_user_values(int crit, int type, int count, ...);
rcode DTLAPI DTL_check_AV_norm_values(int type, int count, ...);

/*** Evaluation commands ***/
rcode DTLAPI DTL_evaluate_frame(int crit, int method, int Ai, int Aj, e_matrix e_result);
rcode DTLAPI DTL_evaluate_full(int crit, int method, int Ai, int Aj, e_matrix e_result);
rcode DTLAPI DTL_evaluate_omega(int Ai, int mode, cr_col o_result, ci_col o_rank);
rcode DTLAPI DTL_evaluate_omega1(int Ai, int mode, cr_col o_result, ci_col o_node);
// belief mass
rcode DTLAPI DTL_get_mass_range(int crit, double lo_level, double up_level, double *mass);
rcode DTLAPI DTL_get_mass_above(int crit, double lo_level, double *mass);
rcode DTLAPI DTL_get_mass_below(int crit, double up_level, double *mass);
rcode DTLAPI DTL_get_mass_density(int crit, double ev_level, double *density);
rcode DTLAPI DTL_get_support_mass(int crit, double belief_level, double *lobo, double *upbo);
rcode DTLAPI DTL_get_support_lower(int crit, double belief_level, double *lobo, double *upbo);
rcode DTLAPI DTL_get_support_upper(int crit, double belief_level, double *lobo, double *upbo);
rcode DTLAPI DTL_get_aversion_value(int crit, double risk_aversion, double *ra_value);
// compound evaluations
rcode DTLAPI DTL_compare_alternatives(int crit, int method, double belief_level, ar_col lo_value, ar_col up_value);
rcode DTLAPI DTL_delta_mass(int crit, int mode, ar_matrix delta_value, ar_matrix delta_mass);
// orderings and rankings
rcode DTLAPI DTL_rank_alternatives(int crit, int mode, double gamma_tolerance, double omega_tolerance, 
		ai_col gamma_rank, ai_col omega_rank, ar_col gamma_value, ar_col omega_value);
rcode DTLAPI DTL_daisy_chain(int crit, ai_col omega_rank, ar_col daisy_value, ar_col omega_value);
rcode DTLAPI DTL_daisy_chain1(int crit, int mode, ai_col omega_rank, ar_col daisy_value, ar_col omega_value);
rcode DTLAPI DTL_daisy_chain2(int crit, int mode, double radius, ai_col omega_rank, ar_col daisy_value, ar_col omega_value);
rcode DTLAPI DTL_pie_chart(int crit, ar_col pie_value);
rcode DTLAPI DTL_pie_chart1(int crit, double moderation, ar_col pie_value);
rcode DTLAPI DTL_pie_chart2(int crit, int mode, double moderation1, double moderation2, ar_col pie_value);
rcode DTLAPI DTL_sec_level(int crit, double v_min, s_matrix s_result);

/*** Dominance commands ***/
rcode DTLAPI DTL_get_dominance(int crit, int Ai, int Aj, double *cd_value, int *d_order);
rcode DTLAPI DTL_get_dominance_matrix(int crit, double threshold, ai_matrix dominance_mx);
rcode DTLAPI DTL_get_dominance_nt_matrix(int crit, double threshold, ai_matrix dominance_mx);
rcode DTLAPI DTL_get_dominance_rank(int crit, int mode, int dmode, double threshold, ai_vector dom_rank);
rcode DTLAPI DTL_get_cardinal_dominance_matrix(int crit, int dmode, double threshold, ar_matrix cardinal_mx);
rcode DTLAPI DTL_get_abs_dominance_matrix(int dmode, double threshold, ai_matrix dominance_mx);

/*** Sensitivity commands ***/
rcode DTLAPI DTL_get_W_tornado(int mode, h_matrix t_lobo, h_matrix t_upbo);
rcode DTLAPI DTL_get_W_tornado_alt(int alt, int mode, h_vector t_lobo, h_vector t_upbo);
rcode DTLAPI DTL_get_P_tornado(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo);
rcode DTLAPI DTL_get_MCP_tornado(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo);
rcode DTLAPI DTL_get_V_tornado(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo);
rcode DTLAPI DTL_get_MCV_tornado(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo);
rcode DTLAPI DTL_get_cons_influence(int crit, int mode, h_matrix result);

/*** Error commands ***/
char* DTLAPI DTL_get_errtxt(rcode drc);
char* DTLAPI DTL_get_errtxt_p(rcode drc);
rcode DTLAPI DTL_get_errtxt_i(rcode drc, char *str, unsigned *len);
rcode DTLAPI DTL_get_errtxt_i16(rcode drc, char *str, unsigned *len, bool LE);
bool  DTLAPI DTL_error(rcode drc);
int   DTLAPI DTL_error2(rcode drc);
bool  DTLAPI DTL_u_error(rcode drc);
int   DTLAPI DTL_u_error2(rcode drc);

/*** Miscellaneous commands ***/
void DTLAPI DTL_get_release(char *relstrg, unsigned c_size);
void DTLAPI DTL_get_release_long(char *relstrg, unsigned c_size);
void DTLAPI DTL_get_capacity(char *capstrg, unsigned c_size);
void DTLAPI DTL_get_J_properties(char *J_strg, unsigned c_size);
// number peeks
int DTLAPI DTL_nbr_of_W_stmts();
int DTLAPI DTL_nbr_of_P_stmts(int crit);
int DTLAPI DTL_nbr_of_V_stmts(int crit);
int DTLAPI DTL_nbr_of_weights();
int DTLAPI DTL_nbr_of_crit();
int DTLAPI DTL_nbr_of_alts();
int DTLAPI DTL_total_cons(int crit);
int DTLAPI DTL_nbr_of_cons(int crit, int alt);
int DTLAPI DTL_total_nodes(int crit);
int DTLAPI DTL_nbr_of_nodes(int crit, int alt);


 /********************************************************
  *
  *  SML: UNEDA stakeholder layer
  *
  ********************************************************/

#ifdef C_SML
#include "SML.h"
#endif
