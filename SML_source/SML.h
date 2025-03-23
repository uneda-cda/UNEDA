/*
 *
 *                   _/_/_/      _/_/_/_/_/   _/
 *                  _/    _/        _/       _/
 *                 _/      _/      _/       _/
 *                _/      _/      _/       _/
 *               _/      _/      _/       _/
 *              _/      _/      _/       _/
 *             _/     _/       _/       _/
 *            _/_/_/_/        _/       _/_/_/_/_/
 *
 *
 *              The Decision Tree Library (DTL)
 *              -------------------------------
 *
 *  +---- o o o --------------------------------------------+
 *  |   o       o            Prof. Mats Danielson           |
 *  |  o  STHLM  o           DECIDE Research Group          |
 *  |  o         o  Dept. of Computer and Systems Sciences  |
 *  |  o   UNI   o           Stockholm University           |
 *  |   o       o    PO Box 7003, SE-164 07 Kista, SWEDEN   |
 *  +---- o o o --------------------------------------------+
 *
 *           Copyright (c) 2022-2024 Mats Danielson
 *                Email: mats.danielson@su.se
 *                    Phone: +46 816 1540
 *
 *
 *   This code is proprietary, NOT free or shared!
 *   It may under NO circumstances be used for any
 *   purpose without a written agreement with the
 *   author.
 *
 */

/*
 *   File: SML.h
 *
 *   Purpose: header file for the SML API
 *
 *
 *   Functions exported outside SML
 *   ------------------------------
 *   SML_abort
 *   SML_new_PM_crit_tree
 *   SML_delete_PM_crit
 *   SML_PM_crit_exists
 *   SML_write_frame
 *   SML_get_W_hull
 *   SML_get_P_hull
 *   SML_set_scale
 *   SML_copy_scale
 *   SML_check_user_values
 *   SML_check_norm_values
 *   SML_error
 *   SML_error2
 *   SML_u_error
 *   SML_u_error2
 *   SML_get_release
 *   SML_nbr_of_weights
 *   SML_nbr_of_crit
 *   SML_nbr_of_alts
 *   SML_total_cons
 *   SML_nbr_of_cons
 *   SML_total_nodes
 *   SML_nbr_of_nodes
 *   SML_open_W_phull
 *   SML_close_W_phull
 *   SML_check_W_phull
 *   SML_prune_W_phull
 *   SML_cut_W_phull
 *   SML_equal_W_phull
 *   SML_mod_W_phull
 *   SML_CAR_open_W_phull
 *   SML_CAR_close_W_phull
 *
 *   Macros exported outside SML
 *   ---------------------------
 *   shc
 *   sc
 *
 *
 *   Version history
 *
 *   Ver.   Date   Main reasons
 *   ----  ------  ------------
 *   1.22  220520  SML introduced
 *   1.24  221010  New SML calls
 *   1.27  230606  New SML CAR calls
 *
 */


 /*****************************************************
  *
  *  Included layers
  *
  *****************************************************/

#include "CAR.h" // sibling at the same layer level


 /*****************************************************
  *
  *  Exported data types
  *
  *****************************************************/

#define MAX_CDF 100
typedef double c_vector[MAX_CDF+1];

#define P_MAX_CRIT 12
#define P_MAX_ALTS 10
/* For correlation */
typedef cr_col cr_matrix[MAX_CRIT+1];
/* For correlation */
typedef int p_ivector[P_MAX_CRIT*(P_MAX_ALTS+1)];
typedef double p_dvector[P_MAX_CRIT*(P_MAX_ALTS+1)];
typedef double p_avector[P_MAX_ALTS+1];


 /*****************************************************
  *
  *  SML-specific interface functions
  *
  *****************************************************/

// rcode interpreters
rcode DTLAPI scall(rcode func);
rcode DTLAPI ucall(rcode func);
// system commands
rcode DTLAPI SML_init();
rcode DTLAPI SML_init2(int mode);
rcode DTLAPI SML_exit();
// frame functions
rcode DTLAPI SML_new_DM_flat_frame(int ufnbr, int n_crit, int n_alts);
rcode DTLAPI SML_new_DM_tree_frame(int ufnbr, int n_alts, int n_wtnodes, ta_tree wtree);
rcode DTLAPI SML_new_PM_flat_frame(int ufnbr, int n_crit, int n_alts);
rcode DTLAPI SML_new_PM_tree_frame(int ufnbr, int n_alts, int n_wtnodes, ta_tree wtree);
rcode DTLAPI SML_new_SM_tree_frame(int ufnbr, int type, int n_alts, int n_sh,	int n_wtnodes, ta_tree wt_tree);
rcode DTLAPI SML_read_frame(int ufnbr, int type, int n_sh, char *fn, char *folder);
rcode DTLAPI SML_load_frame(int ufnbr);
rcode DTLAPI SML_unload_frame();
rcode DTLAPI SML_dispose_frame(int ufnbr);
rcode DTLAPI SML_frame_type(int ufnbr);
// data entry functions
rcode DTLAPI SML_set_W_base(h_vector lobox, h_vector mbox, h_vector upbox);
rcode DTLAPI SML_set_W_base2(h_vector lobox, h_vector mbox, h_vector upbox, int *inc_var);
rcode DTLAPI SML_set_P_base(int crit, h_matrix lobox, h_matrix mbox, h_matrix upbox);
rcode DTLAPI SML_set_P_base2(int crit, h_matrix lobox, h_matrix mbox, h_matrix upbox, int *inc_var);
rcode DTLAPI SML_set_V_base(int crit, bool rev, int renorm, h_matrix lobox, h_matrix mbox, h_matrix upbox);
rcode DTLAPI SML_set_V_base2(int crit, bool rev, int renorm, h_matrix lobox, h_matrix mbox, h_matrix upbox, int *inc_var);
rcode DTLAPI SML_get_V_hull(int crit, h_matrix lobo, h_matrix mid, h_matrix upbo);
rcode DTLAPI SML_set_W_correlations(cr_matrix corr_mx);
// basic evaluation
rcode DTLAPI SML_evaluate_frame(int crit, int method, int Ai, int Aj, e_matrix e_result);
// belief mass functions
rcode DTLAPI SML_get_mass_range(double lo_level, double up_level, double *mass);
rcode DTLAPI SML_get_mass_above(double lo_level, double *mass);
rcode DTLAPI SML_get_mass_below(double up_level, double *mass);
rcode DTLAPI SML_get_support_mass(double belief_level, double *lobo, double *upbo);
rcode DTLAPI SML_get_support_lower(double belief_level, double *lobo, double *upbo);
rcode DTLAPI SML_get_support_upper(double belief_level, double *lobo, double *upbo);
rcode DTLAPI SML_evaluate_cdf(int crit, int Ai, c_vector level, c_vector cdf);
// compound evaluations
rcode DTLAPI SML_compare_alternatives(int crit, int method, double belief_level, ar_col lo_value, ar_col up_value);
rcode DTLAPI SML_evaluate_mid(int Ai, int mode, cr_col o_result, ci_col o_rank);
rcode DTLAPI SML_evaluate_omega(int Ai, cr_col o_result);
rcode DTLAPI SML_evaluate_omega1(int Ai, cr_col o_result, ci_col o_node);
rcode DTLAPI SML_evaluate_omega2(int Ai, int mode, cr_col o_result, ci_col o_node);
rcode DTLAPI SML_delta_mass(int crit, ar_matrix delta_mass, ai_col delta_order);
rcode DTLAPI SML_delta_mass2(int crit, int mode, ar_matrix delta_mass, ai_col delta_order);
rcode DTLAPI SML_rank_alternatives(int crit, int mode, double gamma_tolerance, double omega_tolerance, 
  ai_col gamma_rank, ai_col omega_rank, ar_col gamma_value, ar_col omega_value);
rcode DTLAPI SML_rank_gamma(int crit, ai_col gamma_rank, ar_col gamma_value);
rcode DTLAPI SML_rank_omega(int crit, ai_col omega_rank, ar_col omega_value);
rcode DTLAPI SML_daisy_chain(int crit, ai_col daisy_rank, ar_col daisy_value);
rcode DTLAPI SML_daisy_chain2(int crit, int mode, ai_col omega_rank, ar_col daisy_value, ar_col omega_value);
rcode DTLAPI SML_pie_chart(int crit, ar_col pie_value);
rcode DTLAPI SML_pie_chart1(int crit, ar_col pie_value);
rcode DTLAPI SML_pie_chart2(int crit, double moderation, ar_col pie_value);
// sensitivity analyses
rcode DTLAPI SML_get_W_tornado(h_matrix t_lobo, h_matrix t_upbo);
rcode DTLAPI SML_get_W_tornado2(int mode, h_matrix t_lobo, h_matrix t_upbo);
rcode DTLAPI SML_get_P_tornado(int crit, h_matrix t_lobo, h_matrix t_upbo);
rcode DTLAPI SML_get_P_tornado2(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo);
rcode DTLAPI SML_get_MCP_tornado(int crit, h_matrix t_lobo, h_matrix t_upbo);
rcode DTLAPI SML_get_MCP_tornado2(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo);
rcode DTLAPI SML_get_V_tornado(int crit, h_matrix t_lobo, h_matrix t_upbo);
rcode DTLAPI SML_get_V_tornado2(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo);
rcode DTLAPI SML_get_MCV_tornado(int crit, h_matrix t_lobo, h_matrix t_upbo);
rcode DTLAPI SML_get_MCV_tornado2(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo);
rcode DTLAPI SML_get_cons_influence(int crit, int mode, h_matrix result);
// error text
char* DTLAPI SML_get_errtxt(rcode drc);
char* DTLAPI SML_get_errtxt2(rcode drc, int style);
// service function (only for test/debug)
rcode DTLAPI SML_is_stakeholder(int node);


 /*****************************************************
  *
  *  Stakeholder management
  *
  *****************************************************/

/* MS VC++ 6.0 (pre-C99) lacks variadic macros -> function macros cannot
 * be defined using variable argument lists. Instead, use sh/crit macros
 * for accessing data and results in stakeholder addressing mode.
 *
 * Macros
 * ------
 * shc(sh,cr,nn) find a criterion for a stakeholder
 * sc(sh,cr)     find a criterion w/o supplying tot nbr of crit
 *
 * Macro parameters
 * ----------------
 * sh  the stakeholder number
 * cr  the criterion number within the stakeholder
 * nc  the number of criteria per stakeholder */

int DTLAPI sml_nc();
#define shc(sh,cr,nc) (sh-1)*(nc)+(cr)
#define sc(sh,cr) shc(sh,cr,sml_nc())


/* If the crit number is not verified using a call to any of the other SML
 * functions, use SML_crit_nbr instead which has built-in parameter checks. */

int DTLAPI SML_crit_nbr(int sh, int crit);


 /*****************************************************
  *
  *  Functions piggybacked onto the DTL interface
  *
  *****************************************************/

#define SML_abort DTL_abort
#define SML_new_PM_crit_tree DTL_new_PM_crit_tree
#define SML_delete_PM_crit DTL_delete_PM_crit
#define SML_PM_crit_exists DTL_PM_crit_exists
#define SML_write_frame DTL_write_frame
#define SML_get_W_hull DTL_get_W_hull
#define SML_get_P_hull DTL_get_P_hull
#define SML_set_scale DTL_set_AV_MC_scale
#define SML_copy_scale DTL_copy_AV_MC_scale
#define SML_check_user_values DTL_check_AV_user_values
#define SML_check_norm_values DTL_check_AV_norm_values
#define SML_error DTL_error
#define SML_error2 DTL_error2
#define SML_u_error DTL_u_error
#define SML_u_error2 DTL_u_error2
#define SML_get_release DTL_get_release
#define SML_nbr_of_weights DTL_nbr_of_weights
#define SML_nbr_of_crit DTL_nbr_of_crit
#define SML_nbr_of_alts DTL_nbr_of_alts
#define SML_total_cons DTL_total_cons
#define SML_nbr_of_cons DTL_nbr_of_cons
#define SML_total_nodes DTL_total_nodes
#define SML_nbr_of_nodes DTL_nbr_of_nodes


 /*****************************************************
  *
  *  SML functions piggybacked onto SML CAR
  *
  *****************************************************/

#define SML_open_W_phull SML_CAR_open_W_phull
#define SML_close_W_phull SML_CAR_close_W_phull
#define SML_check_W_phull SML_CAR_check_W_phull
#define SML_prune_W_phull SML_CAR_prune_W_phull
#define SML_cut_W_phull SML_CAR_cut_W_phull
#define SML_equal_W_phull SML_CAR_equal_W_phull
#define SML_mod_W_phull SML_CAR_mod_W_phull


 /*****************************************************
  *
  *  SML interface - CAR functions
  *
  *  Note that the check functions deviate from the
  *  DTL standard of returning a result/return code
  *
  *****************************************************/

bool  DTLAPI SML_CAR_check_W_nodes(int n_nodes, car_vector ord_nodes);
bool  DTLAPI SML_CAR_check_P_nodes(int crit, int alt, int n_nodes, car_vector ord_nodes);
bool  DTLAPI SML_CAR_check_V_nodes(int crit, car_vector ord_alts, car_vector ord_nodes);
rcode DTLAPI SML_CAR_set_W_base(int n_nodes, car_vector ord_nodes, car_vector rel, bool reset);
rcode DTLAPI SML_CAR_set_P_base(int crit, int alt, int n_nodes, car_vector ord_nodes, car_vector rel, bool reset);
rcode DTLAPI SML_CAR_set_V_base(int crit, car_vector ord_alts, car_vector ord_nodes, car_vector rel);
rcode DTLAPI SML_CAR_check_W_phull(struct user_w_stmt_rec* swp, double *tradeoff);
rcode DTLAPI SML_CAR_prune_W_phull(struct user_w_stmt_rec* swp);
rcode DTLAPI SML_CAR_cut_W_phull(struct user_w_stmt_rec* swp);
rcode DTLAPI SML_CAR_equal_W_phull(struct user_w_stmt_rec* swp);
rcode DTLAPI SML_CAR_mod_W_phull(int type, struct user_w_stmt_rec* swp);


 /*****************************************************
  *
  *  Functions piggybacked onto the DTL CAR interface
  *
  *****************************************************/

#define SML_CAR_open_W_phull CAR_open_W_phull
#define SML_CAR_close_W_phull CAR_close_W_phull


 /*****************************************************
  *
  *  SML survey functions interface
  *
  *****************************************************/

rcode DTLAPI SML_open_survey(int crit, int alts, int smin, int smax, double tse);
rcode DTLAPI SML_add_survey_row(p_ivector p_entry);
rcode DTLAPI SML_convert_survey(int ufnbr, p_avector p_result);
rcode DTLAPI SML_close_survey();


 /*****************************************************
  *
  *  Piggybacked return codes
  *
  *****************************************************/

#define SML_OK               DTL_OK
#define SML_KERNEL_ERROR     DTL_KERNEL_ERROR
#define SML_INPUT_ERROR      DTL_INPUT_ERROR
#define SML_TREE_ERROR       DTL_TREE_ERROR
#define SML_OUTPUT_ERROR     DTL_OUTPUT_ERROR
#define SML_FRAME_EXISTS     DTL_FRAME_EXISTS
#define SML_FRAME_UNKNOWN    DTL_FRAME_UNKNOWN
#define SML_FRAME_IN_USE     DTL_FRAME_IN_USE
#define SML_FRAME_NOT_LOADED DTL_FRAME_NOT_LOADED
#define SML_FRAME_CORRUPT    DTL_FRAME_CORRUPT
#define SML_WRONG_FRAME_TYPE DTL_WRONG_FRAME_TYPE
#define SML_WRONG_STMT_TYPE  DTL_WRONG_STMT_TYPE
#define SML_CONS_OVERFLOW    DTL_CONS_OVERFLOW
#define SML_CRIT_OVERFLOW    DTL_CRIT_OVERFLOW
#define SML_LOGFILE_ERROR    DTL_LOGFILE_ERROR
#define SML_INCONSISTENT     DTL_INCONSISTENT
#define SML_DIFFERING_RANKS  DTL_DIFFERING_RANKS
#define SML_STMT_ERROR       DTL_STMT_ERROR
#define SML_SYS_CORRUPT      DTL_SYS_CORRUPT
#define SML_ALT_OVERFLOW     DTL_ALT_OVERFLOW
#define SML_NODE_OVERFLOW    DTL_NODE_OVERFLOW
#define SML_CRIT_MISSING     DTL_CRIT_MISSING
#define SML_TOO_FEW_ALTS     DTL_TOO_FEW_ALTS
#define SML_USER_ABORT       DTL_USER_ABORT
#define SML_STATE_ERROR      DTL_STATE_ERROR
#define SML_CRIT_UNKNOWN     DTL_CRIT_UNKNOWN
#define SML_CRIT_EXISTS      DTL_CRIT_EXISTS
#define SML_ALT_UNKNOWN      DTL_ALT_UNKNOWN
#define SML_ALT_MISMATCH     DTL_ALT_MISMATCH
#define SML_BUSY             DTL_BUSY
#define SML_NAME_MISSING     DTL_NAME_MISSING
#define SML_NAME_TOO_LONG    DTL_NAME_TOO_LONG
#define SML_NAME_EXISTS      DTL_NAME_EXISTS
#define SML_NOT_ALLOWED      DTL_NOT_ALLOWED
#define SML_WRONG_METHOD     DTL_WRONG_METHOD
#define SML_WRONG_TOLERANCE  DTL_WRONG_TOLERANCE
#define SML_FILE_UNKNOWN     DTL_FILE_UNKNOWN
#define SML_SCALE_CHANGE     DTL_SCALE_CHANGE
#define SML_INTERNAL_ERROR   DTL_INTERNAL_ERROR
#define SML_WEAK_MASS_DISTR  DTL_WEAK_MASS_DISTR
#define SML_MEMORY_LEAK      DTL_MEMORY_LEAK
#define SML_BUFFER_OVERRUN   DTL_BUFFER_OVERRUN
#define SML_ASSERT_FAILED    DTL_ASSERT_FAILED
