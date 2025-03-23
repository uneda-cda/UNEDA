/*
 *
 *                 _/_/_/_/       _/        _/_/_/_/
 *               _/             _/ _/      _/     _/
 *              _/            _/    _/    _/      _/
 *             _/           _/      _/   _/      _/
 *            _/           _/_/_/_/_/   _/_/_/_/
 *           _/           _/      _/   _/    _/
 *           _/          _/      _/   _/      _/
 *            _/_/_/    _/      _/   _/       _/
 *
 *
 *        The Cardinal Alternative Ranking Method (CAR)
 *        ---------------------------------------------
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
 */

/*
 *   File: CAR.h
 *
 *   Purpose: header file for the CAR layer API
 *
 */


/* Borrow/piggyback error codes from DTL */

#define CAR_OK               DTL_OK
#define CAR_INPUT_ERROR      DTL_INPUT_ERROR
#define CAR_STATE_ERROR      DTL_STATE_ERROR
#define CAR_NOT_ALLOWED      DTL_NOT_ALLOWED
#define CAR_CRIT_UNKNOWN     DTL_CRIT_UNKNOWN
#define CAR_ALT_UNKNOWN      DTL_ALT_UNKNOWN
#define CAR_WRONG_FRAME_TYPE DTL_WRONG_FRAME_TYPE
#define CAR_FRAME_NOT_LOADED DTL_FRAME_NOT_LOADED
#define CAR_NOT_ACTIVATED    DTL_STATE_ERROR
#define CAR_INCONSISTENT     DTL_INCONSISTENT
#define CAR_ILLEGAL_TREE     DTL_TREE_ERROR
#define CAR_SYS_CORRUPT      DTL_SYS_CORRUPT
#define CAR_SAME_RANKINGS    DTL_DIFFERING_RANKS


/* CAR data structure defintion */

typedef int car_vector[MAX_CONS+1];


/* The CAR method layer calls on top of DTL */

rcode DTLAPI CAR_init(int method, int mode);
rcode DTLAPI CAR_exit();
rcode DTLAPI CAR_set_compat(double w_unc, double v_unc);
rcode DTLAPI CAR_get_W_ordinal(int n_nodes, cr_col ord_wts);
rcode DTLAPI CAR_set_W_base(int n_nodes, car_vector ord_crit, car_vector rel);
rcode DTLAPI CAR_set_P_base(int crit, int alt, int n_nodes, car_vector ord_nodes, car_vector rel);
rcode DTLAPI CAR_set_V_base(int crit, car_vector ord_alts, car_vector ord_nodes, car_vector rel);
/*** Weight partial hull (DURENO-II) ***/
rcode DTLAPI CAR_open_W_phull();
rcode DTLAPI CAR_close_W_phull();
rcode DTLAPI CAR_check_W_phull(struct user_w_stmt_rec* swp, double *tradeoff);
rcode DTLAPI CAR_prune_W_phull(struct user_w_stmt_rec* swp);
rcode DTLAPI CAR_cut_W_phull(struct user_w_stmt_rec* swp);
rcode DTLAPI CAR_equal_W_phull(struct user_w_stmt_rec* swp);
/*** Distance ranking ***/
rcode DTLAPI CAR_rank_W_base(int n_nodes, car_vector ord_crit, double dist);
rcode DTLAPI CAR_rank_P_base(int crit, int alt, int n_nodes, car_vector ord_nodes, double dist);
