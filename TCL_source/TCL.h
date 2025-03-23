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
 *   File: TCL.h
 *
 *   Purpose: the main header file for Tree Core
 *
 */

/* Return codes */
#define TCL_OK              0
#define TCL_INCONSISTENT    1
#define TCL_INPUT_ERROR     2
#define TCL_TREE_ERROR      3
#define TCL_ILLEGAL_NODE    4
#define TCL_TOO_MANY_ALTS   5
#define TCL_TOO_MANY_CONS   6
#define TCL_TOO_MANY_STMTS  7
#define TCL_TOO_NARROW_STMT 8
#define TCL_TOO_FEW_ALTS    9
#define TCL_CORRUPTED      10
#define TCL_ATTACHED       11
#define TCL_DETACHED       12
#define TCL_NO_FILE        13
#define TCL_UNLIMITED      14
#define TCL_OUT_OF_MEMORY  15
#define TCL_MEMORY_LEAK    16
#define MAX_RCODE TCL_MEMORY_LEAK

#ifndef bool
/* Microsoft C/C++ version 14.00.50727.762, which comes with Visual C++ 2005, and
 * version 15.00.30729.01, which comes with Visual C++ 2008, do not define bool.
 * Neither do any of the older ones, even though bool was introduced already in C99.
 * Suddenly, starting with Visual C++ 2019 (_MSC_VER > 1919), it is not possible
 * to declare bool as int any longer (even though that is what it is intentionally
 * defined as by the C standard to which there has been no changes). */
typedef int bool;
#endif
typedef int rcode; // result code

/* Library structures */
#define MAX_TERMS 2

struct stmt_rec {
	int n_terms;
	int alt[MAX_TERMS+1];
	int cons[MAX_TERMS+1];
	int sign[MAX_TERMS+1];
	double lobo;
	double upbo;
	};

/* Max number of user statements */
#define MAX_STMTS 301

/* Solver data types */
#define MAX_ROWS 2*MAX_STMTS+MAX_COPA
#define MAX_COLS 3*MAX_ROWS

typedef double d_row[MAX_NODES+1];
typedef d_row d_matrix[MAX_ROWS+1];
typedef int b_col[MAX_ROWS+1];
typedef char c_col[MAX_ROWS+1];
typedef double a_row[MAX_ALTS+1];
typedef char c_row[MAX_NODES+1];
typedef double a_vector[MAX_ALTS+1];
typedef a_vector a_result[MAX_ERESULT+1];
typedef int i_row[MAX_NODES+1];
typedef int tn_row[MAX_NODES+1];
typedef int t_row[MAX_NOPA+1];
typedef t_row t_matrix[MAX_ALTS+1];

struct base {
	int watermark;
	int n_stmts;
	struct stmt_rec stmt[MAX_STMTS+1];
	d_row lo_midbox;
	d_row up_midbox;
	d_row lo_im_midbox;
	d_row up_im_midbox;
	bool box;
	d_row box_lobo;
	d_row box_upbo;
	d_row im_box_lobo;
	d_row im_box_upbo;
	};

struct d_frame {
	int watermark;
	char name[FNSIZE+6];
	bool tree;
	bool attached;
	int n_alts;
	int n_cons[MAX_ALTS+1];   /* Real cons */
	int im_cons[MAX_ALTS+1];  /* Intermediate cons */
	int tot_cons[MAX_ALTS+1]; /* Both types of cons */
	t_matrix next;
	t_matrix prev;
	t_matrix down;
	t_matrix up;
	/* Bases */
	struct base *P_base;
	struct base *V_base;
	};


/* Link data structures */
struct link_rec {
	int h_alt;
	int head;
	int t_alt;
	int tail;
	char cmp;
	double value;
	} l_rec;
typedef struct link_rec link_set[2*MAX_STMTS+1];


/*** Input parameters ***/

/* Evaluation method */
#define OMEGA 0
#define DELTA 1
#define GAMMA 2
#define PSI 3
#define DIGAMMA 4
#define MAX_EMETHOD DIGAMMA


/*** Library calls ***/

/*** Frame procedures ***/
rcode TCL_create_flat_frame(struct d_frame **dfp, int n_alts, int n_cons[]);
rcode TCL_create_tree_frame(struct d_frame **dfp, int n_alts, int tot_cons[], 
				t_matrix next, t_matrix down);
rcode TCL_dispose_frame(struct d_frame *df);
rcode TCL_attach_frame(struct d_frame *df);
rcode TCL_detach_frame(struct d_frame *df);
int TCL_pure_tree(struct d_frame *df, int alt);
int TCL_different_parents(struct d_frame *df, int alt, int node1, int node2);
int TCL_nbr_of_siblings(struct d_frame *df, int alt, int node);

/*** P-base procedures ***/
rcode TCL_add_P_constraint(struct d_frame *df, struct stmt_rec *P_stmt);
rcode TCL_replace_P_constraint(struct d_frame *df, int stmt_nbr, struct stmt_rec *P_stmt);
rcode TCL_change_P_constraint(struct d_frame *df, int stmt_nbr, double lobo, double upbo);
rcode TCL_delete_P_constraint(struct d_frame *df, int stmt_nbr);
rcode TCL_add_P_mstatement(struct d_frame *df, struct stmt_rec *P_stmt);
rcode TCL_delete_P_mstatement(struct d_frame *df, struct stmt_rec *P_stmt);
rcode TCL_set_P_box(struct d_frame *df, d_row tbox_lobo, d_row tbox_upbo);
rcode TCL_unset_P_box(struct d_frame *df);
rcode TCL_set_P_mbox(struct d_frame *df, d_row tmbox_lobo, d_row tmbox_upbo);
rcode TCL_get_P_box(struct d_frame *df, d_row lo_box, d_row up_box);
rcode TCL_get_P_mbox(struct d_frame *df, d_row lo_mid, d_row up_mid);
rcode TCL_get_P_index(struct d_frame *df, int alt, int node);
rcode TCL_get_P_hull(struct d_frame *df, d_row hlobo, d_row hupbo,d_row llobo, d_row lupbo);
rcode TCL_get_P_masspoint(struct d_frame *df, d_row P_mid, d_row P_lmid);
rcode TCL_reset_P_base(struct d_frame *df);

/*** V-base procedures ***/
rcode TCL_add_V_constraint(struct d_frame *df, struct stmt_rec *V_stmt);
rcode TCL_replace_V_constraint(struct d_frame *df, int stmt_nbr, struct stmt_rec *V_stmt);
rcode TCL_change_V_constraint(struct d_frame *df, int stmt_nbr, double lobo, double upbo);
rcode TCL_delete_V_constraint(struct d_frame *df, int stmt_nbr);
rcode TCL_add_V_mstatement(struct d_frame *df, struct stmt_rec *V_stmt);
rcode TCL_delete_V_mstatement(struct d_frame *df, struct stmt_rec *V_stmt);
rcode TCL_set_V_box(struct d_frame *df, d_row tbox_lobo, d_row tbox_upbo);
rcode TCL_unset_V_box(struct d_frame *df);
rcode TCL_set_V_mbox(struct d_frame *df, d_row tmbox_lobo, d_row tmbox_upbo);
rcode TCL_get_V_box(struct d_frame *df, d_row lo_box, d_row up_box);
rcode TCL_get_V_mbox(struct d_frame *df, d_row lo_mid, d_row up_mid);
rcode TCL_get_V_index(struct d_frame *df, int alt, int node);
rcode TCL_get_V_hull(struct d_frame *df, d_row lobo, d_row upbo);
rcode TCL_get_V_masspoint(struct d_frame *df, d_row V_mid);
rcode TCL_reset_V_base(struct d_frame *df);

/*** Evaluation procedures ***/
rcode TCL_evaluate(struct d_frame *df, int Ai, int Aj, int eval_method, a_result result);
rcode TCL_evaluate_omega(struct d_frame *df, int Ai, double *result);

/*** Security levels ***/
rcode TCL_security_level(struct d_frame *df, double sec_level, 
				a_vector strong, a_vector marked, a_vector weak);

/*** Optimisation ***/
rcode TCL_get_P_max(struct d_frame *df, int alt, int snode, 
				d_row V_pt, d_row P_pt, d_row im_P_pt, bool positive, double *maxval);
rcode TCL_get_P_min(struct d_frame *df, int alt, int snode, 
				d_row V_pt, d_row P_pt, d_row im_P_pt, bool positive, double *minval);

/*** Moment calculus ***/
rcode TCL_get_moments(struct d_frame *df, a_row rm1, a_row cm2, a_row cm3);
rcode TCL_get_mc_moments(struct d_frame *df, int snode, d_row Vx_rm1, d_row Vx_cm2, 
			d_row Vx_cm3, double *rm1, double *cm2, double *cm3);
rcode TCL_get_P_sd(struct d_frame *df, int inx, double *sd);
rcode TCL_get_V_sd(struct d_frame *df, int inx, int im, double *sd);

/*** Error procedure ***/
char *TCL_get_errtxt(rcode rc);

/*** Memory management (common to TCL and DTL) ***/
void *mem_alloc(size_t size, char *type, char *source);
rcode mem_free(void* mem_ptr);
void mem_map();
rcode mem_exit();

/*** Unofficial entry points ***/
int TCL_get_real_index(int alt, int cons);
int TCL_get_tot_index(int alt, int cons);
