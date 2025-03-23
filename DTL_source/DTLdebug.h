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
 *   File: DTLdebug.h
 *
 *   Purpose: header file for the DTI API
 *
 *
 *     ********************************
 *     *** DO NOT RELEASE THIS FILE ***
 *     *** ONLY FOR DEVELOPMENT USE ***
 *     *** NOT PART OF THE USER API ***
 *     ********************************
 *
 */


 /*******************************************************
  *
  *  DTI (Developer's Toolkit Interface) for test/debug.
  *  Most void-declared functions print data to console.
  *  (*) denotes undocumented call only for specialists.
  *
  *******************************************************/

// system calls
rcode DTLAPI DTI_set_folder(char *folder, int style); // (*)
rcode DTLAPI DTI_set_folder16(char *folder16);        // (*)
rcode DTLAPI DTI_reset_folder();                      // (*)
rcode DTLAPI DTI_get_folder(char *folder, unsigned *c_size, int style); // (*)
void  DTLAPI DTI_raw_trace(char *msg);
void  DTLAPI DTI_get_API_type(char *typestrg, unsigned c_size); // (*)

// frame calls
void  DTLAPI DTI_list_all_frames();
rcode DTLAPI DTI_split_DM_frame(int ufnbr); // (*)

// tree structure calls
void  DTLAPI DTI_tree_structure(int crit);
bool  DTLAPI DTI_is_tree(int crit);         // (*)
rcode DTLAPI DTI_crit_exists(int crit);     // (*)

// PW-base calls
int   DTLAPI DTI_nbr_of_P_mbox(int crit);

// V-base calls
int   DTLAPI DTI_real_V_node(int crit, int alt, int node); // (*)
int   DTLAPI DTI_nbr_of_V_mbox(int crit);

// autoscale calls (undocumented)
rcode DTLAPI DTI_set_AV_crit_scale(int crit, double v_min, double v_max);       // (*)
rcode DTLAPI DTI_reset_AV_crit_scale(int crit);                                 // (*)
rcode DTLAPI DTI_AV_scale_ratio(int c_from, int c_to, int mode, double *ratio); // (*)
rcode DTLAPI DTI_check_AV_values(int crit, int type, int count, ...);           // (*)
rcode DTLAPI DTI_is_AV_default_scale(int crit);                                 // (*)

// evaluation call
rcode DTLAPI DTI_get_support_mass(int crit, double belief_level, double *lobo, double *upbo); // (*)

// misc call
int   DTLAPI DTI_prune_alts(int fnr, int new_alts);

// base content print calls
void  DTLAPI DTI_show_W_base();
void  DTLAPIDTI_show_W_base_crit(int crit);
void  DTLAPI DTI_show_W_box();
void  DTLAPI DTI_show_W_mbox();
void  DTLAPI DTI_show_P_base(int crit);
void  DTLAPI DTI_show_P_box(int crit);
void  DTLAPI DTI_show_P_mbox(int crit);
void  DTLAPI DTI_show_V_base(int crit);
void  DTLAPI DTI_show_V_box(int crit);
void  DTLAPI DTI_show_V_mbox(int crit);

// moment calls (only for specialists with mathematical knowledge of the algorithms)
rcode DTLAPI DTI_bn_to_dtl_cdf(int crit, double cdf_bn, double *cdf_dtl);
rcode DTLAPI DTI_get_support_mid(int crit, double *cdf);
rcode DTLAPI DTI_get_mass_moments(int crit, double *rm1, double *cm2, double *cm3);
rcode DTLAPI DTI_get_psi_moments(int crit, int alt, double *rm1, double *cm2, double *cm3);
rcode DTLAPI DTI_get_bn_params(int crit, double *loc, double *scale, double *alpha);


 /*******************************************************
  *
  *  DTI piggybacking onto internal DTL functions
  *
  *******************************************************/

#if (CALL_STACK == 1) // C89/C99 = _cdecl
#define DTI_node2crit dtl_node2crit
#define DTI_crit2node dtl_crit2node
#define DTI_pure_W_tree dtl_pure_W_tree
#define DTI_W_node_parents dtl_W_node_parents
#define DTI_P_node_parents dtl_P_node_parents
#define DTI_real_W_crit dtl_real_W_crit
#define DTI_nbr_W_midpoints dtl_nbr_W_midpoints
#else // other stack calling conventions -> refurbish
int DTLAPI DTI_node2crit(int node);
int DTLAPI DTI_crit2node(int node);
int DTLAPI DTI_pure_W_tree();
int DTLAPI DTI_W_node_parents(int node1, int node2);
int DTLAPI DTI_P_node_parents(int crit, int alt, int node1, int node2);
int DTLAPI DTI_real_W_crit(int node);
int DTLAPI DTI_nbr_W_midpoints();
#endif
