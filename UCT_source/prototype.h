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
 *   File: prototype.h
 *
 *   Purpose: function prototypes for UCT - UNEDA core tester
 *
 */

#define noARTICLE // cmp for article

int call(int rc);
void prep_input();
void set_display_level();
void set_folder();
void new_frame();
void open_frame();
void close_frame();
void save_frame();
void save_frame_as();
void revert_frame();
void gen_frame();
void read_cac(struct d_frame *df, int *alt, int *node, int mode);
void read_P_interval(double *lobo, double *upbo);
double x1v(double x);
double x2v(double x);
double v0x(double v, int n_terms);
double v1x(double v);
double v2x(double v);
void read_V1_interval(double *lobo, double *upbo);
void read_V2_interval(double *lobo, double *upbo);
void make_P_midpoint();
void remove_P_midpoint();
void add_P_constraint();
void remove_P_constraint();
void make_V_midpoint();
void remove_V_midpoint();
void add_V_constraint();
void remove_V_constraint();
void replace_P_constraint();
void change_P_constraint();
void replace_V_constraint();
void change_V_constraint();
void show_index();
void show_frame_info();
void show_P_base();
void show_V_base();
void show_P_midpoints();
void show_V_midpoints();
void show_P_hull();
void show_V_hull();
void show_fp();
void show_moments();
void show_P_core();
void show_V_core();
void show_P_links();
void show_V_links();
void show_all();
void show_memory();
void show_frame_info();
void show_tree_structure();
void compare_delta();
void compare_gamma();
void compare_psi();
void compare_omega();
void security_level();

