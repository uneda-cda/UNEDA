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
 *   File: uct.c
 *
 *   Purpose: main loop for UCT - UNEDA core tester
 */

#include "uct.h"
#include <signal.h>
#include <float.h>

#define N_CMDS 54
#define N_CMDS_PUBLIC N_CMDS-1

#define N_BOTHLOAD 4
#define N_NONLOAD 2

int cmd_nbr=0;

struct f_entry {
	char cmd_name[4];
	void (*func)();
	char *help_text;
	};

struct f_entry f_table[N_CMDS];


static void user_signal(int sig) {

	fprintf(stderr,"\007\n- UCT runtime abort -\n");
	/* Save and close open frame */
	if (uf->loaded)
		close_frame();
	mem_free(uf);
	exit(0);
	}


int init_uct() {
	int rc=TCL_OK;

	/* Create a new global user frame */
	uf = (struct user_frame *) mem_alloc(sizeof(struct user_frame),"struct user_frame","init_user");
	uf->frame_name[0] = '\0';
	uf->df = NULL;
	uf->loaded = FALSE;
	uf->weak = 0;
	signal(SIGINT, user_signal);
#ifdef THINKC
	console_options.title = "\pUCT - UNEDA Core Tester - UCT";
	console_options.nrows = N_CONSOLE_ROWS;
#endif
	return rc;
	}


void quit_uct() {

	/* Save and close open frame */
	if (uf->loaded)
		close_frame();
	mem_free(uf);
	if (user_error)
		printf("%d error%s during this hour session\n",user_error,user_error==1?"":"s");
#ifdef PERFORMANCE
	rc = PerfDump(PerfBlock,"\pdelta.perf",TRUE,10);
	if (rc)
		printf("Performance dump error (%d)\n",rc);
	TermPerf(PerfBlock);
#endif
	/* Exit to shell */
#ifdef THINKC
	if (!user_error && (secs>0L))
		console_options.pause_atexit = 0;
#endif
	exit(0);
	}


void welcome() {

	printf("\n\n\n");
  printf("         _/       _/   _/       _/    _/_/_/_/_/   _/_/_/          _/  \n");
  printf("        _/       _/   _/_/     _/    _/           _/    _/       _/ _/ \n");
  printf("       _/       _/   _/ _/    _/    _/           _/      _/    _/    _/\n");
  printf("      _/       _/   _/  _/   _/    _/_/_/_/     _/      _/   _/      _/\n");
  printf("     _/       _/   _/   _/  _/    _/           _/      _/   _/_/_/_/_/ \n");
  printf("    _/       _/   _/    _/ _/    _/           _/      _/   _/      _/  \n");
  printf("    _/     _/    _/     _/_/    _/           _/     _/    _/      _/   \n");
  printf("     _/_/_/     _/       _/    _/_/_/_/_/   _/_/_/_/     _/      _/\n\n\n");
	printf("\t     (+)\n");
	printf("    +----- o &&& o --------------------------------------------------+\n");
 	printf("    |    o         o               Prof. Mats Danielson              |\n");
 	printf("    |   o  SCIENCE  o              DECIDE Research Group             |\n");
	printf("    |   o    AND    o     Dept. of Computer and Systems Sciences     |\n");
	printf("    |   o    ART    o              Stockholm University              |\n");
	printf("    |    o         o       PO Box 1203, SE-164 25 Kista, SWEDEN      |\n");
	printf("    +------ o x o ---------------------------------------------------+\n\n");
	printf("            UNEDA %d.%d.%d Core Tester  (c) 2025 Mats Danielson\n\n",
				DTL_MAIN,DTL_FUNC,DTL_TECH);
	}

#define SHOW_MAX 2*N_CONSOLE_ROWS

void show_commands() {
	int i;

	for (i=0; i<(uf->loaded?N_CMDS_PUBLIC:N_BOTHLOAD+N_NONLOAD); i+=2) {
		printf(" %s  %s\t",f_table[i].cmd_name,f_table[i].help_text);
		printf(" %s  %s\n",f_table[i+1].cmd_name,f_table[i+1].help_text);
		}
	}


void show_version() {
#ifdef UNIX
	char hostname[HOST_NAME_MAX+1];
#endif

#ifdef UNIX
  gethostname(hostname,HOST_NAME_MAX+1);
	printf("Host name: %s\n\n",hostname);
#endif
	printf("<- UNEDA ->\n");
	printf("Ver. %d.%d.%d\n",DTL_MAIN,DTL_FUNC,DTL_TECH);
	printf("Alts  %5d\n",MAX_ALTS);
	printf("Cons  %5d\n",MAX_CONS);
	printf("Nodes %5d\n",MAX_NODES);
	printf("Stmts %5d\n",MAX_STMTS);
	}


void show_memory() {

	mem_map();
	}


// 10-17 chars = tab, 18-25 chars = no tab, other = out of range
struct f_entry f_table[N_CMDS] = {
	/* Commands for both modes */
	{"cmd",show_commands,"This command list\t"},
	{"bye",quit_uct,"Quit from UCT\t"},
	{"ver",show_version,"Show release version"},
	{"sfo",set_folder,"Set folder name\t"},
	/* Commands for non-loaded mode */
	{"new",new_frame,"Create new frame\t"},
	{"opn",open_frame,"Open existing frame"},
	/* Commands for loaded mode */
	{"mem",show_memory,"Show memory allocation"},
	{"sav",save_frame,"Save current frame"},
	{"sas",save_frame_as,"Save as another name"},
	{"rev",revert_frame,"Revert current frame"},
	{"cls",close_frame,"Close current frame"},
	{"aps",add_P_constraint,"Add prob constraint"},
	{"avs",add_V_constraint,"Add value constraint"},
	{"cps",change_P_constraint,"Change prob constraint"},
	{"cvs",change_V_constraint,"Change value constraint"},
	{"rps",replace_P_constraint,"Replace prob constraint"},
	{"rvs",replace_V_constraint,"Replace value constraint"},
	{"dps",remove_P_constraint,"Delete prob constraint"},
	{"dvs",remove_V_constraint,"Delete value constraint"},
	{"apm",make_P_midpoint,"Add prob midpoint"},
	{"rpm",remove_P_midpoint,"Remove prob midpoint"},
	{"avm",make_V_midpoint,"Add value midpoint"},
	{"rvm",remove_V_midpoint,"Remove value midpoint"},
	{"sal",show_all,"Show all info\t"},
	{"six",show_index,"Show frame indices"},
	{"spb",show_P_base,"Show prob base\t"},
	{"sph",show_P_hull,"Show prob hull\t"},
	{"spm",show_P_midpoints,"Show prob midpoints"},
	{"spc",show_P_core,"Show prob core\t"},
	{"svb",show_V_base,"Show value base\t"},
	{"svh",show_V_hull,"Show value hull\t"},
	{"svm",show_V_midpoints,"Show value midpoints"},
	{"svc",show_V_core,"Show value core\t"},
	{"sfp",show_fp,"Show focal point\t"},
	{"sfi",show_frame_info,"Show frame info\t"},
	{"sts",show_tree_structure,"Show tree structure"},
	{"cad",compare_delta,"Compare DELTA\t"},
	{"cag",compare_gamma,"Compare GAMMA\t"},
	{"cap",compare_psi,"Compare PSI\t"},
	{"cao",compare_omega,"Compare OMEGA\t"},
	{"smo",show_moments,"Show NEMO moments\t"},
	{"sel",security_level,"Security level\t"},
	{"ban",welcome,"Show welcome banner\t"}
	};


void dispatch(char *cmd) {
	int i,launch;

	for (i=0; i<N_CMDS; i++)
		if (!strncmp(cmd,f_table[i].cmd_name,3)) {
			/* Found command, ok to launch ? */
			if (i < N_BOTHLOAD)
				launch = TRUE;
			else if (i < N_BOTHLOAD+N_NONLOAD)
				launch = !uf->loaded;
			else
				launch = uf->loaded;
			if (!launch) {
				printf("Frame %s open\n",(uf->loaded?"already":"not"));
				user_error++;
				return;
				}
			(f_table[i].func)();
			return;
			}
	printf("%s: unknown command\n",cmd);
	user_error++;
	}


static int batch_NI_err=0;
static int batch_err=0;

#define CMD_SIZE 20

d_row P_ccp,P_lccp,V_ccp;
d_row V_ccp2,V_ccp3,V_ccp4;
#ifndef UNIX
#include <sys\timeb.h>
struct timeb run_start,run_end;
int run_time;
#endif

static i_row strong,marked,weak;
static a_result cube1;
static d_row P_min,P_mid,P_max,V_min,V_mid,V_max;
static a_row rm1,cm2;

int main() {
	char cmd[CMD_SIZE];
	int len;
#ifdef BATCH_MODE
	FILE *fp,*zp;
	char f_name[80],batseq[24];
	int i,j,k,Ai,Aj,rc,cont_mode,calc_method;
	struct d_frame *df;
	struct frame_info fi;
	int fexact,texact,flat_method;
	double min_value;
	int n_min,n_mid,n_max,loopies;
	char now[40];
	time_t v_start,v_stop;
#endif

	/* Init environment */
	init_uct();
	init_random();

	/* Display welcome banner */
	welcome();

	/* Loop on user commands */
	for (;;) {
		/* prompt for command (3 chars) */
		printf("UCT> ");
		fflush(stdout);
		fflush(stdin);
		len = 0;
		/* Leaving the nice scanf world for low level I/O.
		   Must work on Mac, Unix, Win. */
		while (len < CMD_SIZE) { 
			/* low level read works on Unix too */
			read(0,cmd+len,1);
			/* scanf and getchar fails on Unix */
			/* Should really do C_RAW mode instead */
			if (!len && (cmd[len] == 0x20))
				continue;
			if (cmd[len] == 0x0A) {
				cmd[len] = '\0';
				break;
				}
			else
				cmd[len] = tolower(cmd[len]);
			len++;
			}
		/* Back to the nice scanf world for high level I/O */
		user_cmd++;
		if (strlen(cmd) == 3) {
			/* dispatch on main command */
			dispatch(cmd);
			}
		else if (len) {
			printf("Try 'cmd' for a command list\n");
			user_error++;
			}
		}
	}
