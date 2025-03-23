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
 *   File: file.c
 *
 *   Purpose: file handling in UCT - UNEDA core tester
 *
 */

#include "uct.h"
#include <signal.h>
#ifndef UNIX
#include <io.h>
#endif
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

static t_matrix tnext,tdown;

void set_folder() {
	char input[32],u_input[32];
#ifndef WINNT
	int fd;
#endif

	printf("Folder [%s]: ",folder[0]?folder:".");
	prep_input();
	scanf("%29s",input);
	if (!strcmp(input,".")) {
		strcpy(folder,"");
		strcpy(u_folder,HOME_FOLDER);
		}
	else {
#ifdef MAC
		*u_input = ':';
		strcpy(u_input+1,input);
		strcat(u_input,":");
#else
#ifdef WINNT
		*u_input = '.';
		*(u_input+1) = '\\';
		strcpy(u_input+2,input);
		strcat(u_input,"\\");
#else	/* Unix */
		*u_input = '.';
		*(u_input+1) = '/';
		strcpy(u_input+2,input);
		strcat(u_input,"/");
		if ((fd=open(u_input,0)) == -1)
			printf("%s: no such folder - creating new\n",input);
		else
			close(fd);
#endif
#ifdef _MSC_VER
		_mkdir(u_input);
#else
		mkdir(u_input,0700);
#endif
#endif
		strcpy(folder,input);
		strcpy(u_folder,u_input);
		}
	}


int fill_tree(int alt, int cnstart, t_row ftnext, t_row ftdown) {
	int cn,step;
	char ijrs;

	for (cn=cnstart;;) {
		printf("N%d.%d (I/J/R/S): ",alt,cn);
		fflush(stdin);
		prep_input();
		scanf(" %c",&ijrs);
		if (tolower(ijrs)=='i') {
			ftdown[cn] = cn+1;
			step = fill_tree(alt,cn+1,ftnext,ftdown);
			ftnext[cn] = cn+step+1;
			cn += step+1;
			}
		else if (tolower(ijrs)=='j') {
			ftdown[cn] = cn+1;
			step = fill_tree(alt,cn+1,ftnext,ftdown);
			ftnext[cn] = 0;
			return cn+step+1-cnstart;
			}
		else if (tolower(ijrs)=='r') {
			ftdown[cn] = 0;
			ftnext[cn] = cn+1;
			cn++;
 			}
		else if (tolower(ijrs)=='s') {
			ftdown[cn] = 0;
			ftnext[cn] = 0;
			return cn+1-cnstart;
			}
		}
	}


void new_frame() {
	int i,j,n_alts,n_nodes[MAX_ALTS+1],rem_nodes,t_nodes;
	char yn;

	if (folder[0])
		printf("Decision frame (folder '%s'): ",folder);
	else
		printf("Decision frame (home folder): ");
	prep_input();
	scanf("%29s",uf->frame_name);
	/* Create a new d_frame with at least two alts */
	do {
		printf("Number of alternatives (2..%d): ",MAX_ALTS);
		prep_input();
		scanf(" %d",&n_alts);
		if (!n_alts)
			return;
		} while ((n_alts < 2) || (n_alts > MAX_ALTS));

	/* Get alternative names */
	for (i=1; i<=n_alts; i++) {
reenter_alt:
		printf("Name of alternative %d: ",i);
		fflush(stdin);
		prep_input();
		scanf("%29s",uf->alt_name[i]);
		for (j=1; j<i; j++)
			if (!strcmp(uf->alt_name[i],uf->alt_name[j])) {
				printf("\007Same name as alternative %d\n",j);
				goto reenter_alt;
				}
		}	

	/* Specify number of levels */
	uf->multilevel = FALSE;
	printf("Multi-level tree (Y/N/A): ");
	fflush(stdin);
	prep_input();
	scanf(" %c",&yn);
	if (tolower(yn)=='a')
		return;
	if (tolower(yn)=='y') {
		uf->multilevel = TRUE;
		/* Tree: specify number of nodes */
		rem_nodes = MAX_NODES;
		for (i=1; i<=n_alts; i++) {
			do {
				printf("Number of nodes in A%d (1..%d): ",i,rem_nodes+i-n_alts);
				prep_input();
				scanf(" %d",n_nodes+i);
				if (!n_nodes[i])
					return;
				} while ((n_nodes[i] < 1) || (n_nodes[i] > rem_nodes+i-n_alts));
			rem_nodes -= n_nodes[i];
			if (rem_nodes < n_alts-i) {
				printf("Too many nodes\n");
				return;
				}
			}
		/* Get tree structure */
		for (i=1; i<=n_alts; i++) {
			t_nodes = fill_tree(i,1,tnext[i],tdown[i]);
			if (t_nodes != n_nodes[i]) {
				printf("Alt %d: expected %d nodes, got %d\n",i,n_nodes[i],t_nodes);
				return;
				}
			}
		}
	else {
		/* Flat: specify number of consequences */
		rem_nodes = MAX_CONS;
		for (i=1; i<=n_alts; i++) {
			do {
				printf("Number of consequences in A%d (1..%d): ",i,rem_nodes+i-n_alts);
				prep_input();
				scanf(" %d",n_nodes+i);
				if (!n_nodes[i])
					return;
				} while ((n_nodes[i] < 1) || (n_nodes[i] > rem_nodes+i-n_alts));
			rem_nodes -= n_nodes[i];
			if (rem_nodes < n_alts-i) {
				printf("Too many consequences\n");
				return;
				}
			}
		}

	/* Range of value statements */
	printf("Lower value limit: ");
	prep_input();
	scanf(" %lf",&v_min);
	do {
		printf("Upper value limit: ");
		prep_input();
		scanf(" %lf",&v_max);
		} while (v_min >= v_max);

	/* Backup twin and create frame */
	if (!backup_ufile(uf->frame_name,u_folder))
		printf("Existing frame '%s' backed up as '%s.bkp'\n",uf->frame_name,uf->frame_name);
	if (uf->multilevel)	{
		if (call(TCL_create_tree_frame(&(uf->df),n_alts,n_nodes,tnext,tdown)))
			return;
		}
	else {
		if (call(TCL_create_flat_frame(&(uf->df),n_alts,n_nodes)))
			return;
		}
	strcpy(uf->df->name,uf->frame_name);
	// UNEDA version
	uf_ver_main = DTL_MAIN;
	uf_ver_func = DTL_FUNC;
	uf_ver_tech = DTL_TECH;
	if (call(TCL_attach_frame(uf->df))) {
		TCL_dispose_frame(uf->df);
		return;
		}
	uf->v_lo = v_min;
	uf->v_up = v_max;
	v_mm = v_max - v_min;
	uf->loaded = TRUE;
	uf->dirty = TRUE;
	ask_to_save = TRUE;
	}


/* Made public to enable performance batch testing */
#ifndef BATCH_MODE
static
#endif
int open_frame2(char *f_name) {

	/* read_ufile creates df */
	if (call(read_ufile(f_name,u_folder,uf))) {
		if (uf->df)
			TCL_dispose_frame(uf->df);
		return -1;
		}
	v_min = uf->v_lo;
	v_max = uf->v_up;
	v_mm = v_max - v_min;
	uf->loaded = TRUE;
	uf->dirty = FALSE;
	return 0;
	}


void open_frame() {
	char f_name[32],fn_ddt[64];
#ifdef _MSC_VER
	struct _stat buf;
#else
	struct stat buf;
#endif
  int rc,fh,result;

	if (folder[0])
		printf("Decision frame (folder '%s'): ",folder);
	else
		printf("Decision frame (home folder): ");
	prep_input();
	scanf(" %29s",f_name);
	if ((rc = open_frame2(f_name))==0) {
		strcpy(fn_ddt,u_folder);
		strcat(fn_ddt,f_name);
		strcat(fn_ddt,".ddt");
#ifdef _MSC_VER
		if((fh=_open(fn_ddt,_O_RDONLY)) < 0)
#else
		if((fh=open(fn_ddt,O_RDONLY)) < 0)
#endif
			printf("No frame file data available\n");
		else {
		 /* Get data associated with "fh": */
#ifdef _MSC_VER
			result = _fstat(fh,&buf);
#else
			result = fstat(fh,&buf);
#endif
		 /* Check if statistics are valid: */
			if(result)
				printf("Bad frame file handle %s\n",fn_ddt);
			else {
				printf("PS-%s '%s' ",uf->multilevel?"tree":"frame",uf->frame_name);
				printf("contains %ld bytes\n",buf.st_size);
				printf("Created  %s",ctime(&buf.st_ctime));
				if (buf.st_mtime > buf.st_ctime+5)
					printf("Modified %s",ctime(&buf.st_mtime));
				}
			}
#ifdef _MSC_VER
//	_close(fh);
#else
//	close(fh);
#endif
		}
	}


void close_frame() {
	char yn; 

	/* Note: also entered from 'quit_delta' and 'signal' */
	if (uf->dirty) {
		/* Ask to save modified frame */
		if (folder[0])
			printf("Save frame '%s' in folder '%s' (y/n): ",uf->frame_name,folder);
		else
			printf("Save frame '%s' in home folder (y/n): ",uf->frame_name);
		prep_input();
		scanf(" %c",&yn);
		if (tolower(yn) == 'y') {
			/* Try first current folder and then home folder */
			if (call(write_ufile(uf->frame_name,u_folder,uf))) {
				if (!folder[0]) // already in home folder
					goto finish_close; // error message printed
				if (call(write_ufile(uf->frame_name,HOME_FOLDER,uf)))
					goto finish_close;
				printf("Folder '%s' corrupt - frame '%s' saved in home folder\n",
								folder,uf->frame_name);
				}
			uf->dirty = FALSE;
			}
		else if (tolower(yn) != 'n')
			printf("Unknown response - not saved\n\007");
		}
finish_close:
	/* Release resources */
	call(TCL_detach_frame(uf->df));
	call(TCL_dispose_frame(uf->df));
	uf->df = NULL;
	uf->frame_name[0] = '\0';
	uf->loaded = FALSE;
	}


void save_frame() {

	if (uf->dirty) {
		if (call(write_ufile(uf->frame_name,u_folder,uf)))
			return;
		if (folder[0])
			printf("Saved frame '%s' in folder '%s'\n",uf->frame_name,folder);
		else
			printf("Saved frame '%s' in home folder\n",uf->frame_name);
		uf->dirty = FALSE;
		}
	else
		printf("No changes need saving\n");
	}


void save_frame_as() {

	if (folder[0])
		printf("New frame name (folder '%s'): ",folder);
	else
		printf("New frame name (home folder): ");
	prep_input();
	scanf("%29s",uf->frame_name);
	if (!backup_ufile(uf->frame_name,u_folder))
		printf("Existing frame '%s' backed up onto file '%s.bkp'\n",uf->frame_name,uf->frame_name);
	if (call(write_ufile(uf->frame_name,u_folder,uf)))
		return;
	if (folder[0])
		printf("Saved frame '%s' in folder '%s'\n",uf->frame_name,folder);
	else
		printf("Saved frame '%s' in home folder\n",uf->frame_name);
	uf->dirty = FALSE;
	}


void revert_frame() {
	char yn;

	if (!uf->dirty)
		/* Nothing to revert */
		return;
	printf("Are you sure (y/n): ");
	prep_input();
	scanf(" %c",&yn);
	if (tolower(yn) != 'y')
		return;
	uf->loaded = FALSE;
	if (call(TCL_detach_frame(uf->df)))
		return;
	if (call(TCL_dispose_frame(uf->df)))
		return;
	/* read_ufile creates df */
	if (call(read_ufile(uf->frame_name,u_folder,uf)))
		return;
	if (call(TCL_attach_frame(uf->df))) {
		if (uf->df)
			TCL_dispose_frame(uf->df);
		return;
		}
	uf->loaded = TRUE;
	uf->dirty = FALSE;
	}


void remove_file(char *rfn, char *rfolder) {
	int rc;
	char fn_ddt[64];

	strcpy(fn_ddt,rfolder);
	strcat(fn_ddt,rfn);
	strcat(fn_ddt,".ddt");
	rc = remove(fn_ddt);
	if (!rc)
		printf("Removed file '%s' from disk\n",rfn);
	else
		printf("Failed to remove file '%s' (%d)\n",rfn,rc);
	}
