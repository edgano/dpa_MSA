/************************************************************************
MSA subroutines.
Version 1.0   June 27, 1989

	Program by David Lipman

The following code is for inputting and the using "fixed points" in the
multiple alignment so that the user can force certain residues to be aligned.

Each line of input file of fixed points must have following format:

1 2 3  	S 22 21 25 	L 10  	S 35 36 41 	L 14
         
seqs.| "S" precedes block start | "L" precedes block length

*************************************************************************/

#include <stdio.h>
#include <string.h>
#include "defs.h"
#define MAXLINE 250

extern void fatal();
void fill(), fix();

typedef struct block {		/* data structure for fixed points  */
	int i,j,len;		/* of each seq. pair -- i,j are starting */
	struct block *next;	/* points of fixed BLOCK, len is length, and */
} BLOCK;			/* next points to next BLOCK for seq. pair */

BLOCK 	*pFIX[NUMBER + 1][NUMBER + 1];  	/* BLOCKS, indexed by seq. pair */

void fixread(fname)
	char *fname;
{
	FILE *fp,*fopen();
	char s[MAXLINE],*p,fw[10];
	int i,j,len,sqnum,sqnames[NUMBER],sqstarts[NUMBER];	

	if ((fp=fopen(fname,"r"))==NULL) fatal(stderr,"Cannot open %s.",fname);
	for (i=1;i<=NUMBER;++i) for (j=1;j<=NUMBER;++j) pFIX[i][j]=NULL;

	while ((fgets(s,MAXLINE,fp)!=NULL)) {	/*  read in entire line */
		p=s;
		sqnum=0;
		p= strtok(p," ");
		sscanf(p,"%s",fw);	
		while (strcmp(fw,"S")!=0) {	/*get Sq names*/
			p=strtok(NULL," ");
			sqnames[sqnum++]=atoi(fw);
			sscanf(p,"%s",fw);	
		} 

		do {			/* parse string through to end */
			i=0;
			p=strtok(NULL," ");
			sscanf(p,"%s",fw);
			while(strcmp(fw,"L")!=0){	/*get Sq starts*/
				p=strtok(NULL," ");
				sqstarts[i++]=atoi(fw);
				sscanf(p,"%s",fw);
			}

			/* check consistency of block file */

			if (i!=sqnum) fatal("Block starts do not equal sequence number.");
			p=strtok(NULL," ");
			sscanf(p,"%d",&len);		/* get Block length */
			fill(sqnum,sqnames,sqstarts,len);
		} while ((p=strtok(NULL," "))!=NULL);  /* parse string to end */
	}			/*  end loop for each line of input  */
	fclose(fp);
}

/**** fill takes information from input file and derives pFIX structure ****/

void fill(sqnum,sqnames,sqstarts,Len)
	int sqnum,*sqnames,*sqstarts,Len;
{
	int I,J,ni,nj;
	/*char *malloc();*/
	BLOCK *bp1,*bp2;

	for (I=sqnum-1;I;--I) for (J=I-1;J>=0;--J) {
		bp1=(BLOCK *)malloc(sizeof(BLOCK));
		bp2=(BLOCK *)malloc(sizeof(BLOCK));
		bp1->i=bp2->j=sqstarts[I];
		bp1->j=bp2->i=sqstarts[J];
		bp1->len=bp2->len=Len;
		ni=sqnames[I];		nj=sqnames[J];
		bp1->next=pFIX[ni][nj];	bp2->next=pFIX[nj][ni];
		pFIX[ni][nj]=bp1;	pFIX[nj][ni]=bp2;
	}
}

/********* fix uses pFIX structure to "fix" points of alignment***********
********** fix is called from primer() and faces()***********************/

void fix(I,J,C,cl)
	int I,J,*C,cl;
{
	int i,j,k;
	BLOCK *p;

	for (i=1;i<=cl;i++) C[i]=0;
	for (p=pFIX[I][J];p!=NULL;p=p->next) {
		i=p->i;
		j=p->j;
		for (k=0;k<p->len;k++) C[i++]=j++;
	}
}
