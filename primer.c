/*****************************************************************************
MSA subroutines.
Version 1.0   June 27, 1989

	Program by JD Kececioglu & SF Altschul

	This program calculates pairwise costs using the current cost file
	and pairwise distances based on number of identities.
	The costs are used to estimate epsilons for the pairs (in ecalc()).
	The distances are used to calculate an evolutionary tree (in bias()).

******************************************************************************/

#include "defs.h"
#include <stdio.h>

extern char	*S[NUMBER+1];
extern int	K,G,GG,fflag,oflag;
extern int	N[NUMBER+1];
extern int	scale[NUMBER+1][NUMBER+1];
extern short	D[SIGMA][SIGMA];
extern int	min3();
extern void	fatal(),fix();
extern int 	costs[NUMBER+1][NUMBER+1];
extern int	Con[NUMBER+1][NUMBER+1][LENGTH+1];
extern int 	*dd [LENGTH+1],   /* forward diagonal   distance */
	    	*hh [LENGTH+1],   /* forward horizontal distance */
    		*vv [LENGTH+1];   /* forward vertical   distance */

static  char	*A, *B;

void primer()
{
	register int 	I,J,i,j,Gi,Gj,n,m,q;
	auto  int	C[LENGTH+1];	/* Forced alignment positions */
    
	if (oflag) {
		for (i=(K*K-K)/2;i;--i) fprintf(stderr,".");
		fprintf(stderr,"\n");
	}
	if (!fflag) for (i=1;i<=LENGTH;++i) C[i]=0;
	for (I=1;I<=K;I++) for (i=N[I];i>=0;i--) Con[I][I][i]=i;
	for (I=1;I<K;I++) for (n=N[I],A=S[I],J=I+1;J<=K;J++) {
		if (fflag) fix(I,J,C,n);
		m=N[J];
		B=S[J];

/* compute distance from <0,0> to <i,j> */

		dd[0][0]=0;
		hh[0][0]=vv[0][0]=GG;
		Con[I][J][0]=Con[J][I][0]=0;
		for (j=1;j<=m;j++) {
			vv[0][j]=dd[0][j]=BIG;
			hh[0][j]=hh[0][j-1]+D[DASH][B[j]];
			Con[J][I][j]= -1;
		}
		for (i=1;i<=n;i++) {
			hh[i][0]=dd[i][0]=BIG;
			vv[i][0]=vv[i-1][0]+(C[i] ? BIG : D[A[i]][DASH]);
			Con[I][J][i]= -1;
		}
		for (i=1;i<=n;i++) for (q=C[i],Gi= i==n?GG:G,j=1;j<=m;j++) {
			Gj= j==m ? GG : G;
			dd[i][j]=min3(dd[i-1][j-1],hh[i-1][j-1],vv[i-1][j-1])+
				(q && q!=j ? BIG : D[A[i]][B[j]]);
			hh[i][j]=min3(dd[i][j-1]+Gi,hh[i][j-1],vv[i][j-1]+Gi)+
				D[DASH][B[j]];
			vv[i][j]=min3(dd[i-1][j]+Gj,hh[i-1][j]+Gj,vv[i-1][j])+
				(q ? BIG : D[A[i]][DASH]);
		}
		costs[I][J] = min3(dd[n][m],hh[n][m],vv[n][m]);
		scale[J][I]=convert(I,J,n,m);
		if (oflag) fprintf(stderr,"*");

#ifdef BIAS2
	if (scale[J][I]<=0) scale[J][I]=1;
#endif

	}
	if (oflag) fprintf(stderr,"\n");
}

/*	Evolutionary distances based on PAM model of molecular evolution */
	
#ifdef PAMDIST
int ED[101]={						0,
	  10,  20,  30,  40,  51,  62,  73,  84,  95, 107,
	 118, 130, 142, 154, 166, 179, 192, 205, 218, 231,
	 245, 259, 273, 287, 302, 317, 332, 348, 364, 380,
	 400, 410, 430, 450, 470, 480, 500, 520, 540, 560,
	 580, 600, 630, 650, 670, 700, 720, 740, 770, 800,
	 820, 850, 880, 910, 940, 970,1010,1040,1080,1120,
	1150,1190,1240,1280,1330,1380,1430,1480,1530,1590,
	1660,1720,1790,1870,1950,2030,2130,2230,2340,2460,
	2590,2730,2890,3080,3280,3520,3800,4200,4700,5300,
	6000,6800,8000,9999,9999,9999,9999,9999,9999,9999};
#endif

/*	Calculate percent mismatch between sequences and convert it to
*	a distance measure for use in calculating an evolutionary tree.	*/

convert(I,J,n,m)
	int I,J,n,m;
{
	int i,j,V,H,M;
	int dir=DIAG;	/* Direction of previous edge in traceback	*/
	int match=0;	/* Number of matches in an optimal alignment	*/
	float f,g;

	for (i=n,j=m; i||j;) {
		V=vv[i][j]-(dir==VERT ? (j==m ? GG:G) : 0);
		H=hh[i][j]-(dir==HORZ ? (i==n ? GG:G) : 0);
		M=min3(V,H,dd[i][j]);
		if	(!j || M==V)	{ dir=VERT; --i; }
		else if	(!i || M==H)	{ dir=HORZ; --j; }
		else {
			dir=DIAG;
			match+= S[I][i]==S[J][j];
			Con[I][J][i]=j;
			Con[J][I][j]=i;
			--i; --j;
		}
	}

#ifdef PAMDIST
	i=f=(100.0*(n-match+m-match))/(n+m);
	g=f-i;
	return((int) (0.5+g*ED[i+1]+(1-g)*ED[i]));
#else
	return((int) (0.5+1000.0*(n-match+m-match)/(n+m)));
#endif

}
