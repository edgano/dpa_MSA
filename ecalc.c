/*****************************************************************************
MSA subroutines.
Version 1.0   June 27, 1989

	Program by SF Altschul

	This program estimates epsilons for all pairs of input sequences.
	It does so by constructing an heuristic multiple alignment, and
	calculating the difference between the imposed pairwise costs and
	the optimal pairwise costs.  The heuristic involves finding all
	positions consistent among all optimal pairwise alignments, and
	forcing them into alignment.  For intervening segments, a progessive
	alignment approach is taken.

******************************************************************************/

#include <stdio.h>
#include "defs.h"
#define MINE 5
#define MAXE 50

extern char	*S[NUMBER+1];
extern int	K,G,GG,gflag;
extern int	N[NUMBER+1];
extern int	epsi[NUMBER+1][NUMBER+1];
extern int	scale[NUMBER+1][NUMBER+1];
extern short	D[SIGMA][SIGMA];
extern int	T[3][3][3][3];
extern int 	costs[NUMBER+1][NUMBER+1];
extern int	Con[NUMBER+1][NUMBER+1][LENGTH+1];
extern int	min3();
extern int 	*dd [LENGTH+1],
   		*hh [LENGTH+1],
    		*vv [LENGTH+1];
extern void	fatal();

int Num[36] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
char	*SS[NUMBER+1];
int	NN[NUMBER+1];
int	AL[NUMBER+1][2*LENGTH];
int	NAL[NUMBER+1][2*LENGTH];
int	OAL[NUMBER+1][2*LENGTH];
int	Ord[NUMBER+1];
void	order();

void ecalc(len)
	int len;
{
	register int	I,J,i,j,c,pos,p1,oldp1,test;

	for (I=1;I<=K;++I) for (J=1;J<=K;++J) Con[I][J][N[I]+1]=N[J]+1;
	for (I=1;I<=K;++I) AL[I][0]=DASH;
	pos=oldp1=0;
	for (p1=1;p1<=N[1]+1;++p1) {

/*	Search for positions consistent among all pairwise alignments and
*	align segments of all sequences from last consistent position to
*	current one.   */

		for (I=test=2;test && I<=K;++I) {
			test=((i=Con[1][I][p1])>=0);
			for (J=I+1;test && J<=K;++J)
				test=((j=Con[1][J][p1])>=0 && Con[I][J][i]==j);
		}
		if (test) {

/*	Pick an order in which to construct a progressive multiple alignment  */

			order(oldp1,p1);

/*	Align all segments in the order chosen	*/

			j=align(oldp1,p1,len);

/*	Add the aligned segments to the heuristic alignment	*/

			for (I=1;I<=K;++I) {
				J=Ord[I];
				for (i=1;i<=j;++i) AL[J][pos+i]=S[J][NAL[I][j-i]];
			}

/*	Record whether a given alignment position is consistent among all
*	optimal alignments, or has arisen by a progressive alignment	*/

			for (i=1;i<=j;++i) AL[0][pos+i]=i>K?' ':Num[Ord[i]];
			AL[0][pos+=i]='*';
			for (I=1;I<=K;++I) AL[I][pos]=S[I][Con[1][I][p1]];
			oldp1=p1;
		}
	}

/*	Print out the heuristic multiple alignment	*/

	printf("\n                 ***  Heuristic Multiple Alignment  ***\n\n");
	for (i=(--pos);i>0;i-=LINE) {
		for (I=0;I<=K;++I) {
			for (j=1;j<=LINE&&j<=i;++j) printf("%c",AL[I][pos-i+j]);
			printf("\n");
		}
		printf("\n");
	}

/*	For each pair, calculate difference between imposed alignment cost
*	and optimal alignment cost	*/

	for (I=1;I<K;++I) for (J=I+1;J<=K;++J) {
		c=pcost(I,J,pos)-costs[I][J];
		epsi[I][J]=(c<MINE)?MINE:(c>MAXE?MAXE:c);
	}
}

/*	Subroutine to choose an order in which to construct a progressive
*	multiple alignment	*/

void order(oldp1,p1)
	int oldp1,p1;
{
	int I,J,i,j,i1,j1,PI,PJ,t,mind,p,M,end;
	int test[NUMBER+1];
	int sum[NUMBER+1];
	int dis[NUMBER+1][NUMBER+1];
	int a1[LENGTH+1];
	int a2[LENGTH+1];

	a1[0]=a2[0]=DASH;
	mind=BIG;

/*	Calculate all pairwise costs for the segments in question	*/

	for (I=1;I<K;++I) for (J=I+1;J<=K;++J) {
		PJ=Con[1][J][oldp1];
		PI=Con[1][I][oldp1];
		end=Con[1][I][p1];
		for (p=0,i=PI+1;i<=end;++i) {
			j=(t=Con[I][J][i-1])>=0 ? Con[I][J][i]-t:Con[I][J][i]>0;
			for (M=j>0?j:1;M;--M) {
				a1[++p]= M==1 ? S[I][++PI] : DASH;
				a2[p]  = M<=j ? S[J][++PJ] : DASH;
			}
		}
		for (t=0,i=1;i<p;++i) t+=D[a1[i]][a2[i]]+T[a1[i-1]!=DASH]
			[a2[i-1]!=DASH] [a1[i]!=DASH] [a2[i]!=DASH];
		if (t<mind) { mind=t; i1=I; j1=J; }
		dis[I][J]=dis[J][I]=t;
	}

/*	Make lowest cost pair the first two segments in the order	*/

	for (I=1;I<=K;++I) test[I]=1;
	Ord[1]=i1; SS[1]=S[i1]; NN[1]=N[i1]; test[i1]=0;
	Ord[2]=j1; SS[2]=S[j1]; NN[2]=N[j1]; test[j1]=0;
	for (I=1;I<=K;++I) sum[I]=dis[I][i1]+dis[I][j1];

/*	Fill out the order using average distances	*/

	for (j=3;j<=K;++j) {
		mind=BIG;
		for (I=1;I<=K;++I) if (test[I] && sum[I]<mind) mind=sum[i=I];
		Ord[j]=i; SS[j]=S[i]; NN[j]=N[i]; test[i]=0;
		for (I=1;I<=K;++I) sum[I]+=dis[I][i];
	}
}

/*	Subroutine to calculate an optimal progressive
*	alignment of the segments in the order chosen.		*/

static  char	*A;

align(oldp1,p1,len)
	int oldp1,p1,len;
{
	int I,J,i,low,m;

	for (I=1;I<=K;++I) OAL[I][0]=0;
	low=Con[1][Ord[1]][oldp1];
	A=SS[1]+low;
	m=Con[1][Ord[1]][p1]-low-1;
	for (i=1;i<=m;++i) NAL[1][m-i]=low+i;
	for (I=2;I<=K;++I) {
		for (J=1;J<I;++J) for (i=1;i<=m;++i) OAL[J][i]=NAL[J][m-i];
		low=Con[1][Ord[I]][oldp1];
		A=SS[I]+low;
		m=newal(I,low,Con[1][Ord[I]][p1]-low-1,m);
		if (m>len && I<K) fatal("Heuristic alignment segment is longer than %d.",len);
	}
	return(m);
}

/*	Subroutine to calculate the optimal alignment of a new segment and
*	a multiple alignment	*/

newal(I,low,n,m)
	int I,low,n,m;
{
	register int 	i,j,k,x,q,qq,dg,vg,hg,sum,sum2;
	int		d[2*LENGTH],v[2*LENGTH],h[2*LENGTH];
	int		dn[2*LENGTH],vn[2*LENGTH],hn[2*LENGTH];
	int		sc[NUMBER+1];
    
	for (sum=sum2=0,k=1;k<I;++k) {
		sum+= sc[k] = Ord[k]<Ord[I] ?
			scale[Ord[k]][Ord[I]] : scale[Ord[I]][Ord[k]];
		if (OAL[k][1]) sum2+= sc[k];
	}
	for (j=0;j<=m;j++) hn[j]=dn[j]=vn[j]=BIG;
	for (i=0;i<=n;i++) {
		for (j=0;j<=m;j++) {
			h[j]=hn[j]; d[j]=dn[j]; v[j]=vn[j];
		}
		if (i==0) {
			dn[0]=0;
			vn[0]=sum* (low?G:GG);
			hn[0]=sum2*(low?G:GG);
		}
		else {
			vn[0]=v[0]+sum*D[A[i]][DASH];
			hn[0]=dn[0]=BIG;
		}
		vv[i][0]=dd[i][0]=hh[i][0]=1;

/*	Calculate optimal alignment in the case that terminal gap costs
*	are different than internal gap costs				*/

		if (gflag) for (j=1;j<=m;j++) {
			dg=d[j-1]; hg=h[j-1];
			for (x=0,k=1;k<I;++k) {
				q=OAL[k][j-1];
				qq=OAL[k][j];
				if (q<NN[k]) {
					dg += sc[k] * T[1][q>0][1][qq>0];
					if (q || i+low>1) hg +=
						sc[k]*T[0][q>0][1][qq>0];
				}
				x+= sc[k] * D [SS[k][qq]] [A[i]];
			}
			dn[j]=x+(k=min3(dg,hg,v[j-1]));
			dd[i][j]= k==dg ? DIAG : (k==hg ? HORZ : VERT);

			dg=d[j]; hg=h[j];
			for (k=1;k<I;++k) if ((qq=OAL[k][j])<NN[k]) {
				dg += sc[k] * T[1][qq>0][1][0];
				if (qq || i+low>1) hg+=sc[k] * T[0][qq>0][1][0];
			}
			vn[j]=sum*D[DASH][A[i]]+(k=min3(dg,hg,v[j]));
                   /* A. A. Schaffer fixed bug reported by A. Ropelewski*/
                        if (k > BIG) {
			  fprintf(stderr, "\n BIG in defs.h is not big enough");
			  exit(1);
			}

			vv[i][j]= k==dg ? DIAG : (k==hg ? HORZ : VERT);

			dg=dn[j-1]; vg=vn[j-1]; hg=hn[j-1];
			for (x=0,k=1;k<I;++k) {
				q=OAL[k][j-1];
				qq=OAL[k][j];
				if (qq>1) hg+=sc[k] * T [0] [q>0] [0] [1];
				if (low+i<NN[I]) {
					dg += sc[k] * T [1] [q>0] [0] [qq>0];
					vg += sc[k] * T [1] [0]   [0] [qq>0];
				}
				x += sc[k] * D[SS[k][qq]] [DASH];
			}
			hn[j]=x+(k=min3(dg,hg,vg));
                   /* A. A. Schaffer fixed bug reported by A. Ropelewski*/
                        if (k > BIG){
			  fprintf(stderr, "\n BIG in defs.h is not big enough");
			  exit(1);
			}
			hh[i][j]= k==dg ? DIAG : (k==hg ? HORZ : VERT);
		}

/*	Calculate optimal alignment in the case that terminal gap costs
*	are the same as internal gap costs				*/

		else for (j=1;j<=m;j++) {
			dg=d[j-1]; hg=h[j-1];
			for (x=0,k=1;k<I;++k) {
				q=OAL[k][j-1];
				qq=OAL[k][j];
				dg += sc[k] * T [1] [q>0] [1] [qq>0];
				hg += sc[k] * T [0] [q>0] [1] [qq>0];
				x  += sc[k] * D [SS[k][qq]] [A[i]];
			}
			dn[j]=x+(k=min3(dg,hg,v[j-1]));
			dd[i][j]= k==dg ? DIAG : (k==hg ? HORZ : VERT);

			dg=d[j]; hg=h[j];
			for (k=1;k<I;++k) {
				qq=OAL[k][j];
				dg += sc[k] * T [1] [qq>0] [1] [0];
				hg += sc[k] * T [0] [qq>0] [1] [0];
			}
			vn[j]=sum*D[DASH][A[i]]+(k=min3(dg,hg,v[j]));
			vv[i][j]= k==dg ? DIAG : (k==hg ? HORZ : VERT);

			dg=dn[j-1]; vg=vn[j-1]; hg=hn[j-1];
			for (x=0,k=1;k<I;++k) {
				q=OAL[k][j-1];
				qq=OAL[k][j];
				dg += sc[k] * T [1] [q>0] [0] [qq>0];
				vg += sc[k] * T [1] [0]   [0] [qq>0];
				hg += sc[k] * T [0] [q>0] [0] [qq>0];
				x  += sc[k] * D [SS[k][qq]] [DASH];
			}
			hn[j]=x+(k=min3(dg,hg,vg));
			hh[i][j]= k==dg ? DIAG : (k==hg ? HORZ : VERT);
		}
	}

/*	Traceback to reconstruct optimal alignment	*/

	j=0;
	k= dn[m]<=vn[m] ? (dn[m]<=hn[m]?DIAG:HORZ) : (vn[m]<=hn[m]?VERT:HORZ);
	while ( n || m)
		if (k==DIAG) {
			for (i=1;i<I;++i) NAL[i][j]=OAL[i][m];
			NAL[I][j++]=low+n;
			k=dd[n--][m--];
		}
		else if (k==VERT) {
			for (i=1;i<I;++i) NAL[i][j]=0;
			NAL[I][j++]=low+n;
			k=vv[n--][m];
		}
		else {
			for (i=1;i<I;++i) NAL[i][j]=OAL[i][m];
			NAL[I][j++]=0;
			k=hh[n][m--];
		}
	return(j);
}

/*	Subroutine to calculate imposed cost for any pair of sequences	*/

pcost(I,J,b)
	int I,J,b;
{
	int i;
	int s=0;

	for (i=1;i<=b;++i) s+=D[AL[I][i]][AL[J][i]]+T[AL[I][i-1]!=DASH]
		[AL[J][i-1]!=DASH] [AL[I][i]!=DASH] [AL[J][i]!=DASH];
	i=1;
	while (AL[I][i]==DASH && AL[J][i]==DASH) ++i;
	if (AL[I][i]==DASH || AL[J][i]==DASH) s+=(GG-G);
	i=b;
	while (AL[I][i]==DASH && AL[J][i]==DASH) --i;
	if (AL[I][i]==DASH || AL[J][i]==DASH) s+=(GG-G);
	return(s);
}
