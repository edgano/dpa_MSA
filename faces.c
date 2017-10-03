/*****************************************************************************
MSA subroutines.
Version 1.0   June 27, 1989

	Program by JD Kececioglu & SF Altschul

	This program calculates pairwise costs using the current cost file.
	Certain alignment positions may be forced by calling fix().
	Terminal gaps may be counted differently than internal gaps.

******************************************************************************/

#include "defs.h"
#include <stdio.h>

extern char	*S [NUMBER+1];
extern int	K,G,GG,Lower,fflag,oflag;
extern int	N [NUMBER+1];
extern int	epsi[NUMBER+1][NUMBER+1];
extern int	scale[NUMBER+1][NUMBER+1];
extern short	D [SIGMA] [SIGMA];
extern int	*vector(), min3();
extern int	*dd [LENGTH+1],   /* forward diagonal   distance */
		*hh [LENGTH+1],   /* forward horizontal distance */
		*vv [LENGTH+1];   /* forward vertical   distance */

typedef struct {
	int	*column;	/* column index of vertex */
	int	width;		/* number of vertices on row in region */
} ROW;

extern ROW	face [NUMBER+1] [NUMBER+1] [LENGTH+1];	/* faces of lattice */
extern int	costs [NUMBER+1] [NUMBER+1]; 	    	/* pairwise costs */
extern void	fix();

static char	*A, *B;

void faces()
{
	register int	i, j, Gi, Gj, q;
	auto  int	d_ [LENGTH+1];      /* reverse diagonal   distance */
        auto  int       h_ [LENGTH+1];	    /* reverse horizontal distance */
	auto int	v_ [LENGTH+1];	    /* reverse vertical   distance */
	auto int	n, m, U, I, J, col [LENGTH+1], w,
			h_l, h_r, d_l, d_r, v_l, v_r;
	register ROW	*row;
	auto     int	C[LENGTH+1];	/* Forced alignment positions	*/
    
	if (oflag) {
		for (i=(K*K-K)/2;i;--i) fprintf(stderr,".");
		fprintf(stderr,"\n");
	}
	if (!fflag) for (i=1;i<=LENGTH;++i) C[i]=0;
	for (Lower=0,I=1;I<K;I++) for (n=N[I],A=S[I],J=I+1;J<=K;J++) {
		if (fflag) fix(I,J,C,n);
		m = N[J];
		B = S[J];

	/* compute distance from <0,0> to <i,j> */

		dd[0][0] = 0; hh[0][0] = vv[0][0] = GG;
		for (j=1;j<=m;j++) {
			vv[0][j] = dd[0][j] = BIG;
			hh[0][j] = hh[0][j-1] + D[DASH][B[j]];
		}
		for (i=1;i<=n;i++) {
			hh[i][0] = dd[i][0] = BIG;
			vv[i][0] = vv[i-1][0] + (C[i] ? BIG : D[A[i]][DASH]);
		}
		for (i=1;i<=n;i++) for (q=C[i], Gi= i==n ? GG:G,j=1;j<=m;j++) {
			Gj= j==m ? GG : G;
			dd[i][j] = min3(dd[i-1][j-1],hh[i-1][j-1],vv[i-1][j-1])
				+ (q && q!=j ? BIG : D[A[i]][B[j]]);
			hh[i][j] = min3(dd[i][j-1]+Gi,hh[i][j-1],vv[i][j-1]+Gi)
				+ D[DASH][B[j]];
			vv[i][j] = min3(dd[i-1][j]+Gj,hh[i-1][j]+Gj,vv[i-1][j])
				+ (q ? BIG : D[A[i]][DASH]);
    		}
		U= (costs[I][J]= min3(dd[n][m],hh[n][m],vv[n][m])) + epsi[I][J];
		Lower += scale[I][J] * costs[I][J];

	/* compute distance from <n,m> to <i,j> */

		d_[m] = 0; h_[m] = v_[m] = GG;
		for ( j = m-1; j >= 0; j-- )
        		v_[j]= (d_[j]= h_[j]= h_[j+1] + D[DASH][B[j+1]]) + G;
    		for (j=w=0;j<m;j++)
			if (min3(hh[n][j]-GG,dd[n][j],vv[n][j])+h_[j]<= U) {
				col[w++] = j;
#if 0
				dist[w++] = h_[j] - GG;
#endif
			}
    		col[w++] = m;
#if 0
		dist[w++] = 0;
#endif
		row = face[I][J] + n;
		row->width = w;
	    	row->column = vector(col,w);
            	for (i=n-1;i>=0;i--) {
			q=C[i+1];
			Gi= i==0 ? GG : G;
			h_r = h_[m]; d_r = d_[m]; v_r = v_[m];
			h_[m] = (d_[m] = v_[m] = v_r + D[A[i+1]][DASH]) + G;
			for ( j = m-1; j >= 0; j-- ) {
	    			Gj= j==0 ? GG : G;
	    			h_l = h_[j]; d_l = d_[j]; v_l = v_[j];
	    			d_[j]= min3(d_r, h_r, v_r) +
					(q && q!=j+1 ? BIG : D[A[i+1]][B[j+1]]);
				h_[j]= min3(d_[j+1]+Gi,h_[j+1],v_[j+1]+Gi) +
					D[DASH][B[j+1]];
	    			v_[j]= min3(d_l + Gj, h_l + Gj, v_l) +
					(q ? BIG : D[A[i+1]][DASH]);
	    			h_r = h_l; d_r = d_l; v_r = v_l;
			}
	       		for (j=w=0;j<=m;j++) {
				Gj= j==0 || j==m ? GG : G;
				if (min3(
				hh[i][j]+min3(h_[j]-Gi,d_[j],v_[j]   ),
				dd[i][j]+min3(h_[j]   ,d_[j],v_[j]   ),
				vv[i][j]+min3(h_[j]   ,d_[j],v_[j]-Gj) ) <=U) {
					col[w++]=j;
#if 0
					dist[w++]=min3(h_[j]-Gi,d_[j],v_[j]-Gj);
#endif
				}
			}
			row = face[I][J] + i; row->width = w;
			row->column = vector(col,w);
		}
		if (oflag) fprintf(stderr,"*");
	}
	if (oflag) fprintf(stderr,"\n");
}
