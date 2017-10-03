/************************ Multiple Sequence Alignment **************************

	Version 2.0     November 22, 1994
	Program (version 1.0, June 1989)
        by  JD Kececioglu, SF Altschul, DJ Lipman & R Miner

        Improvements (from 1.0 to 2.0) by SK Gupta, AA Schaffer,
        with some guidance from JD Kececioglu

        Please cite papers 6 and 7 below if you wish to cite the software
        package MSA.

See:
    1. Carrillo & Lipman, "The Multiple Sequence Alignment Problem in Biology",
		SIAM J. Appl. Math. 48 (1988) 1073-1082;
    2. Altschul & Lipman, "Trees, Stars, and Multiple Biological Sequence
		Alignment", SIAM J. Appl. Math. 49 (1989) 197-209;
    3. Altschul, "Gap Costs for Multiple Sequence Alignment",
		J. Theor. Biol. 138 (1989) 297-309;
    4. Altschul, Carroll & Lipman, "Weights for Data Related by a Tree",
		J. Molec. Biol. 207 (1989) 647-653;
    5. Altschul, "Leaf Pairs and Tree Dissections",
		SIAM J. Discrete Math. 2 (1989) 293-299;
    6. Lipman, Altschul & Kececioglu, "A Tool for Multiple Sequence Alignment",
		Proc. Natl. Acad. Sci. USA 86 (1989) 4412-4415.
    7. Gupta, Kececioglu & Schaffer, "Improving the Practical Time
                and Space Efficiency of the Shortest-Paths Approach to
                Sum-of-Pairs Multiple Sequence Alignment", J. Computational
                Biology 2(1995) 459-472.

	Computes an optimal multiple alignment within a lattice defined by
	the join of vertices from two-dimensional path graphs.  There
	is one such path graph for each pair of sequences, and all vertices
	contained in paths whose cost is within an epsilon of the minimal
	are included.  Features include:

	The program computes a minimal SP alignment (Sum of Pairwise costs).
	Each pairwise cost may be given a different weight.  These weights
	may be calculated from an evolutionary tree using either of two
	rationales (Altschul et al., "Weights for Data Related by a Tree",
	J. Molec. Biol. 208) or equal weights may be specified.  The
	evolutionary tree is estimated from pairwise distances using the
	Neighbor Joining method (Saitou & Nei, Mol. Biol. Evol. 4:406-425).

	Epsilons for each sequence pair may be input by the user or
	estimated.  An heuristic multiple alignment is computed and
	the epsilons are taken to be the difference between the
	imposed and minimal pairwise costs.

	Gap costs are of the affine type generalized for multiple sequence
	alignments (Altschul, "Gap Costs for Multiple Sequence Alignment",
	J. Theor. Biol. 138:297-309).  The user may specify whether terminal
	gaps are counted.

	The user may specify residues in any of the sequences to be
	forced into alignment.

*******************************************************************************/

#include <stdio.h>
#include "defs.h"
#include <sys/time.h>

#define TRUE        1	    	/* boolean truth */
#define FALSE       0	    	/* boolean falsity */
#define	FREE(x)		free((char *) x)

typedef struct coordinate COORDINATE;
typedef struct coordinate_values COORDINATE_VALUES;
typedef struct vertex VERTEX;
typedef struct edge EDGE;
typedef struct heap HEAP;

/* array for accessing vertices by lattice coordinate */
struct coordinate {
  int	lo, hi;		/* lower and upper limits on array indices */
  COORDINATE_VALUES *coord_vals;	/* next     coordinate array indices */
  COORDINATE_VALUES *prev_coord_val;/* previous coordinate array index */
  COORDINATE *next_on_free_list; /*maintains available records*/
  int	refer;		/* how many valid coordinate values do I have */
};

/* index in array */
struct coordinate_values {
  COORDINATE	*next_coord;	/* next    coordinate array */
  COORDINATE	*curr_coord;	/* current coordinate array */
/*  short value;          value of this coordinate in absolute terms*/
};

/* lattice vertex */
struct vertex {
  EDGE *out;		/* outgoing edge adjacency list */
  COORDINATE_VALUES	*prev_coord_val;	/* father in array */
  EDGE *nonextracted_inedges;/* incoming edges still not extracted from heap */
};

/* lattice edge */
struct edge {
	VERTEX	*tail, *head;	/* edge tail and head vertices */
	int	dist;		/* distance to head from source along edge */
        int     refer;          /* how many backtrack edges point to me */
	EDGE *next, *prev;	/* edge adjacency list links */
	EDGE *heap_succ, *heap_pred;	/* heap bucket links */
	EDGE *nonextracted_next,
             *nonextracted_prev;   /* nonextracted_inedges links */
	EDGE *backtrack;	/* edge to previous edge in path */
};

/* discrete heap of edges ranked by distance */
struct heap {
	EDGE	**min, **max;	/* minimum and maximum buckets of heap */
	EDGE	*bucket[1];	/* buckets of edges */
};

char	*S [NUMBER+1];			/* sequences */
int	K;				/* number of sequences */
int	N [NUMBER+1];			/* lengths of sequences */
int	delta;				/* upper and lower bound difference */
int	epsi[NUMBER+1][NUMBER+1];	/* projected & pairwise cost diff. */
int	scale[NUMBER+1][NUMBER+1];	/* pairwise cost weight scale */
int	Con[NUMBER+1][NUMBER+1][LENGTH+1];	/* consistency check */
int	proj [NUMBER+1] [NUMBER+1];		/* projected costs */
char	cname[20]="pam250";			/* name of cost file */
int	Upper, Lower;		/* upper and lower bounds on alignment distance */
int	bflag,gflag,fflag,oflag;
int 	*dd [LENGTH+1],   		/* forward diagonal   distance */
	*hh [LENGTH+1],   		/* forward horizontal distance */
	*vv [LENGTH+1];   		/* forward vertical   distance */
VERTEX  *presource;              /* Vertex before source; tail of first edge */

void	main(), data(), faces(), display(), column(), output(),
coord(), project(), free_vertex(), free_edge(), fatal(), primer(),
safe_coord(),
ecalc(), bias(), fixread(), adjacent();
char	*alloc();
int	intersect(), cost(), min3(), *vector(), free_coordinate();
FILE	*must_open();
VERTEX	*source(), *sink(), *create_vertex();
EDGE	*msa(), *extract(), *create_edge();
COORDINATE	*create_coordinate();
HEAP	*heap();

/*************************   FLAGS   *******************************************

	-m		suppresses computation of optimal multiple alignment
	-g		penalizes terminal gaps
	-b		sets all pairwise weights to 1 (does not compute tree)
	-o		suppresses status reports to stderr
-d<delta>	 	user specified delta (upper bound for total alignment
			cost)
-e <epsilon file>	user specified epsilons for each sequence pair
			(does not compute heuristic alignment)
-c <cost file>		user specified cost file (default is pam250)
-f <fixer file>		allows user to force residues in alignment (see
			comment in fixer.c for file format)
<input filename>	sequences to be aligned

**************************   FLAGS   *****************************************/

#define	USAGE   "Usage:  msa [-m] [-g] [-b] [-o] [-d<delta>] [-e <epsilon file>] [-c <cost file>] [-f <fixer file>] <input filename>"

struct timeval starttime, endtime;
double deltatime;

void main(argc, argv)
int	argc;
char	*argv[];
{
	int i,j,len,size;
	int eflag=1;
	int mflag=1;
	char ename[FILE_NAME_SIZE];
	char fname[FILE_NAME_SIZE],Fname[FILE_NAME_SIZE];
	FILE    *stream, *efile, *fopen();

	gettimeofday(&starttime, NULL);  /* Get start time */
	/* process arguments */
	oflag=bflag=gflag=1;
	fflag=0;
	delta = -1;
	stream = NULL;
	if ( argc > MAX_ARGS )
          fatal(USAGE);
	for ( ++argv; --argc; argv++ )
          if ( (*argv)[0] == '-' )
	    if ( (*argv)[1] == 'm' )
              mflag=0;
	    else if ( (*argv)[1] == 'g' ) gflag = 0;
	    else if ( (*argv)[1] == 'b' ) bflag = 0;
	    else if ( (*argv)[1] == 'o' ) oflag = 0;
	    else if ( (*argv)[1] == 'd' ) delta = atoi(&(*argv)[2]);
	    else if ( (*argv)[1] == 'e' ) {
	      ++argv;
	      --argc;
	      sscanf(*argv,"%s",ename);
	      eflag=0;
	    }
	    else if ( (*argv)[1] == 'c' ) {
	      ++argv;
	      --argc;
	      sscanf(*argv,"%s",cname);
	      stream = must_open(cname,"r");
	    }
	    else if ( (*argv)[1] == 'f' ) {
	      ++argv;
	      --argc;
	      sscanf(*argv,"%s",fname);
	      fflag=1;
	    }
	    else fatal(USAGE);
	  else
            sscanf(*argv,"%s",Fname);

	data(stream,must_open(Fname,"r"));    /* get input sequences */
	if (fflag)
          fixread(fname);
	if (!eflag) {
		efile = must_open(ename,"r");
		for (i=1;i<K;++i)
                  for (j=i+1;j<=K;++j) {
			fscanf(efile,"%d",&epsi[i][j]);
			if (epsi[i][j]<0)
                          fatal("Epsilon must be positive.");
		  }
	}

	/* allocate memory */

	for (len=i=1;i<=K;++i)
          if (N[i]>len)
            len=N[i];
	++len;
	size=sizeof(int)*len*len*(1+eflag);
	dd[0] = (int *) alloc(size);
	hh[0] = (int *) alloc(size);
	vv[0] = (int *) alloc(size);
	for (i=1;i<len;++i) {
		dd[i]=dd[i-1]+(1+eflag)*len;
		hh[i]=hh[i-1]+(1+eflag)*len;
		vv[i]=vv[i-1]+(1+eflag)*len;
	}

	/* compute multiple sequence alignment */

	if (K==2)
          bflag=0;
	if (!bflag)
          for (i=1;i<K;++i)
            for (j=i+1;j<=K;++j)
              scale[i][j]=1;
	if (bflag || eflag) {
		if (oflag)
                  fprintf(stderr,"Calculating pairwise alignments.\n");
		primer();
	}
	if (bflag) {
		if (oflag)
                  fprintf(stderr,"Calculating weights.\n");
		bias();
	}
	if (eflag) {
		if (oflag)
                  fprintf(stderr,"Calculating epsilons.\n");
		ecalc(2*len-2);
		if (oflag) {

			fprintf(stderr,"----Estimated epsilons----\n");
			for (i=1;i<K;++i)
                for (j=i+1;j<=K;++j)
                    fprintf(stderr,"I =%2d  J =%2d  epsilon = %2d\n",i,j,epsi[i][j]);

		}
	}
	if (mflag) {
		if (oflag)
                  fprintf(stderr,"Calculating pairwise projection costs.\n");
		faces();
		FREE(vv[0]);
                FREE(hh[0]);
                FREE(dd[0]);
		if (delta<0)
                  for (delta=0,i=1;i<K;++i)
                    for (j=i+1;j<=K;++j)
			delta+=(scale[i][j]*epsi[i][j]);
		Upper=Lower+delta;
		for (i=1;i<K;++i)
                  for (j=i+1;j<=K;++j) {
                    proj[i][j]=0;
		    scale[j][i] = scale[i][j]; /*make symmetric for msa() */
		  }
		if (oflag)
                  fprintf(stderr,"Calculating multiple alignment.\n");
		display(msa());
	}
	gettimeofday(&endtime, NULL);  /* Get end time */
	deltatime = (endtime.tv_sec - starttime.tv_sec)
	     + (endtime.tv_usec - starttime.tv_usec)/1000000.0;
	printf("Elapsed time = %7.3f\n",deltatime);
}

/*Define strings in two different ways depending on whether ANSI*/
/*standards are used or not*/
#if defined(__STDC__)
#define strng(a) #a
#else
#define strng(a) "a"
#endif

#define SUB(a,b)  D[strng(a)[0]][strng(b)[0]] = D[strng(b)[0]][strng(a)[0]]
#define DAG(a)    D[strng(a)[0]][strng(a)[0]]

short	D [SIGMA] [SIGMA];	/* symbol distance */
int	T[3][3][3][3];		/* Altschul gap counts */
int     *Tpointer[NUMBER+1][NUMBER+1][3]; /* cuts off first two dimensions of T*/
int	G,GG;			/* gap cost */

/* data :   read input data */
void data (stream,Ffile)
FILE	*stream,*Ffile;
{
	static	 	char	buffer [LENGTH+1];
	auto		char	a, b;
	auto		int	c, n;
        auto            int     i, j;
	/* read sequences */
	K=0;
	while((n=getseq(buffer,LENGTH+1,Ffile))){
		if (++K>NUMBER)
                  fatal("Cannot exceed %d sequences.",NUMBER);
		N[K] = n;
		S[K] = alloc(n + 2); /* initial DASH and terminal '\0' */
		S[K][0] = DASH;
		strcpy(&S[K][1],buffer);
	}

       /* initialize gap cost, gap count table,
	and positive symmetric integer distance table */

	if ( stream ) {

		/* distance file contains gap cost
		followed by triples of the form <a, b, d(a,b)> */

		fscanf(stream," %d\n",&G);
		while (fscanf(stream," %c %c %d\n",&a,&b,&c)==3)
			D[a][b] = D[b][a] = c;
		fclose(stream);
	}
	else {           /* default is dayhoff matrix and gap cost 8	*/
		G = 8;
		DAG(-) = 0;
		DAG(W)=0;
		DAG(Y)=7;
		DAG(F)=8;
		DAG(V)=13;
		DAG(L)=11;
		DAG(I)=12;
		DAG(M)=11;
		DAG(K)=12;
		DAG(R)=11;
		DAG(H)=11;
		DAG(Q)=13;
		DAG(E)=13;
		DAG(D)=13;
		DAG(N)=15;
		DAG(G)=12;
		DAG(A)=15;
		DAG(P)=11;
		DAG(T)=14;
		DAG(S)=15;
		DAG(C)=5;
		SUB(-,C) = 12;
		SUB(-,S) = 12;
		SUB(-,T) = 12;
		SUB(-,P) = 12;
		SUB(-,A) = 12;
		SUB(-,G) = 12;
		SUB(-,N) = 12;
		SUB(-,D) = 12;
		SUB(-,E) = 12;
		SUB(-,Q) = 12;
		SUB(-,H) = 12;
		SUB(-,R) = 12;
		SUB(-,K) = 12;
		SUB(-,M) = 12;
		SUB(-,I) = 12;
		SUB(-,L) = 12;
		SUB(-,V) = 12;
		SUB(-,F) = 12;
		SUB(-,Y) = 12;
		SUB(-,W) = 12;
		SUB(W,C) = 25;
		SUB(W,S) = 19;
		SUB(W,T) = 22;
		SUB(W,P) = 23;
		SUB(W,A) = 23;
		SUB(W,G) = 24;
		SUB(W,N) = 21;
		SUB(W,D) = 24;
		SUB(W,E) = 24;
		SUB(W,Q) = 22;
		SUB(W,H) = 20;
		SUB(W,R) = 15;
		SUB(W,K) = 20;
		SUB(W,M) = 21;
		SUB(W,I) = 22;
		SUB(W,L) = 19;
		SUB(W,V) = 23;
		SUB(W,F) = 17;
		SUB(W,Y) = 17;
		SUB(Y,C) = 17;
		SUB(Y,S) = 20;
		SUB(Y,T) = 20;
		SUB(Y,P) = 22;
		SUB(Y,A) = 20;
		SUB(Y,G) = 22;
		SUB(Y,N) = 19;
		SUB(Y,D) = 21;
		SUB(Y,E) = 21;
		SUB(Y,Q) = 21;
		SUB(Y,H) = 17;
		SUB(Y,R) = 21;
		SUB(Y,K) = 21;
		SUB(Y,M) = 19;
		SUB(Y,I) = 18;
		SUB(Y,L) = 18;
		SUB(Y,V) = 19;
		SUB(Y,F) = 10;
		SUB(F,C) = 21;
		SUB(F,S) = 20;
		SUB(F,T) = 20;
		SUB(F,P) = 22;
		SUB(F,A) = 21;
		SUB(F,G) = 22;
		SUB(F,N) = 21;
		SUB(F,D) = 23;
		SUB(F,E) = 22;
		SUB(F,Q) = 22;
		SUB(F,H) = 19;
		SUB(F,R) = 21;
		SUB(F,K) = 22;
		SUB(F,M) = 17;
		SUB(F,I) = 16;
		SUB(F,L) = 15;
		SUB(F,V) = 18;
		SUB(V,C) = 19;
		SUB(V,S) = 18;
		SUB(V,T) = 17;
		SUB(V,P) = 18;
		SUB(V,A) = 17;
		SUB(V,G) = 18;
		SUB(V,N) = 19;
		SUB(V,D) = 19;
		SUB(V,E) = 19;
		SUB(V,Q) = 19;
		SUB(V,H) = 19;
		SUB(V,R) = 19;
		SUB(V,K) = 19;
		SUB(V,M) = 15;
		SUB(V,I) = 13;
		SUB(V,L) = 15;
		SUB(L,C) = 23;
		SUB(L,S) = 20;
		SUB(L,T) = 19;
		SUB(L,P) = 20;
		SUB(L,A) = 19;
		SUB(L,G) = 21;
		SUB(L,N) = 20;
		SUB(L,D) = 21;
		SUB(L,E) = 20;
		SUB(L,Q) = 19;
		SUB(L,H) = 19;
		SUB(L,R) = 20;
		SUB(L,K) = 20;
		SUB(L,M) = 13;
		SUB(L,I) = 15;
		SUB(I,C) = 19;
		SUB(I,S) = 18;
		SUB(I,T) = 17;
		SUB(I,P) = 19;
		SUB(I,A) = 18;
		SUB(I,G) = 20;
		SUB(I,N) = 19;
		SUB(I,D) = 19;
		SUB(I,E) = 19;
		SUB(I,Q) = 19;
		SUB(I,H) = 19;
		SUB(I,R) = 19;
		SUB(I,K) = 19;
		SUB(I,M) = 15;
		SUB(M,C) = 22;
		SUB(M,S) = 19;
		SUB(M,T) = 18;
		SUB(M,P) = 19;
		SUB(M,A) = 18;
		SUB(M,G) = 20;
		SUB(M,N) = 19;
		SUB(M,D) = 20;
		SUB(M,E) = 19;
		SUB(M,Q) = 18;
		SUB(M,H) = 19;
		SUB(M,R) = 17;
		SUB(M,K) = 17;
		SUB(K,C) = 22;
		SUB(K,S) = 17;
		SUB(K,T) = 17;
		SUB(K,P) = 18;
		SUB(K,A) = 18;
		SUB(K,G) = 19;
		SUB(K,N) = 16;
		SUB(K,D) = 17;
		SUB(K,E) = 17;
		SUB(K,Q) = 16;
		SUB(K,H) = 17;
		SUB(K,R) = 14;
		SUB(R,C) = 21;
		SUB(R,S) = 17;
		SUB(R,T) = 18;
		SUB(R,P) = 17;
		SUB(R,A) = 19;
		SUB(R,G) = 20;
		SUB(R,N) = 17;
		SUB(R,D) = 18;
		SUB(R,E) = 18;
		SUB(R,Q) = 16;
		SUB(R,H) = 15;
		SUB(H,C) = 20;
		SUB(H,S) = 18;
		SUB(H,T) = 18;
		SUB(H,P) = 17;
		SUB(H,A) = 18;
		SUB(H,G) = 19;
		SUB(H,N) = 15;
		SUB(H,D) = 16;
		SUB(H,E) = 16;
		SUB(H,Q) = 14;
		SUB(Q,C) = 22;
		SUB(Q,S) = 18;
		SUB(Q,T) = 18;
		SUB(Q,P) = 17;
		SUB(Q,A) = 17;
		SUB(Q,G) = 18;
		SUB(Q,N) = 16;
		SUB(Q,D) = 15;
		SUB(Q,E) = 15;
		SUB(E,C) = 22;
		SUB(E,S) = 17;
		SUB(E,T) = 17;
		SUB(E,P) = 18;
		SUB(E,A) = 17;
		SUB(E,G) = 17;
		SUB(E,N) = 16;
		SUB(E,D) = 14;
		SUB(D,C) = 22;
		SUB(D,S) = 17;
		SUB(D,T) = 17;
		SUB(D,P) = 18;
		SUB(D,A) = 17;
		SUB(D,G) = 16;
		SUB(D,N) = 15;
		SUB(N,C) = 21;
		SUB(N,S) = 16;
		SUB(N,T) = 17;
		SUB(N,P) = 18;
		SUB(N,A) = 17;
		SUB(N,G) = 17;
		SUB(G,C) = 20;
		SUB(G,S) = 16;
		SUB(G,T) = 17;
		SUB(G,P) = 18;
		SUB(G,A) = 16;
		SUB(A,C) = 19;
		SUB(A,S) = 16;
		SUB(A,T) = 16;
		SUB(A,P) = 16;
		SUB(P,C) = 20;
		SUB(P,S) = 16;
		SUB(P,T) = 17;
		SUB(T,C) = 19;
		SUB(T,S) = 16;
		SUB(S,C) = 17;
	}
	GG= gflag ? 0 : G;
	T[0][0][0][0] = 0; /* -- : -- */
	T[0][0][0][1] = G; /* -- : -x */
	T[0][0][1][0] = G; /* -x : -- */
	T[0][0][1][1] = 0; /* -x : -x */
	T[0][1][0][0] = 0; /* -- : x- */
	T[0][1][0][1] = 0; /* -- : xx */
	T[0][1][1][0] = G; /* -x : x- */
	T[0][1][1][1] = 0; /* -x : xx */
	T[1][0][0][0] = 0; /* x- : -- */
	T[1][0][0][1] = G; /* x- : -x */
	T[1][0][1][0] = 0; /* xx : -- */
	T[1][0][1][1] = 0; /* xx : -x */
	T[1][1][0][0] = 0; /* x- : x- */
	T[1][1][0][1] = G; /* x- : xx */
	T[1][1][1][0] = G; /* xx : x- */
	T[1][1][1][1] = 0; /* xx : xx */
	T[2][0][2][0] = 0;
	T[2][0][2][1] = 0;
	T[2][1][2][0] = 0;
	T[2][1][2][1] = 0;
	T[0][2][0][2] = 0;
	T[0][2][1][2] = 0;
	T[1][2][0][2] = 0;
	T[1][2][1][2] = 0;
	T[2][2][2][2] = 0;
}

getseq(seq,maxl,fp)
char *seq;
int  maxl;
FILE *fp;
{
	char buf[256];
	int n,j,l;

	n=0;
	fgets(buf,256,fp);
	while (fgets(buf,256,fp) && buf[0]!='>') {
		l=strlen(buf);
		for (j=0;j<l;j++)
                  if (buf[j] >= 'A' && buf[j] <= 'Z') {
			seq[n]=buf[j];
			if (++n > maxl)
                          fatal("Cannot exceed %d characters per sequence.",maxl-1);
		  }
	}
	seq[n]='\0';
	if (buf[0]=='>')
          fseek(fp,-strlen(buf),1);
	return(n);
}

/* heap insertion and deletion */

#define INSERT(e,h)\
	{ register EDGE **b = (h)->bucket + (e)->dist;\
	  if (*b != NULL)\
            (*b)->heap_pred = (e);\
	  (e)->heap_succ = *b; (e)->heap_pred = NULL;\
          *b = (e); }
#define DELETE(e,h)\
	{ register EDGE **b = (h)->bucket + (e)->dist;\
	  if ((e)->heap_pred != NULL)\
            (e)->heap_pred->heap_succ = (e)->heap_succ;\
	  else\
            *b = (e)->heap_succ;\
	  if ((e)->heap_succ != NULL)\
            (e)->heap_succ->heap_pred = (e)->heap_pred; }

HEAP	*h;

/* msa :    compute multiple sequence alignment */

EDGE *msa ()
{
	static		int	p [NUMBER+1], q [NUMBER+1], r [NUMBER+1];
	register	int	d,inc;
	register	int	ccost=0;
	register	VERTEX	*v, *w, *t;
	auto		VERTEX	*s;
	register	EDGE	*e, *f;
	auto		char	C [NUMBER+1];
	register	int	I, J;
	auto		int	delta0[NUMBER+1],delta1[NUMBER+1];
        register        int     difference;
        auto            int     ends[NUMBER];
        auto            int     endcount, endindex;

	/* compute shortest paths to vertices in intersected region of lattice */

	s = source();
	t = sink();
	h = heap(Upper);
        presource = create_vertex(NULL);
	e = create_edge(presource,s);
	e->dist = 0;
        e->refer++;          /* Make sure edge does not get freed */
        e->backtrack = NULL;
	INSERT(e,h);
	if (oflag) {
		fprintf(stderr,"....1....2....3....4....5....6....7....8....9....0\n");
		inc=1+Upper/50;
	}
	while ((e=extract()) != NULL && (v=e->head) != t) {
		if (oflag && e->dist>ccost) {
			fprintf(stderr,"*");
			ccost+=inc;
		}
		if (e->dist <= Upper) {

              /*put coordiantes of tail into p array
                and coordinates of head of edge into q*/
	          coord(e->tail,p);
	          safe_coord(e->head,q);

             /*next loop is from cost function
                    difference between p and q*/
		  for (I=1;I<=K;I++)
		    delta0[I] = q[I] - p[I];
		  if(gflag) {
		    endcount = 0;
		    for (I=1;I<=K;I++)
		      if (q[I]==0||q[I]==N[I]) {
			delta0[I] = 2;
			ends[endcount++] = I;
		      }
                  }

		  for (I=2;I<=K;I++)
		    for (J=1;J<I;J++){
		      Tpointer[I][J][0] = T[delta0[I]][delta0[J]][0];
		      Tpointer[I][J][1] = T[delta0[I]][delta0[J]][1];
		      Tpointer[I][J][2] = T[delta0[I]][delta0[J]][2];
		    }

                  if(v->out == NULL)  /* If first time visiting v */
                    adjacent(e,q);
                  for(f=v->out; f!=NULL; f=f->next) {
                    difference = f->dist - e->dist;
                    if (difference > 0) {
		  /*get coordinates of next vertex into r*/
                      safe_coord(f->head,r);

	for (I=1;I<=K;I++) {
		C[I] = (r[I]>q[I] ? S[I][r[I]] : DASH);
		delta1[I] = r[I] - q[I];
	}
	if(gflag)
/*          for (I=1;I<=K;I++)
            if (q[I]==0||q[I]==N[I])
              delta1[I] = 2;*/
          for(endindex = 0; endindex < endcount; endindex++)
	    delta1[ends[endindex]] = 2;
        d = 0;
        d+=scale[1][2]*(D[C[1]][C[2]]+Tpointer[2][1][delta1[2]][delta1[1]]);
	for (I=K;I>=3;I--){
          if (d >= difference)
	    goto nextedge;
          for (J=1;J<I;J++)
	    d+=scale[I][J]*(D[C[I]][C[J]]+Tpointer[I][J][delta1[I]][delta1[J]]);
	  }


		      if (d < difference) {
			DELETE(f,h);
                        e->refer++;
                        if(f->backtrack != NULL)
                          if( -- f->backtrack->refer == 0)
                            free_edge(f->backtrack);
			f->dist = d + e->dist;
                        f->backtrack = e;
			INSERT(f,h);
                      }
		    }
		  nextedge: continue;
                  }
                }
		if (e->refer==0)
                  free_edge(e);
	}
	if (oflag)
          fprintf(stderr,"\n");
	return e;
}

/* row in region on face of lattice */
typedef struct {
	int	*column;	/* column index of vertex */
	int	width;		/* number of vertices on row in region */
} ROW;

/*first index of face is first sequence, second index of face is second
  sequence. Third index of face is position in first sequence.
  Value of face item is set of positions in second sequence that are
  consistent with the position in the first sequence based on the
  pairwise alignment of the two sequences. "Consistent" means that
  this pair of positions was of use in the dynamic programming graph
  for the pairwise alignment*/
ROW	face [NUMBER+1] [NUMBER+1] [LENGTH+1];	/* faces of lattice */
int	costs [NUMBER+1] [NUMBER+1];		/* pairwise costs */

COORDINATE	*A [NUMBER+2];  /* array to access vertices by lattice coordinate */

/* source : create source vertex of lattice */

VERTEX *source()
{
	auto		int	p[NUMBER+1],i,index[LENGTH+1];
	register	COORDINATE	*a;

	for (i=1;i<=K;i++)
          p[i] = 0;
	for (i=N[1];i>=0;i--) {
	  index[i] = i;
	}
	for (i=2, a=A[1]=create_coordinate(index,N[1]+1,NULL);i<=K;i++)
		a=a->coord_vals->next_coord
		    = create_coordinate(index,intersect(p,i,index),a->coord_vals);
        A[1]->refer++;    /* Make sure coordinate does not get freed */
	return (VERTEX *) (a->coord_vals->next_coord = (COORDINATE *)create_vertex(a->coord_vals));
}

/* sink :   create sink vertex of lattice */

VERTEX *sink ()
{
	register	int	i;
	auto		int	p[NUMBER+1],index[LENGTH+1];
	register	COORDINATE	*a;
	register	COORDINATE_VALUES	*f;

	for (i=1;i<=K;i++)
          p[i] = N[i];
	for (i=2,a=A[1];i<=K;i++) {
		f = a->coord_vals + p[i-1] - a->lo;
		a = f->next_coord = create_coordinate(index,intersect(p,i,index),f);
	}
	f = a->coord_vals + p[K] - a->lo;
	return (VERTEX *) (f->next_coord = (COORDINATE *)create_vertex(f));
}

/* intersect :  intersect regions on rows of faces */

/* Finds the possible values for the "seqnum"th coordinate of point,
   given point[1],point[2],...,point[seqnum-1].  The corresponding
   bound values are copied to bound[]. */
int intersect (point, seqnum, possible_values)
register	int	possible_values[];
int	point[], seqnum;
{
	register	int	i, j, k, *c;
	auto		int	J, m, n;
	register	ROW	*r;

       /*retrieve values that are consistent with the pairwise alignment of
         sequenecs 1 and seqnum*/
	r = face[1][seqnum] + point[1];
	c = r->column;
	m = r->width;
	for (i=0;i<m;i++) {
		possible_values[i] = c[i];
	}
	for (J=2;J<seqnum;J++) {
           /*Get values that are consistent with the pairwise
             alignement of sequences J and seqnum*/
		r = face[J][seqnum] + point[J];
		c = r->column;
		n = r->width;
             /*compute intersection of possible values array and
               c array for sequence j, keeping result in possible
               values*/
		for (i=j=k=0; i<m && j<n;)
			if (possible_values[i]<c[j])
                          i++;
			else if (possible_values[i] > c[j])
                          j++;
			else {              /* possible_values[i] == c[j] */
			  possible_values[k++] = possible_values[i++];
			  j++;
			}
		if ((m=k) <= 0) break;
	}
	return m;
}

/* init_adjacent :  initialize adjacent vertex generator */

/*void init_adjacent (e, p, q)
EDGE			*e;
register	int	p[], q[];
{
	register    int     i;
	register    COORDINATE_VALUES   *a;

	coord(e->tail,p);
	safe_coord(e->head,q);
	for (i=K, a=e->head->prev_coord_val; i>=1; i--, a=a->curr_coord->prev_coord_val) {
		A[i] = a->curr_coord;
		B[i] = (P[i] = q[i]) + 1;
	}
}*/

/* adjacent :   generate adjacent vertices within region */

void adjacent (e, q)
EDGE			*e;
register	int	q[];
{
	register	int	I, i;
	auto		int     index [LENGTH+1];
	auto		int	n;
	register	COORDINATE	*a;
	register	COORDINATE_VALUES	*f;
        int     P [NUMBER+1];   /* coordinates of generated vertex */
        int     B [NUMBER+1];   /* bound on coordinates of generated vertex */


	for (i=K, f=e->head->prev_coord_val; i>=1; i--, f=f->curr_coord->prev_coord_val) {
		A[i] = f->curr_coord;
		B[i] = (P[i] = q[i]) + 1;
	}

	for (I=K;I<=K;) {
		for (;I>=1;I--)
		  if (P[I]<B[I])
		    break;
		if (I<1)
		  return;
     /*first test is to see whether i is in the proper range,
       second test is to see whether i is a valid coordinate*/
		if ( (i = P[I] = B[I]) >= (a = A[I])->lo && i <= a->hi
		    && (f = a->coord_vals + i - a->lo)->curr_coord != NULL ) {
     /*if child in trie does not exist*/
		  if (f->next_coord==NULL)
		    if ( I < K )
			if ( (n = intersect(P,I+1,index)) > 0 )
				f->next_coord = create_coordinate(index,n,f);
			else {
			  f->curr_coord=NULL;
			  if (free_coordinate(a))
			    I--;
			  continue;
			}
		    else
		      f->next_coord = (COORDINATE *)create_vertex(f);
		  for ( a = A[++I] = f->next_coord; I <= K; a = A[++I] = f->next_coord )
		    if ((i=P[I]=B[I]-1) >= a->lo && i<=a->hi &&
			(f=a->coord_vals+i-a->lo)->curr_coord != NULL) {
                      if (f->next_coord==NULL)
			if ( I < K )
			  if ( (n = intersect(P,I+1,index)) > 0 )
				f->next_coord = create_coordinate(index,n,f);
			  else {
			    f->curr_coord=NULL;
			    if (free_coordinate(a))
			      I--;
			    break;
			  }
			else f->next_coord = (COORDINATE *)create_vertex(f);
		    }
		    else break;
		}
            if (I == (K+1)) {
              create_edge(e->head,a);
              I = K;
            }
	}
}

/* safe_coord :  compute lattice coordinate of vertex */

void safe_coord (v, p)
VERTEX			*v;
register	int	p[];
{
	register	int	i;
	register	COORDINATE_VALUES	*a;
	register	COORDINATE	*b;

          for (i=K, a=v->prev_coord_val; i>=1; i--, a=b->prev_coord_val) {
		b = a->curr_coord;
		p[i] = a - b->coord_vals + b->lo;
                /*p[i] = a->value;*/
	  }
}

/* coord :  compute lattice coordinate of vertex */

void coord (v, p)
VERTEX			*v;
register	int	p[];
{
	register	int	i;
	register	COORDINATE_VALUES	*a;
	register	COORDINATE	*b;

	if (v!=presource)
          for (i=K, a=v->prev_coord_val; i>=1; i--, a=b->prev_coord_val) {
		b = a->curr_coord;
		p[i] = a - b->coord_vals + b->lo;
                /*p[i] = a->value;*/
	  }
	else
          for (i=K;i>=1;i--)
            p[i] = -1;
}

/* cost : compute cost of edge <q,r> preceded by edge <p,q> */
/*int cost (p, q, r)
int	p[], q[], r[];
{
	register	char	C [NUMBER+1];
	register	int	w, I, J, t [2] [NUMBER+1];

	for (I=1;I<=K;I++) {
		C[I] = (r[I]>q[I] ? S[I][r[I]] : DASH);
		t[0][I] = q[I] - p[I];
		t[1][I] = r[I] - q[I];
	}
	if(gflag)
          for (I=1;I<=K;I++)
            if (q[I]==0||q[I]==N[I])
              t[0][I]=t[1][I]=2;
	for (I=2,w=0;I<=K;I++)
          for (J=1;J<I;J++)
	    w+=scale[J][I]*(D[C[J]][C[I]]+T[t[0][I]][t[0][J]][t[1][I]][t[1][J]]);
	return w;
}*/

/* project :	compute projected cost of edge <q,r> preceded by <p,q> */

void project (p, q, r)
int	    p[], q[], r[];
{
	char    C [NUMBER+1];
	register    int     I, J;
	int     t [2] [NUMBER+1];

	for (I=1;I<=K;I++) {
		C[I] = (r[I]>q[I] ? S[I][r[I]] : DASH);
		t[0][I] = q[I] - p[I];
		t[1][I] = r[I] - q[I];
	}
	if(gflag)
          for (I=1;I<=K;I++)
            if (q[I]==0||q[I]==N[I])
              t[0][I]=t[1][I]=2;
	for (I=1;I<K;I++)
          for (J=I+1;J<=K;J++)
	    proj[I][J]+=D[C[I]][C[J]]+T[t[0][I]][t[0][J]][t[1][I]][t[1][J]];
}

/* display :    display multiple sequence alignment */

void display (e)
register	EDGE	*e;
{
	static  	int     p [NUMBER+1], q [NUMBER+1], r [NUMBER+1];
	auto    	int 	d, I, J;
	register	EDGE	*f;

	/* shortest sink to source path in lattice within bound? */

	if (e==NULL || (d=e->dist)>Upper)
		fatal("Multiple alignment within bound does not exist.");

	/* recover shortest path to source tracing backward from sink */

	for (e->next=NULL;e->tail!=presource;e=f) {
                f=e->backtrack;
		coord(f->tail,p);
		safe_coord(f->head,q);
		safe_coord(e->head,r);
		f->next = e;
		project(p,q,r);
	}

	/* display alignment corresponding to path */

	printf("\n                  ***  Optimal Multiple Alignment  ***\n\n");
	for (e=e->next; e!=NULL; e=e->next) {
		coord(e->tail,p);
		safe_coord(e->head,q);
		column(p,q);
	}
	// this output() is the align WE WANT
	if (outputAlignment()==1){
	    printf("well done \n");
	    }

	/* output statistics */
/* // REMOVE STADISTICS TO MAKE IT CLEAR
***********************
	if (gflag)
          printf( "End gaps not penalized.");
	printf( "\nCostfile:                   %s\n",cname);
	printf( "Alignment cost: %7d     Lower bound: %7d\n",d,Lower);
	printf( "Delta:          %7d     Max. Delta:  %7d\n\n",d-Lower,delta);
	printf("Sequences  Proj. Cost  Pair. Cost  Epsilon");
	if (bflag) {
		printf("  Max. Epsi.  Weight  Weight*Cost\n");
		for (I=1;I<K;I++)
                  for (J=I+1;J<=K;J++)
		    printf("%3d%4d%12d%12d%9d%10d%10d%12d\n",
			    I,J,proj[I][J],costs[I][J],proj[I][J]-costs[I][J],
			    epsi[I][J],scale[I][J],scale[I][J]*proj[I][J]);
	}
	else {
		printf("  Max. Epsi.\n");
		for (I=1;I<K;I++)
                  for (J=I+1;J<=K;J++)
		    printf("%3d%4d%12d%12d%9d%10d\n",I,J,proj[I][J],
			    costs[I][J],proj[I][J]-costs[I][J],epsi[I][J]);
	}
***********************
*/
}

char    M [NUMBER+1] [LINE];    /* multiple sequence alignment matrix */
int     C = 0;                  /* current column in alignment matrix */

/* column : compute multiple sequence alignment column for edge <p,q> */

void column (p, q)
register    int     p[], q[];
{
	register    int     k;

	for (k=1;k<=K;k++)
          M[k][C] = S[k][q[k] * (q[k]-p[k])];
	if (++C >= LINE)
          output();
}

/* output : output multiple sequence alignment rows */

void output()
{
	register    int     k, c;

	if (C==0)
          return;
	for (k=1;k<=K;k++) {
		for (c=0;c<C;c++)
                  printf("%c",M[k][c]);
		printf("\n");
	}
	C=0;
}
// outputfunction for the optimal alignment
int outputAlignment(){
    	register    int     k, c;

    	if (C==0){
    	    fatal("Not able to generate Optimistic alignment");
            return 0;
        }
    	for (k=1;k<=K;k++) {
    		for (c=0;c<C;c++)
                      printf("%c",M[k][c]);
    		printf("\n");
    	}
    	C=0;
    	return 1;

}
/* heap :   create a discrete heap */

HEAP *heap (max)
int     max;
{
	register    EDGE	**b;
	register    HEAP	*h;

	h = (HEAP *) alloc(sizeof( HEAP ) + max*sizeof( EDGE * ));
	h->min = h->bucket;
	h->max = h->bucket + max;
	for (b=h->min; b<=h->max; b++)
          *b=NULL;
	return h;
}

/* extract :    extract minimum distance edge from heap */

EDGE *extract ()
{
	register    EDGE    **b, *e;

	for (b=h->min; b<=h->max; b++)
          if (*b != NULL)
            break;
	if ((h->min=b) > h->max)
          return NULL;
	e = *b;
	DELETE(e,h);
      /*remove e from the list of non-extracted in-edges incident to
        the head of e*/
	if (e->nonextracted_prev != NULL)
          e->nonextracted_prev->nonextracted_next = e->nonextracted_next;
	else
          e->head->nonextracted_inedges = e->nonextracted_next;
	if (e->nonextracted_next != NULL)
          e->nonextracted_next->nonextracted_prev = e->nonextracted_prev;

        e->nonextracted_prev = NULL;
        e->nonextracted_next = NULL;
	return e;
}

VERTEX  *avail_vertex = NULL;    /* available vertex list */

/* create_vertex : create a vertex */

VERTEX *create_vertex (prev_coord_val)
COORDINATE_VALUES   *prev_coord_val;
{
	register    VERTEX  *v;

	if (avail_vertex!=NULL) {
		v = avail_vertex;
		avail_vertex = (VERTEX *) v->out;
	}
	else
          v= (VERTEX *) alloc(sizeof( VERTEX ));
	v->prev_coord_val = prev_coord_val;
	v->out = NULL;
        v->nonextracted_inedges = NULL;
	return v;
}

/* free_vertex :    free a vertex */

void free_vertex (v)
register    VERTEX  *v;
{
	register    COORDINATE   *a;
        register    EDGE         *e, *e2;

	a = v->prev_coord_val->curr_coord;
	v->prev_coord_val->curr_coord = NULL;
	free_coordinate(a);
    /*remove all the remaining edges coming into v from the heap*/
        e = v->nonextracted_inedges;
        while(e != NULL) {
          DELETE(e,h);
          e2 = e->nonextracted_next;
          e->refer = -1;   /* kludge telling 
			      free_edge() not to call free_vertex() */
          free_edge(e);
          e = e2;
        }
	v->out = (EDGE *) avail_vertex;
	avail_vertex = v;
}

EDGE	*avail_edge = NULL;  	/* available edge list */

/* edge :  create edge v -> w */

EDGE *create_edge (v, w)
VERTEX  *v, *w;
{
	register	EDGE    *e;

	if (avail_edge!=NULL) {
		e = avail_edge;
		avail_edge = e->next;
	}
	else
          e = (EDGE *) alloc(sizeof( EDGE ));
    /*set endpoints of e*/
	e->tail = v;
	e->head = w;
    /*insert e into v's list of outgoing edges*/
	e->prev = NULL;
	e->next = v->out;
        e->backtrack = NULL;
        e->refer = 0;
	v->out = e;
	if ( e->next != NULL )
          e->next->prev = e;
	e->heap_succ = e->heap_pred = e->nonextracted_next
                     = e->nonextracted_prev = e;
	e->dist = Upper+1;   /*set e's distance to upper bound
			       for lack of more information */
      /*insert e into w's list of nonextracted incoming edges*/
	if (e->head->nonextracted_inedges != NULL)
          e->head->nonextracted_inedges->nonextracted_prev = e;
	e->nonextracted_next = e->head->nonextracted_inedges;
        e->nonextracted_prev = NULL;
        e->head->nonextracted_inedges = e;
	return e;
}

/* free_edge :	free an edge */

void free_edge (e)
register	EDGE    *e;
{
    if (e->heap_succ != NULL && e->heap_pred != NULL) {
        DELETE(e, h);
    }

    if (e->nonextracted_prev != NULL && e->nonextracted_next != NULL) {
	if (e->nonextracted_prev != NULL)
          e->nonextracted_prev->nonextracted_next = e->nonextracted_next;
	else
          e->head->nonextracted_inedges = e->nonextracted_next;
	if (e->nonextracted_next != NULL)
          e->nonextracted_next->nonextracted_prev = e->nonextracted_prev;
    }

	if (e->prev!=NULL)
          e->prev->next = e->next; 
	else
          e->tail->out = e->next;
	if (e->next!=NULL)
          e->next->prev = e->prev;
        if ((e->backtrack) && (-- e->backtrack->refer == 0))
          free_edge(e->backtrack);
	if ((e->head->out == NULL) && (e->refer != -1))
          free_vertex(e->head);
	e->next = avail_edge;
	avail_edge = e;
}

COORDINATE   *avail_coordinate = NULL;     /* available coordinate list */

/* coordinate :  create an coordinate */

COORDINATE *create_coordinate (index, n, prev_coord_val)
register    int     index[], n;
COORDINATE_VALUES       *prev_coord_val;
{
	auto        int     b;
	register    COORDINATE_VALUES   *p;
	register    COORDINATE   *a;
        register    short     counter;

      /*look for available COORDINATE on free list */
	if ( avail_coordinate != NULL ) {
		a = avail_coordinate;
		avail_coordinate = a->next_on_free_list;
	}
	else
          a = (COORDINATE *) alloc(sizeof( COORDINATE ));
	a->lo = index[0];
	a->hi = index[n-1];
	a->prev_coord_val = prev_coord_val;
	a->refer = n;
	a->coord_vals = (COORDINATE_VALUES *) alloc((a->hi - a->lo + 1) * sizeof( COORDINATE_VALUES ));
	for (p= a->coord_vals + a->hi - a->lo, counter = a->hi;
	     p>=a->coord_vals;p--,counter--){
          p->next_coord=p->curr_coord=NULL;
	  /*p->value = counter;*/
        }
  /*link each new COORDINATE_VALUE to its parent COORDINATE and adjust
    bounds; note that this loop works only on the possible values which
    are stored in the array index*/
	for (n--;n>=0;n--) {
		(p = a->coord_vals + index[n] - a->lo)->curr_coord = a;
	}
	return a;
}

/* free_coordinate :	free an coordinate */

int free_coordinate (a)
register	COORDINATE	*a;
{
	register	COORDINATE	*b;

	if (-- a->refer <= 0) {
                b = a->prev_coord_val->curr_coord;
		a->prev_coord_val->curr_coord = NULL;
		free((char *)a->coord_vals);
		a->next_on_free_list =  avail_coordinate;
		avail_coordinate =  a;
		if (b != NULL)
		  free_coordinate(b);
		return TRUE;
	}
	else return FALSE;
}

/* vector : create a vector of integers */

int *vector (a, n)
register	int	a[], n;
{
	register	int	*v, i;

	v = (int *) alloc(sizeof( int ) * n);
	for (i=0;i<n;i++)
          v[i] = a[i];
	return v;
}

/* min3 :   three argument minimum */

int min3 (x,y,z)
register	int	x,y,z;
{
	if (x<y)
		if (x<z)
                  return x;
		else
                  return z;
	else
		if (y<z)
                  return y;
		else
                  return z;
}

/* must_open :  force open file */

FILE *must_open (name, mode)
char	*name, *mode;
{
	FILE	*stream;

	if ((stream=fopen(name,mode))==NULL)
          fatal("Cannot open `%s'.",name);
	return stream;
}

/* alloc :  allocate memory */

char *alloc (bytes)
register	int	bytes;
{
	register	char	*memory;
/*	auto		char	*malloc();*/
 char	*malloc();

	if ((memory=malloc((unsigned)bytes))== NULL)
		fatal("Cannot allocate memory.");
	return memory;
}

/* fatal :  fatal error */

void fatal (format, message1, message2)
char	*format, *message1, *message2;
{
	fprintf(stderr,format,message1,message2);
	fprintf(stderr,"\n");
	gettimeofday(&endtime, NULL);  /* Get end time */
	deltatime = (endtime.tv_sec - starttime.tv_sec)
	     + (endtime.tv_usec - starttime.tv_usec)/1000000.0;
	printf("Elapsed time = %7.3f\n",deltatime);
	exit(1);
}
