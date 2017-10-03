/*************************************************************************
MSA subroutines.
Version 1.0   June 27, 1989

Routines for computing pairwise weights for SP alignment.
One first estimates an unrooted tree and then uses the topology
(Rationale-1) or the topology and branch lengths (Rationale-2)
to determine the weights.

**************************************************************************

Tree Program written by R Miner, with assistance by DJ Lipman

Program to make unrooted evolutionary tree from pairwise distance
matrix using the Neighbor Joining method.

See:	Saitou & Nei, Mol. Biol. Evol. 4 (1987) 406-425

Start with array (node_index) of NODES, where each NODE is an operational
taxonomic unit (OTU).  Array indicates starting assumption of tree
topology.  Find best neighbors and coalesce -> remove from array, and
replace with NODE which is ancestor of the neighbors.  Continue this
process until there are only 3 NODES in node_index (only one topology
possible with 3 OTU's, so problem is solved).  Kludge slightly to 
form into binary tree with appropriate root.  NOTE: Root has no biological
significance, but is for convenience in dealing with structure.
Nodes are labelled with 1) positive numbers corresponding to original
distance matrix 2) Negative one's  corresponding to internal nodes
3) Negative 9 corresponding to ancestral node.

Output is preorder traversal of tree (depth first).  
See associated program "read.c" for input routine of tree structure
*******************************************************************/

#include <stdio.h>
#include "defs.h"

#define ROOT -9
#define INOD -1

extern int	K;
extern int	scale[NUMBER+1][NUMBER+1];

typedef struct tn{
	int sqn;
	float ltop;
	float w;
	float W;
	float v;
	float V;
	struct tn *lt;
	struct tn *rt;
	struct tn *par;
	struct tn *bro;
	} NODE;

/****************************** Global Variables ************************/

NODE	**node_index;
int	*Dij[NUMBER];
int	vcount, indexlen, count, count2;
NODE	*vlist[2*NUMBER];
NODE	*makenode();
float	dist(),subdist(),br_len(),compute_S();
int	rdist();
void	whts(), coalesce(), minimize_Sij(), init_data(), rpt();

/***********************************************************************/

void bias()
{
	NODE *ancestor;
	float len;

/* read in dist. matrix, initialize starting data structure	*/

	init_data(); 	

/* determine optimum node pairs, and make them neighbors	*/

	minimize_Sij();

/* finish up turning initial data structure into binary tree	*/
/* by forming root node						*/

	ancestor = makenode(ROOT,0.0,node_index[0],node_index[1]);
	count2 = 1;
        len = dist(0,1);
	len -= subdist(node_index[0],0.0) / count2;
	count2 = 1;
	len -= subdist(node_index[1],0.0) / count2;
	node_index[0]->ltop = len;

/*	Print out rooted tree	*/

        rpt(ancestor);
	printf("\n");

/*	Calculate weights	*/

	whts();
}

/************************************************************************/
/* Routine:    makenode							*/
/* Arguments:  sqn  - sequence number or -1 for internal nodes		*/
/*	       dtop - distance to parent node				*/
/*	       lt,rt - pointers to child nodes				*/
/* purpose:    creates nodes for the tree				*/
/************************************************************************/
NODE *makenode(sqn, ltop, lt, rt)
	int sqn;
	float ltop;
	NODE *lt, *rt;
{
	NODE *ptr;

     	ptr= (NODE *) malloc(sizeof(NODE));
	vlist[vcount++] = ptr;
	ptr->ltop = ltop;
	ptr->lt = lt;
	ptr->rt = rt;
	ptr->sqn = sqn;
	if (sqn<0) {
		lt->bro = rt;
		rt->bro = lt;
		lt->par = rt->par = ptr;
	} 
	return(ptr);
}

/************************************************************************/
/* Function:  dist							*/
/* Arguments: i,j - indices to node_index				*/
/* Purpose:   computes the average distance between the nodes pointed   */
/* 	      to by node_index[i] and node_index[j]			*/
/************************************************************************/
float dist(i,j)
	int i,j;
{
	count=1;
	return((float)rdist(node_index[i], node_index[j])/count);
}

/************************************************************************/
/* Function:  rdist							*/
/* Arguments: A,B - pointers to nodes					*/
/* Purpose:   recursively traverse the tree to compute the sum of path  */
/* 	      lengths between *A and *B					*/
/************************************************************************/
int rdist(A, B)
	NODE *A, *B;
{
	if (A->sqn < 0) {
		++count;
        	return(rdist(A->lt, B) + rdist(A->rt, B));
	}
	else if (B->sqn < 0) {
		++count;
		return(rdist(A, B->lt) + rdist(A, B->rt));
	}
	return(Dij[A->sqn][B->sqn]);
}

/************************************************************************/
/* Function:  subdist							*/
/* Arguments: A - a pointer to a node, total - cumulative path length   */
/* Purpose:   recursively traverses the tree to compute the sum of path */
/*	      lengths from a given node out to the leaves		*/
/************************************************************************/
float subdist(A,total)
	NODE *A;
	float total;
{
	if (A->sqn >= 0) return(total + A->ltop);
	count2++;
	return(subdist(A->lt,A->ltop+total) + subdist(A->rt,A->ltop+total));
}

/************************************************************************/
/* Function:  br_len							*/
/* Arguments: i,j - indices to node_index				*/
/* Purpose:   computes the least squares estimate for the length of the */
/*	      edge between *node_indes[i] and *node_index[j]		*/
/************************************************************************/
float br_len(i,j)
	int i,j;
{
	float diz=0,djz=0;
        int t;

        count2 = 1;
	for (t=0; t<indexlen; ++t) if (t!=i && t!=j) {
		diz += dist(i,t);
		djz += dist(j,t);
	}
	diz = diz / (indexlen - 2);
	djz = djz / (indexlen - 2);
	return((dist(i,j) + diz - djz)/2 - subdist(node_index[i],0.0)/count2);
}

/************************************************************************/
/* Function:  compute_S							*/
/* Arguments: i,j - indices to node_index                               */
/* Purpose:   Computes the sum of path lengths resulting from making    */
/*	      nodes *node_index[i] and *node_index[j] neighbors         */
/************************************************************************/
float compute_S(i,j)
	int i, j;
{
	int t, tt;
	float s1=0, s2=0;

	for (t=0;t<indexlen;t++) if (t!=i && t!=j) s1 += dist(i,t)+dist(j,t);
	s1 = s1 / (2 * (indexlen- 2));
	for (t=0; t<indexlen; t++) for (tt=t+1; tt<indexlen; tt++)
	     if (t!=i && t!=j && tt!=i && tt!=j) s2 += dist(t,tt);
        s2 = s2 / (indexlen- 2);
	return(s1 + s2 + dist(i,j) / 2);
}

/************************************************************************/
/* Function:  coalesce							*/
/* Arguments: i,j - indices to node_index				*/
/* Purpose:   creates an internal node with *node_index[i] and 		*/
/*	      *node_indes[j] as children, and replaces node_index[i]    */
/*	      and node_index[j] with a single pointer to the new node   */
/************************************************************************/
void coalesce(i,j)
	int i, j;
{
	NODE *par;

	node_index[i]->ltop = br_len(i,j);
	node_index[j]->ltop = br_len(j,i);
	par = makenode(INOD, 0.0, node_index[i], node_index[j] );
	node_index[i] = par;
	node_index[j] = node_index[indexlen-1];
	--indexlen;
}

/************************************************************************/
/* Function:  minimize_Sij						*/
/* Arguments: none                                                      */
/* Purpose:   selects the optimum pair of nodes from amongst the nodes  */
/*	      currently pointed to by node_index to make neighbors,     */
/*	      joins them through and internal node, and repeats the     */
/*	      process until only two nodes remain			*/
/************************************************************************/
void minimize_Sij()
{
	int i, j, min_i=0, min_j=0;
	float tmp, min = BIG +1;

	for (i=0; i<indexlen; i++) for (j=i+1; j<indexlen; j++) {
		tmp = compute_S(i,j);
 	 	if (tmp < min) {
	     		min_i = i;
			min_j = j;
			min = tmp;
		}
        }
	coalesce(min_i,min_j);
	if (indexlen > 2) minimize_Sij();
}

/************************************************************************/
/* Function:  init_data							*/
/* Arguments: none                                                      */
/* Purpose:   reads in a pairwise distance matrix, creates leaf nodes,  */
/* 	      and initializes node_index to point at them.		*/
/************************************************************************/
void init_data()
{
	int i,j,t;

	indexlen = K;

	for (t=0;t<K;t++) Dij[t] = (int *) calloc(K, sizeof(int));
	node_index = (NODE **) malloc(K * sizeof(NODE *));
	for (i=0;i<K-1;++i) for (j=i+1;j<K;++j) Dij[j][i]=Dij[i][j]=scale[j+1][i+1];
	for (i=vcount=0;i<K;i++) node_index[i]=makenode(i,0.0,NULL,NULL);
}

/************************************************************************/
/* Function:  rpt							*/
/* Arguments: A - pointer to the root node                              */
/* Purpose:   traverses the tree in preorder form, writes node out to   */
/*	      a file and optionally the screen				*/
/************************************************************************/
void rpt(A)
	NODE *A;
{
	if (A->sqn >= 0) printf("Leaf #%d        Distance to parent = %7.2f\n",1+A->sqn,A->ltop);
	else {
		if (A->sqn == ROOT) printf("---------------- Tree given from ancestor ----------------\n");
		else printf("Internal Node  Distance to parent = %7.2f\n",A->ltop);
		printf("On the left:   ");
		rpt(A->lt);
		printf("On the right:  ");
		rpt(A->rt);
	}
}

/*************************************************************************

Program written by SF Altschul
Version 1.0   June 27, 1989

Program to calculate pair weights given an evolutionary tree.

See:	Altschul, "Leaf Pairs and Tree Dissections",
		SIAM J. Discrete Math. 2 (1989).
	Altschul, Carrol & Lipman, "Weights for Data Related by a Tree",
		J. Molec. Biol. 208 (1989). 

*************************************************************************/

NODE	**pN;
static  float	B[NUMBER][NUMBER];
void	trace();

#ifdef BIAS2			/* Rationale-2 weights */

void whts()
{
	int i,j;
	NODE *no;
	float sm;

/* Calculate the weights of all trees hanging from all internal nodes */

	for (pN=vlist; (no= *pN)->sqn > INOD; ++pN) {
		no->w = 1.0;
		no->W = no->ltop;
	}
	for (; (no= *pN)->sqn > ROOT; ++pN) {
		no->w = no->lt->w * no->rt->W + no->lt->W * no->rt->w;
		no->W = no->ltop  * no->w     + no->lt->W * no->rt->W;
	} 
	no->V = 1;
	no->v = 0;
	do {
		no= *(--pN);
		no->v = no->par->v * no->bro->W + no->par->V * no->bro->w;
		no->V = no->ltop   * no->v      + no->par->V * no->bro->W;
	}
	while (pN != vlist);

/* Calculate weights for leaf pairs using precomputed subtree weights */

	for (; (no= *pN)->sqn >INOD; ++pN) trace(1.0,no->ltop,no->par,no->bro);

/* Scale pair weights so that smallest is about 8	*/

	sm=1.0E+30;
	for (j=1;j<K;++j) for (i=0;i<j;++i) if (B[i][j]<sm) sm=B[i][j];
	sm /= 7.9;
	for (i=0;i<K-1;++i) for (j=i+1;j<K;++j) scale[i+1][j+1]=B[i][j]/sm+0.5;
}

/* trace is a recursive function that traverses tree to find all leaf pairs
whose first member is a given leaf  (Rationale-2)			*/

void trace(prod,sum,no,sis)
	float prod,sum;
	NODE *no, *sis;
{
	if (no->sqn > INOD) B[(*pN)->sqn][no->sqn] = sum*prod;
	else if (sis==NULL) {
		trace(prod * no->lt->W, sum + no->rt->ltop, no->rt, NULL);	
		trace(prod * no->rt->W, sum + no->lt->ltop, no->lt, NULL);	
	}
	else {
		trace(prod * no->V, sum + sis->ltop, sis, NULL);
		if (no->sqn != ROOT)
			trace(prod * sis->W, sum + no->ltop, no->par, no->bro);
	}
}

#else			/* Rationale-1 weights */

void whts()
{
	int i,j;
	NODE *no;
	float sm;

	for (pN=vlist; (no= *pN)->sqn >INOD; ++pN) trace(1.0, no->par, no->bro);
	for (sm=BIG,j=1;j<K;++j) for (i=0;i<j;++i) if (B[i][j]<sm) sm=B[i][j];
	for (i=0;i<K-1;++i) for (j=i+1;j<K;++j) scale[i+1][j+1]=B[i][j]/sm+0.5;
}

/* trace is a recursive function that traverses tree to find all leaf pairs
whose first member is a given leaf  (Rationale-1)			*/

void trace(prod,no,sis)
	float prod;
	NODE *no, *sis;
{
	if (no->sqn > INOD) B[(*pN)->sqn][no->sqn] = prod;
	else if (sis==NULL) {
		trace(prod/2, no->rt, NULL);
		trace(prod/2, no->lt, NULL);
	}
	else if (no->sqn != ROOT) {
		trace(prod/2, sis, NULL);
		trace(prod/2, no->par, no->bro);
	}
	else trace(prod, sis, NULL);
}

#endif








