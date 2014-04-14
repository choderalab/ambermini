/*
 *  AVL routines by Gregory Tseytin of St. Petersburg, Russia,
 *	Adapted by him to use the same interface of the BPLUS
 *	B-tree shareware package. Error checking added by Bill
 *	Ross.
 */

#include "avl.h"

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>


/*
 *  'inner' avl stuff
 */
/* way3.h */

typedef char	way3; /* -1, 0, 1 */

#define way3stop  ((way3)0)
#define way3left ((way3)-1)
#define way3right ((way3)1)

#define way3sum(x,y) ((x)+(y))
/* assume x!=y */

#define way3opp(x) (-(x))


/* node.h */

typedef struct _node {
	struct _node 	*ptr[2]; /* left, right */
	way3 		balance, trace;
	rectype 	data;
} node;

#define stepway(n,x) (((n)->ptr)[way3ix(x)])
#define stepopp(n,x) (((n)->ptr)[way3ix(way3opp(x))])


/* tree.h */

#define SRF_FINDEQUAL 1
#define SRF_FINDLESS  2
#define SRF_FINDGREAT 4
#define SRF_SETMARK   8
#define SRF_FROMMARK 16

#define avltree_init(x) (*(x)=NULL)



static int ix_keylength,ix_dupkeys; /* set from IX_DESC */
static int rec_keylength;           /* set from actual key */
static int node_overhead=sizeof(node)-DEFAULTKEYLEN;

/******************************************************************************
								WAY3
 ******************************************************************************/
static way3
makeway3(int n)
{
	return n>0 ? way3right : n<0 ? way3left : way3stop;
}

static way3
way3opp2(way3 x, way3 y)
{
	return x==y ? way3opp(x) : way3stop;
}

#if 0
static way3
way3random()
{
	return rand()>rand() ? way3left : way3right;
}
#endif


/*****************************************************************************/

static void
freenode(node *n)
{
	free(n);
}


static int
compkey(rectype *r1, rectype *r2)
{
	int n= ix_keylength ?
       		memcmp(r1->key,r2->key,ix_keylength) :
       		strcmp(r1->key,r2->key);
	if (n || !ix_dupkeys) return n;
	return memcmp(&(r1->recptr),&(r2->recptr),sizeof(RECPOS));
}

static void
copydata(rectype *r1, rectype *r2)
{
	r1->recptr = r2->recptr;
	r1->count = r2->count;
	if (ix_keylength)
		memcpy(r1->key, r2->key, ix_keylength);
	else 
		strcpy(r1->key, r2->key);
}

static void
duprec(rectype *r)
{
	if (r->count++==UINT_MAX) {
		fprintf(stderr,"avltrees: repeat count exceeded\n");
		exit(1);
	}
}

static node *
allocnode()
{
	int size=(ix_keylength ? ix_keylength : rec_keylength);
	node *n=(node *)malloc(size+node_overhead);
	if (n==NULL) {
   		fprintf(stderr,"avltrees: out of memory\n");
   		exit(1);
   	}
	if (ix_dupkeys)
		n->data.count=1;
	return n;
}



/******************************************************************************
								NODE
 ******************************************************************************/
static node *
swapptr(node **ptrptr, node *new)
{
	node *old=*ptrptr;
	*ptrptr=new;
	return old;
}

static int
way3ix(way3 x) /* assume x!=0 */
{
	return x==way3right ? 1 : 0;
}

/******************************************************************************
								TREE
 ******************************************************************************/

typedef int bool;

static node **t;
static node *r,*s;
static way3 wayhand;

static bool
restruct(bool op_del)
{
	way3 n=r->balance,c;
	node *p;
	bool g= n==way3stop ? op_del : n==wayhand;
	if (g) p=r;
	else {
   		p=stepopp(r,wayhand);
   		stepopp(r,wayhand)=swapptr(&stepway(p,wayhand),r);
   		c=p->balance;
   		s->balance= way3opp2(c,wayhand);
   		r->balance= way3opp2(c,way3opp(wayhand));
   		p->balance= way3stop;
   	}
	stepway(s,wayhand)=swapptr(&stepopp(p,wayhand),s);
	*t=p;
#ifdef TESTING
	if (op_del)
   		{if (g) rstd1++; else rstd2++;}
else
   	{if (g) rsti1++; else rsti2++;}
#endif
	return g;
}

static rectype *
avltree_search(node **tt, rectype *key, unsigned short searchflags)
{
	way3 aa;
	node *p,*q,*pp;
	way3 waydir,wayopp;
	if (!(~searchflags & (SRF_FINDGREAT|SRF_FINDLESS))) return NULL;
	if (!(searchflags & (SRF_FINDGREAT|SRF_FINDEQUAL|SRF_FINDLESS)))
		return NULL;
	waydir=searchflags & SRF_FINDGREAT ? way3right :
       		searchflags & SRF_FINDLESS ? way3left : way3stop;
	wayopp=way3opp(waydir);
	p=q=NULL;
	while ((pp=*tt)!=NULL) {
   		aa= searchflags & SRF_FROMMARK ? pp->trace :
                                    makeway3(compkey(key,&(pp->data)));
   		if (searchflags & SRF_SETMARK) pp->trace=aa;
   		if (aa==way3stop) {
      			if (searchflags & SRF_FINDEQUAL) return &(pp->data);
      			if ((q=stepway(pp,waydir))==NULL) break;
      			if (searchflags & SRF_SETMARK) pp->trace=waydir;
      			for (;;) {
         			if ((pp=stepway(q,wayopp))==NULL) {
            				if (searchflags & SRF_SETMARK)
						q->trace=way3stop;
            				return &(q->data);
            			}
         			if (searchflags & SRF_SETMARK) q->trace=wayopp;
         			q=pp;
         		}
      		}
   		/* remember the point where we can change direction to waydir */
   		if (aa==wayopp) p=pp;
   		tt=&stepway(pp,aa);
   	}
	if (p==NULL || !(searchflags & (SRF_FINDLESS|SRF_FINDGREAT)))
		return NULL;
	if (searchflags & SRF_SETMARK) p->trace=way3stop;
	return &(p->data);
}

static void
avltree_first(node **tt)
{
	node *pp;
	while ((pp=*tt)!=NULL) {
   		pp->trace=way3left;
   		tt=&stepway(pp,way3left);
   	}
}

static void
avltree_last(node **tt)
{
	node *pp;
	while ((pp=*tt)!=NULL) {
   		pp->trace=way3right;
   		tt=&stepway(pp,way3right);
   	}
}

static rectype *
avltree_insert(node **tt, rectype *key)
{
	way3 aa,b;
	node *p,*q,*pp;

	t=tt;
	p=*tt;
	while ((pp=*tt)!=NULL) {
		aa= makeway3(compkey(key,&(pp->data)));
		if (aa==way3stop) {
			if (ix_dupkeys == 2)
				duprec(&(pp->data));
			return NULL;
		}
   		if (pp->balance!=way3stop)
			t=tt; /* t-> the last disbalanced node */
   		pp->trace=aa;
   		tt=&stepway(pp,aa);
   	}
	*tt=q=allocnode();
	q->balance=q->trace=way3stop;
	stepway(q,way3left)=stepway(q,way3right)=NULL;
	key->count = 1;
	copydata(&(q->data),key);
	/* balancing */
	s=*t; wayhand=s->trace;
	if (wayhand!=way3stop) {
   		r=stepway(s,wayhand);
   		for (p=r; p!=NULL; p=stepway(p,b))
      			b=p->balance=p->trace;
   		b=s->balance;
   		if (b!=wayhand) s->balance=way3sum(wayhand,b);
   		else if (restruct(0)) s->balance=r->balance=way3stop;
   	}
	return &(q->data);
}

static rectype *
avltree_delete(node **tt, rectype *key, unsigned short searchflags)
{
	way3 aa,aaa,b,bb;
	node *p,*q,*pp,*p1;
	node **t1,**tt1,**qq1,**rr=tt;
	t=t1=tt1=qq1=tt;
	p=*tt; q=NULL;
	aaa=way3stop;
	while ((pp=*tt)!=NULL) {
   		aa= aaa!=way3stop ? aaa :
            		searchflags & SRF_FROMMARK ? pp->trace :
            		makeway3(compkey(key,&(pp->data)));
   		b=pp->balance;
   		if (aa==way3stop) {
      			qq1=tt; q=pp; rr=t1;
      			aa= b!=way3stop ? b : way3left;
      			aaa=way3opp(aa); /* will move opposite to aa */
      		}
   		t=t1;
   		if (b==way3stop || (b!=aa && stepopp(pp,aa)->balance==way3stop))
			t1=tt;
   		tt1=tt;
   		tt=&stepway(pp,aa);
   		pp->trace=aa;
   	}
	if (aaa==way3stop) return NULL;
	copydata(key,&(q->data));
	p=*tt1;
	*tt1=p1=stepopp(p,p->trace);
	if (p!=q) {
   		*qq1=p; memcpy(p->ptr,q->ptr,sizeof(p->ptr));
   		p->balance=q->balance;
   		wayhand=p->trace=q->trace;
   		if (t==&stepway(q,wayhand)) t=&stepway(p,wayhand);
   	}
	while ((s=*t)!=p1) {
   		wayhand=way3opp(s->trace);
   		b=s->balance;
   		if (b!=wayhand) s->balance=way3sum(wayhand,b);
   		else {
      			r=stepway(s,wayhand);
      			if (restruct(1)) {
         			if ((bb=r->balance)!=way3stop)
					s->balance=way3stop;
         			r->balance=way3sum(way3opp(wayhand),bb);
         		}
      		}
   		t=&stepopp(s,wayhand);
   	}
	while ((p=*rr)!=NULL) {
		/* adjusting trace */
  		aa= makeway3(compkey(&(q->data),&(p->data)));
  		p->trace=aa; rr=&stepway(p,aa);
  	}
	freenode(q);
	return key;
}

static long
avltree_clear(node **tt)
{
	long nodecount=0L;
	node *p=*tt,*q=NULL,*x,**xx;
	if (p != NULL) {
		for (;;) {
   			if ((x=stepway(p,way3left))!=NULL ||
      					(x=stepway(p,way3right))!=NULL) {
      				stepway(p,way3left)=q;
      				q=p; p=x; continue;
      			}
   			freenode(p); nodecount++;
   			if (q==NULL) break;
   			if (*(xx=&stepway(q,way3right))==p) *xx=NULL;
   			p=q; q=*(xx=&stepway(p,way3left)); *xx=NULL;
   		}
		*tt=NULL;
	}
	return nodecount;
}


/******************************************************************************
						'PLUS' interface style
 ******************************************************************************/

int
create_index(IX_DESC *pix, int dup, int keylength)
{
	if (dup < 0  ||  dup > 2) {
		fprintf(stderr,
			"create_index 'dup'=%d: programming error\n", dup);
		exit(1);
	}
	if (keylength < 0) {
		fprintf(stderr,
			"create_index 'keylength'=%d: programming error\n",
			keylength);
		exit(1);
	}
	pix->root=NULL;
	pix->keylength=keylength;
	pix->dup_keys=dup;
	return IX_OK;
}

int
destroy_index(IX_DESC *pix)
{
	ix_keylength=pix->keylength;
	avltree_clear((node **)&(pix->root));
	pix->root=NULL;
	return IX_OK;
}

int
find_key(IX_REC *pe, IX_DESC *pix)
{
	rectype *ptr;

	ix_keylength=pix->keylength; ix_dupkeys=pix->dup_keys;
	memset((void *)&(pe->recptr), 0, sizeof(RECPOS));
	ptr=avltree_search((node **)&(pix->root),pe,
			SRF_FINDEQUAL|SRF_SETMARK|SRF_FINDGREAT);
	if (ptr == NULL) return IX_FAIL;
	pe->recptr=ptr->recptr;
	pe->count = ptr->count;
	return compkey(pe,ptr) ? IX_FAIL : IX_OK;
}

int
locate_key(IX_REC *pe, IX_DESC *pix)
{
	rectype *ptr; int ret;
	ix_keylength=pix->keylength; ix_dupkeys=pix->dup_keys;
	memset((void *)&(pe->recptr),0,sizeof(RECPOS));
	ptr=avltree_search((node **)&(pix->root),pe,
			SRF_FINDEQUAL|SRF_SETMARK|SRF_FINDGREAT);
	if (ptr==NULL) return EOIX;
	ret= compkey(pe,ptr) ? IX_FAIL : IX_OK;
	copydata(pe,ptr); return ret;
}

int
add_key(IX_REC *pe, IX_DESC *pix)
{
	ix_keylength=pix->keylength; ix_dupkeys=pix->dup_keys;
	if (ix_keylength == 0)
		rec_keylength = strlen(pe->key) + 1;
	if (avltree_insert((node **)&(pix->root),pe)==NULL && ix_dupkeys < 2)
		return IX_FAIL;
	return IX_OK;
}

int
delete_key(IX_REC *pe, IX_DESC *pix)
{
	rectype *ptr;

	ix_keylength=pix->keylength; ix_dupkeys=pix->dup_keys;
	ptr=avltree_search((node **)&(pix->root),pe,SRF_FINDEQUAL|SRF_SETMARK);
	if (ptr==NULL) 
		return IX_FAIL;
	if (ix_dupkeys==2 && --pe->count) 
		return IX_OK;
	avltree_delete((node **)&(pix->root),pe,SRF_FROMMARK);
	return IX_OK;
}

int
first_key(IX_DESC *pix)
{
	avltree_first((node **)&(pix->root)); return IX_OK;
}

int
last_key(IX_DESC *pix)
{
	avltree_last((node **)&(pix->root)); return IX_OK;
}

int
next_key(IX_REC *pe, IX_DESC *pix)
{
	rectype *ptr;
	ix_keylength=pix->keylength; ix_dupkeys=pix->dup_keys;
	if ((ptr=avltree_search((node **)&(pix->root),pe, /* pe not used */
                       SRF_FROMMARK|SRF_SETMARK|SRF_FINDGREAT))==NULL)
	   	return EOIX;
	copydata(pe,ptr); return IX_OK;
}

int
prev_key(IX_REC *pe, IX_DESC *pix)
{
	rectype *ptr;
	ix_keylength=pix->keylength; ix_dupkeys=pix->dup_keys;
	if ((ptr=avltree_search((node **)&(pix->root),pe, /* pe not used */
                       SRF_FROMMARK|SRF_SETMARK|SRF_FINDLESS))==NULL)
	   	return EOIX;
	copydata(pe,ptr); return IX_OK;
}

int
find_exact(IX_REC *pe, IX_DESC *pix)
{
	rectype *ptr;
	ix_keylength=pix->keylength; ix_dupkeys=pix->dup_keys;
	ptr=avltree_search((node **)&(pix->root),pe,
                  SRF_FINDEQUAL|SRF_FINDGREAT|SRF_SETMARK);
	if (ptr==NULL) return IX_FAIL;
	if (!ix_dupkeys && pe->recptr!=ptr->recptr) return IX_FAIL;
	return IX_OK;
}
