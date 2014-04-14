
/*
 *  'inner' avl stuff
 */
/* way3.h */

typedef int way3; /* -1, 0, 1 */

#define way3stop  ((way3)0)
#define way3left ((way3)-1)
#define way3right ((way3)1)

way3 makeway3();

#define way3sum(x,y) ((x)+(y))
/* assume x!=y */

#define way3opp(x) (-(x))

way3 way3opp2();

way3 way3random();

/* node.h */

typedef struct _node {
   struct _node *ptr[2]; /* left, right */
   way3 balance:2, trace:2;
   rectype data;
} node;

#define stepway(n,x) (((n)->ptr)[way3ix(x)])
#define stepopp(n,x) (((n)->ptr)[way3ix(way3opp(x))])

node *allocnode();
void freenode();
node *swapptr();
int way3ix();


/* tree.h */

#define SRF_FINDEQUAL 1
#define SRF_FINDLESS  2
#define SRF_FINDGREAT 4
#define SRF_SETMARK   8
#define SRF_FROMMARK 16

rectype *avltree_search();
rectype *avltree_insert();
rectype *avltree_delete();
void avltree_first();
void avltree_last();
long avltree_clear();

#define avltree_init(x) (*(x)=NULL)

int compkey();
void copydata();

/* 'PLUS' interface */

#define IX_OK   1
#define IX_FAIL 0
#define EOIX (-2)

int create_index();
int destroy_index();
int find_key();
int locate_key();
int add_key();
int delete_key();
int first_key();
int last_key();
int next_key();
int prev_key();
int find_exact();

