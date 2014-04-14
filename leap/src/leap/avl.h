/*
 *  'application' avl stuff (BPLUS-style interface)
 */

#define DEFAULTKEYLEN (4 * sizeof(int))	/* size of default key */

typedef void *RECPOS;

typedef struct {
	RECPOS recptr;
	unsigned int count;
	char key[DEFAULTKEYLEN]; /* actually can be of any length */
} rectype;

typedef rectype	IX_REC;

typedef struct {
	void *root;
	int keylength; /* zero for null-terminated strings */
	int dup_keys;
		/*
		 *  0 -- repeated key causes an error message;
		 *  1 -- repeated key & rec cause an error message;
		 *  2 -- complete repetitions allowed, use repetition count.
		 */
} IX_DESC;

#define IX_OK   1
#define IX_FAIL 0
#define EOIX (-2)

extern int	create_index(IX_DESC *pix, int dup, int keylength);
extern int	destroy_index(IX_DESC *pix);
extern int	find_key(IX_REC *pe, IX_DESC *pix);
extern int	locate_key(IX_REC *pe, IX_DESC *pix); 
extern int	add_key(IX_REC *pe, IX_DESC *pix);
extern int	delete_key(IX_REC *pe, IX_DESC *pix);
extern int	first_key(IX_DESC *pix);
extern int	last_key(IX_DESC *pix); 
extern int	next_key(IX_REC *pe, IX_DESC *pix);
extern int	prev_key(IX_REC *pe, IX_DESC *pix);
extern int	find_exact(IX_REC *pe, IX_DESC *pix);
