typedef struct CMNT_t {
    char *record;
    struct CMNT_t *next;
} CMNT;

typedef char WRD[8];

typedef struct MAPANM_t {
    WRD atmname[5];
    struct MAPANM_t *next;
} MAPANM;

typedef struct MAPRDX_t {
    int residx[5];
    struct MAPRDX_t *next;
} MAPRDX;

typedef struct {
    char title[256];
    WRD *reslist;
    WRD *creslist;
    WRD *nreslist;
    int nres;
    CMNT *cmnt;
    int resolution;
    double *map;
    int active;
    int mapidx;
    int termmap;  // applicable to terminal
    WRD atmname[5];
    WRD catmname[5];
    WRD natmname[5];
    int residx[5];
    int cresidx[5];
    int nresidx[5];
} CMAP;

typedef struct {
    WRD res;
    int atoms[5];
    int mapid;
} PHIPSI;

typedef struct CMAPLST_t {
    CMAP *cmap;
    struct CMAPLST_t *next;
} CMAPLST;

extern CMAP *cmap;
extern CMAPLST *cmaplst;
extern int mapnum;
