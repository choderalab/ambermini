/*
 *   MxNameLen applies to block, category, item and keyword names...
 *   
 */
#define MxNameLen    		 80 
/*
 *   This is the size of the buffer length allocated to individual
 *   values.
 */

#define CifDefaultSpace           2

#undef   YYLMAX
#define  YYLMAX       1024
#if 0
#undef   YYMAXDEPTH  
#define  YYMAXDEPTH  20000
#undef   YYINITDEPTH 
#define  YYINITDEPTH  1000
#define  YYPRINT(file, type, value) cifpprint(file, type, value)
#endif

#define MAXVALUELENGTH   8182

#define TRUE                   1
#define FALSE                  0

typedef struct NdbCifRowFormat {
  char **columns;
} NdbCifRowFormat;

typedef struct NdbCifCategoryFormat {
  int numCol;
  int allCol;
  int curCol;
  int allRow;
  int numRow;
  int curRow;
  char categoryName[MxNameLen];
  char **colNames;
  NdbCifRowFormat *rows;
} NdbCifCategoryFormat;

typedef struct NdbCifDatablockFormat {
  int numCategory; /* Number of categories in this datablock */
  int allCategory; /* Allocated space */
  int curCategory; /* index of the current category */
  char datablockName[MxNameLen]; 
  NdbCifCategoryFormat *categories;
} NdbCifDatablockFormat;  

typedef struct NdbCifDatablocksFormat {
  int numDatablock; /* Number of datablocks in this structure */
  int curDatablock; /* index of the current datablock */
  int allDatablock; /* Allocated space */
  NdbCifDatablockFormat *datablocks;
} NdbCifDatablocksFormat;


int cifpparse();
char *ndb_cif_copy_item_value();
char ndb_cif_set_null_char();
int ndb_cif_init();
int ndb_cif_read_file(FILE *fp);
int get_category_index(int blockIndex, char *categoryName);
int get_column_index(int blockIndex, int categoryIndex, char *columnName);
void ndb_cif_process_item_name_list();
void ndb_cif_process_value_list();
void ndb_cif_process_item_name_value_pair();
void ndb_cif_process_loop_declaration();
int ndb_cif_current_category();
int ndb_cif_get_category_name_from_item_name(char*, char*);
int ndb_cif_get_item_keyword_from_item_name(char*, char*);
int ndb_cif_put_item_keyword(char*);
int ndb_cif_new_row();
int ndb_cif_put_item_value(int, char*);
int ndb_cif_new_category(char*);
int ndb_cif_current_col();
int ndb_cif_rewind_column();
int ndb_cif_get_datablock_id(char*);
int ndb_cif_new_datablock(char*);
int ndb_cif_move_datablock(char*);
int ndb_cif_reset_datablocks();
int ndb_cif_reset_datablock();
int ndb_cif_rewind_datablock();
int ndb_cif_current_datablock_name(char*);
int ndb_cif_current_datablock();
int ndb_cif_count_datablock();
int ndb_cif_count_column();
int ndb_cif_count_row();
int ndb_cif_get_item_name(int, char*);
int ndb_cif_get_item_value(int, char*, int);
int ndb_cif_next_row();
int ndb_cif_next_category();
int ndb_cif_next_datablock();
int ndb_cif_rewind_row();
int ndb_cif_rewind_category();
int ndb_cif_remove_category();
int ndb_cif_reset_category_by_id(int, int);
int ndb_cif_remove_row_by_id(int, int, int);
int ndb_cif_get_category_id(char*, char*);


#ifdef CIF_GLOBAL
	FILE *cifpin;
	char TempKeyword[MxNameLen+1], TempValue[MAXVALUELENGTH+1];
	NdbCifDatablocksFormat cifFiles;
	int  lineNo;      
#else
	extern char TempKeyword[MxNameLen+1], TempValue[MAXVALUELENGTH+1];
	extern FILE *cifpin;
	extern int  lineNo;
	extern NdbCifDatablocksFormat cifFiles;
#endif
