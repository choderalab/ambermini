# define COLORTEXT "YES" 
# define REDUCEFACTOR 1
# define PI 3.1415926
# define Bohr 0.529177249
# define MAXCHAR 2048
# define MAXATOM 2560
# define MAXRES  512 
# define MAXBOND 512
# define MAXRING 500 
# define MAXGAS  500 /*maximum gasiteger parameters*/
/*
For MAXRING, no dynamic memory applied since the actuall number is determined 
using a recursive function. However, for small and middle-sized molecules, 
it is unlikely that the ring num is larger than 1000
*/
# define MAXCYCLE 100
# define OUTPUTSTEP 10
# define MAXTWIST 10
# define ECSLONG 2
# define COSCUT 120
# define DEGRAD 3.1415926/180
# define VDWIDIST 10
# define ESIDIST 14
# define THETACUT 15
# define CUBE 2.0
# define MAXWILDATOM 20
# define MAXSCHAIN 100 
# define MAXCES 20
# define MAXBEED 20
# define MAXATOMTYPE 250
# ifdef NCSU_PENALTIES_H
#  define PSCUTOFF 10
# else
#  define PSCUTOFF 10
#  define MAXVASTATE 8192
# endif /* NCSU_PENALTIES_H */
# define MAX_CES_BOND 100
# define MAXAPUNIT 20 
# define MAXAPCONS 20

