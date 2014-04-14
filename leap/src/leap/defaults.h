#define MAXRESID 50000

typedef struct {
	double	dDSearchDistance;
	double	dDESPGridSpace;
	double	dDESPBoxSize;
	float	fDESPDielectric;
	int	iDESPConstant;
	int	pdbwritecharges;
	int nocenter;
	double	dGridSpace;
	double	dShellExtent;
	int	iDielectricFlag;
	int	iGBparm;
	int 	iOldPrmtopFormat;
	int		iGibbs;
	int 	iCharmm; 
	int	iResidueImpropers;
	int	iDeleteExtraPointAngles;
	int	iPdbLoadSequential;
	int	iHaveResIds;
	int iUseResIds;
	int iFlexibleWater;
	char	sResidueId[MAXRESID][7];
	double  dDipoleDampFactor;
	double  dSceeScaleFactor;
	double  dScnbScaleFactor;
	int     iCMAP;
	int     iIPOL;    
	int     iIPOLset;    /* indicate IPOL set in frcmod/parm.dat */
} defaultstruct ;

extern defaultstruct GDefaults; 

