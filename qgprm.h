/******************************************************************************
 *
 * File:        qgprm.h        
 *
 * Created:     15/04/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: - Header for qgprm.c.
 *              - Defines default settings for the model.
 *
 * Revisions:   
 *
 *****************************************************************************/

#if !defined(_QGPRM_H)

#define MREFIN_DEF 7
#define NX1_DEF 3
#define NY1_DEF 3
#define LX_DEF 1.0
#define RKB_DEF 3.0e-4
#define RKH_DEF -5.0e-8
#define RKH2_DEF 1.0e-12
#define F_DEF 1600.0
#define R_DEF 1.0e-5
#define A_DEF 6.29
#define K_DEF 1.0
#define SCHEME_DEF SCHEME_RK4
#define SAVE_Q_DEF 0
#define DT_DEF 1.0
#define TEND_DEF 50000.0
#define DTOUT_DEF 20.0
#define OUTFNAME_DEF "qg.nc"
#define OUTFNAMEAVE_DEF "qgave.nc"

typedef struct {
    char* prmfname;
    int mrefin;
    int nx1, ny1;
    double lx;
    double dt;
    double tend;
    double dtout;
    double dtoutave;
    double rkb;
    double rkh;
    double rkh2;
    double f;
    double r;
    double a;
    double k;
    int scheme;
    int rstart;
    char* infname;
    char* outfname;
    char* outfnameave;
    int save_q;

    int m, n;
} qgprm;

qgprm* qgprm_read(char* fname);
void qgprm_destroy(qgprm* prm);
void qgprm_print(qgprm* prm);
void qgprm_describe(void);

#define _QGPRM_H
#endif
