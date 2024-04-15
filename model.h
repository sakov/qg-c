/******************************************************************************
 *
 * File:        model.h        
 *
 * Created:     15/04/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: header for model.c
 *
 * Revisions:   
 *
 *****************************************************************************/

#if !defined(_MODEL_H)

typedef struct {
    char* prmfname;
    int mrefin;
    int nx1, ny1;

    int m, n, mn;
    double lx;
    double dt;
    double tend;
    double dtout;
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
    int save_q;

    double t;
    double** psi;
    double** q;
    double** psiguess;
    double** J;
    double** zeta;
    double** zeta2;
    double** zeta4;
    double** qflux1;
    double** qflux2;
    double** qflux3;
    double** qflux4;
    double** qflux5;
    double** q2;
    double** q3;
    double** q4;

    double* curlt;
} model;

model* model_create(qgprm* prm);
void model_destroy(model* m);
void model_readinput(model* m);
void model_createoutput(model* m);
void model_writedump(model* m);

#define _MODEL_H
#endif
