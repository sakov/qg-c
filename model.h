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
    qgprm* prm;
    int mn;

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
    double tave;
    double* psiave;
    double* qave;

    double* curlt;
} model;

model* model_create(qgprm* prm);
void model_destroy(model* m);
void model_readinput(model* m);
void model_createoutput(model* m);
void model_writedump(model* m, int ave);
void model_writeobs(model*qg);

#define _MODEL_H
#endif
