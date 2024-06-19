/******************************************************************************
 *
 * File:        model.h        
 *
 * Created:     15/04/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: handles structure "model"
 *
 * Revisions:   
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ncw.h"
#include "qg.h"
#include "qgprm.h"
#include "utils.h"
#include "model.h"

#define EPS_DOUBLE 1.0e-8
#define EPS_DOUBLE2 1.0e-15

/**
 */
model* model_create(qgprm* prm)
{
    model* qg = malloc(sizeof(model));
    int m = prm->m;
    int n = prm->n;
    int scheme = prm->scheme;
    int i;

    qg->prm = prm;
    qg->mn = m * n;

    qg->t = NAN;
    qg->psi = alloc2d(m, n, sizeof(double));
    qg->q = alloc2d(m, n, sizeof(double));
    qg->psiguess = alloc2d(m, n, sizeof(double));
    qg->J = alloc2d(m, n, sizeof(double));
    qg->zeta = alloc2d(m, n, sizeof(double));
    qg->zeta2 = alloc2d(m, n, sizeof(double));
    qg->zeta4 = alloc2d(m, n, sizeof(double));
    qg->qflux1 = alloc2d(m, n, sizeof(double));
    qg->qflux2 = (scheme >= SCHEME_RK4) ? alloc2d(m, n, sizeof(double)) : NULL;
    qg->qflux3 = (scheme >= SCHEME_RK4) ? alloc2d(m, n, sizeof(double)) : NULL;
    qg->qflux4 = (scheme >= SCHEME_RK4) ? alloc2d(m, n, sizeof(double)) : NULL;
    qg->qflux5 = (scheme == SCHEME_DP5) ? alloc2d(m, n, sizeof(double)) : NULL;
    qg->q2 = (scheme >= SCHEME_ORDER2) ? alloc2d(m, n, sizeof(double)) : NULL;
    qg->q3 = (scheme == SCHEME_ORDER2 || scheme == SCHEME_RK4) ? alloc2d(m, n, sizeof(double)) : NULL;
    qg->q4 = (scheme == SCHEME_RK4) ? alloc2d(m, n, sizeof(double)) : NULL;
    if (isfinite(prm->dtoutave)) {
        qg->psiave = calloc(qg->mn, sizeof(double));
        if (prm->save_q)
            qg->qave = calloc(qg->mn, sizeof(double));
    } else {
        qg->psiave = NULL;
        qg->qave = NULL;
    }

    qg->curlt = malloc(n * sizeof(double));
    for (i = 0; i < n; ++i) {
        double tmp = sin(prm->k * 2.0 * M_PI * (double) i / (double) n);

        qg->curlt[i] = -prm->a * tmp * fabs(tmp);
    }

    return qg;
}

/**
 */
void model_destroy(model* qg)
{
    free(qg->psi);
    free(qg->q);
    free(qg->psiguess);
    free(qg->J);
    free(qg->zeta);
    free(qg->zeta2);
    free(qg->zeta4);
    free(qg->qflux1);
    if (qg->qflux2 != NULL)
        free(qg->qflux2);
    if (qg->qflux3 != NULL)
        free(qg->qflux3);
    if (qg->qflux4 != NULL)
        free(qg->qflux4);
    if (qg->qflux5 != NULL)
        free(qg->qflux5);
    free(qg->q2);
    if (qg->q3 != NULL)
        free(qg->q3);
    if (qg->q4 != NULL)
        free(qg->q4);

    if (qg->psiave != NULL)
        free(qg->psiave);
    if (qg->qave != NULL)
        free(qg->qave);
    free(qg->curlt);
    qgprm_destroy(qg->prm);
    free(qg);
}

/**
 */
void model_readinput(model* qg)
{
    qgprm* prm = qg->prm;
    char* fname = prm->infname;
    int ncid, varid_psi, varid_t;
    size_t dimlen[3];

    if (verbose)
        printf("  reading restart from \"%s\":\n", fname);

    ncw_open(fname, NC_NOWRITE, &ncid);

    if (prm->rstart < 0)
        prm->rstart = ncw_inq_nrecords(ncid) - 1;
    printf("    start record = %d\n", prm->rstart);

    /*
     * check that saved and current model settings are the same
     */
    {
        double tmp;

        ncw_check_dimlen(ncid, "j", prm->m);
        ncw_check_dimlen(ncid, "i", prm->n);
        ncw_get_att_double(ncid, NC_GLOBAL, "lx", &tmp);
        if (fabs(tmp - prm->lx) > EPS_DOUBLE)
            quit("%s: lx = %f, %s: lx = %f\n", fname, tmp, prm->prmfname, prm->lx);
        ncw_get_att_double(ncid, NC_GLOBAL, "rkb", &tmp);
        if (fabs(tmp - prm->rkb) > EPS_DOUBLE)
            quit("%s: rkb = %f, %s: rkb = %f\n", fname, tmp, prm->prmfname, prm->rkb);
        ncw_get_att_double(ncid, NC_GLOBAL, "rkh", &tmp);
        if (fabs(tmp - prm->rkh) > EPS_DOUBLE)
            quit("%s: rkh = %f, %s: rkh = %f\n", fname, tmp, prm->prmfname, prm->rkh);
        ncw_get_att_double(ncid, NC_GLOBAL, "rkh2", &tmp);
        if (fabs(tmp - prm->rkh2) > EPS_DOUBLE2)
            quit("%s: rkh2 = %f, %s: rkh2 = %f\n", fname, tmp, prm->prmfname, prm->rkh2);
        ncw_get_att_double(ncid, NC_GLOBAL, "F", &tmp);
        if (fabs(tmp - prm->f) > EPS_DOUBLE)
            quit("%s: F = %f, %s: F = %f\n", fname, tmp, prm->prmfname, prm->f);
        ncw_get_att_double(ncid, NC_GLOBAL, "r", &tmp);
        if (fabs(tmp - prm->r) > EPS_DOUBLE)
            quit("%s: r = %f, %s: r = %f\n", fname, tmp, prm->prmfname, prm->r);
    }

    ncw_inq_varid(ncid, "psi", &varid_psi);
    ncw_check_varndims(ncid, varid_psi, 3);
    ncw_inq_vardims(ncid, varid_psi, 3, NULL, dimlen);
    if (dimlen[0] <= prm->rstart)
        quit("%s: can not open dump %d: # nr = %d\n", fname, prm->rstart, dimlen[0]);
    if (dimlen[1] != prm->m)
        quit("%s: psi: dimlen[1] = %zu != model->m = %d\n", fname, dimlen[1], prm->m);
    if (dimlen[2] != prm->n)
        quit("%s: psi: dimlen[2] = %zu != model->n = %d\n", fname, dimlen[2], prm->n);
    ncw_get_var_double_record(ncid, varid_psi, prm->rstart, qg->psi[0]);

    ncw_inq_varid(ncid, "t", &varid_t);
    ncw_check_varndims(ncid, varid_t, 1);
    ncw_get_var_double_record(ncid, varid_t, prm->rstart, &qg->t);
    if (verbose)
        printf("    start time = %.1f\n", qg->t);

    ncw_close(ncid);
}

/**
 */
static void _model_createoutput(model* qg, int ave)
{
    qgprm* prm = qg->prm;
    char* fname = (!ave) ? prm->outfname : prm->outfnameave;
    int ncid;
    int dimids[3];

    if (file_exists(fname) & !force)
        quit("\"%s\" exists; use \"-f\" to overwrite", fname);
    ncw_create(fname, NC_CLOBBER | NETCDF_FORMAT, &ncid);
    ncw_put_att_double(ncid, NC_GLOBAL, "dt", 1, &prm->dt);
    ncw_put_att_double(ncid, NC_GLOBAL, "lx", 1, &prm->lx);
    ncw_put_att_double(ncid, NC_GLOBAL, "rkb", 1, &prm->rkb);
    ncw_put_att_double(ncid, NC_GLOBAL, "rkh", 1, &prm->rkh);
    ncw_put_att_double(ncid, NC_GLOBAL, "rkh2", 1, &prm->rkh2);
    ncw_put_att_double(ncid, NC_GLOBAL, "F", 1, &prm->f);
    ncw_put_att_double(ncid, NC_GLOBAL, "r", 1, &prm->r);
    ncw_put_att_double(ncid, NC_GLOBAL, "A", 1, &prm->a);
    ncw_put_att_double(ncid, NC_GLOBAL, "k", 1, &prm->k);
    ncw_put_att_text(ncid, NC_GLOBAL, "scheme", SCHEME_STR[prm->scheme]);
    if (prm->rstart >= 0) {
        ncw_put_att_text(ncid, NC_GLOBAL, "restart", prm->infname);
        ncw_put_att_int(ncid, NC_GLOBAL, "rstart", 1, &prm->rstart);
    } else
        ncw_put_att_text(ncid, NC_GLOBAL, "comment", "spun up from zero");
    ncw_def_dim(ncid, "record", NC_UNLIMITED, &dimids[0]);
    ncw_def_dim(ncid, "j", prm->m, &dimids[1]);
    ncw_def_dim(ncid, "i", prm->n, &dimids[2]);
    ncw_def_var(ncid, "t", NC_DOUBLE, 1, &dimids[0], NULL);
    ncw_def_var(ncid, "psi", NC_FLOAT, 3, dimids, NULL);
    if (prm->save_q)
        ncw_def_var(ncid, "q", NC_FLOAT, 3, dimids, NULL);
    ncw_close(ncid);
    if (!ave)
        qg->t = 0.0;
    else
        qg->tave = 0.0;
}

/**
 */
static void model_createobsoutput(model* qg)
{
    qgprm* prm = qg->prm;
    char* fname = prm->obsfname;
    int ncid, varid, dimids[2];

    if (file_exists(fname) & !force)
        quit("\"%s\" exists; use \"-f\" to overwrite", fname);
    
    ncw_create(fname, NC_CLOBBER | NETCDF_FORMAT, &ncid);
    ncw_def_dim(ncid, "record", NC_UNLIMITED, &dimids[0]);
    ncw_def_dim(ncid, "nobs", prm->nobs, &dimids[1]);
    ncw_def_var(ncid, "t", NC_DOUBLE, 1, dimids, NULL);
    ncw_def_var(ncid, "ij", NC_INT, 2, dimids, NULL);
    ncw_def_var(ncid, "psi", NC_FLOAT, 2, dimids, &varid);
    ncw_put_att_double(ncid, varid, "estd", 1, &prm->estd);
    ncw_close(ncid);
}

/**
 */
void model_createoutput(model* qg)
{
    _model_createoutput(qg, 0);
    if (isfinite(qg->prm->dtoutave))
        _model_createoutput(qg, 1);
    model_createobsoutput(qg);
}

/**
 */
void model_writedump(model* qg, int ave)
{
    qgprm* prm = qg->prm;
    char* fname = (!ave) ? prm->outfname : prm->outfnameave;
    int ncid, varid;
    int nr;

    ncw_open(fname, NC_WRITE, &ncid);
    nr = ncw_inq_nrecords(ncid);
    ncw_inq_varid(ncid, "t", &varid);
    ncw_put_var_double_record(ncid, varid, nr, (!ave) ? &qg->t : &qg->tave);
    ncw_inq_varid(ncid, "psi", &varid);
    ncw_put_var_double_record(ncid, varid, nr, (!ave) ? qg->psi[0] : qg->psiave);
    if (prm->save_q) {
        ncw_inq_varid(ncid, "q", &varid);
        ncw_put_var_double_record(ncid, varid, nr, (!ave) ? qg->q[0] : qg->qave);
    }
    ncw_close(ncid);

    if (verbose) {
        printf((!ave) ? ":" : ".");
        fflush(stdout);
    }
}

/**
 */
void model_writeobs(model* qg)
{
    qgprm* prm = qg->prm;
    int* pos = malloc(prm->nobs * sizeof(int));
    float* obs = malloc(prm->nobs * sizeof(float));
    float* error = NULL;
    int ncid, varid;
    int nr, i;

    if (prm->estd > 0.0)
        error = malloc(prm->nobs * sizeof(float));
    get_obspos(prm->nobs, qg->mn, (int) qg->t, pos, prm->estd, error);
    if (error == NULL)
        for (i = 0; i < prm->nobs; ++i)
            obs[i] = (float) qg->psi[0][pos[i]];
    else
        for (i = 0; i < prm->nobs; ++i)
            obs[i] = (float) qg->psi[0][pos[i]] + error[i];
    
    ncw_open(prm->obsfname, NC_WRITE, &ncid);
    nr = ncw_inq_nrecords(ncid);
    ncw_inq_varid(ncid, "t", &varid);
    ncw_put_var_double_record(ncid, varid, nr, &qg->t);
    ncw_inq_varid(ncid, "ij", &varid);
    ncw_put_var_int_record(ncid, varid, nr, pos);
    ncw_inq_varid(ncid, "psi", &varid);
    ncw_put_var_float_record(ncid, varid, nr, obs);
    ncw_close(ncid);
    
    free(pos);
    free(obs);
}
