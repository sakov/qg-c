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
    model* m = malloc(sizeof(model));
    int i;

    m->prmfname = strdup(prm->prmfname);
    m->mrefin = prm->mrefin;
    m->nx1 = prm->nx1;
    m->ny1 = prm->ny1;

    m->m = prm->m;
    m->n = prm->n;
    m->mn = m->m * m->n;
    m->lx = prm->lx;
    m->dt = prm->dt;
    m->tend = prm->tend;
    m->dtout = prm->dtout;
    m->rkb = prm->rkb;
    m->rkh = prm->rkh;
    m->rkh2 = prm->rkh2;
    m->f = prm->f;
    m->r = prm->r;
    m->scheme = prm->scheme;
    m->rstart = prm->rstart;
    m->infname = (prm->infname == NULL) ? NULL : strdup(prm->infname);
    m->outfname = strdup(prm->outfname);
    m->save_q = prm->save_q;
    m->a = prm->a;
    m->k = prm->k;

    m->t = NAN;
    m->psi = alloc2d(m->m, m->n, sizeof(double));
    m->q = alloc2d(m->m, m->n, sizeof(double));
    m->psiguess = alloc2d(m->m, m->n, sizeof(double));
    m->J = alloc2d(m->m, m->n, sizeof(double));
    m->zeta = alloc2d(m->m, m->n, sizeof(double));
    m->zeta2 = alloc2d(m->m, m->n, sizeof(double));
    m->zeta4 = alloc2d(m->m, m->n, sizeof(double));
    m->qflux1 = alloc2d(m->m, m->n, sizeof(double));
    m->qflux2 = (m->scheme >= SCHEME_RK4) ? alloc2d(m->m, m->n, sizeof(double)) : NULL;
    m->qflux3 = (m->scheme >= SCHEME_RK4) ? alloc2d(m->m, m->n, sizeof(double)) : NULL;
    m->qflux4 = (m->scheme >= SCHEME_RK4) ? alloc2d(m->m, m->n, sizeof(double)) : NULL;
    m->qflux5 = (m->scheme == SCHEME_DP5) ? alloc2d(m->m, m->n, sizeof(double)) : NULL;
    m->q2 = (m->scheme >= SCHEME_ORDER2) ? alloc2d(m->m, m->n, sizeof(double)) : NULL;
    m->q3 = (m->scheme == SCHEME_ORDER2 || m->scheme == SCHEME_RK4) ? alloc2d(m->m, m->n, sizeof(double)) : NULL;
    m->q4 = (m->scheme == SCHEME_RK4) ? alloc2d(m->m, m->n, sizeof(double)) : NULL;

    m->curlt = malloc(m->n * sizeof(double));
    for (i = 0; i < m->n; ++i) {
        double tmp = sin(m->k * 2.0 * M_PI * (double) i / (double) m->n);

        m->curlt[i] = -m->a * tmp * fabs(tmp);
    }

    return m;
}

/**
 */
void model_destroy(model* m)
{
    free(m->prmfname);
    free(m->infname);
    free(m->outfname);
    free(m->psi);
    free(m->q);
    free(m->psiguess);
    free(m->J);
    free(m->zeta);
    free(m->zeta2);
    free(m->zeta4);
    free(m->qflux1);
    if (m->qflux2 != NULL)
        free(m->qflux2);
    if (m->qflux3 != NULL)
        free(m->qflux3);
    if (m->qflux4 != NULL)
        free(m->qflux4);
    if (m->qflux5 != NULL)
        free(m->qflux5);
    free(m->q2);
    if (m->q3 != NULL)
        free(m->q3);
    if (m->q4 != NULL)
        free(m->q4);
    free(m->curlt);
    free(m);
}

/**
 */
void model_readinput(model* m)
{
    char* fname = m->infname;
    int ncid, varid_psi, varid_t;
    size_t dimlen[3];

    if (verbose)
        printf("  reading restart from \"%s\":\n", fname);

    ncw_open(fname, NC_NOWRITE, &ncid);

    if (m->rstart < 0)
        m->rstart = ncw_inq_nrecords(ncid) - 1;
    printf("    start record = %d\n", m->rstart);

    /*
     * check that saved and current model settings are the same
     */
    {
        double tmp;

        ncw_check_dimlen(ncid, "j", m->m);
        ncw_check_dimlen(ncid, "i", m->n);
        ncw_get_att_double(ncid, NC_GLOBAL, "lx", &tmp);
        if (fabs(tmp - m->lx) > EPS_DOUBLE)
            quit("%s: lx = %f, %s: lx = %f\n", fname, tmp, m->prmfname, m->lx);
        ncw_get_att_double(ncid, NC_GLOBAL, "rkb", &tmp);
        if (fabs(tmp - m->rkb) > EPS_DOUBLE)
            quit("%s: rkb = %f, %s: rkb = %f\n", fname, tmp, m->prmfname, m->rkb);
        ncw_get_att_double(ncid, NC_GLOBAL, "rkh", &tmp);
        if (fabs(tmp - m->rkh) > EPS_DOUBLE)
            quit("%s: rkh = %f, %s: rkh = %f\n", fname, tmp, m->prmfname, m->rkh);
        ncw_get_att_double(ncid, NC_GLOBAL, "rkh2", &tmp);
        if (fabs(tmp - m->rkh2) > EPS_DOUBLE2)
            quit("%s: rkh2 = %f, %s: rkh2 = %f\n", fname, tmp, m->prmfname, m->rkh2);
        ncw_get_att_double(ncid, NC_GLOBAL, "F", &tmp);
        if (fabs(tmp - m->f) > EPS_DOUBLE)
            quit("%s: F = %f, %s: F = %f\n", fname, tmp, m->prmfname, m->f);
        ncw_get_att_double(ncid, NC_GLOBAL, "r", &tmp);
        if (fabs(tmp - m->r) > EPS_DOUBLE)
            quit("%s: r = %f, %s: r = %f\n", fname, tmp, m->prmfname, m->r);
    }

    ncw_inq_varid(ncid, "psi", &varid_psi);
    ncw_check_varndims(ncid, varid_psi, 3);
    ncw_inq_vardims(ncid, varid_psi, 3, NULL, dimlen);
    if (dimlen[0] <= m->rstart)
        quit("%s: can not open dump %d: # nr = %d\n", fname, m->rstart, dimlen[0]);
    if (dimlen[1] != m->m)
        quit("%s: psi: dimlen[1] = %zu != model->m = %d\n", fname, dimlen[1], m->m);
    if (dimlen[2] != m->n)
        quit("%s: psi: dimlen[2] = %zu != model->n = %d\n", fname, dimlen[2], m->n);
    ncw_get_var_double_record(ncid, varid_psi, m->rstart, m->psi[0]);

    ncw_inq_varid(ncid, "t", &varid_t);
    ncw_check_varndims(ncid, varid_t, 1);
    ncw_get_var_double_record(ncid, varid_t, m->rstart, &m->t);
    if (verbose)
        printf("    start time = %.1f\n", m->t);

    ncw_close(ncid);
}

/**
 */
void model_createoutput(model* m)
{
    char* fname = m->outfname;
    int ncid;
    int dimids[3];

    if (file_exists(fname) & !force)
        quit("\"%s\" exists; use \"-f\" to overwrite", fname);
    ncw_create(fname, NC_CLOBBER | NETCDF_FORMAT, &ncid);
    ncw_put_att_double(ncid, NC_GLOBAL, "dt", 1, &m->dt);
    ncw_put_att_double(ncid, NC_GLOBAL, "lx", 1, &m->lx);
    ncw_put_att_double(ncid, NC_GLOBAL, "rkb", 1, &m->rkb);
    ncw_put_att_double(ncid, NC_GLOBAL, "rkh", 1, &m->rkh);
    ncw_put_att_double(ncid, NC_GLOBAL, "rkh2", 1, &m->rkh2);
    ncw_put_att_double(ncid, NC_GLOBAL, "F", 1, &m->f);
    ncw_put_att_double(ncid, NC_GLOBAL, "r", 1, &m->r);
    ncw_put_att_double(ncid, NC_GLOBAL, "A", 1, &m->a);
    ncw_put_att_double(ncid, NC_GLOBAL, "k", 1, &m->k);
    ncw_put_att_text(ncid, NC_GLOBAL, "scheme", SCHEME_STR[m->scheme]);
    if (m->rstart >= 0) {
        ncw_put_att_text(ncid, NC_GLOBAL, "restart", m->infname);
        ncw_put_att_int(ncid, NC_GLOBAL, "rstart", 1, &m->rstart);
    } else
        ncw_put_att_text(ncid, NC_GLOBAL, "comment", "spun up from zero");
    ncw_def_dim(ncid, "record", NC_UNLIMITED, &dimids[0]);
    ncw_def_dim(ncid, "j", m->m, &dimids[1]);
    ncw_def_dim(ncid, "i", m->n, &dimids[2]);
    ncw_def_var(ncid, "t", NC_DOUBLE, 1, &dimids[0], NULL);
    ncw_def_var(ncid, "psi", NC_FLOAT, 3, dimids, NULL);
    if (m->save_q)
        ncw_def_var(ncid, "q", NC_FLOAT, 3, dimids, NULL);
    ncw_close(ncid);

    m->t = 0;
}

/**
 */
void model_writedump(model* m)
{
    char* fname = m->outfname;
    int ncid, varid;
    int nr;

    ncw_open(fname, NC_WRITE, &ncid);
    nr = ncw_inq_nrecords(ncid);
    ncw_inq_varid(ncid, "t", &varid);
    ncw_put_var_double_record(ncid, varid, nr, &m->t);
    ncw_inq_varid(ncid, "psi", &varid);
    ncw_put_var_double_record(ncid, varid, nr, m->psi[0]);
    if (m->save_q) {
        ncw_inq_varid(ncid, "q", &varid);
        ncw_put_var_double_record(ncid, varid, nr, m->q[0]);
    }
    ncw_close(ncid);

    if (verbose) {
        printf(".");
        fflush(stdout);
    }
}
