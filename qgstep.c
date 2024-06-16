/******************************************************************************
 *
 * File:        qgstep.c
 *
 * Created:     15/04/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Implements quasi-geostrophic equation for potential vorticity.
 *              Contains 4 time-stepping schemes:
 *                1-order
 *                2-order
 *                rk4
 *                dp5
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "qg.h"
#include "qgprm.h"
#include "model.h"
#include "calcs.h"
#include "helmholtz.h"

/**
 */
void calc_psi(model* qg, double** psiin, double** q, double** psiout)
{
    qgprm* prm = qg->prm;
    int one = 1;
    int two = 2;
    int three = 3;
    double h1 = prm->lx / (double) prm->nx1;

    __helmholtz_mod_MOD_helmholtz(psiin[0], q[0], &prm->f, &prm->nx1, &prm->ny1, &prm->mrefin, &h1, &one, &two, &three, &prm->n, &prm->m, psiout[0]);
}

/**
 */
static void qgflux(model* qg, double** q, double** psiin, double** psiout, double** qflux)
{
    qgprm* prm = qg->prm;
    double h = prm->lx / (double) (prm->n - 1);
    int i, j;

    calc_psi(qg, psiin, q, psiout);
    arakawa(h, prm->m, prm->n, psiout, q, qg->J);
    laplacian(h, prm->m, prm->n, psiout, qg->zeta);
    laplacian(h, prm->m, prm->n, qg->zeta, qg->zeta2);
    laplacian(h, prm->m, prm->n, qg->zeta2, qg->zeta4);

    /*
     * (assume that qgflux was initialised to zero and therefore stays zero at
     *  the boundaries)
     */
    for (j = 1; j < prm->m - 1; ++j) {
        double* qfluxj = qflux[j];
        double* Jj = qg->J[j];
        double* zetaj = qg->zeta[j];
        double* zeta2j = qg->zeta2[j];
        double* zeta4j = qg->zeta4[j];
        double* psijp1 = psiout[j + 1];
        double* psijm1 = psiout[j - 1];

        for (i = 1; i < prm->n - 1; ++i)
            qfluxj[i] = -prm->r * Jj[i] - prm->rkb * zetaj[i] + prm->rkh * zeta2j[i] - prm->rkh2 * zeta4j[i] + qg->curlt[i] - 0.5 / h * (psijp1[i] - psijm1[i]);
    }
}

/**
 */
void qg_step_order1(model* qg)
{
    double* q0 = qg->q[0];
    double* qflux0 = qg->qflux1[0];
    double dt = qg->prm->dt;
    int i;

    qgflux(qg, qg->q, qg->psi, qg->psiguess, qg->qflux1);
    memcpy(qg->psi[0], qg->psiguess[0], qg->mn * sizeof(double));
    for (i = 0; i < qg->mn; ++i)
        q0[i] += dt * qflux0[i];

    qg->t += dt;
}

/**
 */
void qg_step_order2(model* qg)
{
    double* q0 = qg->q[0];
    double* q20 = qg->q2[0];
    double* q30 = qg->q3[0];
    double* qflux0 = qg->qflux1[0];
    double dt = qg->prm->dt;
    double dt2 = 0.5 * dt;
    int i;

    qgflux(qg, qg->q, qg->psiguess, qg->psi, qg->qflux1);
    for (i = 0; i < qg->mn; ++i)
        q20[i] = q0[i] + dt2 * qflux0[i];
    for (i = 0; i < qg->mn; ++i)
        q30[i] = q0[i] + dt * qflux0[i];

    qgflux(qg, qg->q2, qg->psi, qg->psiguess, qg->qflux1);
    for (i = 0; i < qg->mn; ++i)
        q20[i] += dt2 * qflux0[i];
    for (i = 0; i < qg->mn; ++i)
        q0[i] = 2.0 * q20[i] - q30[i];

    qg->t += dt;
}

/**
 */
void qg_step_rk4(model* qg)
{
    double* q0 = qg->q[0];
    double* q20 = qg->q2[0];
    double* q30 = qg->q3[0];
    double* q40 = qg->q4[0];
    double* qflux10 = qg->qflux1[0];
    double* qflux20 = qg->qflux2[0];
    double* qflux30 = qg->qflux3[0];
    double* qflux40 = qg->qflux4[0];
    double dt = qg->prm->dt;
    double dt2 = 0.5 * dt;
    double dt6 = dt / 6.0;
    int i;

    qgflux(qg, qg->q, qg->psiguess, qg->psi, qg->qflux1);
    for (i = 0; i < qg->mn; ++i)
        q20[i] = q0[i] + dt2 * qflux10[i];

    qgflux(qg, qg->q2, qg->psi, qg->psiguess, qg->qflux2);
    for (i = 0; i < qg->mn; ++i)
        q30[i] = q0[i] + dt2 * qflux20[i];

    qgflux(qg, qg->q3, qg->psiguess, qg->psi, qg->qflux3);
    for (i = 0; i < qg->mn; ++i)
        q40[i] = q0[i] + dt * qflux30[i];

    qgflux(qg, qg->q4, qg->psi, qg->psiguess, qg->qflux4);
    for (i = 0; i < qg->mn; ++i)
        q0[i] += (qflux10[i] + 2.0 * (qflux20[i] + qflux30[i]) + qflux40[i]) * dt6;

    qg->t += dt;
}

/**
 */
void qg_step_dp5(model* qg)
{
    double* q0 = qg->q[0];
    double* q20 = qg->q2[0];
    double* qflux10 = qg->qflux1[0];
    double* qflux20 = qg->qflux2[0];
    double* qflux30 = qg->qflux3[0];
    double* qflux40 = qg->qflux4[0];
    double* qflux50 = qg->qflux5[0];
    double dt = qg->prm->dt;
    int i;

    qgflux(qg, qg->q, qg->psiguess, qg->psi, qg->qflux1);
    for (i = 0; i < qg->mn; ++i)
        q20[i] = q0[i] + dt * 0.2 * qflux10[i];

    qgflux(qg, qg->q2, qg->psi, qg->psiguess, qg->qflux2);
    for (i = 0; i < qg->mn; ++i)
        q20[i] = q0[i] + dt * ((3.0 / 40.0) * qflux10[i] + (9.0 / 40.0) * qflux20[i]);

    qgflux(qg, qg->q2, qg->psiguess, qg->psi, qg->qflux3);
    for (i = 0; i < qg->mn; ++i)
        q20[i] = q0[i] + dt * ((44.0 / 45.0) * qflux10[i] - (56.0 / 15.0) * qflux20[i] + (32.0 / 9.0) * qflux30[i]);

    qgflux(qg, qg->q2, qg->psi, qg->psiguess, qg->qflux4);
    for (i = 0; i < qg->mn; ++i)
        q20[i] = q0[i] + dt * ((19372.0 / 6561.0) * qflux10[i] - (25360.0 / 2187.0) * qflux20[i] + (64448.0 / 6561.0) * qflux30[i] - (212.0 / 729.0) * qflux40[i]);

    qgflux(qg, qg->q2, qg->psiguess, qg->psi, qg->qflux5);
    for (i = 0; i < qg->mn; ++i)
        q20[i] = q0[i] + dt * ((9017.0 / 3168.0) * qflux10[i] - (355.0 / 33.0) * qflux20[i] + (46732.0 / 5247.0) * qflux30[i] + (49.0 / 176.0) * qflux40[i] - (5103.0 / 18656.0) * qflux50[i]);

    qgflux(qg, qg->q2, qg->psi, qg->psiguess, qg->qflux2);
    for (i = 0; i < qg->mn; ++i)
        q0[i] += dt * ((35.0 / 384.0) * qflux10[i] + (500.0 / 1113.0) * qflux30[i] + (125.0 / 192.0) * qflux40[i] - (2187.0 / 6784.0) * qflux50[i] + (11.0 / 84.0) * qflux20[i]);

    qg->t += dt;
}
