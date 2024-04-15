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
void calc_psi(model* m, double** psiin, double** q, double** psiout)
{
    int one = 1;
    int two = 2;
    int three = 3;
    double h1 = m->lx / (double) m->nx1;

    __helmholtz_mod_MOD_helmholtz(psiin[0], q[0], &m->f, &m->nx1, &m->ny1, &m->mrefin, &h1, &one, &two, &three, &m->n, &m->m, psiout[0]);
}

/**
 */
static void qgflux(model* m, double** q, double** psiin, double** psiout, double** qflux)
{
    double h = m->lx / (double) (m->n - 1);
    int i, j;

    calc_psi(m, psiin, q, psiout);
    arakawa(h, m->m, m->n, psiout, q, m->J);
    laplacian(h, m->m, m->n, psiout, m->zeta);
    laplacian(h, m->m, m->n, m->zeta, m->zeta2);
    laplacian(h, m->m, m->n, m->zeta2, m->zeta4);

    /*
     * (assume that qgflux was initialised to zero and therefore stays zero at
     *  the boundaries)
     */
    for (j = 1; j < m->m - 1; ++j) {
        double* qfluxj = qflux[j];
        double* Jj = m->J[j];
        double* zetaj = m->zeta[j];
        double* zeta2j = m->zeta2[j];
        double* zeta4j = m->zeta4[j];
        double* psijp1 = psiout[j + 1];
        double* psijm1 = psiout[j - 1];

        for (i = 1; i < m->n - 1; ++i)
            qfluxj[i] = -m->r * Jj[i] - m->rkb * zetaj[i] + m->rkh * zeta2j[i] - m->rkh2 * zeta4j[i] + m->curlt[i] - 0.5 / h * (psijp1[i] - psijm1[i]);
    }
}

/**
 */
void qg_step_order1(model* m)
{
    double* q0 = m->q[0];
    double* qflux0 = m->qflux1[0];
    int i;

    qgflux(m, m->q, m->psi, m->psiguess, m->qflux1);
    memcpy(m->psi[0], m->psiguess[0], m->m * m->n * sizeof(double));
    for (i = 0; i < m->mn; ++i)
        q0[i] += m->dt * qflux0[i];

    m->t += m->dt;
}

/**
 */
void qg_step_order2(model* m)
{
    double* q0 = m->q[0];
    double* q20 = m->q2[0];
    double* q30 = m->q3[0];
    double* qflux0 = m->qflux1[0];
    double dt2 = 0.5 * m->dt;
    int i;

    qgflux(m, m->q, m->psiguess, m->psi, m->qflux1);
    for (i = 0; i < m->mn; ++i)
        q20[i] = q0[i] + dt2 * qflux0[i];
    for (i = 0; i < m->mn; ++i)
        q30[i] = q0[i] + m->dt * qflux0[i];

    qgflux(m, m->q2, m->psi, m->psiguess, m->qflux1);
    for (i = 0; i < m->mn; ++i)
        q20[i] += dt2 * qflux0[i];
    for (i = 0; i < m->mn; ++i)
        q0[i] = 2.0 * q20[i] - q30[i];

    m->t += m->dt;
}

/**
 */
void qg_step_rk4(model* m)
{
    double* q0 = m->q[0];
    double* q20 = m->q2[0];
    double* q30 = m->q3[0];
    double* q40 = m->q4[0];
    double* qflux10 = m->qflux1[0];
    double* qflux20 = m->qflux2[0];
    double* qflux30 = m->qflux3[0];
    double* qflux40 = m->qflux4[0];
    double dt2 = 0.5 * m->dt;
    double dt6 = m->dt / 6.0;
    int i;

    qgflux(m, m->q, m->psiguess, m->psi, m->qflux1);
    for (i = 0; i < m->mn; ++i)
        q20[i] = q0[i] + dt2 * qflux10[i];

    qgflux(m, m->q2, m->psi, m->psiguess, m->qflux2);
    for (i = 0; i < m->mn; ++i)
        q30[i] = q0[i] + dt2 * qflux20[i];

    qgflux(m, m->q3, m->psiguess, m->psi, m->qflux3);
    for (i = 0; i < m->mn; ++i)
        q40[i] = q0[i] + m->dt * qflux30[i];

    qgflux(m, m->q4, m->psi, m->psiguess, m->qflux4);
    for (i = 0; i < m->mn; ++i)
        q0[i] += (qflux10[i] + 2.0 * (qflux20[i] + qflux30[i]) + qflux40[i]) * dt6;

    m->t += m->dt;
}

/**
 */
void qg_step_dp5(model* m)
{
    double* q0 = m->q[0];
    double* q20 = m->q2[0];
    double* qflux10 = m->qflux1[0];
    double* qflux20 = m->qflux2[0];
    double* qflux30 = m->qflux3[0];
    double* qflux40 = m->qflux4[0];
    double* qflux50 = m->qflux5[0];
    int i;

    qgflux(m, m->q, m->psiguess, m->psi, m->qflux1);
    for (i = 0; i < m->mn; ++i)
        q20[i] = q0[i] + m->dt * 0.2 * qflux10[i];

    qgflux(m, m->q2, m->psi, m->psiguess, m->qflux2);
    for (i = 0; i < m->mn; ++i)
        q20[i] = q0[i] + m->dt * ((3.0 / 40.0) * qflux10[i] + (9.0 / 40.0) * qflux20[i]);

    qgflux(m, m->q2, m->psiguess, m->psi, m->qflux3);
    for (i = 0; i < m->mn; ++i)
        q20[i] = q0[i] + m->dt * ((44.0 / 45.0) * qflux10[i] - (56.0 / 15.0) * qflux20[i] + (32.0 / 9.0) * qflux30[i]);

    qgflux(m, m->q2, m->psi, m->psiguess, m->qflux4);
    for (i = 0; i < m->mn; ++i)
        q20[i] = q0[i] + m->dt * ((19372.0 / 6561.0) * qflux10[i] - (25360.0 / 2187.0) * qflux20[i] + (64448.0 / 6561.0) * qflux30[i] - (212.0 / 729.0) * qflux40[i]);

    qgflux(m, m->q2, m->psiguess, m->psi, m->qflux5);
    for (i = 0; i < m->mn; ++i)
        q20[i] = q0[i] + m->dt * ((9017.0 / 3168.0) * qflux10[i] - (355.0 / 33.0) * qflux20[i] + (46732.0 / 5247.0) * qflux30[i] + (49.0 / 176.0) * qflux40[i] - (5103.0 / 18656.0) * qflux50[i]);

    qgflux(m, m->q2, m->psi, m->psiguess, m->qflux2);
    for (i = 0; i < m->mn; ++i)
        q0[i] += m->dt * ((35.0 / 384.0) * qflux10[i] + (500.0 / 1113.0) * qflux30[i] + (125.0 / 192.0) * qflux40[i] - (2187.0 / 6784.0) * qflux50[i] + (11.0 / 84.0) * qflux20[i]);

    m->t += m->dt;
}
