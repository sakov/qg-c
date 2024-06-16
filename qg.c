/******************************************************************************
 *
 * File:        qg.c
 *
 * Created:     15/04/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: main() for the QG model
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qg.h"
#include "utils.h"
#include "qgprm.h"
#include "model.h"
#include "calcs.h"
#include "qgstep.h"

int verbose = VERBOSE_DEF;
int force = 0;

#define EPS 1.0e-8

/**
 */
static void usage()
{
    printf("  Usage: %s { - | <prm file> } [{-f | --force}]\n", PROGRAM_NAME);
    printf("         %s --describe-prm\n", PROGRAM_NAME);
    printf("         %s {-v | --version}\n", PROGRAM_NAME);
    printf("  Options:\n");
    printf("         -f, --force - clobber the output file\n");
    exit(0);
}

/**
 */
static void parse_commandline(int argc, char* argv[], char** fname)
{
    int usedefaults = 0;
    int i;

    if (argc < 2)
        usage();

    i = 1;
    while (i < argc) {
        if (argv[i][0] != '-') {
            if (*fname == NULL) {
                *fname = argv[i];
                i++;
                continue;
            } else
                usage();
        } else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--force") == 0) {
            force = 1;
            i++;
        } else if (strcmp(argv[i], "--describe-prm") == 0) {
            qgprm_describe();
            exit(0);
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0) {
            printf("  %s version %s\n", PACKAGE_NAME, QG_VERSION);
            exit(0);
        } else if (strcmp(argv[i], "-") == 0) {
            usedefaults = 1;
            i++;
            continue;
        } else
            quit("command line: option \"%s\" not recognised", argv[i]);
    }

    if (*fname == NULL && !usedefaults)
        quit("command line: parameter file not specified");
}

/**
 */
int main(int argc, char* argv[])
{
    char* prmfname = NULL;
    qgprm* prm = NULL;
    model* qg;
    int nstep, step, nr, dnout, dnoutave;

    parse_commandline(argc, argv, &prmfname);
    if (verbose)
        printf("  parameter file = %s\n", (prmfname == NULL) ? "-" : prmfname);

    prm = qgprm_read(prmfname);
    if (verbose)
        qgprm_print(prm);

    qg = model_create(prm);

    if (prm->infname == NULL || strcmp(prm->infname, prm->outfname) != 0)
        model_createoutput(qg);

    if (prm->infname == NULL) {
        qg->t = 0;
        if (verbose)
            printf("  starting model from zero state\n");
    } else {
        if (strcmp(prm->infname, prm->outfname) == 0 && prm->rstart != -1)
            quit("%s: can't specify \"rstart\" when starting from and writing to the same file", prm->infname);
        model_readinput(qg);
        if (verbose)
            printf("  starting the model\n");
    }

    if (prm->tend < qg->t) {
        prm->tend = qg->t + prm->tend;
        if (verbose)
            printf("    readjusted the end time: tend = %f\n", prm->tend);
    }

    laplacian(prm->lx / (double) (prm->n - 1), prm->m, prm->n, qg->psi, qg->q);
    {
        double* psi0 = qg->psi[0];
        double* q0 = qg->q[0];
        double* psiguess0 = qg->psiguess[0];
        int i;

        for (i = 0; i < qg->mn; ++i)
            q0[i] = q0[i] - prm->f * psi0[i];
        for (i = 0; i < qg->mn; ++i)
            psiguess0[i] = psi0[i];
    }

    nstep = (prm->tend - qg->t) / prm->dt;
    dnout = (int) (prm->dtout / prm->dt);
    dnoutave = (isfinite(prm->dtoutave)) ? (int) (prm->dtoutave / prm->dt) : 0;

    if (verbose) {
        printf("    nstep = %d\n", nstep);
        printf("    nrecord = %d\n", (int) ((prm->tend - qg->t) / prm->dtout + EPS));
    }
    if (verbose) {
        printtime("  ");
        printf("  main cycle:");
    }
    for (step = 0, nr = 0; step < nstep; ++step) {
        double* q0 = qg->q[0];
        double weight;
        int i;

        for (i = 0; i < qg->mn; ++i)
            if (fabs(q0[i]) > QMAX)
                quit("|qmax| > %.3g at t = %.6f, step = %d\n", QMAX, qg->t, step);
        if (prm->scheme == SCHEME_ORDER1)
            qg_step_order1(qg);
        else if (prm->scheme == SCHEME_ORDER2)
            qg_step_order2(qg);
        else if (prm->scheme == SCHEME_RK4)
            qg_step_rk4(qg);
        else if (prm->scheme == SCHEME_DP5)
            qg_step_dp5(qg);

        if (step > 0 && step % dnout == 0) {
            calc_psi(qg, qg->psiguess, qg->q, qg->psi);
            model_writedump(qg, 0);
            nr++;
        }

        if (!isfinite(prm->dtoutave))
            continue;

        weight = (step % dnoutave == 0) ? 0.5 : 1.0;
        for (i = 0; i < qg->mn; ++i)
            qg->psiave[i] += qg->psi[0][i] * weight;
        if (prm->save_q)
            for (i = 0; i < qg->mn; ++i)
                qg->qave[i] += qg->q[0][i] * weight;
        qg->tave += qg->t * weight;

        if (step > 0 && step % dnoutave == 0) {
            for (i = 0; i < qg->mn; ++i)
                qg->psiave[i] /= (double) dnoutave;
            if (prm->save_q)
                for (i = 0; i < qg->mn; ++i)
                    qg->qave[i] /= (double) dnoutave;
            qg->tave /= (double) dnoutave;
            model_writedump(qg, 1);
            for (i = 0; i < qg->mn; ++i)
                qg->psiave[i] = qg->psi[0][i] * 0.5;
            if (prm->save_q)
                for (i = 0; i < qg->mn; ++i)
                    qg->qave[i] = qg->q[0][i] * 0.5;
            qg->tave = qg->t * 0.5;
        }
    }
    if (verbose && nr > 0) {
        printf("\n");
        printtime("  ");
    }

    model_destroy(qg);

    return 0;
}
