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
    printf("         %s --version\n", PROGRAM_NAME);
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
        } else if (strcmp(argv[i], "--version") == 0) {
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
    model* m;
    int nstep, step, nr;

    parse_commandline(argc, argv, &prmfname);
    if (verbose)
        printf("  parameter file = %s\n", (prmfname == NULL) ? "-" : prmfname);

    prm = qgprm_read(prmfname);
    if (verbose)
        qgprm_print(prm);

    m = model_create(prm);
    qgprm_destroy(prm);

    if (m->infname == NULL || strcmp(m->infname, m->outfname) != 0)
        model_createoutput(m);

    if (m->infname == NULL) {
        m->t = 0;
        if (verbose)
            printf("  starting model from zero state\n");
    } else {
        if (strcmp(m->infname, m->outfname) == 0 && m->rstart != -1)
            quit("%s: can't specify \"rstart\" when starting from and writing to the same file", m->infname);
        model_readinput(m);
        if (verbose)
            printf("  starting the model\n");
    }

    if (m->tend < m->t) {
        m->tend = m->t + m->tend;
        if (verbose)
            printf("    readjusted the end time: tend = %f\n", m->tend);
    }

    laplacian(m->lx / (double) (m->n - 1), m->m, m->n, m->psi, m->q);
    {
        double* psi0 = m->psi[0];
        double* q0 = m->q[0];
        double* psiguess0 = m->psiguess[0];
        int i;

        for (i = 0; i < m->mn; ++i)
            q0[i] = q0[i] - m->f * psi0[i];
        for (i = 0; i < m->mn; ++i)
            psiguess0[i] = psi0[i];
    }

    nstep = (m->tend - m->t) / m->dt;
    if (verbose) {
        printf("    nstep = %d\n", nstep);
        printf("    nrecord = %d\n", (int) ((m->tend - m->t) / m->dtout + EPS));
    }
    if (verbose) {
        printtime("  ");
        printf("  main cycle:");
    }
    for (step = 0, nr = 0; step < nstep; ++step) {
        double* q0 = m->q[0];
        int i;

        for (i = 0; i < m->mn; ++i)
            if (fabs(q0[i]) > QMAX)
                quit("|qmax| > %.3g at t = %.6f, step = %d\n", QMAX, m->t, step);
        if (m->scheme == SCHEME_ORDER1)
            qg_step_order1(m);
        else if (m->scheme == SCHEME_ORDER2)
            qg_step_order2(m);
        else if (m->scheme == SCHEME_RK4)
            qg_step_rk4(m);
        else if (m->scheme == SCHEME_DP5)
            qg_step_dp5(m);

        if (floor(m->t / m->dtout + EPS) != floor((m->t - m->dt) / m->dtout + EPS)) {
            calc_psi(m, m->psiguess, m->q, m->psi);
            model_writedump(m);
            nr++;
        }
    }
    if (verbose && nr > 0) {
        printf("\n");
        printtime("  ");
    }

    model_destroy(m);

    return 0;
}
