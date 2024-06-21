/******************************************************************************
 *
 * File:        qgprm.c        
 *
 * Created:     15/04/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: handles parameter file for the QG model
 *
 * Revisions:   
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include "qg.h"
#include "qgprm.h"
#include "utils.h"

const char* SCHEME_STR[] = {
    "order1",
    "order2",
    "rk4",
    "dp5"
};

static int ipow(int base, unsigned int exp)
{
    int tmp;

    if (exp == 0)
        return 1;
    tmp = ipow(base, exp / 2);
    if (exp % 2 == 0)
        return tmp * tmp;
    else
        return base * tmp * tmp;
}

static qgprm* qgprm_create(void)
{
    qgprm* prm = malloc(sizeof(qgprm));

    prm->prmfname = NULL;
    prm->mrefin = MREFIN_DEF;
    prm->nx1 = NX1_DEF;
    prm->ny1 = NY1_DEF;
    prm->lx = LX_DEF;
    prm->m = NX1_DEF * ipow(2, MREFIN_DEF - 1) + 1;
    prm->n = NY1_DEF * ipow(2, MREFIN_DEF - 1) + 1;
    prm->dt = DT_DEF;
    prm->tend = TEND_DEF;
    prm->dtout = DTOUT_DEF;
    prm->dtoutave = NAN;
    prm->rkb = RKB_DEF;
    prm->rkh = RKH_DEF;
    prm->rkh2 = RKH2_DEF;
    prm->f = F_DEF;
    prm->r = R_DEF;
    prm->a = A_DEF;
    prm->k = K_DEF;
    prm->scheme = SCHEME_DEF;
    prm->rstart = -1;
    prm->infname = NULL;
    prm->outfname = strdup(OUTFNAME_DEF);
    prm->outfnameave = strdup(OUTFNAMEAVE_DEF);
    prm->nobs = 0;
    prm->estd = 0.0;
    prm->dtobs = NAN;
    prm->obsfname = NULL;
    prm->save_q = SAVE_Q_DEF;

    return prm;
}

qgprm* qgprm_read(char* fname)
{
    qgprm* prm = qgprm_create();
    FILE* f = NULL;
    char buf[MAXSTRLEN];
    int line;

    if (fname == NULL) {
        prm->prmfname = strdup("-");
        return prm;
    }

    f = fopen(fname, "r");
    if (f == NULL) {
        int errno_saved = errno;

        quit("could not open \"%s\": %s", fname, strerror(errno_saved));
    }
    prm->prmfname = strdup(fname);

    line = 0;
    while (fgets(buf, MAXSTRLEN, f) != NULL) {
        char seps[] = " =\t\n";
        char* token;

        line++;
        if (buf[0] == '#')
            continue;
        if ((token = strtok(buf, seps)) == NULL)
            continue;
        if (strcasecmp(token, "MREFIN") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: MREFIN not specified", fname, line);
            if (!str2int(token, &prm->mrefin))
                quit("%s, l.%d: could not convert \"%s\" entry to int", fname, line, token);
        } else if (strcasecmp(token, "NX1") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: NX1 not specified", fname, line);
            if (!str2int(token, &prm->nx1))
                quit("%s, l.%d: could not convert \"%s\" entry to int", fname, line, token);
        } else if (strcasecmp(token, "NY1") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: NY1 not specified", fname, line);
            if (!str2int(token, &prm->ny1))
                quit("%s, l.%d: could not convert \"%s\" entry to int", fname, line, token);
        } else if (strcasecmp(token, "LX") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: LX not specified", fname, line);
            if (!str2double(token, &prm->lx))
                quit("%s, l.%d: could not convert \"%s\" entry to double", fname, line, token);
        } else if (strcasecmp(token, "TEND") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: TEND not specified", fname, line);
            if (!str2double(token, &prm->tend))
                quit("%s, l.%d: could not convert \"%s\" entry to double", fname, line, token);
        } else if (strcasecmp(token, "DTOUT") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: DTOUT not specified", fname, line);
            if (!str2double(token, &prm->dtout))
                quit("%s, l.%d: could not convert \"%s\" entry to double", fname, line, token);
        } else if (strcasecmp(token, "DTOUTAVE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: DTOUTAVE not specified", fname, line);
            if (!str2double(token, &prm->dtoutave))
                quit("%s, l.%d: could not convert \"%s\" entry to double", fname, line, token);
        } else if (strcasecmp(token, "DT") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: DT not specified", fname, line);
            if (!str2double(token, &prm->dt))
                quit("%s, l.%d: could not convert \"%s\" entry to double", fname, line, token);
        } else if (strcasecmp(token, "RKB") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: RKB not specified", fname, line);
            if (!str2double(token, &prm->rkb))
                quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "RKH") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: RKH not specified", fname, line);
            if (!str2double(token, &prm->rkh))
                quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "RKH2") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: RKH2 not specified", fname, line);
            if (!str2double(token, &prm->rkh2))
                quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "F") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: F not specified", fname, line);
            if (!str2double(token, &prm->f))
                quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "R") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: R not specified", fname, line);
            if (!str2double(token, &prm->r))
                quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "A") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: A not specified", fname, line);
            if (!str2double(token, &prm->a))
                quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "VERBOSE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: VERBOSE not specified", fname, line);
            if (!str2int(token, &verbose))
                quit("%s, l.%d: could not convert \"%s\" to int", fname, line, token);
        } else if (strcasecmp(token, "SCHEME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: SCHEME not specified", fname, line);
            if (strcasecmp(token, "ORDER1") == 0)
                prm->scheme = SCHEME_ORDER1;
            else if (strcasecmp(token, "ORDER2") == 0)
                prm->scheme = SCHEME_ORDER2;
            else if (strcasecmp(token, "RK4") == 0)
                prm->scheme = SCHEME_RK4;
            else if (strcasecmp(token, "DP5") == 0)
                prm->scheme = SCHEME_DP5;
            else
                quit("%s, l.%d: unexpected scheme \"%s\"", fname, line, token);
        } else if (strcasecmp(token, "INFNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: INFNAME not specified", fname, line);
            prm->infname = strdup(token);
        } else if (strcasecmp(token, "OUTFNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: OUTFNAME not specified", fname, line);
            if (prm->outfname != NULL)
                free(prm->outfname);
            prm->outfname = strdup(token);
        } else if (strcasecmp(token, "OUTFNAMEAVE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: OUTFNAMEAVE not specified", fname, line);
            if (prm->outfnameave != NULL)
                free(prm->outfnameave);
            prm->outfnameave = strdup(token);
        } else if (strcasecmp(token, "RSTART") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: RSTART not specified", fname, line);
            if (!str2int(token, &prm->rstart))
                quit("%s, l.%d: could not convert \"%s\" to int", fname, line, token);
        } else if (strcasecmp(token, "SAVE_Q") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: SAVE_Q not specified", fname, line);
            if (!str2bool(token, &prm->save_q))
                quit("%s, l.%d: could not convert \"%s\" to bool", fname, line, token);
        } else if (strcasecmp(token, "NOBS") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: NOBS not specified", fname, line);
            if (!str2int(token, &prm->nobs))
                quit("%s, l.%d: could not convert \"%s\" entry to int", fname, line, token);
        } else if (strcasecmp(token, "ESTD") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: ESTD not specified", fname, line);
            if (!str2double(token, &prm->estd))
                quit("%s, l.%d: could not convert \"%s\" entry to double", fname, line, token);
        } else if (strcasecmp(token, "DTOBS") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: DTOBS not specified", fname, line);
            if (!str2double(token, &prm->dtobs))
                quit("%s, l.%d: could not convert \"%s\" entry to double", fname, line, token);
        } else if (strcasecmp(token, "OBSFNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                quit("%s, l.%d: OBSFNAME not specified", fname, line);
            if (prm->obsfname != NULL)
                free(prm->obsfname);
            prm->obsfname = strdup(token);
        } else
            quit("%s, l.%d: unexpected entry \"%s\"", fname, line, token);
    }
    fclose(f);

    prm->n = prm->nx1 * ipow(2, prm->mrefin - 1) + 1;
    prm->m = prm->ny1 * ipow(2, prm->mrefin - 1) + 1;

    if (!isfinite(prm->dtout) && !isfinite(prm->dtoutave) && !isfinite(prm->dtobs))
        quit("no output specified (need to set some of dtout, dtoutave, dtobs)");
    if (prm->dtout / prm->dt - floor(prm->dtout / prm->dt) != 0.0)
        quit("dtout has to be integer multiple of dt (dtout / dt = %f)", prm->dtout / prm->dt);
    if (isfinite(prm->dtoutave) && prm->dtoutave / prm->dt - floor(prm->dtoutave / prm->dt) != 0.0)
        quit("dtoutave has to be integer multiple of dt (dtoutave / dt = %f)", prm->dtoutave / prm->dt);
    if (isfinite(prm->dtobs) && prm->dtobs / prm->dt - floor(prm->dtobs / prm->dt) != 0.0)
        quit("dtobs has to be integer multiple of dt (dtobs / dt = %f)", prm->dtobs / prm->dt);
    if (prm->nobs != 0 && prm->dtobs == 0)
        quit("dtobs must be defined when nobs != 0");
    if (prm->nobs != 0 && prm->obsfname == NULL)
        prm->obsfname = strdup(OBSFNAME_DEF);

    return prm;
}

/**
 */
void qgprm_destroy(qgprm* prm)
{
    if (prm->prmfname != NULL)
        free(prm->prmfname);
    if (prm->infname != NULL)
        free(prm->infname);
    if (prm->outfname != NULL)
        free(prm->outfname);
    if (prm->outfnameave != NULL)
        free(prm->outfnameave);
    free(prm);
}

/**
 */
void qgprm_print(qgprm* prm)
{
    printf("  QG parameters:\n");
    if (prm->infname != NULL)
        printf("    infile   = %s\n", prm->infname);
    if (prm->rstart >= 0)
        printf("    rstart   = %d\n", prm->rstart);
    printf("    outfile  = %s\n", prm->outfname);
    if (isfinite(prm->dtoutave))
        printf("    outfileave  = %s\n", prm->outfnameave);
    printf("    save_q   = %s\n", (prm->save_q) ? "yes" : "no");
    printf("    scheme   = %s\n", SCHEME_STR[prm->scheme]);
    printf("    dt       = %.2f\n", prm->dt);
    printf("    tend     = %.0f\n", prm->tend);
    printf("    dtout    = %.1f\n", prm->dtout);
    if (isfinite(prm->dtoutave))
        printf("    dtoutave = %.1f\n", prm->dtoutave);
    if (verbose > 1) {
        printf("    MREFIN   = %d\n", prm->mrefin);
        printf("    NX1      = %d\n", prm->nx1);
        printf("    NY1      = %d\n", prm->ny1);
    }
    printf("    ny       = %d\n", prm->m);
    printf("    nx       = %d\n", prm->n);
    printf("    lx       = %.2f\n", prm->lx);
    if (prm->rkb != 0.0)
        printf("    rkb      = %.1e\n", prm->rkb);
    if (prm->rkh != 0.0)
        printf("    rkh      = %.1e\n", prm->rkh);
    if (prm->rkh2 != 0.0)
        printf("    rkh2     = %.1e\n", prm->rkh2);
    printf("    f        = %.1f\n", prm->f);
    printf("    r        = %.1e\n", prm->r);
    printf("    A        = %.2f\n", prm->a);
    printf("    k        = %.1f\n", prm->k);
    if (prm->nobs != 0) {
        printf("    nobs     = %d\n", prm->nobs);
        printf("    estd     = %.3f\n", prm->estd);
        printf("    dtobs    = %.2f\n", prm->dtobs);
        printf("    obsfile  = %s\n", prm->obsfname);
    }
}

/**
 */
void qgprm_describe(void)
{
    printf("  QG parameter file format:\n");
    printf("  # run parameters:\n");
    printf("    scheme   = {order1 | order2 | rk4 | dp5}");
    if (SCHEME_DEF == SCHEME_ORDER1)
        printf(" (order1)\n");
    else if (SCHEME_DEF == SCHEME_ORDER2)
        printf(" (order2)\n");
    else if (SCHEME_DEF == SCHEME_RK4)
        printf(" (rk4)\n");
    else if (SCHEME_DEF == SCHEME_DP5)
        printf(" (dp5)\n");
    printf("    dt       = <time step>                    (%.2f)\n", DT_DEF);
    printf("    tend     = <end time>                     (%.0f)\n", TEND_DEF);
    printf("    dtout    = <output time interval>         (%.1f)\n", DTOUT_DEF);
    printf("    dtoutave = <output time interval for av.> (NaN)\n");
    printf("    infname  = <input file name>              (<none>)\n");
    printf("    rstart   = <start dump>                   (<last dump>)\n");
    printf("    outfname = <output file name>             (%s)\n", OUTFNAME_DEF);
    printf("    outfnameave = <output file name for av.>  (%s)\n", OUTFNAMEAVE_DEF);
    printf("    save_q   = {yes | no}                     (%s)\n", (SAVE_Q_DEF) ? "yes" : "no");
    printf("  # model parameters:\n");
    printf("    MREFIN   = <small number>                 (%d)\n", MREFIN_DEF);
    printf("    NX1      = <small number>                 (%d)\n", NX1_DEF);
    printf("    NY1      = <small number>                 (%d)\n", NY1_DEF);
    printf("    LX       = <domain size in X direction>   (%.2f)\n", LX_DEF);
    printf("    rkb      = <value>                        (%.1e)\n", RKB_DEF);
    printf("    rkh      = <value>                        (%.1e)\n", RKH_DEF);
    printf("    rkh2     = <value>                        (%.1e)\n", RKH2_DEF);
    printf("    F        = <value>                        (%.1f)\n", F_DEF);
    printf("    r        = <value>                        (%.1e)\n", R_DEF);
    printf("    A        = <value>                        (%.2f)\n", A_DEF);
    printf("    k        = <value>                        (%.1f)\n", K_DEF);
    printf("  # other parameters:\n");
    printf("    verbose  = {0 | 1 | 2}                    (%d)\n", VERBOSE_DEF);
    printf("    nobs     = <number>                       (0)\n");
    printf("    estd     = <value>                        (0.0)\n");
    printf("    dtobs    = <observation time interval>    (NaN)\n");
    printf("    obsfname = <file name>                    (obs_psi.nc)\n");
    printf("\n");
    printf("  Notes:\n");
    printf("    1. (...) denotes a default value\n");
    printf("    2. {...} lists posible choices\n");
    printf("    3. <...> describes the entry\n");
    printf("\n");
    printf("    The model grid has dimensions:\n");
    printf("      (ny, nx) = (NY1 * 2^(MREFIN - 1) + 1, NX1 * 2^(MREFIN - 1) + 1)\n");
    printf("    The cells are of square shape, width = height = LX / (NX - 1)\n");
    printf("    For the model equation see qg.pdf.\n");
}
