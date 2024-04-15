/******************************************************************************
 *
 * File:        utils.c
 *
 * Created:     15/04/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: utilities for QG-C package
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include "qg.h"
#include "utils.h"

/**
 */
void quit(char* format, ...)
{
    va_list args;

    fflush(stdout);

    fprintf(stderr, "\n  error: %s: ", PROGRAM_NAME);
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, "\n\n");
    abort();                    /* raise SIGABRT for debugging */
    exit(1);
}

/**
 */
int str2double(char* token, double* value)
{
    char* end = NULL;

    if (token == NULL) {
        *value = NAN;
        return 0;
    }

    *value = strtod(token, &end);

    if (end == token) {
        *value = NAN;
        return 0;
    }

    return 1;
}

/**
 */
int str2int(char* token, int* value)
{
    long int tmp;
    char* end = NULL;

    if (token == NULL) {
        *value = INT_MAX;
        return 0;
    }

    tmp = strtol(token, &end, 10);

    if (end == token || tmp > INT_MAX || tmp < INT_MIN) {
        *value = INT_MAX;
        return 0;
    }

    *value = (int) tmp;
    return 1;
}

/**
 */
int str2bool(char* token, int* value)
{
    if (token[0] == 'y' || token[0] == 'Y' || token[0] == 't' || token[0] == 'T')
        *value = 1;
    else if (token[0] == 'n' || token[0] == 'N' || token[0] == 'f' || token[0] == 'F')
        *value = 0;
    else if (!str2int(token, value))
        return 0;

    if (*value == 0 || *value == 1)
        return 1;
    return 0;
}

/**
 */
void* alloc2d(size_t nj, size_t ni, size_t unitsize)
{
    size_t size;
    void* p;
    void** pp;
    int i;

    if (ni <= 0 || nj <= 0)
        quit("alloc2d(): invalid size (nj = %d, ni = %d)", nj, ni);

    size = nj * sizeof(void*) + nj * ni * unitsize;
    if ((p = malloc(size)) == NULL) {
        int errno_saved = errno;

        quit("alloc2d(): %s", strerror(errno_saved));
    }
    memset(p, 0, size);

    pp = p;
    p = &((size_t*) p)[nj];
    for (i = 0; i < nj; ++i)
        pp[i] = &((char*) p)[i * ni * unitsize];

    return pp;
}

/**
 */
void printtime(const char offset[])
{
    time_t t;
    struct tm tm;

    t = time(NULL);
    tm = *localtime(&t);

    printf("%s%d-%02d-%02d %02d:%02d:%02d\n", offset, tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
}

/**
 */
int file_exists(char* fname)
{
    FILE* f;

    f = fopen(fname, "r");
    if (f == NULL)
        return 0;
    fclose(f);
    return 1;
}