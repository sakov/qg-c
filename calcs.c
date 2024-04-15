/******************************************************************************
 *
 * File:        calcs.c
 *
 * Created:     15/04/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: contains mathematical utilities
 *
 * Revisions:
 *
 *****************************************************************************/

#include "calcs.h"

/**
 */
void arakawa(double h, int m, int n, double** A, double** B, double** J)
{
    double c = 12.0 * h * h;
    int i, j;

    for (j = 1; j < m - 1; ++j) {
        double* am1 = A[j - 1];
        double* a = A[j];
        double* ap1 = A[j + 1];
        double* bm1 = B[j - 1];
        double* b = B[j];
        double* bp1 = B[j + 1];

        for (i = 1; i < n - 1; ++i)
            J[j][i] = ((am1[i] - a[i - 1]) * bm1[i - 1]
                       + (am1[i - 1] + am1[i] - ap1[i - 1] - ap1[i]) * b[i - 1]
                       + (a[i - 1] - ap1[i]) * bp1[i - 1]
                       + (am1[i + 1] + a[i + 1] - am1[i - 1] - a[i - 1]) * bm1[i]
                       + (a[i - 1] + ap1[i - 1] - a[i + 1] - ap1[i + 1]) * bp1[i]
                       + (a[i + 1] - am1[i]) * bm1[i + 1]
                       + (ap1[i] + ap1[i + 1] - am1[i] - am1[i + 1]) * b[i + 1]
                       + (ap1[i] - a[i + 1]) * bp1[i + 1]) / c;
    }
}

/**
 */
void laplacian(double h, int m, int n, double** A, double** L)
{
    double h2 = 1.0 / h / h;
    double h4 = 4.0 * h2;
    int i, j;

    for (j = 1; j < m - 1; ++j) {
        double* am1 = A[j - 1];
        double* a = A[j];
        double* ap1 = A[j + 1];
        double* lj = L[j];

        for (i = 1; i < n - 1; ++i)
            lj[i] = (am1[i] + ap1[i]) * h2 + (a[i - 1] + a[i + 1]) * h2 - a[i] * h4;
    }
}
