/******************************************************************************
 *
 * File:        calcs.h        
 *
 * Created:     15/04/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: header for calcs.c
 *
 * Revisions:   
 *
 *****************************************************************************/

#if !defined(_CALCS_H)

void arakawa(double h, int m, int n, double** A, double** B, double** J);
void laplacian(double h, int m, int n, double** A, double** L);

#define _CALCS_H
#endif
