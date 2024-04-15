/******************************************************************************
 *
 * File:        helmholtz.h        
 *
 * Created:     15/04/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: header for helmholtz.f90
 *
 * Revisions:   
 *
 *****************************************************************************/

#if !defined(_HELMHOLTZ_H)

void __helmholtz_mod_MOD_helmholtz(double* G, double* F, double* zk2, int* nx1, int* ny1, int* mrefin, double* dx1, int* nu1, int* nu2, int* ncyc, int* nx, int* my, double* P);

#define _HELMHOLTZ_H
#endif
