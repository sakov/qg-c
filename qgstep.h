/******************************************************************************
 *
 * File:        qgstep.h        
 *
 * Created:     15/04/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: header for qgstep.c
 *
 * Revisions:   
 *
 *****************************************************************************/

#if !defined(_QGSTEP_H)

void calc_psi(model* m, double** psiguess, double** q, double** psi);
void qg_step_order1(model* m);
void qg_step_order2(model* m);
void qg_step_rk4(model* m);
void qg_step_dp5(model* m);

#define _QGSTEP_H
#endif
