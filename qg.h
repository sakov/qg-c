/******************************************************************************
 *
 * File:        qg.h        
 *
 * Created:     15/04/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: header for qg.c
 *
 * Revisions:   
 *
 *****************************************************************************/

#if !defined(_QG_H)

#define PROGRAM_NAME "qg"
#define PACKAGE_NAME "QG-C"
#define MAXSTRLEN 4096
#define VERBOSE_DEF 1
#define SCHEME_ORDER1 0
#define SCHEME_ORDER2 1
#define SCHEME_RK4 2
#define SCHEME_DP5 3
#define NETCDF_FORMAT NC_NETCDF4
#define QMAX 1.0e+6

extern int verbose;
extern int force;
extern const char* QG_VERSION;
extern const char* SCHEME_STR[];

#define _QG_H
#endif
