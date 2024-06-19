/******************************************************************************
 *
 * File:        utils.h        
 *
 * Created:     15/04/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: header for utils.c
 *
 * Revisions:   
 *
 *****************************************************************************/

#if !defined(_UTILS_H)

void quit(char* format, ...);
int str2double(char* token, double* value);
int str2int(char* token, int* value);
int str2bool(char* token, int* value);
void* alloc2d(size_t nj, size_t ni, size_t unitsize);
void printtime(const char offset[]);
int file_exists(char* fname);
void get_obspos(int n, int mn, long int seed, int pos[], double estd, float error[]);

#define _UTILS_H
#endif
