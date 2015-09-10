#include <stdlib.h>
#include <math.h>
#include <stdint.h>

//#include <R.h>

#ifndef CAFR_UTIL_H_
#define CAFR_UTIL_H_

typedef int bool;
#define true 1
#define false 0

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

void QuickSortR(double const*, int*, int, int);
void Order(double const*, int, int*);
void QuickSortRL(double const*, uint32_t*, uint32_t, uint32_t);
void OrderL(double const*, uint32_t, uint32_t*, bool);
int RowIndex(
    long, //i
    int); //m
int ColIndex(
    long, //i
    int,  //m
    int); //row

long TriangularIndex(
    int, //i
    int, //j
    int); //m

#ifdef __cplusplus
} //extern C
#endif //__cplusplus

#endif //CAFR_UTIL_H_
