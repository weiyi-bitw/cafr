#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <R.h>

#ifndef CAFR_SPLINE_MI_H_
#define CAFR_SPLINE_MI_H_

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

float log2f(float);
double log2d(double);
void xToZ(const double*, double*, int, int, int, double, double);
double mean(double*, int);
double std(double*,int);
double meani(int*, int);
double stdi(int*,int);
/* void SplineKnots(int*,int,int); */
void knotVector(double*, int, int);
void findWeights(const double *, const double *, double *, int, int, int, double, double);
double entropy1(const double*, int, int);
double mi2(const double*, const double*, int, int, int, int, int);
double mid2(const double*, const double*, int, int, int);

// export R function
void mi2R(const double *, const double *, int *, int *, int *, double *, int *, int *);
void mid2R(const double *, const double *, int *, int *, int *, double *);
void mi2DiffBins(const double *, const double *, int *, int *, int *, int *, int *, double *, int *, int *);
void mi2vs1(const double *, const double *, const double *, int *, int *, int *, double *, int *);
void getAllMIWz(const double *, const double* , double *, int , int , int , int , int , int );
void getAllMIWz_R(const double *, const double* , double *, int *, int *, int *, int *, int *, int *);

#ifdef __cplusplus
} //extern C
#endif //__cplusplus

#endif //CAFR_SPLINE_MI_H_
