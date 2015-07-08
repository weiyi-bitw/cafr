#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "util.h"

#ifndef CAFR_PEARSON_CORR_H_
#define CAFR_PEARSON_CORR_H_

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

double PearsonCorr(const double*, const double*, int);
void corR(const double* , const double* , int* , double *);
void getAllCorWz(const double *, const double * , int *, int *, double *);


SEXP AllPairwiseCorCC(
    SEXP); //arr
/*
SEXP PairwiseCor(
    SEXP, //_data
    SEXP, //_idx_start
    SEXP, //_all_tassk 
    SEXP, //_m
    SEXP, //_n 
    SEXP); // buffer_exp
*/
#ifdef __cplusplus
} //extern C
#endif //__cplusplus

#endif //CAFR_PEARSON_CORR_H_
