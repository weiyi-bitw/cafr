#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>

double cor(const double*, const double*, int);
void corR(const double* , const double* , int* , double *);
void getAllCorWz(const double *, const double * , int *, int *, double *);
SEXP PairwiseCor(
  SEXP, //_data
  SEXP, //_idx_start
  SEXP, //_all_tassk 
  SEXP, //_m
  SEXP, //_n 
  SEXP //_buffer_exp
);
