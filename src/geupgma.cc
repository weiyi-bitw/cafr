#include "geupgma.h"

void GEUPGMA (const double *psim, const double *px, 
              uint32_t k, int m, double *ptrace) {
  uint32_t *o = new uint32_t[k];
  Rprintf("Sorting similarity matrix ... \n");
  OrderL(psim, k, o);
  double maxcor = psim[o[k-1]];
  int cnt = 0;
  int *cluster_sizes = new int[m];
  for(int i = 0; i < m; i++) cluster_sizes[i] = 1;
  while (maxcor > 0) {
    int i = RowIndex(o[k-1], m);
    int j = ColIndex(o[k-1], m, i);
    Rprintf("%d\t%d\t%f\n", i, j, maxcor);
    ptrace[cnt + 0*m] = i;
    ptrace[cnt + 1*m] = j;
    ptrace[cnt + 2*m] = maxcor;
    
    break;
  }

  delete [] o;
  delete [] cluster_sizes;
}

SEXP GEUPGMACC(SEXP sim, SEXP x){
  double *psim = REAL(sim);
  double *px = REAL(x);
  SEXP dim = getAttrib(sim, R_DimSymbol);
  uint32_t k = (uint32_t) length(sim);
  int m = (int)(0.5 + sqrt(1+2*k));
  SEXP trace = PROTECT(allocMatrix(REALSXP, m, 3));
  double *ptrace = REAL(trace);
  UPGMA(psim, px, k, m, ptrace);
  
  UNPROTECT(1);
  return trace;
}


