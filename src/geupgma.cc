#include "geupgma.h"

void PrintTriangularMatrix(const double *px, int m, int n){
  uint32_t c = 0;
  for(int i = 0; i < m; i++){
    for(int j = 0; j < m; j++){
      if (i < j) {
        Rprintf("%f,", px[c++]);
      } else {
        Rprintf("X,");
      }
    }
    Rprintf("\n");
  }
}

void UpdateSimMatrix(
    double *psim,
    int *cluster_sizes,
    double *px,
    const int m,
    const int n,
    int i,      // row index of sim in psim matrix
    int j       // col index of sim in psim matrix
    ){
  //Rprintf("**UpdateSimMatrix\n");
  for (int k = 0; k < j; k++){ 
    psim[TriangularIndex(k, j, m)] = 0;
  }
  for (int k = j+1; k < m; k++){ 
    psim[TriangularIndex(j, k, m)] = 0;
  }
  for (int k = 0; k < n; k++) {
    px[i + k*m] = (px[i + k*m] * cluster_sizes[i] + 
                   px[j + k*m] * cluster_sizes[j]) / 
                  (cluster_sizes[i] + cluster_sizes[j]);
    px[j + k*m] = 0;
  }
  cluster_sizes[i] += cluster_sizes[j];
  cluster_sizes[j] = 0;
  
  double *rs = new double[m];
  double *vec = new double[n];
  for (int k = 0; k < n; k++) vec[k] = px[i + k*m];
  GetAllCorWz(px, vec, m, n, rs);
  for (int k = 0; k < i; k++){ 
    psim[TriangularIndex(k, i, m)] = rs[k];
  }
  for (int k = i+1; k < m; k++){ 
    psim[TriangularIndex(i, k, m)] = rs[k];
  }
  delete [] vec;
  delete [] rs;
  //Rprintf("**UpdateSimMatrix Done\n");
}

void GEUPGMA (
    double *psim, 
    const double *px, 
    const uint32_t k, 
    const int m,
    const int n,
    double *ptrace) {
  uint32_t *o = new uint32_t[k];
  int *cluster_sizes = new int[m];
  for(int i = 0; i < m; i++) cluster_sizes[i] = 1;
  double *px_copy = new double[m*n];
  std::copy (px, px + m*n, px_copy);
  //PrintTriangularMatrix(psim, m, n);
  Rprintf("Sorting similarity matrix ... \n");
  OrderL(psim, k, o, true);
  double maxcor = psim[o[k-1]];
  
  int cnt = 0;
  while (fabs(maxcor) > 1E-10) {
    if (cnt >= m-1) {
      Rprintf("ERROR: There is something wrong!\n");
      break;
    }
    int i = RowIndex(o[k-1], m);
    int j = ColIndex(o[k-1], m, i);
    Rprintf("(%d) %d\t%d\t%f\n", cnt+1, i, j, maxcor);
    ptrace[cnt + 0*m] = i;
    ptrace[cnt + 1*m] = j;
    ptrace[cnt + 2*m] = maxcor;
    UpdateSimMatrix(psim, cluster_sizes, px_copy, m, n, i, j);
    //PrintTriangularMatrix(psim, m, n);
    //Rprintf("**OrderL\n");
    OrderL(psim, k, o, false);
    //Rprintf("**OrderL Done\n");
    maxcor = psim[o[k-1]];
    cnt++;
  }
  ptrace[cnt + 0*m] = -999;
  delete [] px_copy;
  delete [] o;
  delete [] cluster_sizes;
}

SEXP GEUPGMACC(SEXP sim, SEXP x){
  double *psim = REAL(sim);
  double *px = REAL(x);
  SEXP dim = getAttrib(sim, R_DimSymbol);
  uint32_t k = (uint32_t) length(sim);
  int m = (int)(0.5 + sqrt(1+2*k));
  int n = (int)(length(x) / m);
  SEXP trace = PROTECT(allocMatrix(REALSXP, m, 3));
  double *ptrace = REAL(trace);
  GEUPGMA(psim, px, k, m, n, ptrace);
  
  UNPROTECT(1);
  return trace;
}

