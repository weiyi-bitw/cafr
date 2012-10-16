#include <R.h>
#include <Rdefines.h>

SEXP testNA(SEXP i){
	R_len_t n;
	int j;
	double *xi;
	int *xo;
	SEXP o;
	
	n = LENGTH(i);
	PROTECT(i = AS_NUMERIC(i));
	PROTECT(o = NEW_NUMERIC(n));
	xi = NUMERIC_POINTER(i);
	xo = NUMERIC_POINTER(o);

	for(j = 0; j < n; j++){
		if(ISNA(xi[j])){
			xo[j] = NA_REAL;
		}else{
			xo[j] = j;
		}
	}
	UNPROTECT(2);
	return(o);

}
