#include "splineMI.h"
#include "pearsonCorr.h"
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>


// filtering out the probes that has low correlation with others
void probe_filter(double const *x, int m, int n, int *grp, int mg, int useCorr, double th){ 
	int i, j;
	double r;
	double *xx = (double*) calloc(n, sizeof(double));
	double *yy = (double*) calloc(n, sizeof(double));
	double *xsum = (double*) calloc(mg*n, sizeof(double));

	for(i = 0; i < m; i++){
		if(grp[i] < 0) continue;
		for(j = 0; j < n; j++){
			xsum[grp[i] + j*mg] += x[i + j*m];
		}
	}

	for(i = 0; i < m; i++){
		if(grp[i] < 0) continue;
		for(j = 0; j < n; j++){
			xx[j] = x[i + j*m];
			yy[j] = xsum[grp[i] + j*mg];
		}
		r = useCorr? cor(xx, yy, n) : mi2(xx, yy, n, 6, 3, 1, 1);
		Rprintf("row: %d\t group: %d\t corr: %f\n", i, grp[i], r);
		if(r < th) grp[i] = -1;
	}

	free(xx);
	free(yy);
	free(xsum);
}


void probe_summarization(double const *x, int m, int n, int *grp, double *xs, int mg, int useCorr, double th){
	int i, j;
	int *npbs = (int*) calloc(mg, sizeof(int));

	probe_filter(x, m, n, grp, mg, useCorr, th);

	for(i = 0; i < mg; i++) npbs[i] = 0;

	for(i = 0; i < m; i++){
		if(grp[i] < 0) continue;
		npbs[grp[i]]++;
		for(j = 0; j < n; j++){
			xs[grp[i] + j*mg] += x[i + j*m];
		}
	}

	//printf("npbs: \n");
	//for(i = 0; i < mg; i++) printf("%d\t", npbs[i]);
	//printf("\n");

	for(i = 0; i < mg; i++){
		for(j = 0; j < n; j++){
			xs[i + j*mg] /= npbs[i];
		}
	}

	free(npbs);
}

SEXP probe_summarizationR2C(SEXP xIn, SEXP mIn, SEXP nIn, SEXP grpIn, SEXP mgIn, SEXP useCorrIn, SEXP thIn){
	SEXP out;
	double *x, *o, th;
	int *grp, m, n, mg, useCorr; 

	
	PROTECT(xIn = AS_NUMERIC(xIn));
	PROTECT(mIn = AS_INTEGER(mIn));
	PROTECT(mgIn = AS_INTEGER(mgIn));
	PROTECT(nIn = AS_INTEGER(nIn));
	PROTECT(grpIn = AS_NUMERIC(grpIn));
	PROTECT(useCorrIn = AS_INTEGER(useCorrIn));
	PROTECT(thIn = AS_NUMERIC(thIn));
	PROTECT(out = NEW_NUMERIC(mg * n));	

	x = NUMERIC_POINTER(xIn);
	m = INTEGER_POINTER(mIn)[0];
	mg = INTEGER_POINTER(mgIn)[0];
	n = INTEGER_POINTER(nIn)[0];
	grp = INTEGER_POINTER(grpIn);
	useCorr = INTEGER_POINTER(useCorrIn)[0];
	thIn = AS_NUMERIC(th)[0];
	o = NUMERIC_POINTER(out);

	probe_summarization(x, m, n grp, o, mg, useCorr, th);

	UNPROTECT(8);
	return (out);

}


/*
int main(){
	double x[16] = { 1, 1, 2, 3, 2, 2, 2, 4, 3, 3, 2, 5, 4, 4, 2, 6};
	int grp[4] = {0, 0, 0, 1};
	int m=4, n=4, mg=2, i, j;

	double *xs = (double*) calloc(16, sizeof(double));
	for(i = 0; i < mg; i++) for(j = 0; j < n; j++) xs[i + j*m] = 0;

	printf("Before summarization: \n");
	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			printf("%f\t",x[i + j*m]);
		}printf("\n");
	}
	probe_summarization(x, m, n, grp, xs, mg, 1, 0.6);
	printf("After summarization: \n");
	for(i = 0; i < mg; i++){
		for(j = 0; j < n; j++){
			printf("%f\t",xs[i + j*mg]);
		}printf("\n");
	}
	
	free(xs);
	return 0;
}
*/

