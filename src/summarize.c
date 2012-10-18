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
		//Rprintf("row: %d\t group: %d\t corr: %f\n", i, grp[i], r);
		if(r < th) grp[i] = -1;
	}

	free(xx);
	free(yy);
	free(xsum);
}


void probe_summarization(double const *x, int m, int n, int *grp, double *xs, int mg, int useCorr, double th, int verbose){
	int i, j;
	int *npbs = (int*) calloc(mg, sizeof(int));

	if(verbose) Rprintf("Filtering out uncorrelated probes...\n");
	probe_filter(x, m, n, grp, mg, useCorr, th);

	for(i = 0; i < mg; i++){
		 npbs[i] = 0;
		for(j = 0; j < n; j++){
			xs[i + j*mg] = 0;
		}
	}

	if(verbose) Rprintf("Summarize gene-level expression values...\n");
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
		if(npbs[i] <= 0) continue;
		for(j = 0; j < n; j++) xs[i + j*mg] /= npbs[i];
		//Rprintf("\n");
	}
	if(verbose) Rprintf("Done.\n");
	free(npbs);
}

SEXP probe_summarizationR2C(SEXP xIn, SEXP mIn, SEXP nIn, SEXP grpIn, SEXP mgIn, SEXP useCorrIn, SEXP thIn, SEXP vIn){
	SEXP out;
	double *x, *o, th;
	int *grp, m, n, mg, useCorr, verbose; 
	R_len_t no;
	

	PROTECT(xIn = AS_NUMERIC(xIn));
	PROTECT(mIn = AS_INTEGER(mIn));
	PROTECT(mgIn = AS_INTEGER(mgIn));
	PROTECT(nIn = AS_INTEGER(nIn));
	PROTECT(grpIn = AS_INTEGER(grpIn));
	PROTECT(useCorrIn = AS_INTEGER(useCorrIn));
	PROTECT(thIn = AS_NUMERIC(thIn));
	PROTECT(vIn = AS_INTEGER(vIn));
	

	x = NUMERIC_POINTER(xIn);
	m = INTEGER_POINTER(mIn)[0];
	mg = INTEGER_POINTER(mgIn)[0];
	n = INTEGER_POINTER(nIn)[0];
	grp = INTEGER_POINTER(grpIn);
	useCorr = INTEGER_POINTER(useCorrIn)[0];
	th = NUMERIC_POINTER(thIn)[0];
	verbose = INTEGER_POINTER(vIn)[0];


	no = (R_len_t) mg * n;
	//Rprintf("numO: %d\n", no);
	PROTECT(out = NEW_NUMERIC(no));	
	o = NUMERIC_POINTER(out);

	probe_summarization(x, m, n, grp, o, mg, useCorr, th, verbose);

	UNPROTECT(9);
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

