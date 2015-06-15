#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <R.h>
#include <Rdefines.h>
#include "spline_mi.h"
#include "fuf.h"

void getWeightedMetagene(const double *data, const double *mi, double *o, double a, int m, int n){
	int i, j;
	double wsum = 0, wi;
	for(j = 0; j < n; j++) o[j] = 0;
	

	for(i = 0; i < m; i++){
		if(mi[i] > 0){
			wi = exp( a * log(mi[i]) );
			wsum += wi;
			for(j = 0; j < n; j++) o[j] += data[i + j*m] * wi;
		}
	}
	for(j = 0; j < n; j++) o[j] /= wsum;
}

double calcMSE(const double *mi, const double *premi, int m){
	int i;
	double err = 0;
	for(i = 0; i < m; i++) err += ( (mi[i] - premi[i]) * (mi[i] - premi[i]) );
	return err;
}

void caf(const double *data, const double *vec, double *mi, double a, int maxIter, double epsilon, int m, int n, int bin, int so, int negateMI, int verbose){
	int c = 0, converge = 0, i, nm;
	double err;
	double *premi = (double*) calloc(m, sizeof(double));
	double *metag = (double*) calloc(n, sizeof(double));
	int *oi = (int*) calloc(m, sizeof(int));

	getAllMIWz(data, vec, mi, m, n, bin, so, 1, negateMI);
	memcpy(premi, mi, m*sizeof(double));

	while(c < maxIter){
		getWeightedMetagene(data, mi, metag, a, m, n);
		getAllMIWz(data, metag, mi, m, n, bin, so, 1, negateMI);
		if(verbose){
			order(mi, m, oi);
			Rprintf("Iteration %d\n", c);
			Rprintf("Gene Index\tMI\n");
			nm = 20 > m? m : 20;
			for(i = 0; i < nm; i++){
				Rprintf("%d \t %f \n", oi[m-i-1], mi[oi[m-i-1]]);
			}
		}


		err = calcMSE(mi, premi, m);
		if(verbose) Rprintf("Delta: %e\n\n", err);
		if(err < epsilon) {
			converge = 1;
			break;
		}
		memcpy(premi, mi, m*sizeof(double));
		c++;
	}
	if(converge != 1){
		mi[0] = -999;
	}

	free(oi);
	free(metag);
	free(premi);
}

SEXP cafR2C(SEXP dataIn, SEXP vecIn, SEXP aIn, SEXP maxIterIn, SEXP epsilonIn, SEXP mIn, SEXP nIn, SEXP binIn, SEXP soIn, SEXP negateMIIn, SEXP verboseIn){
	SEXP miOut;
	double *data, *vec, *mi, a, epsilon;
	int maxIter, m, n, bin, so, negateMI, verbose;
	R_len_t mo;
	
	PROTECT(dataIn = AS_NUMERIC(dataIn));
	PROTECT(vecIn = AS_NUMERIC(vecIn));
	PROTECT(aIn = AS_NUMERIC(aIn));
	PROTECT(maxIterIn = AS_INTEGER(maxIterIn));
	PROTECT(epsilonIn = AS_NUMERIC(epsilonIn));
	PROTECT(mIn = AS_INTEGER(mIn));
	PROTECT(nIn = AS_INTEGER(nIn));
	PROTECT(binIn = AS_INTEGER(binIn));
	PROTECT(soIn = AS_INTEGER(soIn));
	PROTECT(negateMIIn = AS_INTEGER(negateMIIn));
	PROTECT(verboseIn = AS_INTEGER(verboseIn));

	//Rprintf("finish PROTECT\n");

	data = NUMERIC_POINTER(dataIn);
	vec = NUMERIC_POINTER(vecIn);
	a = NUMERIC_POINTER(aIn)[0];
	maxIter = INTEGER_POINTER(maxIterIn)[0];
	epsilon = NUMERIC_POINTER(epsilonIn)[0];
	m = INTEGER_POINTER(mIn)[0];
	n = INTEGER_POINTER(nIn)[0];
	bin = INTEGER_POINTER(binIn)[0];
	so = INTEGER_POINTER(soIn)[0];
	negateMI = INTEGER_POINTER(negateMIIn)[0];
	verbose = INTEGER_POINTER(verboseIn)[0];

	//Rprintf("finish POINTER\n");

	mo = (R_len_t) m;
	PROTECT(miOut = NEW_NUMERIC(mo));
	mi = NUMERIC_POINTER(miOut);

	//Rprintf("finish output\n");

	caf(data, vec, mi, a, maxIter, epsilon, m, n, bin, so, negateMI, verbose);

	UNPROTECT(12);
	return(miOut);
}

/*
int main(){
	int i, j, m = 6, n = 15;
	double a[15] = {};
	double b[90] = {};
	double *mi = (double*) calloc(m, sizeof(double));

	srand(time(NULL));
	
	for(i = 0; i < 15; i++) a[i] = (double) (rand() % 100) / 100;
	for(i = 0; i < 90; i++) b[i] = (double) (rand() % 100) / 100;

	printf("a:\n");
	for(i = 0; i < 15; i++) printf("%f\t", a[i]);
	printf("\n");

	printf("b:\n");
	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			printf("%f\t", b[i + j*m]);	
		}
		printf("\n");
	}
	//getAllMIWz(b, a, mi, m, n, 6, 3, 1, 1);
	caf(b, a, mi, 5, 100, 1E-14, m, n, 6, 3, 1, 1);	

	printf("mi:\n");
	for(i = 0; i < m; i++) printf("%f\t", mi[i]);
	printf("\n");

	free(mi);
	return 0;
}

*/

