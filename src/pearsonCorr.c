#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#ifndef NaN
#define NaN -sqrt(-1)
#endif

double cor(const double* x, const double* y, int* n, double *cOut){
	double xMean = 0, yMean = 0, xSd = 0, ySd = 0, rho = 0;
	int i;
	for(i = 0; i < *n; i++){
		xMean += x[i];
		yMean += y[i];
		xSd += x[i] * x[i];
		ySd += y[i] * y[i];
		rho += x[i] * y[i];
	}
	xMean /= *n;
	xSd = sqrt(xSd - *n * xMean * xMean);
	yMean /= *n;
	ySd = sqrt(ySd - *n * yMean * yMean);
	rho = (rho - *n * xMean * yMean) / xSd / ySd;
	*cOut = rho;
}
double getAllCorWz(const double *data, const double *vec , int *m, int *n, double *rs){
	double vMean = 0, yMean, vSd = 0, ySd;
	int i, j;
	double y;

	for(j = 0; j < *n; j++){
		vMean += vec[j];
		vSd += vec[j] * vec[j];
	}
	vMean /= *n;
	vSd = sqrt(vSd - *n * vMean * vMean);

	for(i = 0; i < *m; i++){
		yMean = 0;
		ySd = 0;
		rs[i] = 0;
		for(j = 0; j < *n; j++){
			y = data[i + j * (*m)];
			yMean += y;
			ySd += y * y;
			rs[i] += vec[j] * y;
		}
		yMean /= *n;
		ySd = sqrt(ySd - *n * yMean * yMean);
		rs[i] = (rs[i] - *n * vMean * yMean) / vSd / ySd;
	}
}
