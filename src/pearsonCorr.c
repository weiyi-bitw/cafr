#include "pearsonCorr.h"


double cor(const double* x, const double* y, int n){
	double xMean = 0, yMean = 0, xSd = 0, ySd = 0, rho = 0;
	int i;
	for(i = 0; i < n; i++){
		xMean += x[i];
		yMean += y[i];
		xSd += x[i] * x[i];
		ySd += y[i] * y[i];
		rho += x[i] * y[i];
	}
	xMean /= n;
	xSd = xSd - n * xMean * xMean;
	xSd = xSd <= 0? 1 : sqrt(xSd);
	yMean /= n;
	ySd = ySd - n * yMean * yMean;
	ySd = ySd <= 0? 1 : sqrt(ySd);
	rho = (rho - n * xMean * yMean) / xSd / ySd;
	return rho;
}
void corR(const double* x, const double* y, int* n, double *cOut){
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

int RowIndex(long i, int m){
  double m_d = m-1;
  double row = (-2*m_d - 1 + sqrt((4*m_d*(m_d+1) - 8*(double)i - 7))) / -2;
  if(row == (double)(int) row) row -= 1;
  return (int)row;
}

int ColIndex(long i, int m, int row){
  return i - (m-2)*row + row*(row-1)/2 + 1;
}

void PairwiseCor(
  const double *data, 
  long *idx_start,
  long *all_tasks,
  int *m, 
  int *n,
  int *buffer_exp, 
  double *out
){
  const int kBufferSize = 1 << *(buffer_exp);
  const int kOutRowNum = 3;
  unsigned int i, j;
  double x_mean, y_mean, x_sd, y_sd, rho;
  int x, y, r1, r2;

  i = *(idx_start) & (kBufferSize-1);
  do {
    r1 = RowIndex(*(idx_start), *m);
    r2 = ColIndex(*(idx_start), *m, r1);
    x_mean = 0;
    y_mean = 0;
    x_sd = 0;
    y_sd = 0;
    rho = 0;
    for(j = 0; j < *n; ++j) {
      x = data[r1 + j*(*m)];
      y = data[r2 + j*(*m)];
      x_mean += x;
      y_mean += y;
      x_sd += x*x;
      y_sd += y*y;
      rho += x*y;
    }
    x_mean /= *n;
    y_mean /= *n;
    x_sd = sqrt(x_sd - *n * x_mean * x_mean);
    y_sd = sqrt(y_sd - *n * y_mean * y_mean);
    rho = (rho - *n * x_mean * y_mean) / x_sd / y_sd;
    out[kOutRowNum * i] = r1;
    out[kOutRowNum * i + 1] = r2;
    out[kOutRowNum * i + 2] = rho;
    ++*(idx_start);
    ++i;

  } while ( (i < kBufferSize) && (*(idx_start) < *(all_tasks)) );
}

void getAllCorWz(const double *data, const double *vec , int *m, int *n, double *rs){
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
