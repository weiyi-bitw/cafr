#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include "fuf.h"

void cwknn(double const *db, double const *x, int k, int m, int n, int const *target, R_len_t nt, double* o){
	int i, j, mm;
	R_len_t t;
	double y, xsum = 0, xmean, ymean, xsqr = 0, xsd, ysd, wsum; 
	double *rs = (double*) calloc(n, sizeof(double));
	int *sgn = (int*) calloc(n, sizeof(int));	
	int *oo = (int*) calloc(n, sizeof(int));	


	for(i = 0; i < m; i++){
		if(!ISNA(x[i])){
			xsum += x[i]; 
			xsqr += x[i] * x[i];
		}
	}
	//Rprintf("%f \t %f \n", xsum, xsqr);
	//Rprintf("Calculating correlations...\n");
	for(j = 0; j < n; j++){
		rs[j] = 0;
		mm = m;
		ymean = 0;
		xmean = xsum;
		ysd = 0;
		xsd = xsqr;
		for(i = 0; i < m; i++){
			y = db[i + j * m];
			if( ISNA(x[i]) || ISNA(y)){
				mm--;
				if(!ISNA(x[i])){
					xmean -= x[i];
					xsd -= x[i]*x[i];
				}
			}else{
				ymean += y;
				ysd += y*y;
				rs[j] += x[i] * y;
			}
		}
		xmean /= mm;
		ymean /= mm;
		xsd = sqrt(xsd - mm * xmean * xmean);
		ysd = sqrt(ysd - mm * ymean * ymean);
		rs[j] = (rs[j] - mm * xmean * ymean) / xsd / ysd;
		//Rprintf("%d \t %f \t %f \t %f \t %f \n", mm, xmean, ymean, xsd, ysd);
	}
	/*
	for(j = 0; j < n; j++){
		if(rs[j] < 0){
			rs[j] = -rs[j];
			sgn[j] = -1;
		}else sgn[j] = 1;
	}
	*/
	//Rprintf("Calculating weighted avg...\n");
	order(rs, n, oo);
	for(t = 0; t < nt; t++){
		o[t] = 0;
		wsum = 0;
		for(j = 0; j < k; j++){
			//Rprintf("%f*%f\t", rs[oo[n-1-j]], db[target[t] + oo[n-1-j] * m]);
			o[t] += db[target[t] + oo[n-1-j] * m] * rs[oo[n-1-j]];
			wsum += rs[oo[n-1-j]];
		}//Rprintf("\n");
		o[t] /= wsum;
	}
	//Rprintf("Done!\n");
	free(rs);
	free(sgn);
	free(oo);
}



SEXP cwknnR2C(SEXP dbIn, SEXP xIn, SEXP kIn, SEXP mIn, SEXP nIn, SEXP tgtIn){
	SEXP out;
	double *db, *x, *o;
	int m, n, k, *target;
	R_len_t nt;

	nt = LENGTH(tgtIn);

	PROTECT(dbIn = AS_NUMERIC(dbIn));
	PROTECT(xIn = AS_NUMERIC(xIn));
	PROTECT(mIn = AS_NUMERIC(mIn));
	PROTECT(nIn = AS_NUMERIC(nIn));
	PROTECT(kIn = AS_NUMERIC(kIn));
	PROTECT(tgtIn = AS_INTEGER(tgtIn));
	PROTECT(out = NEW_NUMERIC(nt));

	db = NUMERIC_POINTER(dbIn);
	x = NUMERIC_POINTER(xIn);
	m = NUMERIC_POINTER(mIn)[0];
	n = NUMERIC_POINTER(nIn)[0];
	k = NUMERIC_POINTER(kIn)[0];
	target = INTEGER_POINTER(tgtIn);
	o = NUMERIC_POINTER(out);	
	
	cwknn(db, x, k, m, n, target, nt, o);

	UNPROTECT(7);
	return(out);
}


void ewknn(double const *db, double const *x, double const *wvec, int k, int m, int n, int const *target, R_len_t nt, double* o){
	int i, j, mm;
	R_len_t t;
	double y, wsum; 
	double *ds = (double*) calloc(n, sizeof(double));
	double minimumDec = 0.0001; // minimum delta to be added
	int *oo = (int*) calloc(n, sizeof(int));	
	
	for(j = 0; j < n; j++){
		ds[j] = 0;
		mm = m;
		wsum = 0;
		for(i = 0; i < m; i++){
			y = db[i + j*m];
			if(ISNA(x[i]) || ISNA(y)) mm--;
			else {
				ds[j] += (wvec[i]*wvec[i]*(x[i] - y)*(x[i] - y) );
				wsum += wvec[i] * wvec[i];
			}
		}
		//ds[j] = mm/sqrt(ds[j] + minimumDec);
		ds[j] = wsum / sqrt(ds[j] + minimumDec );
	}
	order(ds, n, oo);
	for(t = 0; t < nt; t++){
		o[t] = 0;
		wsum = 0;
		for(j = 0; j < k; j++){
			//Rprintf("%f*%f\t", ds[oo[n-1-j]], db[target[t] + oo[n-1-j] * m]);
			o[t] += ( db[target[t] + oo[n-1-j] * m] * ds[oo[n-1-j]] ) ;
			wsum += ds[oo[n-1-j]] ;
		}//Rprintf("\n");
		//Rprintf("o: %f\nwsum: %f\n",o[t],wsum );
		o[t] /= wsum;
	}
	//Rprintf("Done!\n");
	free(ds);
	free(oo);
}


SEXP ewknnR2C(SEXP dbIn, SEXP xIn, SEXP wvecIn, SEXP kIn, SEXP mIn, SEXP nIn, SEXP tgtIn){
	SEXP out;
	double *db, *x, *wvec, *o;
	int m, n, k, *target;
	R_len_t nt;

	nt = LENGTH(tgtIn);

	PROTECT(dbIn = AS_NUMERIC(dbIn));
	PROTECT(xIn = AS_NUMERIC(xIn));
	PROTECT(wvecIn = AS_NUMERIC(wvecIn));
	PROTECT(mIn = AS_NUMERIC(mIn));
	PROTECT(nIn = AS_NUMERIC(nIn));
	PROTECT(kIn = AS_NUMERIC(kIn));
	PROTECT(tgtIn = AS_INTEGER(tgtIn));
	PROTECT(out = NEW_NUMERIC(nt));

	db = NUMERIC_POINTER(dbIn);
	x = NUMERIC_POINTER(xIn);
	wvec = NUMERIC_POINTER(wvecIn);
	m = NUMERIC_POINTER(mIn)[0];
	n = NUMERIC_POINTER(nIn)[0];
	k = NUMERIC_POINTER(kIn)[0];
	target = INTEGER_POINTER(tgtIn);
	o = NUMERIC_POINTER(out);	
	
	ewknn(db, x, wvec,  k, m, n, target, nt, o);

	UNPROTECT(8);
	return(out);
}

void ewknn_predict(double const *x, double const *y, double const *q, double const *wvec, int k, int m, int n, int nq, double* o){
	int i, j, qj, qj_offset,  mm;
	double xx, wsum; 
	double *ds = (double*) calloc(n, sizeof(double));
	double minimumDec = 0.0001; // minimum delta to be added
	int *oo = (int*) calloc(n, sizeof(int));	

	for(qj = 0; qj < nq; qj++){
	  qj_offset = qj * m;
	  for(j = 0; j < n; j++){
		ds[j] = 0;
		mm = m;
		wsum = 0;
		for(i = 0; i < m; i++){
			xx = x[i + j*m];
			if(ISNA(q[i + qj_offset]) || ISNA(xx)) mm--;
			else {
				ds[j] += (wvec[i]*wvec[i]*(xx - q[i+qj_offset])*(xx - q[i+qj_offset]) );
				wsum += wvec[i] * wvec[i];
			}
		}
		//ds[j] = mm/sqrt(ds[j] + minimumDec);
		ds[j] = wsum / sqrt(ds[j] + minimumDec );
	  }
	  order(ds, n, oo);
	  o[qj] = 0;
	  wsum = 0;
	  for(j = 0; j < k; j++){
		//Rprintf("%f*%f\t", ds[oo[n-1-j]], y[oo[n-1-j]]);
		o[qj] += ( y[ oo[n-1-j] ] * ds[oo[n-1-j]] ) ;
		wsum += ds[oo[n-1-j]] ;
	  }//Rprintf("\n");
	  //Rprintf("o: %f\nwsum: %f\n",o[qj],wsum );
	  o[qj] /= wsum;
	  }
	  //Rprintf("Done!\n");
	free(ds);
	free(oo);
}
SEXP ewknnPredictR2C(SEXP xIn, SEXP yIn, SEXP qIn, SEXP wvecIn, SEXP kIn, SEXP mIn, SEXP nIn, SEXP nqIn){
	SEXP out;
	double *x, *y, *q, *wvec, *o;
	int m, n, k, nq;

	PROTECT(xIn = AS_NUMERIC(xIn));
	PROTECT(yIn = AS_NUMERIC(yIn));
	PROTECT(qIn = AS_NUMERIC(qIn));
	PROTECT(wvecIn = AS_NUMERIC(wvecIn));
	PROTECT(mIn = AS_NUMERIC(mIn));
	PROTECT(nIn = AS_NUMERIC(nIn));
	PROTECT(kIn = AS_NUMERIC(kIn));
	PROTECT(nqIn = AS_NUMERIC(nqIn));

	x = NUMERIC_POINTER(xIn);
	y = NUMERIC_POINTER(yIn);
	q = NUMERIC_POINTER(qIn);
	wvec = NUMERIC_POINTER(wvecIn);
	m = NUMERIC_POINTER(mIn)[0];
	n = NUMERIC_POINTER(nIn)[0];
	k = NUMERIC_POINTER(kIn)[0];
	nq = NUMERIC_POINTER(nqIn)[0];
	PROTECT(out = NEW_NUMERIC(nq));
	o = NUMERIC_POINTER(out);
	
	ewknn_predict(x, y, q, wvec,  k, m, n, nq, o);

	UNPROTECT(9);
	return(out);
}
