#include <stdio.h>
#include <stdlib.h>
//#include <R.h>
//#include <Rdefines.h>


// filtering out the probes that has low correlation with others
void probe_filter(double const *x, int m, int n, int *grp, double *xsum, int mg){ 

}


void probe_summarization(double const *x, int m, int n, int *grp, double *xs, int mg){
	int i, j;
	int *npbs = (int*) calloc(mg, sizeof(int));

	







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

int main(){
	double x[16] = { 1, 1, 2, 3, 2, 2, 2, 4, 3, 3, 2, 5, 4, 4, 2, 6};
	int grp[4] = {0, 0, 1, -1};
	int m=4, n=4, mg=2, i, j;

	double *xs = (double*) calloc(16, sizeof(double));
	for(i = 0; i < mg; i++) for(j = 0; j < n; j++) xs[i + j*m] = 0;

	printf("Before summarization: \n");
	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			printf("%f\t",x[i + j*m]);
		}printf("\n");
	}
	probe_summarization(x, m, n, grp, xs, mg);
	printf("After summarization: \n");
	for(i = 0; i < mg; i++){
		for(j = 0; j < n; j++){
			printf("%f\t",xs[i + j*mg]);
		}printf("\n");
	}
	


	free(xs);
	return 0;
}


