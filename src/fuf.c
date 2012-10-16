
void quickSortR(double const * arr, int *idx, int left, int right) {
      int i = left, j = right;
      int tmp;
      double pivot = arr[ idx[(left + right)/2] ];

      /* partition */
      while (i <= j) {
            while (arr[idx[i]] < pivot)
                  i++;
            while (arr[idx[j]] > pivot)
                  j--;
            if (i <= j) {
                  tmp = idx[i];
                  idx[i] = idx[j];
                  idx[j] = tmp;
                  i++;
                  j--;
            }
      };

      if (left < j)
            quickSortR(arr, idx, left, j);
      if (i < right)
            quickSortR(arr, idx, i, right);
}

void order(double const *x, int n, int *o){
  int i;
  for(i = 0; i < n; i++){
    o[i] = i;
  }
  quickSortR(x, o, 0, n-1);
}
