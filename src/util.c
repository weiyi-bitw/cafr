#include "util.h"

void QuickSortR(double const * arr, int *idx, int left, int right) {
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
            QuickSortR(arr, idx, left, j);
      if (i < right)
            QuickSortR(arr, idx, i, right);
}

void Order(double const *x, int n, int *o){
  int i;
  for(i = 0; i < n; i++){
    o[i] = i;
  }
  QuickSortR(x, o, 0, n-1);
}

void QuickSortRL(double const * arr, uint32_t *idx, 
                 uint32_t left, uint32_t right) {
  uint32_t i = left, j = right;
  uint32_t tmp;
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
    QuickSortRL(arr, idx, left, j);
  if (i < right)
    QuickSortRL(arr, idx, i, right);
}

void OrderL(double const *x, uint32_t n, uint32_t *o){
  uint32_t i;
  for(i = 0; i < n; i++){
    o[i] = i;
  }
  QuickSortRL(x, o, 0, n-1);
}


int RowIndex(long i, int m){
  double m_d = m-1;
  double row = (-2*m_d - 1 + sqrt((4*m_d*(m_d+1) - 8*(double)i - 7))) / -2;
  if(row == (double)(int) row) row -= 1;
  return (int)row;
}

int ColIndex(long i, int m, int row){
  return i - (m-1)*row + row*(row+1)/2 + 1;
}

long TriangularIndex(int i, int j, int m){
  long idx = ( -3*i - (long)i*i ) / 2 + j + m*i-1;
  return idx;
}
