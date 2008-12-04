#include <math.h>

void nullsum(float *Y, int *cind, int* len, int *k, float *sum) 
{
  int i;

/* Compute sums for each cluster */
  for(i=0; i<*len; i++) {
    switch(cind[i]) {
      case 1: 
		sum[0]+=Y[i];
		break;
      case 2: 
		sum[1]+=Y[i];
		break;
      case 3: 
		sum[2]+=Y[i];
		break;
      case 4: 
		sum[3]+=Y[i];
    }

  }
/*  for(i=0; i<*k; i++) sum[i] = log(sum[i]);*/
}
