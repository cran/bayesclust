#include <math.h>

void nullsum(float *Y, int *cind, int* len, int *k, float *sum) 
{
  int i;

/* Compute sums for each cluster */
  for(i=0; i<*len; i++)
	*(sum + cind[i] - 1)+=Y[i];
}
