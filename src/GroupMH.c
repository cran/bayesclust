#include<math.h>

void GroupMH(float *Y, int *cind, int *k, int *p, int *n, 
		float *means, float *variances, int *counts) 
{
  int i, j, m;

/* Compute sub-counts for each cluster */
  for(i=0; i<*n; i++) {
	*(counts + (cind[i] - 1)) += 1;
  }
	
  for(j=0; j<*p; j++) {
    for(i=0; i<*n; i++) {
	*(means + j*(*k) + cind[i] - 1) += *(Y + j*(*n) + i) ;
	*(variances + j*(*k) + cind[i] - 1) += (float) pow(*(Y + j*(*n) + i), 2.0) ;

    }
  }

  for(j=0; j<(*p); j++) {
    for(i=0;i<(*k);i++) { 
      *(means + j*(*k) + i) /=  *(counts + i);
      if(*(counts+i) > 1) { 
      *(variances + j*(*k) + i) =  (1.0/(*(counts+i)  -1.0)) * (*(variances + j*(*k) + i) - 
		      (*(counts + i) * pow(*(means + j*(*k) + i) ,2.0)));
      } 
    }
  }

}

