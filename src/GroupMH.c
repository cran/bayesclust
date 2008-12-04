#include<math.h>

void GroupMH(float *Y, int *cind, int *k, int *p, int *n, 
		float *means, float *variances, int *counts) 
{
  int i, j;

/* Compute sub-counts for each cluster */
  for(i=0; i<*n; i++) {
    switch(cind[i]) {
      case 1: 
		counts[0]+=1;
		break;
      case 2: 
		counts[1]+=1;
		break;
      case 3: 
		counts[2]+=1;
		break;
      case 4: 
		counts[3]+=1;
    }

  }
	
 
  for(j=0; j<*p; j++) {
    for(i=0; i<*n; i++) {
      switch(cind[i]) {
        case 1: 
		*(means + j*(*k)) += *(Y + j*(*n) + i) ;
		*(variances + j*(*k)) += (float) pow(*(Y + j*(*n) + i), 2.0) ;
		break;
        case 2: 
		*(means + j*(*k) + 1) += *(Y + j*(*n) + i) ;
		*(variances + j*(*k) + 1) +=(float)  pow(*(Y + j*(*n) + i), 2.0) ;
		break;
        case 3: 
		*(means + j*(*k) + 2) += *(Y + j*(*n) + i) ;
		*(variances + j*(*k) + 2) +=(float)  pow(*(Y + j*(*n) + i), 2.0) ;
		break;
        case 4: 
		*(means + j*(*k) + 3) += *(Y + j*(*n) + i) ;
		*(variances + j*(*k) + 3) += (float) pow(*(Y + j*(*n) + i), 2.0) ;
      }

    }
  }

  for(j=0; j<(*p); j++) {
    for(i=0;i<(*k);i++) { 
      *(means + j*(*k) + i) /=  *(counts + i);
      if(*(counts+i) > 1) { 
      *(variances + j*(*k) + i) =  (1.0/(*(counts+i)  -1.0)) * (*(variances + j*(*k) + i) - (*(counts + i) * pow(*(means + j*(*k) + i) ,2.0)));
      } 
    }
  }

}

