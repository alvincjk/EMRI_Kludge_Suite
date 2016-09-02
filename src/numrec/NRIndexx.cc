#include "Globals.h"
#include "NRUtil.h"

#define M 7
#define NSTACK 50

void indexx(const unsigned long n, Real arr[], unsigned long indx[])
{
  unsigned long i, indxt, ir=n, j, k, l=1;
  long *istack;
  int jstack=0;
  Real a;
  
  istack = lvector(1, NSTACK);
  for (j=1;j<=n;j++) indx[j]=j;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	indxt=indx[j];
	a=arr[indxt];
	for (i=j-1;i>=1;i--) {
	  if (arr[indx[i]] <= a) break;
	  indx[i+1]=indx[i];
	}
	indx[i+1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(indx[k], indx[l+1]);
      if (arr[indx[l+1]] > arr[indx[ir]]) {
	SWAP(indx[l+1], indx[ir]);
      }
      if (arr[indx[l]] > arr[indx[ir]]) {
	SWAP(indx[l], indx[ir]);
      }
      if (arr[indx[l+1]] > arr[indx[l]]) {
	SWAP(indx[l+1], indx[l]);
      }
      i=l+1;
      j=ir;
      indxt=indx[l];
      a=arr[indxt];
      for (;;) {
	do i++; while (arr[indx[i]] < a);
	do j--; while (arr[indx[j]] > a);
	if (j < i) break;
	SWAP(indx[i], indx[j]);
      }
      indx[l]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack > NSTACK) Die("NSTACK too small in indexx.");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
  free_lvector(istack, 1, NSTACK);
}
#undef M
#undef NSTACK
