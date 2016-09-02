#include <stdio.h>
#include <stdlib.h>
#include "Globals.h"

#define NR_END 1
#define FREE_ARG char*

void Die(const char error_text[])
{
  cerr << error_text << endl;
  exit(0);
}

long *lvector(const long nl, const long nh)
{
  long *v;

  v = (long *)malloc((size_t)((nh - nl + 1 + NR_END)*sizeof(long)));
  if (!v) Die("allocation failure in lvector()");
  return v - nl + NR_END;
}

int *ivector(const long nl, const long nh)
{
  int *v;

  v = (int *)malloc((size_t)((nh - nl + 1 + NR_END)*sizeof(int)));
  if (!v) Die("allocation failure in ivector()");
  return v - nl + NR_END;
}

int **imatrix(const long nrl, const long nrh,
	      const long ncl, const long nch)
{
  long i, nrow = nrh - nrl + 1, ncol= nch - ncl + 1;
  int **m;
  
  m = (int **)malloc((size_t)((nrow + NR_END)*sizeof(int*)));
  if (!m) Die("allocation failure 1 in imatrix()");
  m += NR_END;
  m -= nrl;
  
  m[nrl] = (int *)malloc((size_t)((nrow*ncol + NR_END)*sizeof(int)));
  if (!m[nrl]) Die("allocation failure 2 in imatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;
  
  return m;
}

Real *Realvector(const long nl, const long nh)
{
  Real *v;

  v = (Real *)malloc((size_t)((nh - nl + 1 + NR_END)*sizeof(Real)));
  if (!v) Die("allocation failure in Realvector()");
  return v - nl + NR_END;
}

Real **Realmatrix(const long nrl, const long nrh,
		  const long ncl, const long nch)
{
  long i, nrow = nrh - nrl + 1, ncol= nch - ncl + 1;
  Real **m;
  
  m = (Real **)malloc((size_t)((nrow + NR_END)*sizeof(Real*)));
  if (!m) Die("allocation failure 1 in Realmatrix()");
  m += NR_END;
  m -= nrl;
  
  m[nrl] = (Real *)malloc((size_t)((nrow*ncol + NR_END)*sizeof(Real)));
  if (!m[nrl]) Die("allocation failure 2 in Realmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;
  
  return m;
}

Real ***Real3tensor(const long nrl, const long nrh,
		    const long ncl, const long nch,
		    const long ndl, const long ndh)
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1,
    ndep = ndh - ndl + 1;
  Real ***t;
  
  t = (Real ***)malloc((size_t)((nrow + NR_END)*sizeof(Real**)));
  if (!t) Die("allocation failure 1 in Real3tensor()");
  t += NR_END;
  t -= nrl;
  
  t[nrl] = (Real **)malloc((size_t)((nrow*ncol +
				     NR_END)*sizeof(Real*)));
  if (!t[nrl]) Die("allocation failure 2 in Real3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  t[nrl][ncl] = (Real *)malloc((size_t)((nrow*ncol*ndep +
					 NR_END)*sizeof(Real)));
  if (!t[nrl][ncl]) Die("allocation failure 3 in Real3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  
  for (j = ncl + 1; j <= nch; j++)
    t[nrl][j] = t[nrl][j - 1] + ndep;

  for (i = nrl + 1; i <= nrh; i++) {
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol*ndep;

    for (j = ncl + 1; j <= nch; j++)
      t[i][j] = t[i][j - 1] + ndep;
  }
  
  return t;
}

Real ****Real4tensor(const long nrl, const long nrh,
		     const long ncl, const long nch,
		     const long ndl, const long ndh,
		     const long nhrl, const long nhrh)
{
  long i, j, k, nrow = nrh - nrl + 1, ncol = nch - ncl + 1,
    ndep = ndh - ndl + 1, nhro = nhrh - nhrl + 1;
  Real ****t;
  
  t = (Real ****)malloc((size_t)((nrow + NR_END)*sizeof(Real***)));
  if (!t) Die("allocation failure 1 in Real4tensor()");
  t += NR_END;
  t -= nrl;
  
  t[nrl] = (Real ***)malloc((size_t)((nrow*ncol +
				      NR_END)*sizeof(Real**)));
  if (!t[nrl]) Die("allocation failure 2 in Real4tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  t[nrl][ncl] = (Real **)malloc((size_t)((nrow*ncol*ndep +
					  NR_END)*sizeof(Real*)));
  if (!t[nrl][ncl]) Die("allocation failure 3 in Real4tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  
  t[nrl][ncl][ndl] = (Real *)malloc((size_t)((nrow*ncol*ndep*nhro +
					       NR_END)*sizeof(Real)));
  if (!t[nrl][ncl][ndl]) Die("allocation failure 4 in Real4tensor()");
  t[nrl][ncl][ndl] += NR_END;
  t[nrl][ncl][ndl] -= nhrl;

  for (k = ndl + 1; k <= ndh; k++)
    t[nrl][ncl][k] = t[nrl][ncl][k - 1] + nhro;

  for (j = ncl + 1; j <= nch; j++) {
    t[nrl][j] = t[nrl][j - 1] + ndep;
    t[nrl][j][ndl] = t[nrl][j - 1][ndl] + ndep*nhro;
    for (k = ndl + 1; k <= ndh; k++)
      t[nrl][j][k] = t[nrl][j][k - 1] + nhro;
  }

  for (i = nrl + 1; i <= nrh; i++) {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol*ndep;
    t[i][ncl][ndl] = t[i - 1][ncl][ndl] + ncol*ndep*nhro;
    for (j = ncl + 1; j <= nch; j++) {
      t[i][j] = t[i][j - 1] + ndep;
      t[i][j][ndl] = t[i][j - 1][ndl] + ndep*nhro;
      for (k = ndl + 1; k <= ndh; k++) {
	t[i][ncl][k] = t[i][ncl][k - 1] + nhro;
	t[i][j][k] = t[i][j][k - 1] + nhro;
      }
    }
  }

  return t;
}

Real *****Real5tensor(const long nrl, const long nrh,
		      const long ncl, const long nch,
		      const long ndl, const long ndh,
		      const long nhrl, const long nhrh,
		      const long nhcl, const long nhch)
{
  long i, j, k, l, nrow = nrh - nrl + 1, ncol = nch - ncl + 1,
    ndep = ndh - ndl + 1, nhro = nhrh - nhrl + 1,
    nhco = nhch - nhcl + 1;
  Real *****t;
  
  t = (Real *****)malloc((size_t)((nrow + NR_END)*sizeof(Real****)));
  if (!t) Die("allocation failure 1 in Real5tensor()");
  t += NR_END;
  t -= nrl;
  
  t[nrl] = (Real ****)malloc((size_t)((nrow*ncol +
				       NR_END)*sizeof(Real***)));
  if (!t[nrl]) Die("allocation failure 2 in Real5tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  t[nrl][ncl] = (Real ***)malloc((size_t)((nrow*ncol*ndep +
					   NR_END)*sizeof(Real**)));
  if (!t[nrl][ncl]) Die("allocation failure 3 in Real5tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  
  t[nrl][ncl][ndl] = (Real **)malloc((size_t)((nrow*ncol*ndep*nhro +
					       NR_END)*sizeof(Real*)));
  if (!t[nrl][ncl][ndl]) Die("allocation failure 4 in Real5tensor()");
  t[nrl][ncl][ndl] += NR_END;
  t[nrl][ncl][ndl] -= nhrl;

  t[nrl][ncl][ndl][nhrl] = (Real *)malloc((size_t)((nrow*ncol*ndep*nhro*nhco
						    + NR_END)*sizeof(Real)));
  if (!t[nrl][ncl][ndl][nhrl])
    Die("allocation failure 5 in Real5tensor()");
  t[nrl][ncl][ndl][nhrl] += NR_END;
  t[nrl][ncl][ndl][nhrl] -= nhcl;

  for (l = nhrl + 1; l <= nhrh; l++)
    t[nrl][ncl][ndl][l] = t[nrl][ncl][ndl][l - 1] + nhco;

  for (k = ndl + 1; k <= ndh; k++) {
    t[nrl][ncl][k] = t[nrl][ncl][k - 1] + nhro;
    t[nrl][ncl][k][nhrl] = t[nrl][ncl][k - 1][nhrl]
      + nhro*nhco;
    for (l = nhrl + 1; l <= nhrh; l++)
      t[nrl][ncl][k][l] = t[nrl][ncl][k][l - 1] + nhco;
  }

  for (j = ncl + 1; j <= nch; j++) {
    t[nrl][j] = t[nrl][j - 1] + ndep;
    t[nrl][j][ndl] = t[nrl][j - 1][ndl] + ndep*nhro;
    t[nrl][j][ndl][nhrl] = t[nrl][j - 1][ndl][nhrl] + ndep*nhro*nhco;
    for (k = ndl + 1; k <= ndh; k++) {
      t[nrl][j][k] = t[nrl][j][k - 1] + nhro;
      t[nrl][j][k][nhrl] = t[nrl][j][k - 1][nhrl] + nhro*nhco;
      for (l = nhrl + 1; l <= nhrh; l++) {
 	t[nrl][j][ndl][l] = t[nrl][j][ndl][l - 1] + nhco;
 	t[nrl][j][k][l] = t[nrl][j][k][l - 1] + nhco;
      }
    }
  }

  for (i = nrl + 1; i <= nrh; i++) {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol*ndep;
    t[i][ncl][ndl] = t[i - 1][ncl][ndl] + ncol*ndep*nhro;
    t[i][ncl][ndl][nhrl] = t[i - 1][ncl][ndl][nhrl] + ncol*ndep*nhro*nhco;
    for (j = ncl + 1; j <= nch; j++) {
      t[i][j] = t[i][j - 1] + ndep;
      t[i][j][ndl] = t[i][j - 1][ndl] + ndep*nhro;
      t[i][j][ndl][nhrl] = t[i][j - 1][ndl][nhrl] + ndep*nhro*nhco;
      for (k = ndl + 1; k <= ndh; k++) {
	t[i][j][k] = t[i][j][k - 1] + nhro;
	t[i][ncl][k] = t[i][ncl][k - 1] + nhro;
	t[i][ncl][k][nhrl] = t[i][ncl][k - 1][nhrl] + nhro*nhco;
	t[i][j][k][nhrl] = t[i][j][k - 1][nhrl] + nhro*nhco;
	for (l = nhrl + 1; l <= nhrh; l++) {
	  t[i][ncl][ndl][l] = t[i][ncl][ndl][l - 1] + nhco;
	  t[i][ncl][k][l] = t[i][ncl][k][l - 1] + nhco;
	  t[i][j][ndl][l] = t[i][j][ndl][l - 1] + nhco;
	  t[i][j][k][l] = t[i][j][k][l - 1] + nhco;
	}
      }
    }
  }

  return t;
}

Complex *Complexvector(const long nl, const long nh)
{
  Complex *v;

  v = (Complex *)malloc((size_t)((nh - nl + 1 + NR_END)*sizeof(Complex)));
  if (!v) Die("allocation failure in Complexvector()");
  return v - nl + NR_END;
}

Complex **Complexmatrix(const long nrl, const long nrh,
			const long ncl, const long nch)
{
  long i, nrow = nrh - nrl + 1, ncol= nch - ncl + 1;
  Complex **m;
  
  m = (Complex **)malloc((size_t)((nrow + NR_END)*sizeof(Complex*)));
  if (!m) Die("allocation failure 1 in Complexmatrix()");
  m += NR_END;
  m -= nrl;
  
  m[nrl] = (Complex *)malloc((size_t)((nrow*ncol + NR_END)*sizeof(Complex)));
  if (!m[nrl]) Die("allocation failure 2 in Complexmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;
  
  return m;
}

Complex ***Complex3tensor(const long nrl, const long nrh,
			  const long ncl, const long nch,
			  const long ndl, const long ndh)
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1,
    ndep = ndh - ndl + 1;
  Complex ***t;
  
  t = (Complex ***)malloc((size_t)((nrow + NR_END)*sizeof(Complex**)));
  if (!t) Die("allocation failure 1 in Complex3tensor()");
  t += NR_END;
  t -= nrl;
  
  t[nrl] = (Complex **)malloc((size_t)((nrow*ncol +
					NR_END)*sizeof(Complex*)));
  if (!t[nrl]) Die("allocation failure 2 in Complex3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  t[nrl][ncl] = (Complex *)malloc((size_t)((nrow*ncol*ndep +
					    NR_END)*sizeof(Complex)));
  if (!t[nrl][ncl]) Die("allocation failure 3 in Complex3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  
  for (j = ncl + 1; j <= nch; j++)
    t[nrl][j] = t[nrl][j - 1] + ndep;

  for (i = nrl + 1; i <= nrh; i++) {
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol*ndep;

    for (j = ncl + 1; j <= nch; j++)
      t[i][j] = t[i][j - 1] + ndep;
  }
  
  return t;
}

Complex ****Complex4tensor(const long nrl, const long nrh,
			   const long ncl, const long nch,
			   const long ndl, const long ndh,
			   const long nhrl, const long nhrh)
{
  long i, j, k, nrow = nrh - nrl + 1, ncol = nch - ncl + 1,
    ndep = ndh - ndl + 1, nhro = nhrh - nhrl + 1;
  Complex ****t;
  
  t = (Complex ****)malloc((size_t)((nrow + NR_END)*sizeof(Complex***)));
  if (!t) Die("allocation failure 1 in Complex4tensor()");
  t += NR_END;
  t -= nrl;
  
  t[nrl] = (Complex ***)malloc((size_t)((nrow*ncol +
					 NR_END)*sizeof(Complex**)));
  if (!t[nrl]) Die("allocation failure 2 in Complex4tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  t[nrl][ncl] = (Complex **)malloc((size_t)((nrow*ncol*ndep +
					     NR_END)*sizeof(Complex*)));
  if (!t[nrl][ncl]) Die("allocation failure 3 in Complex4tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  
  t[nrl][ncl][ndl] = (Complex *)malloc((size_t)((nrow*ncol*ndep*nhro +
						 NR_END)*sizeof(Complex)));
  if (!t[nrl][ncl][ndl]) Die("allocation failure 4 in Complex4tensor()");
  t[nrl][ncl][ndl] += NR_END;
  t[nrl][ncl][ndl] -= nhrl;

  for (k = ndl + 1; k <= ndh; k++)
    t[nrl][ncl][k] = t[nrl][ncl][k - 1] + nhro;

  for (j = ncl + 1; j <= nch; j++) {
    t[nrl][j] = t[nrl][j - 1] + ndep;
    t[nrl][j][ndl] = t[nrl][j - 1][ndl] + ndep*nhro;
    for (k = ndl + 1; k <= ndh; k++)
      t[nrl][j][k] = t[nrl][j][k - 1] + nhro;
  }

  for (i = nrl + 1; i <= nrh; i++) {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol*ndep;
    t[i][ncl][ndl] = t[i - 1][ncl][ndl] + ncol*ndep*nhro;
    for (j = ncl + 1; j <= nch; j++) {
      t[i][j] = t[i][j - 1] + ndep;
      t[i][j][ndl] = t[i][j - 1][ndl] + ndep*nhro;
      for (k = ndl + 1; k <= ndh; k++) {
	t[i][ncl][k] = t[i][ncl][k - 1] + nhro;
	t[i][j][k] = t[i][j][k - 1] + nhro;
      }
    }
  }

  return t;
}

Complex *****Complex5tensor(const long nrl, const long nrh,
			    const long ncl, const long nch,
			    const long ndl, const long ndh,
			    const long nhrl, const long nhrh,
			    const long nhcl, const long nhch)
{
  long i, j, k, l, nrow = nrh - nrl + 1, ncol = nch - ncl + 1,
    ndep = ndh - ndl + 1, nhro = nhrh - nhrl + 1,
    nhco = nhch - nhcl + 1;
  Complex *****t;
  
  t = (Complex *****)malloc((size_t)((nrow + NR_END)*sizeof(Complex****)));
  if (!t) Die("allocation failure 1 in Complex5tensor()");
  t += NR_END;
  t -= nrl;
  
  t[nrl] = (Complex ****)malloc((size_t)((nrow*ncol +
					  NR_END)*sizeof(Complex***)));
  if (!t[nrl]) Die("allocation failure 2 in Complex5tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  t[nrl][ncl] = (Complex ***)malloc((size_t)((nrow*ncol*ndep +
					      NR_END)*sizeof(Complex**)));
  if (!t[nrl][ncl]) Die("allocation failure 3 in Complex5tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  
  t[nrl][ncl][ndl] = (Complex **)malloc((size_t)((nrow*ncol*ndep*nhro +
						  NR_END)*sizeof(Complex*)));
  if (!t[nrl][ncl][ndl]) Die("allocation failure 4 in Complex5tensor()");
  t[nrl][ncl][ndl] += NR_END;
  t[nrl][ncl][ndl] -= nhrl;

  t[nrl][ncl][ndl][nhrl] = (Complex *)malloc((size_t)
					     ((nrow*ncol*ndep*nhro*nhco
					       + NR_END)*sizeof(Complex)));
  if (!t[nrl][ncl][ndl][nhrl])
    Die("allocation failure 5 in Complex5tensor()");
  t[nrl][ncl][ndl][nhrl] += NR_END;
  t[nrl][ncl][ndl][nhrl] -= nhcl;

  for (l = nhrl + 1; l <= nhrh; l++)
    t[nrl][ncl][ndl][l] = t[nrl][ncl][ndl][l - 1] + nhco;

  for (k = ndl + 1; k <= ndh; k++) {
    t[nrl][ncl][k] = t[nrl][ncl][k - 1] + nhro;
    t[nrl][ncl][k][nhrl] = t[nrl][ncl][k - 1][nhrl]
      + nhro*nhco;
    for (l = nhrl + 1; l <= nhrh; l++)
      t[nrl][ncl][k][l] = t[nrl][ncl][k][l - 1] + nhco;
  }

  for (j = ncl + 1; j <= nch; j++) {
    t[nrl][j] = t[nrl][j - 1] + ndep;
    t[nrl][j][ndl] = t[nrl][j - 1][ndl] + ndep*nhro;
    t[nrl][j][ndl][nhrl] = t[nrl][j - 1][ndl][nhrl] + ndep*nhro*nhco;
    for (k = ndl + 1; k <= ndh; k++) {
      t[nrl][j][k] = t[nrl][j][k - 1] + nhro;
      t[nrl][j][k][nhrl] = t[nrl][j][k - 1][nhrl] + nhro*nhco;
      for (l = nhrl + 1; l <= nhrh; l++) {
 	t[nrl][j][ndl][l] = t[nrl][j][ndl][l - 1] + nhco;
 	t[nrl][j][k][l] = t[nrl][j][k][l - 1] + nhco;
      }
    }
  }

  for (i = nrl + 1; i <= nrh; i++) {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol*ndep;
    t[i][ncl][ndl] = t[i - 1][ncl][ndl] + ncol*ndep*nhro;
    t[i][ncl][ndl][nhrl] = t[i - 1][ncl][ndl][nhrl] + ncol*ndep*nhro*nhco;
    for (j = ncl + 1; j <= nch; j++) {
      t[i][j] = t[i][j - 1] + ndep;
      t[i][j][ndl] = t[i][j - 1][ndl] + ndep*nhro;
      t[i][j][ndl][nhrl] = t[i][j - 1][ndl][nhrl] + ndep*nhro*nhco;
      for (k = ndl + 1; k <= ndh; k++) {
	t[i][j][k] = t[i][j][k - 1] + nhro;
	t[i][ncl][k] = t[i][ncl][k - 1] + nhro;
	t[i][ncl][k][nhrl] = t[i][ncl][k - 1][nhrl] + nhro*nhco;
	t[i][j][k][nhrl] = t[i][j][k - 1][nhrl] + nhro*nhco;
	for (l = nhrl + 1; l <= nhrh; l++) {
	  t[i][ncl][ndl][l] = t[i][ncl][ndl][l - 1] + nhco;
	  t[i][ncl][k][l] = t[i][ncl][k][l - 1] + nhco;
	  t[i][j][ndl][l] = t[i][j][ndl][l - 1] + nhco;
	  t[i][j][k][l] = t[i][j][k][l - 1] + nhco;
	}
      }
    }
  }

  return t;
}

FILE ****FILE3tensor(const long nrl, const long nrh,
			const long ncl, const long nch,
			const long ndl, const long ndh)
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1,
    ndep = ndh - ndl + 1;
  FILE ****t;
  
  t = (FILE ****)malloc((size_t)((nrow + NR_END)*sizeof(FILE***)));
  if (!t) Die("allocation failure 1 in FILE3tensor()");
  t += NR_END;
  t -= nrl;
  
  t[nrl] = (FILE ***)malloc((size_t)((nrow*ncol +
				      NR_END)*sizeof(FILE**)));
  if (!t[nrl]) Die("allocation failure 2 in FILE3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  t[nrl][ncl] = (FILE **)malloc((size_t)((nrow*ncol*ndep +
					  NR_END)*sizeof(FILE*)));
  if (!t[nrl][ncl]) Die("allocation failure 3 in FILE3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  
  for (j = ncl + 1; j <= nch; j++)
    t[nrl][j] = t[nrl][j - 1] + ndep;

  for (i = nrl + 1; i <= nrh; i++) {
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol*ndep;

    for (j = ncl + 1; j <= nch; j++)
      t[i][j] = t[i][j - 1] + ndep;
  }
  
  return t;
}

char ****char4tensor(const long nrl, const long nrh,
		     const long ncl, const long nch,
		     const long ndl, const long ndh,
		     const long nhrl, const long nhrh)
{
  long i, j, k, nrow = nrh - nrl + 1, ncol = nch - ncl + 1,
    ndep = ndh - ndl + 1, nhro = nhrh - nhrl + 1;
  char ****t;
  
  t = (char ****)malloc((size_t)((nrow + NR_END)*sizeof(char***)));
  if (!t) Die("allocation failure 1 in char4tensor()");
  t += NR_END;
  t -= nrl;
  
  t[nrl] = (char ***)malloc((size_t)((nrow*ncol +
				      NR_END)*sizeof(char**)));
  if (!t[nrl]) Die("allocation failure 2 in char4tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  t[nrl][ncl] = (char **)malloc((size_t)((nrow*ncol*ndep +
					  NR_END)*sizeof(char*)));
  if (!t[nrl][ncl]) Die("allocation failure 3 in char4tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  
  t[nrl][ncl][ndl] = (char *)malloc((size_t)((nrow*ncol*ndep*nhro +
					      NR_END)*sizeof(char)));
  if (!t[nrl][ncl][ndl]) Die("allocation failure 4 in char4tensor()");
  t[nrl][ncl][ndl] += NR_END;
  t[nrl][ncl][ndl] -= nhrl;
  
  for (k = ndl + 1; k <= ndh; k++)
    t[nrl][ncl][k] = t[nrl][ncl][k - 1] + nhro;
  
  for (j = ncl + 1; j <= nch; j++) {
    t[nrl][j] = t[nrl][j - 1] + ndep;
    t[nrl][j][ndl] = t[nrl][j - 1][ndl] + ndep*nhro;
    for (k = ndl + 1; k <= ndh; k++)
      t[nrl][j][k] = t[nrl][j][k - 1] + nhro;
  }
  
  for (i = nrl + 1; i <= nrh; i++) {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol*ndep;
    t[i][ncl][ndl] = t[i - 1][ncl][ndl] + ncol*ndep*nhro;
    for (j = ncl + 1; j <= nch; j++) {
      t[i][j] = t[i][j - 1] + ndep;
      t[i][j][ndl] = t[i][j - 1][ndl] + ndep*nhro;
      for (k = ndl + 1; k <= ndh; k++) {
	t[i][ncl][k] = t[i][ncl][k - 1] + nhro;
	t[i][j][k] = t[i][j][k - 1] + nhro;
      }
    }
  }

  return t;
}

void free_lvector(long *v, const long nl, const long nh)
{
  if (nh < nl) Die("nh < nl in free_lvector!");
  free((FREE_ARG)(v + nl - NR_END));
}

void free_ivector(int *v, const long nl, const long nh)
{
  if (nh < nl) Die("nh < nl in free_ivector!");
  free((FREE_ARG)(v + nl - NR_END));
}

void free_imatrix(int **m,
		  const long nrl, const long nrh,
		  const long ncl, const long nch)
{
  if (nrh < nrl || nch < ncl)
    Die("Arguments out of whack in free_Realmatrix!");
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

void free_Realvector(Real *v, const long nl, const long nh)
{
  if (nh < nl) Die("nh < nl in free_Realvector!");
  free((FREE_ARG)(v + nl - NR_END));
}

void free_Realmatrix(Real **m,
		     const long nrl, const long nrh,
		     const long ncl, const long nch)
{
  if (nrh < nrl || nch < ncl)
    Die("Arguments out of whack in free_Realmatrix!");
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

void free_Real3tensor(Real ***t,
		      long nrl, long nrh,
		      long ncl, long nch,
		      long ndl, long ndh)
{
  if (nrh < nrl || nch < ncl || ndh < ndl)
    Die("Arguments out of whack in free_Real3tensor!");
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
}

void free_Real4tensor(Real ****t,
		      long nrl, long nrh,
		      long ncl, long nch,
		      long ndl, long ndh,
		      long nhrl, long nhrh)
{
  if (nrh < nrl || nch < ncl || ndh < ndl || nhrh < nhrl)
    Die("Arguments out of whack in free_Real4tensor!");
  free((FREE_ARG)(t[nrl][ncl][ndl] + nhrl - NR_END));
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
}

void free_Real5tensor(Real *****t,
		      long nrl, long nrh,
		      long ncl, long nch,
		      long ndl, long ndh,
		      long nhrl, long nhrh,
		      long nhcl, long nhch)
{
  if (nrh < nrl || nch < ncl || ndh < ndl || nhrh < nhrl || nhch < nhcl)
    Die("Arguments out of whack in free_Real5tensor!");
  free((FREE_ARG)(t[nrl][ncl][ndl][nhrl] + nhcl - NR_END));
  free((FREE_ARG)(t[nrl][ncl][ndl] + nhrl - NR_END));
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
}

void free_Complexvector(Complex *v, const long nl, const long nh)
{
  if (nh < nl) Die("nh < nl in free_Complexvector!");
  free((FREE_ARG)(v + nl - NR_END));
}

void free_Complexmatrix(Complex **m,
			const long nrl, const long nrh,
			const long ncl, const long nch)
{
  if (nrh < nrl || nch < ncl)
    Die("Arguments out of whack in free_Complexmatrix!");
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

void free_Complex3tensor(Complex ***t,
			 long nrl, long nrh,
			 long ncl, long nch,
			 long ndl, long ndh)
{
  if (nrh < nrl || nch < ncl || ndh < ndl)
    Die("Arguments out of whack in free_Complex3tensor!");
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
}

void free_Complex4tensor(Complex ****t,
			 long nrl, long nrh,
			 long ncl, long nch,
			 long ndl, long ndh,
			 long nhrl, long nhrh)
{
  if (nrh < nrl || nch < ncl || ndh < ndl || nhrh < nhrl)
    Die("Arguments out of whack in free_Complex4tensor!");
  free((FREE_ARG)(t[nrl][ncl][ndl] + nhrl - NR_END));
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
}

void free_Complex5tensor(Complex *****t,
			 long nrl, long nrh,
			 long ncl, long nch,
			 long ndl, long ndh,
			 long nhrl, long nhrh,
			 long nhcl, long nhch)
{
  if (nrh < nrl || nch < ncl || ndh < ndl || nhrh < nhrl || nhch < nhcl)
    Die("Arguments out of whack in free_Complex5tensor!");
  free((FREE_ARG)(t[nrl][ncl][ndl][nhrl] + nhcl - NR_END));
  free((FREE_ARG)(t[nrl][ncl][ndl] + nhrl - NR_END));
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
}

void free_FILE3tensor(FILE ****t,
		      long nrl, long nrh,
		      long ncl, long nch,
		      long ndl, long ndh)
{
  if (nrh < nrl || nch < ncl || ndh < ndl)
    Die("Arguments out of whack in free_FILE3tensor!");
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
}

void free_char4tensor(char ****t,
		      long nrl, long nrh,
		      long ncl, long nch,
		      long ndl, long ndh,
		      long nhrl, long nhrh)
{
  if (nrh < nrl || nch < ncl || ndh < ndl || nhrh < nhrl)
    Die("Arguments out of whack in free_char4tensor!");
  free((FREE_ARG)(t[nrl][ncl][ndl] + nhrl - NR_END));
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
}
