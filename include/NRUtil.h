// NK: Numerical Recipes

#ifndef _NRUTIL_H
#define _NRUTIL_H

#include "Globals.h"
#include <stdio.h>

// A reasonable C++ square function
template<class NUMBER>
inline NUMBER Pow2(const NUMBER a) { return a*a; };

// Crutches for Numerical Recipes legacy code
#define SQR(a) Pow2(a)
#define DSQR(a) Pow2(a)

// Reasonable maximum
template<class NUMBER>
inline NUMBER Max(const NUMBER a, const NUMBER b) {
  return (a > b ? a : b);
};

#define DMAX(a,b) Max(a,b)
#define FMAX(a,b) Max(a,b)
#define LMAX(a,b) Max(a,b)
#define IMAX(a,b) Max(a,b)

// Reasonable minimum
template<class NUMBER>
inline NUMBER Min(const NUMBER a, const NUMBER b) {
  return (a < b ? a : b);
};

#define DMIN(a,b) Min(a,b)
#define FMIN(a,b) Min(a,b)
#define LMIN(a,b) Min(a,b)
#define IMIN(a,b) Min(a,b)

// Reasonable sign ... chokes on integers.
template <class NUMBER>
inline NUMBER SIGN(const NUMBER a, const NUMBER b) {
  return(b >= 0.0 ? fabs(a) : -fabs(a));
};

// Reasonable swap, used in indexx.
template <class NUMBER>
inline void SWAP(NUMBER & a, NUMBER & b) {
  const NUMBER tmp = a;
  a = b;
  b = tmp;
};

// overload abs for type int ... for some reason, g++ doesn't
// do this.  Stupid.
inline int abs(const int a) {
  return(a >= 0 ? a : -a);
}

// Prototypes of stuff in NRUtil.cc
void Die(const char error_text[]);
long *lvector(const long nl, const long nh);
int *ivector(const long nl, const long nh);
int **imatrix(const long nrl, const long nrh,
	      const long ncl, const long nch);
//
void dfour1(Real *data, unsigned long nn, int isign);
Real *Realvector(const long nl, const long nh);
Real **Realmatrix(const long nrl, const long nrh,
		  const long ncl, const long nch);
Real ***Real3tensor(const long nrl, const long nrh,
		    const long ncl, const long nch,
		    const long ndl, const long ndh);
Real ****Real4tensor(const long nrl, const long nrh,
		     const long ncl, const long nch,
		     const long ndl, const long ndh,
		     const long nhrl, const long nhrh);
Real *****Real5tensor(const long nrl, const long nrh,
		      const long ncl, const long nch,
		      const long ndl, const long ndh,
		      const long nhrl, const long nhrh,
		      const long nhcl, const long nhch);
//
Complex *Complexvector(const long nl, const long nh);
Complex **Complexmatrix(const long nrl, const long nrh,
			const long ncl, const long nch);
Complex ***Complex3tensor(const long nrl, const long nrh,
			  const long ncl, const long nch,
			  const long ndl, const long ndh);
Complex ****Complex4tensor(const long nrl, const long nrh,
			   const long ncl, const long nch,
			   const long ndl, const long ndh,
			   const long nhrl, const long nhrh);
Complex *****Complex5tensor(const long nrl, const long nrh,
			    const long ncl, const long nch,
			    const long ndl, const long ndh,
			    const long nhrl, const long nhrh,
			    const long nhcl, const long nhch);
//
FILE ****FILE3tensor(const long nrl, const long nrh,
		     const long ncl, const long nch,
		     const long ndl, const long ndh);
char ****char4tensor(const long nrl, const long nrh,
		     const long ncl, const long nch,
		     const long ndl, const long ndh,
		     const long nhrl, const long nhrh);
//
void free_lvector(long *v,
		  const long nl, const long nh);
void free_ivector(int *v,
		  const long nl, const long nh);
void free_imatrix(int **v,
		  const long nrl, const long nrh,
		  const long ncl, const long nch);
//
void free_Realvector(Real *v,
		     const long nl, const long nh);
void free_Realmatrix(Real **v,
		     const long nrl, const long nrh,
		     const long ncl, const long nch);
void free_Real3tensor(Real ***t,
		      const long nrl, const long nrh,
		      const long ncl, const long nch,
		      const long ndl, const long ndh);
void free_Real4tensor(Real ****t,
		      const long nrl, const long nrh,
		      const long ncl, const long nch,
		      const long ndl, const long ndh,
		      const long nhrl, const long nhrh);
void free_Real5tensor(Real *****t,
		      const long nrl, const long nrh,
		      const long ncl, const long nch,
		      const long ndl, const long ndh,
		      const long nhrl, const long nhrh,
		      const long nhcl, const long nhch);
//
void free_Complexvector(Complex *v,
			const long nl, const long nh);
void free_Complexmatrix(Complex **v,
			const long nrl, const long nrh,
			const long ncl, const long nch);
void free_Complex3tensor(Complex ***t,
			 const long nrl, const long nrh,
			 const long ncl, const long nch,
			 const long ndl, const long ndh);
void free_Complex4tensor(Complex ****t,
			 const long nrl, const long nrh,
			 const long ncl, const long nch,
			 const long ndl, const long ndh,
			 const long nhrl, const long nhrh);
void free_Complex5tensor(Complex *****t,
			 const long nrl, const long nrh,
			 const long ncl, const long nch,
			 const long ndl, const long ndh,
			 const long nhrl, const long nhrh,
			 const long nhcl, const long nhch);
//
void free_FILE3tensor(FILE ****t,
		      const long nrl, const long nrh,
		      const long ncl, const long nch,
		      const long ndl, const long ndh);
void free_char4tensor(char ****t,
		      const long nrl, const long nrh,
		      const long ncl, const long nch,
		      const long ndl, const long ndh,
		      const long nhrl, const long nhrh);
#endif /* _NRUTIL_H */
