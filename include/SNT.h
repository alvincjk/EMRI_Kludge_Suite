// NK: Sasaki-Nakamura and Teukolsky solutions

#ifndef _SNT_H
#define _SNT_H
#include "Globals.h"

#define EPSILON_Int 1.e-14  // Epsilon used in all integration.

class SNT {
public:
  SNT(const int emm, const Real rad,
      const Real spin, const Real Omega,
      const Real Lambda, const Real epsilon);
  ~SNT();
  //
  // Leaving almost everything public now for debugging purposes.
  // A lot of this will be "privatized" later.
  //
  Complex Ain, Bin;
  Complex TeukRH, TeukRInf;
  Complex dr_TeukRH, dr_TeukRInf;
  Complex ddr_TeukRH, ddr_TeukRInf;
  Complex SNXH, SNXHerr, SNXInf, SNXInferr;
  Complex dr_SNXH, dr_SNXHerr, dr_SNXInf, dr_SNXInferr;
  Complex ddr_SNXH, ddr_SNXInf;

  //
  // These must be left public for lots of stuff in the
  // Green's function.
  //
  inline Real K(const Real rad)
    {
      return((rad*rad + a*a)*omega - m*a);
    };

  inline Real dr_K(const Real rad)
    {
      return(2.*rad*omega);
    };

  inline Real ddr_K()
    {
      return(2.*omega);
    };

private:
  int m;
  Real r, a, omega, lambda, pmomega;
  Real EPSILON, EPSILON_Ext;
  Complex dlmw;

  Complex c0, c1, c2, c3, c4;

  inline Complex eta(const Real rad)
    {
      const Real x = 1./rad;

      return(c0 + x*(c1 + x*(c2 + x*(c3 + x*c4))));
    };

  inline Complex dr_eta(const Real rad)
    {
      const Real x = 1./rad;

      return(-x*x*(c1 + x*(2.*c2 + x*(3.*c3 + x*4.*c4))));
    };

  //
  // These two functions are given in the simplified form found by
  // Mathematica --- quite a bit simpler than their "canonical" form.
  //
  inline Real GG(const Real rad)
    {
      const Real r2pa2 = rad*rad + a*a;
      
      const Real numer = a*a*(2. - rad) - rad*rad*rad;
      const Real denom = r2pa2*r2pa2;
      
      return(numer/denom);
    };

  inline Real dr_G(const Real rad)
    {
      const Real r2pa2 = rad*rad + a*a;
      
      const Real numer = rad*rad*rad*rad - a*a*(a*a + 8.*rad);
      const Real denom = r2pa2*r2pa2*r2pa2;
      
      return(numer/denom);
    };

  inline Complex V(const Real rad)
    {
      const Real Kay = K(rad);
      
      const Complex V_1 = -Kay*(Kay + 4.*II*(rad - 1.))/
	Kerr::Delta(rad,a);
      const Complex V_2 = 8.*II*omega*rad + lambda;
      
      return(V_1 + V_2);
    };
  
  inline Complex beta(Real const rad)
    {
      const Real Delta = Kerr::Delta(rad,a);
      
      return(2.*Delta*(-II*K(rad) + rad - 1. - 2.*Delta/rad));
    };

  inline Complex dr_beta(const Real rad)
    {
      const Real Delta = Kerr::Delta(rad,a);
      const Real dr_Delta = Kerr::dr_Delta(rad);
      const Real Kay = K(rad);
      const Real dr_Kay = dr_K(rad);
      
      const Complex dr_beta_1 = 2.*(rad - 1. - 2.*Delta/rad -
				    II*Kay)*dr_Delta;
      const Complex dr_beta_2 = 2.*(1. + 2.*Delta/(rad*rad) -
				    2.*dr_Delta/rad - II*dr_Kay)*Delta;
      
      return(dr_beta_1 + dr_beta_2);
    };

  inline Complex ddr_beta(const Real rad)
    {
      const Real x = 1./rad;
      const Real Delta = Kerr::Delta(rad,a);
      const Real dr_Delta = Kerr::dr_Delta(rad);
      const Real ddr_Delta = Kerr::ddr_Delta();
      const Real Kay = K(rad);
      const Real dr_Kay = dr_K(rad);
      const Real ddr_Kay = ddr_K();
      
      const Complex ddr_beta_1 = 4.*dr_Delta*(1. - 2.*x*(dr_Delta -
							 x*Delta) -
					      II*dr_Kay);
      const Complex ddr_beta_2 = 2.*(rad - 1. - 2.*x*Delta -
				     II*Kay)*ddr_Delta;
      const Complex ddr_beta_3 = 2.*Delta*(x*(-2*ddr_Delta +
					      x*(4.*dr_Delta +
						 x*(-4*Delta))) -
					   II*ddr_Kay);
      
      return(ddr_beta_1 + ddr_beta_2 + ddr_beta_3);
    };

  inline Complex alpha(const Real rad)
    {
      const Real Kay = K(rad);
      const Real dr_Kay = dr_K(rad);
      const Real Delta = Kerr::Delta(rad,a);
      const Complex bbeta = beta(rad);
      
      return(-II*Kay*bbeta/(Delta*Delta) + 3.*II*dr_Kay + lambda +
	     6.*Delta/(rad*rad));
    };

  Complex dr_alpha(const Real rad)
    {
      const Real x = 1./rad;
      const Real Kay = K(rad);
      const Real dr_Kay = dr_K(rad);
      const Real ddr_Kay = ddr_K();
      const Real Delta = Kerr::Delta(rad,a);
      const Real dr_Delta = Kerr::dr_Delta(rad);
      const Complex bbeta = beta(rad);
      const Complex dr_bbeta = dr_beta(rad);
      
      const Real dr_alpha_1 = 6.*x*x*(dr_Delta - 2.*x*Delta);
      const Complex dr_alpha_2 = -II*(Kay*dr_bbeta +
				      dr_Kay*bbeta)/(Delta*Delta);
      const Complex dr_alpha_3 = II*(2.*bbeta*Kay*dr_Delta/
				     (Delta*Delta*Delta) + 3.*ddr_Kay);
      
      return(dr_alpha_1 + dr_alpha_2 + dr_alpha_3);
    };
  
  inline Complex U1(const Real rad)
    {
      const Real Delta = Kerr::Delta(rad,a);
      const Real dr_Delta = Kerr::dr_Delta(rad);
      const Complex Vee = V(rad);
      const Complex Kay = K(rad);
      const Complex aalpha = alpha(rad);
      const Complex dr_aalpha = dr_alpha(rad);
      const Complex bbeta = beta(rad);
      const Complex dr_bbeta = dr_beta(rad);
      const Complex ddr_bbeta = ddr_beta(rad);
      const Complex eeta = eta(rad);
      const Complex dr_eeta = dr_eta(rad);
      
      const Complex U1_1 = (2.*Delta*Delta*dr_aalpha -
			    dr_bbeta*dr_Delta)/bbeta;
      const Complex U1_2 = -Delta*dr_eeta*(aalpha*Delta +
					   dr_bbeta)/(bbeta*eeta);
      const Complex U1_3 = Delta*ddr_bbeta/bbeta;
      
      return(Vee + U1_1 + U1_2 + U1_3);
    };

  inline Complex F(const Real rad)
    {
      const Real Delta = Kerr::Delta(rad,a);
      const Real r2pa2 = rad*rad + a*a;
      const Complex eeta = eta(rad);
      const Complex dr_eeta = dr_eta(rad);
      
      const Complex tmp = (dr_eeta/eeta)*(Delta/r2pa2);
      
      return(tmp);
    };

  inline Complex U(const Real rad)
    {
      const Real Delta = Kerr::Delta(rad,a);
      const Real r2pa2 = rad*rad + a*a;
      const Real Gee = GG(rad);
      const Real dr_Gee = dr_G(rad);
      const Complex Yoo1 = U1(rad);
      const Complex Eff = F(rad);
      
      const Complex U_1 = Delta*Yoo1/(r2pa2*r2pa2) + Gee*Gee;
      const Complex U_2 = Delta*dr_Gee/r2pa2 - Eff*Gee;
      
      const Complex tmp = U_1 + U_2;
      
      return(tmp);
    };

  void SNT_derivs(const Real rad, Complex SNX[], Complex dr_SNX[]);

  void XInfInitCond(Complex XInfInit[], const Real rout);
  void XHInitCond(Complex XHInit[], const Real rin);

  void CalcXInf(Complex XInf[], Complex XInferr[], const Real rad);
  void CalcXH(Complex XH[], Complex XHerr[], const Real rad);
  Complex CalcAin();

  Complex CalcR(const Complex X, const Complex dr_X);
  Complex Calcdr_R(const Complex X, const Complex dr_X,
		   const Complex ddr_X);
  Complex Calcddr_R(Complex R, Complex dr_R);

  void IntegrateODE(Complex X[], const Real r1, const Real r2,
		    const Real tol);

  Complex P(const Real rad);
  Complex dr_P(const Real rad);

  Real drsdr(const Real r);

  void odeint(Complex ystart[], const int nvar, const Real x1,
	      const Real x2, const Real eps, const Real h1,
	      const Real hmin);
  void bsstep(Complex y[], Complex dydx[], const int nv,
	      Real *xx, Real htry, const Real eps,
	      Real yscal[], Real *hdid, Real *hnext);
  void mmid(Complex y[], Complex dydx[], int nvar, Real xs,
	    Real htot, int nstep, Complex yout[]);
};
#endif

