#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_trig.h>

#include "KSParMap.h"
#include "AKTDI.h"

using namespace std;

void PNevolutionAKTDI(complex<double> **Xp, complex<double> **Xc, complex<double> **Yp, complex<double> **Yc, complex<double> **Zp, complex<double> **Zc,\
 double *Ampvec, double *tvec, double *evec, double *gimvec, double *Phivec, double *alpvec, double *nuvec, double *gimdotvec, double *alpdotvec,\
  double *gimddotvec, double *Phiddotvec, double *alpddotvec, double dt_ph, int Nph, double *par, double nu0, double zeta, int nmodes, double arm, int *i_plunge){

	double mu=par[0];
	double M=par[1];
  	double S=par[2];
  	double e0=par[3];
  	double clam=cos(par[4]);
  	double slam=sin(par[4]);
  	double gim0=par[5];
  	double Phi0=par[6];
  	double qS=par[7];
  	double phiS=par[8];
  	double qK=par[9];
  	double phiK=par[10];
  	double alp0=par[11];

   	double k[3];
   	double uhat[3];
   	double vhat[3];
   	double u[3];
   	double v[3];
   	double kn[3];
   	double kq[3];

   	double ctS=cos(qS);
   	double stS=sin(qS);
   	double cpS=cos(phiS);
   	double spS=sin(phiS);
	double up = (ctS*sin(qK)*cos(phiS-phiK) - cos(qK)*stS);
   	double dw = (sin(qK)*sin(phiS-phiK));
   	double psi;
   	if (dw != 0.0) {
   		psi = -atan2(up, dw);
   	}
   	else {
    	psi = 0.5*M_PI;
   	}
   	double c2psi=cos(2.*psi);
   	double s2psi=sin(2.*psi);

	// note that k = -n, where k is propagation vector and n is sky location
   	k[0] = -stS*cpS;
   	k[1] = -stS*spS;
   	k[2] = -ctS;

	uhat[0] = ctS*cpS;
   	uhat[1] = ctS*spS;
   	uhat[2] = -stS;

   	vhat[0] = spS;
   	vhat[1] = -cpS;
   	vhat[2] = 0.0;
   	double nU, nV;

   	double *R;
   	double **q;
   	double **n;
   	R = (double*)malloc(3*sizeof(double));
   	q = (double**)malloc(3*sizeof(double*));
   	n = (double**)malloc(3*sizeof(double*));
   	for (int i=0; i<3; i++){
       	q[i] = (double*)malloc(3*sizeof(double));
       	n[i] = (double*)malloc(3*sizeof(double));
   	}

   	double e = e0;
   	double e2 = e*e;
   	double nu = nu0;
   	double phi = Phi0;
   	double gam = gim0;
   	double alp = alp0;
   	double Y, Z;
   	double edotm, nudotm, phidotm, alpdotm, gamdotm;
   	double edot, nudot, phidot, alpdot, gamdot;
   	double dalpdnu, dalpde, dgamdnu, dgamde, alpddot, gamddot, phiddot;
   	double de, dnu, dphi, dgam, dalp, rhs;
   	double T=0.0;
   	int ind;
   	double ampfct = zeta;

   	double nn, mm, om;
   	complex<double> chi0;
   	complex<double> chi1;
   	complex<double> chi2;
   	complex<double> img(0.0, 1.0);
   	double x, x2;

   	for (int i=0; i<Nph; i++){

          Y=1./(1.-e2);
          Z=pow(2.*M_PI*M*nu,1./3.);
          if (i != 0){
             edotm=edot;
             nudotm=nudot;
             phidotm=phidot;
             alpdotm=alpdot;
             gamdotm=gamdot;
          }
          edot = e*mu/M/M*(-1./15.*pow(Y,3.5)*pow(Z,8.)*((304.+121.*e2)/Y+Z*Z*(70648.-231960.*e2-56101.*e2*e2)/56.)+\
                   S*clam*pow(Z,11.)*pow(Y,4.)*(8184.+10064.*e2+789.*e2*e2)/30.);
          nudot = 96./(10.*M_PI)*mu/pow(M,3)*(pow(Z,11.)*pow(Y,4.5)*((96.+292.*e2+37.*e2*e2)/Y/96.\
       			       +Z*Z*(20368.-61464.*e2-163170.*e2*e2-13147.*e2*e2*e2)/5376.) - \
       			      pow(Z,14.)*pow(Y,5.)*S*clam*(1168.+9688.*e2+6286.*e2*e2 +195.*e2*e2*e2)/192.);
          phidot = 2.*M_PI*nu;
          alpdot = 8.*M_PI*M_PI*nu*nu*S*M*pow(Y,1.5);
          gamdot = 6.*M_PI*nu*Z*Z*Y*(1.+.25*Z*Z*Y*(26.-15.*e2))-3.*clam*alpdot;

          dalpdnu=16.*S*M_PI*nu*M*pow(Y,1.5);
          dalpde=12.*S*M_PI*nu*nu*M*sqrt(Y)*2.*e*Y*Y;
          dgamdnu=6.*Z*Z*Y*(1.+.25*Z*Z*Y*(26.-15.*e2))-3.*clam*dalpdnu+12.*Z*Y*(1.+.5*Z*Z*Y*(25.-15.*e2))*Z/(3.);
          dgamde=(6.*nu*Z*Z*(1.+.5*Z*Z*Y*(26.-15.*e2)))*2.*e*Y*Y-45.*nu*Z*Z*Z*Z*Y*Y*e-3.*clam*dalpde;

          alpddot=M_PI*(dalpdnu*nudot+dalpde*edot);
          gamddot=M_PI*(dgamdnu*nudot+dgamde*edot);
          phiddot=2.*M_PI*nudot;

          if (i == 0) {
             edotm=edot;
             nudotm=nudot;
             alpdotm=alpdot;
             gamdotm=gamdot;
             phidotm=phidot;
          }

          de=(1.5*edot-.5*edotm)*dt_ph;
          dnu=(1.5*nudot-.5*nudotm)*dt_ph;
          dphi=(1.5*phidot-.5*phidotm)*dt_ph;
          dgam=(1.5*gamdot-.5*gamdotm)*dt_ph;
          dalp=(1.5*alpdot-.5*alpdotm)*dt_ph;

          Phivec[i]=phi;
          alpvec[i]=alp;
          gimvec[i]=gam;
          evec[i]=e;
          nuvec[i]=nu;
          gimdotvec[i]=gamdot;
          alpdotvec[i]=alpdot;
          Ampvec[i]=ampfct*pow(2.*M_PI*M*nu, 2./3.);
          Phiddotvec[i]=phiddot;
          gimddotvec[i]=gamddot;
          alpddotvec[i]=alpddot;
          tvec[i] = T;
          e+=de;
          phi+=dphi;
          gam+=dgam;
          alp+=dalp;
          e2=e*e;
          nu+=dnu;
          T+=dt_ph;

          rhs = pow( (1.0-e2)/(6.0+2.0*e), 1.5 )/(2.*M_PI*M);
          if(rhs - nu <= 0.0){
            *i_plunge = i;
       	    break;
          }

          EccentricLISAMotion(0.0, 0.0, T, R, q, n);
          for(int j=0; j<3; j++){
       	      kn[j] = 0.0;
               kq[j] = 0.0;
       		   nU = 0.0;
       		   nV = 0.0;
       		   for(int ii=0; ii<3; ii++){
       		   	 kn[j] += k[ii]*n[j][ii];
       		 	 kq[j] += k[ii]*q[j][ii];
       			 nU += uhat[ii]*n[j][ii];
       			 nV += vhat[ii]*n[j][ii];
       		   }
       		   u[j] = 0.5*(nU*nU - nV*nV);
       		   v[j] = nU*nV;
       	  }

          ind = 0;
          for (int ii=0; ii<nmodes; ii++){ // nn harmonic
              nn = (double)ii+1.;
              for (int jj=0; jj<5; jj++){
                 mm = (double)jj-2.;
                 om = nn*2.*M_PI*nu + 2.*gamdot + mm*alpdot;
                 x = om*arm;
                 x2 = 0.5*x;
                 chi1 = -x*sin(x)*( gsl_sf_sinc(x2*(1.-kn[1]))*exp(-img*x) \
                                    + gsl_sf_sinc(x2*(1.+kn[1])) )*exp(-img*x2*(3.0 + kq[0] + kq[2]));
                 chi2 = x*sin(x)*( gsl_sf_sinc(x2*(1.-kn[2])) + exp(-img*x)*\
                 							gsl_sf_sinc(x2*(1.+kn[2])) )*exp(-img*x2*(3.0 + kq[1] + kq[0]));

                 Xp[ind][i] = (u[1]*c2psi - v[1]*s2psi)*chi1 + (u[2]*c2psi - v[2]*s2psi)*chi2;
                 Xc[ind][i] = (v[1]*c2psi + u[1]*s2psi)*chi1 + (v[2]*c2psi + u[2]*s2psi)*chi2;


                 chi2 = -x*sin(x)*( gsl_sf_sinc(x2*(1.-kn[2]))*exp(-img*x) \
                                      + gsl_sf_sinc(x2*(1.+kn[2])) )*exp(-img*x2*(3.0 + kq[1] + kq[0]));
                 chi0 = x*sin(x)*( gsl_sf_sinc(x2*(1.-kn[0])) + exp(-img*x)*\
                   							gsl_sf_sinc(x2*(1.+kn[0])) )*exp(-img*x2*(3.0 + kq[2] + kq[1]));

                 Yp[ind][i] = (u[2]*c2psi - v[2]*s2psi)*chi2 + (u[0]*c2psi - v[0]*s2psi)*chi0;
                 Yc[ind][i] = (v[2]*c2psi + u[2]*s2psi)*chi2 + (v[0]*c2psi + u[0]*s2psi)*chi0;


                 chi0 = -x*sin(x)*( gsl_sf_sinc(x2*(1.-kn[0]))*exp(-img*x) \
                                        + gsl_sf_sinc(x2*(1.+kn[0])) )*exp(-img*x2*(3.0 + kq[2] + kq[1]));
                 chi1 = x*sin(x)*( gsl_sf_sinc(x2*(1.-kn[1])) + exp(-img*x)*\
                     							gsl_sf_sinc(x2*(1.+kn[1])) )*exp(-img*x2*(3.0 + kq[0] + kq[2]));

                 Zp[ind][i] = (u[0]*c2psi - v[0]*s2psi)*chi0 + (u[1]*c2psi - v[1]*s2psi)*chi1;
                 Zc[ind][i] = (v[0]*c2psi + u[0]*s2psi)*chi0 + (v[1]*c2psi + u[1]*s2psi)*chi1;

                 ind++;

              }
          }

	 if(i==Nph-1) *i_plunge=i+1;

     }

  for (int i=0; i<3; i++){
    free(q[i]);
    free(n[i]);
  }
  free(R);
  free(q);
  free(n);
  return;

}


void waveformAKTDI(double *Xf_r, double *Xf_im, double *Yf_r, double *Yf_im, double *Zf_r, double *Zf_im, double timestep, int vlength, int nmodes, double zeta, double *par, double nu0){

	double mu=par[0];
	double M=par[1];
  	double S=par[2];
  	double e0=par[3];
  	double clam=cos(par[4]);
  	double slam=sin(par[4]);
  	double gim0=par[5];
  	double Phi0=par[6];
  	double qS=par[7];
  	double phiS=par[8];
  	double qK=par[9];
  	double phiK=par[10];
  	double alp0=par[11];

	double dt_w = timestep;
    int Nps = vlength;
    double Tobs = (Nps-1)*timestep;
    double df = 1./Tobs;
    double dt_ph = 0.01*M*M/mu;
    int Nph = (int)(Tobs/dt_ph);
    dt_ph = Tobs/(Nph-1);
    double arm =  2.5e9/C; // LISA's arm in sec

	complex<double> **Xp,**Xc,**Yp,**Yc,**Zp,**Zc;
  	Xp=(complex<double>**)malloc(5*nmodes*sizeof(complex<double>*));
  	Xc=(complex<double>**)malloc(5*nmodes*sizeof(complex<double>*));
  	Yp=(complex<double>**)malloc(5*nmodes*sizeof(complex<double>*));
  	Yc=(complex<double>**)malloc(5*nmodes*sizeof(complex<double>*));
  	Zp=(complex<double>**)malloc(5*nmodes*sizeof(complex<double>*));
  	Zc=(complex<double>**)malloc(5*nmodes*sizeof(complex<double>*));
	for(int i=0;i<5*nmodes;i++){
    	Xp[i]=(complex<double>*)malloc(Nph*sizeof(complex<double>));
    	Xc[i]=(complex<double>*)malloc(Nph*sizeof(complex<double>));
    	Yp[i]=(complex<double>*)malloc(Nph*sizeof(complex<double>));
    	Yc[i]=(complex<double>*)malloc(Nph*sizeof(complex<double>));
    	Zp[i]=(complex<double>*)malloc(Nph*sizeof(complex<double>));
    	Zc[i]=(complex<double>*)malloc(Nph*sizeof(complex<double>));
    }

  	double *Ampvec,*tvec,*evec,*gimvec,*Phivec,*alpvec,*nuvec,*gimdotvec,*alpdotvec,*gimddotvec,*Phiddotvec,*alpddotvec;
  	Ampvec=(double *)malloc(Nph*sizeof(double));
  	tvec=(double *)malloc(Nph*sizeof(double));
  	evec=(double *)malloc(Nph*sizeof(double));
  	gimvec=(double *)malloc(Nph*sizeof(double));
  	Phivec=(double *)malloc(Nph*sizeof(double));
  	alpvec=(double *)malloc(Nph*sizeof(double));
  	nuvec=(double *)malloc(Nph*sizeof(double));
  	gimdotvec=(double *)malloc(Nph*sizeof(double));
  	alpdotvec=(double *)malloc(Nph*sizeof(double));
  	gimddotvec=(double *)malloc(Nph*sizeof(double));
  	Phiddotvec=(double *)malloc(Nph*sizeof(double));
  	alpddotvec=(double *)malloc(Nph*sizeof(double));

  	int i_plunge;
  	PNevolutionAKTDI(Xp,Xc,Yp,Yc,Zp,Zc,Ampvec,tvec,evec,gimvec,Phivec,alpvec,nuvec,gimdotvec,alpdotvec,gimddotvec,Phiddotvec,alpddotvec,dt_ph,Nph,par,nu0,zeta,nmodes,arm,&i_plunge);

    double cS= cos(qS);
    double sS= sin(qS);
    double cK= cos(qK);
    double sK= sin(qK);
    double cSp= cos(phiS);
    double sSp= sin(phiS);
    double cKp= cos(phiK);
    double sKp= sin(phiK);

    double cX = cS*cK + sS*sK*cos(phiS  - phiK);
    double sX = sqrt( sS*sS*cK*cK - 2.*sS*sSp*cK*cS*sK*sKp + \
                      cS*cS*sK*sK - 2.*cS*sK*cKp*sS*cSp*cK + \
                      sS*sS*cSp*cSp*sK*sK*sKp*sKp - 2.*sS*sS*cSp*sK*sK*sKp*sSp*cKp +\
                      sS*sS*sSp*sSp*sK*sK*cKp*cKp);

    double Apc1, Aps1, Apcn1, Apc2, Aps2, Apcn2, Aqc1, Aqs1, Aqcn1, Aqc2, Aqs2, Aqcn2;

       Apc1 = ( -cK*cKp*sS*cSp - cK*sKp*sS*sSp + sK*cS )/(sX);
       Aps1 = ( sKp*sS*cSp - cKp*sS*sSp )/(sX);
       Apcn1 = ( cK*cS + sK*cKp*sS*cSp + sK*sKp*sS*sSp - cX)*clam/(sX*slam);

       Apc2 = (sS*cSp*sKp - sS*sSp*cKp )*clam/(sX);
       Aps2 = ( cK*cKp*sS*cSp + cK*sKp*sS*sSp - cS*sK )*clam/(sX);
       Apcn2 = 0.0;


       Aqc1 = ( sS*cSp*sKp - sS*sSp*cKp  )*cX/(sX);
       Aqs1 = ( cK*cKp*sS*cSp + cK*sKp*sS*sSp - cS*sK )*cX/(sX);
       Aqcn1 = 0.0;


       Aqc2 = cX*clam*( cK*cKp*sS*cSp + cK*sKp*sS*sSp - sK*cS)/(sX);
       Aqs2 = -cX*clam*sS*( sKp*cSp - cKp*sSp )/(sX);
       Aqcn2 = -( cX*clam*clam*( cK*cS + sK*cKp*sS*cSp + sK*sKp*sS*sSp ) + \
                                1.- cX*cX - clam*clam )/(sX*slam);

    double Bp1c1 = 2.0*(Apc1*Apcn1 - Aqc1*Aqcn1 + Aqc2*Aqcn2 - Apc2*Apcn2);
    double Bp1c2 =  0.5*(Aps2*Aps2 - Aqc1*Aqc1  + Apc1*Apc1  - Aps1*Aps1 + \
                                Aqc2*Aqc2 + Aqs1*Aqs1 - Apc2*Apc2 - Aqs2*Aqs2);
    double Bp1s1 = 2.0*(Aqs2*Aqcn2 - Aps2*Apcn2 - Aqs1*Aqcn1 + Aps1*Apcn1);
    double Bp1s2 = (Apc1*Aps1 + Aqc2*Aqs2 - Apc2*Aps2 - Aqc1*Aqs1);
    double Bp1cn = 0.5*(Apc1*Apc1 + Aps1*Aps1 - Aqc1*Aqc1 - Aqs1*Aqs1 - Apc2*Apc2 \
                                + Aqc2*Aqc2 + Aqs2*Aqs2 - Aps2*Aps2) + Aqcn2*Aqcn2 - Aqcn1*Aqcn1 \
                                + Apcn1*Apcn1 - Apcn2*Apcn2;

    double Bp2c1 = (Apcn1*Apc2 + Apc1*Apcn2 - Aqcn1*Aqc2 - Aqc1*Aqcn2);
    double Bp2c2 = 0.5*(Aqs1*Aqs2 - Aps1*Aps2 + Apc1*Apc2 - Aqc1*Aqc2);
    double Bp2s1 = (Aps1*Apcn2 + Apcn1*Aps2 - Aqcn1*Aqs2 - Aqs1*Aqcn2);
    double Bp2s2 = 0.5*( Apc1*Aps2 - Aqc1*Aqs2 + Aps1*Apc2 - Aqs1*Aqc2);
    double Bp2cn = 0.5*(Aps1*Aps2 - Aqs1*Aqs2 - Aqc1*Aqc2 + Apc1*Apc2) -Aqcn1*Aqcn2 + Apcn1*Apcn2;

    double Bc1c1 = (-Apc2*Aqcn2 - Apcn2*Aqc2 + Apc1*Aqcn1 + Apcn1*Aqc1);
    double Bc1c2 = 0.5*( Apc1*Aqc1 - Aps1*Aqs1 - Apc2*Aqc2 + Aps2*Aqs2);
    double Bc1s1 = (Apcn1*Aqs1 - Aps2*Aqcn2 + Aps1*Aqcn1 - Apcn2*Aqs2);
    double Bc1s2 = 0.5*(-Apc2*Aqs2 + Apc1*Aqs1 - Aps2*Aqc2 + Aps1*Aqc1);
    double Bc1cn = -Apcn2*Aqcn2 + Apcn1*Aqcn1 + 0.5*(Apc1*Aqc1 - Aps2*Aqs2 + Aps1*Aqs1 - Apc2*Aqc2);

    double Bc2c1 = (Aqc1*Apcn2 + Aqcn1*Apc2 + Apc1*Aqcn2 + Apcn1*Aqc2);
    double Bc2c2 = 0.5*( Apc1*Aqc2 - Aps1*Aqs2 + Aqc1*Apc2 - Aqs1*Aps2);
    double Bc2s1 = (Apcn1*Aqs2 + Aqs1*Apcn2 + Aps1*Aqcn2 + Aqcn1*Aps2);
    double Bc2s2 = 0.5*(Aqc1*Aps2 + Apc1*Aqs2 + Aqs1*Apc2 + Aps1*Aqc2);
    double Bc2cn = Aqcn1*Apcn2 + Apcn1*Aqcn2 + 0.5*(Apc1*Aqc2 + Aqs1*Aps2 +Aps1*Aqs2 + Aqc1*Apc2);

    double AApcos[5],AApsin[5],AAccos[5],AAcsin[5];
       AApcos[0]=0.5*(Bp1c2+Bp2s2);
       AApsin[0]=0.5*(Bp2c2-Bp1s2);
       AAccos[0]=0.5*(Bc1c2+Bc2s2);
       AAcsin[0]=0.5*(Bc2c2-Bc1s2);
       AApcos[1]=0.5*(Bp1c1+Bp2s1);
       AApsin[1]=0.5*(Bp2c1-Bp1s1);
       AAccos[1]=0.5*(Bc1c1+Bc2s1);
       AAcsin[1]=0.5*(Bc2c1-Bc1s1);
       AApcos[2]=Bp1cn;
       AApsin[2]=Bp2cn;
       AAccos[2]=Bc1cn;
       AAcsin[2]=Bc2cn;
       AApcos[3]=0.5*(Bp1c1-Bp2s1);
       AApsin[3]=0.5*(Bp2c1+Bp1s1);
       AAccos[3]=0.5*(Bc1c1-Bc2s1);
       AAcsin[3]=0.5*(Bc2c1+Bc1s1);
       AApcos[4]=0.5*(Bp1c2-Bp2s2);
       AApsin[4]=0.5*(Bp2c2+Bp1s2);
       AAccos[4]=0.5*(Bc1c2-Bc2s2);
       AAcsin[4]=0.5*(Bc2c2+Bc1s2);

    double T, ec, hamp;
    int ind_low, ind_up;
    double delta, eps, amp;
    double  xi, Apc, Aps, Acc, Acs;
    int harmms[]={-2,-1,0,1,2};
    double fact;
    complex<double> cFp;
    complex<double> cFc;
    complex<double> x, xp, xc;
    complex<double> img(0.0, 1.0);
    double om_in;
    double om_fin;
    double dOm = 2.*M_PI*df;
    double Om, dom, delom;
    double xi_in, xi_fin, dxi;
    double orbOm = 2.*M_PI/year;
    double DM = AUsec*sS;
    double faza, sinph, cosph;
    int ind;
    double nn,mm;

    for (int i=1; i<i_plunge; i++){
        xi_in = tvec[i-1] - DM*cos(orbOm*tvec[i-1]-phiS);
        xi_fin = tvec[i] - DM*cos(orbOm*tvec[i]-phiS);
        dxi = xi_fin - xi_in;
        ind = 0;
        // loops over harmonics
        for (int j=0;j<nmodes;j++) {
        	nn=(double)(j+1);
        	for(int jj=0; jj<5; jj++){
        	    mm = (double)jj-2.;
        	    om_in = nn*2.*M_PI*nuvec[i-1] + 2.*gimdotvec[i-1]  + mm*alpdotvec[i-1];
        	    om_fin = nn*2.*M_PI*nuvec[i] + 2.*gimdotvec[i]  + mm*alpdotvec[i];
              	delom = om_fin - om_in;
              	ind_low = (int) ceil(om_in/dOm);
              	ind_up = (int) floor(om_fin/dOm);
              	if (om_fin == (double)ind_up*dOm && om_fin != 0.){
                	ind_up = ind_up-1; // otherwise we count twice this bin
              	}
              	if (om_in != 0. && om_fin != 0.){
                 	// loop over fourier bins between two values of harmonic
                 	int ind_max = min(ind_up+1,(Nps+1)/2);
                 	for (int ii=ind_low; ii<ind_max; ii++){
                    	Om = (double)ii * dOm;
                    	delta = (Om - om_in)/delom;
                    	eps = 1.-delta;
                		T =  tvec[i]*delta + tvec[i-1]*eps;
                    	xi = T - DM*cos(orbOm*T - phiS);
                    	delta = (xi - xi_in)/dxi;
                    	eps = 1.-delta;
                    	dom = (nn*Phiddotvec[i-1] + 2.*gimddotvec[i-1]  + mm*alpddotvec[i-1])*eps +
                           (nn*Phiddotvec[i] + 2.*gimddotvec[i]  + mm*alpddotvec[i])*delta;
                    	faza = (nn*Phivec[i-1] + 2.*gimvec[i-1]  + mm*alpvec[i-1])*eps +
                               (nn*Phivec[i] + 2.*gimvec[i]  + mm*alpvec[i])*delta;
                    	amp = Ampvec[i-1]*eps + Ampvec[i]*delta;
                    	ec = evec[i-1]*eps +  evec[i]*delta;
                    	double ne=nn*ec,J[5];
      					if(nn==1){
        					J[0]=-gsl_sf_bessel_J1(ne);
        					J[1]=gsl_sf_bessel_J0(ne);
        					J[2]=gsl_sf_bessel_J1(ne);
        					J[3]=gsl_sf_bessel_Jn(2,ne);
        					J[4]=gsl_sf_bessel_Jn(3,ne);
      					}
      					else gsl_sf_bessel_Jn_array(nn-2,nn+2,ne,J);
                    	hamp = -nn*((J[0]-2.*ec*J[1]+2./nn*J[2]+2.*ec*J[3]-J[4])+sqrt(1-ec*ec)*(J[0]-2.*J[2]+J[4]));
               			delta = (T - tvec[i-1])/dt_ph;
                    	eps = 1.-delta;
                    	cFp = Xp[ind][i-1]*eps + Xp[ind][i]*delta;
                    	cFc = Xc[ind][i-1]*eps + Xc[ind][i]*delta;
                    	Apc = hamp*AApcos[harmms[jj]+2];
         		    	Aps = hamp*AApsin[harmms[jj]+2]; // should be "-" for strain
         		    	Acc = 2.*hamp*AAccos[harmms[jj]+2];
         		    	Acs = 2.*hamp*AAcsin[harmms[jj]+2];

         		    	sinph = sin(faza - Om*xi + M_PI*0.25);
                    	cosph = cos(faza - Om*xi + M_PI*0.25);

                    	fact = 0.5*amp*sqrt(2.*M_PI/dom);
                    	xp = fact*(Apc*cosph + Aps*sinph + img*(Apc*sinph - Aps*cosph));
                    	xc = fact*(Acc*cosph + Acs*sinph + img*(Acc*sinph - Acs*cosph));
                    	Xf_r[ii] += (cFp*xp + cFc*xc).real();
                    	Xf_im[ii] += (cFp*xp + cFc*xc).imag();
                    	cFp = Yp[ind][i-1]*eps + Yp[ind][i]*delta;
                    	cFc = Yc[ind][i-1]*eps + Yc[ind][i]*delta;
                    	Yf_r[ii] += (cFp*xp + cFc*xc).real();
                    	Yf_im[ii] += (cFp*xp + cFc*xc).imag();
                    	cFp = Zp[ind][i-1]*eps + Zp[ind][i]*delta;
                    	cFc = Zc[ind][i-1]*eps + Zc[ind][i]*delta;
                    	Zf_r[ii] += (cFp*xp + cFc*xc).real();
                    	Zf_im[ii] += (cFp*xp + cFc*xc).imag();
                	}
            	}
            	ind ++;
        	}// jj loop
    	} //j loop
	}

	for(int i=0;i<5*nmodes;i++){
    	free(Xp[i]);
    	free(Xc[i]);
    	free(Yp[i]);
    	free(Yc[i]);
    	free(Zp[i]);
    	free(Zc[i]);
    }
  	free(Xp);
  	free(Xc);
  	free(Yp);
  	free(Yc);
  	free(Zp);
  	free(Zc);
  	free(Ampvec);
  	free(tvec);
  	free(evec);
  	free(gimvec);
  	free(Phivec);
  	free(alpvec);
  	free(nuvec);
  	free(gimdotvec);
  	free(alpdotvec);
  	free(gimddotvec);
  	free(Phiddotvec);
  	free(alpddotvec);
  	return;

}


void GenAKTDI(double *f, double *Xf_r, double *Xf_im, double *Yf_r, double *Yf_im, double *Zf_r, double *Zf_im, double timestep, int vlength, double e0, double nu0, double M, double mu, double S, double dist, double inc, double gam0, double Phi0, double qS, double phiS, double alp0, double qK, double phiK, bool mich, bool traj){

  double par[12];
  par[0]=mu*SOLARMASSINSEC;
  par[1]=M*SOLARMASSINSEC;
  par[2]=S;
  par[3]=e0;
  par[4]=inc;
  par[5]=gam0;
  par[6]=Phi0;
  par[7]=qS;
  par[8]=phiS;
  par[9]=qK;
  par[10]=phiK;
  par[11]=alp0;

  // ----- number of modes summed -----
  int nmodes=(int)(30*par[3]);
  if (par[3]<0.135) nmodes=4;
  // ----------

  double zeta=par[0]/dist/Gpc; // M/D

  double Tobs=(vlength-1)*timestep;
  double df=1./Tobs;
  for(int i=0;i<(vlength+1)/2;i++) f[i]=i*df;

  waveformAKTDI(Xf_r,Xf_im,Yf_r,Yf_im,Zf_r,Zf_im,timestep,vlength,nmodes,zeta,par,nu0);
  return;

}
