#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_trig.h>

#include "IEKG.h"
#include "KSParMap.h"
#include "AAKTDI.h"

using namespace std;


// ----- magnitude of azimuthal angular frequency for prograde/retrograde orbits -----
double OmegaPhiTDI(double v, double e, double cosiota, double s, double M){

  double omegaphi;
  if(cosiota>0) omegaphi=dphidm(v,e,cosiota,s)/dtdm(v,e,cosiota,s)/M;
  else omegaphi=dphidm(v,e,-cosiota,-s)/dtdm(v,e,-cosiota,-s)/M;

  return omegaphi;

}
// ----------

void PNevolutionTDI(complex<double> **Xp, complex<double> **Xc, complex<double> **Yp, complex<double> **Yc, complex<double> **Zp, complex<double> **Zc,\
 double *Ampvec, double *tvec, double *evec, double *vvec, double *Mvec, double *Svec, double *gimvec, double *Phivec, double *alpvec,\
  double *nuvec, double *gimdotvec, double *alpdotvec, double *gimddotvec, double *Phiddotvec, double *alpddotvec, double timestep, int vlength,\
   double *par, double e_traj[], double v_map[], double M_phys, double M_map[], double S_phys, double S_map[], double dt_map[], int steps,\
    double zeta, int nmodes, double arm, int *i_plunge, bool backint){

	double mu=par[0];
	double M0=par[1];
  	double S0=par[2];
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
	double v0=v_map[0];

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

   	int ind;
   	double nn, mm, om;
   	complex<double> chi0;
   	complex<double> chi1;
   	complex<double> chi2;
   	complex<double> img(0.0, 1.0);
   	double x, x2;

  	double edot,vdot,Phidot,gimdot,alpdot;
  	double edotp,vdotp,Phidotp,gimdotp,alpdotp;
  	double edotm,vdotm,Phidotm,gimdotm,alpdotm;
  	double Phiddot,gimddot,alpddot;

  	double dt_large=dt_map[1];
  	double dt_small=dt_map[steps-1]-dt_map[steps-2];
  	int j0=(int)((vlength*timestep-dt_map[steps-1])/dt_large+steps+1);
  	int j_start=j0;
  	int j_min=j_start;
  	int j_end=2*j0-1;
  	int j_max=j_end;
  	int j_max_temp=j_max;
  	double t_end=vlength*timestep;

  	double *e_AK,*v_AK,*t_fit,*e_fit,*v_fit,*M_fit,*S_fit;
  	e_AK=(double*)malloc(2*j0*sizeof(double));
  	v_AK=(double*)malloc(2*j0*sizeof(double));
  	t_fit=(double*)malloc(2*j0*sizeof(double));
  	e_fit=(double*)malloc(2*j0*sizeof(double));
  	v_fit=(double*)malloc(2*j0*sizeof(double));
  	M_fit=(double*)malloc(2*j0*sizeof(double));
  	S_fit=(double*)malloc(2*j0*sizeof(double));

  	// ----- evolve AK from t_0 to T_fit -----
  	e_AK[j0]=e0;
  	v_AK[j0]=v0;
  	for(int j=j0;j<j0+steps-1;j++){
    	edotm=edot;
    	vdotm=vdot;
    	edot=dedt(v_AK[j],e_AK[j],clam,mu,M0,S0);
    	vdot=dvdt(v_AK[j],e_AK[j],clam,mu,M0,S0);
    	if(j==j0){
      		edotm=edot;
      		vdotm=vdot;
    	}
    	e_AK[j+1]=e_AK[j]+(1.5*edot-.5*edotm)*(dt_map[j-j0+1]-dt_map[j-j0]);
    	v_AK[j+1]=v_AK[j]+(1.5*vdot-.5*vdotm)*(dt_map[j-j0+1]-dt_map[j-j0]);
    	if(e_AK[j+1]<0. || v_AK[j+1]<v_AK[j] || isnan(e_AK[j+1]) || isnan(v_AK[j+1])){
      		if(j_max==j_max_temp) j_max=j;
      		e_AK[j+1]=e_AK[j_max];
      		v_AK[j+1]=v_AK[j_max];
    	}
  	}
  	j_max_temp=j_max;
  	// ----------

  	// ----- compute fit coefficients -----
  	int points=min(steps,j_max-j0+1);
  	while(e_traj[points-1]>e_traj[points-2]) points--;
  	if(points<3){fprintf(stderr,"Fitting error: t_0 too close to plunge\n"); exit(EXIT_FAILURE);}

  	double *e_diff,*v_diff,*e_coeff,*v_coeff,*M_coeff,*S_coeff;
  	e_diff=(double*)malloc(points*sizeof(double));
  	v_diff=(double*)malloc(points*sizeof(double));
  	e_coeff=(double*)malloc(2*sizeof(double));
  	v_coeff=(double*)malloc(2*sizeof(double));
  	M_coeff=(double*)malloc(2*sizeof(double));
  	S_coeff=(double*)malloc(2*sizeof(double));

  	for(int i=0;i<points;i++){
    	e_diff[i]=e_traj[i]-e_AK[i+j0];
    	v_diff[i]=v_map[i]-v_AK[i+j0];
  	}

  	PolyFit(e_coeff,dt_map,e_diff,points);
  	PolyFit(v_coeff,dt_map,v_diff,points);
  	PolyFit(M_coeff,dt_map,M_map,points);
  	PolyFit(S_coeff,dt_map,S_map,points);
  	// ----------

  	// ----- evolve AK from T_fit to t_end -----
  	for(int j=j0+steps-1;j<j_end;j++){
    	edotm=edot;
    	vdotm=vdot;
    	edot=dedt(v_AK[j],e_AK[j],clam,mu,M0,S0);
    	vdot=dvdt(v_AK[j],e_AK[j],clam,mu,M0,S0);
    	e_AK[j+1]=e_AK[j]+(1.5*edot-.5*edotm)*dt_large;
    	v_AK[j+1]=v_AK[j]+(1.5*vdot-.5*vdotm)*dt_large;
    	if(e_AK[j+1]<0. || v_AK[j+1]<v_AK[j] || isnan(e_AK[j+1]) || isnan(v_AK[j+1])){
      		if(j_max==j_max_temp) j_max=j;
      		e_AK[j+1]=e_AK[j_max];
      		v_AK[j+1]=v_AK[j_max];
    	}
  	}
  	j_max_temp=j_max;
  	// ----------

  	// ----- fit AK from t_0 to t_end -----
  	for(int j=j0;j<=j_end;j++){
  		double dt;
  		if(j<j0+steps) dt=dt_map[j-j0];
    	else dt=dt_map[steps-1]+dt_large*(j-j0-steps+1);
  		double dt2=dt*dt;
    	t_fit[j]=dt;
    	e_fit[j]=e_AK[j]+e_coeff[0]*dt+e_coeff[1]*dt2;
    	v_fit[j]=v_AK[j]+v_coeff[0]*dt+v_coeff[1]*dt2;
    	M_fit[j]=M0+(M_coeff[0]*dt+M_coeff[1]*dt2)*SOLARMASSINSEC;
    	S_fit[j]=S0+S_coeff[0]*dt+S_coeff[1]*dt2;
    	if(e_fit[j]<0. || v_fit[j]<v_fit[max(j0,j-1)]){
      		if(j_max==j_max_temp) j_max=j-1;
      		e_fit[j]=e_fit[j_max];
      		v_fit[j]=v_fit[j_max];
      		M_fit[j]=M_fit[j_max];
      		S_fit[j]=S_fit[j_max];
    	}
  	}
  	j_max_temp=j_max;
  	// ----------

  	// ----- check for plunge -----
  	int j_RR=(int)((SOLARMASSINSEC*M_phys*SOLARMASSINSEC*M_phys/mu)/dt_large)+1;
  	double z1=1.+pow(1.-S_phys*S_phys,1./3.)*(pow(1.+S_phys,1./3.)+pow(1.-S_phys,1./3.));
  	double z2=sqrt(3.*S_phys*S_phys+z1*z1);
  	double LSO_min=3.+z2-sqrt((3.-z1)*(3.+z1+2.*z2));
  	double LSO_max=3.+z2+sqrt((3.-z1)*(3.+z1+2.*z2));

  	for(int j=j0;j<=j_max_temp;j++){
    	if(j>j_max) break;
    	if((1.-e_fit[j]*e_fit[j])*pow(OmegaPhiTDI(v_fit[j],e_fit[j],clam,S_fit[j],M_fit[j])*M_phys*SOLARMASSINSEC,-2./3.)<LSO_min && j_max==j_max_temp){
      		j_min=max(j-j_RR,j0);
      		j_max=min(j+1,j_max_temp);
    	}
    	if(j>j0 && (j-j0)%j_RR==0 && (1.-e_fit[j]*e_fit[j])*pow(OmegaPhiTDI(v_fit[j],e_fit[j],clam,S_fit[j],M_fit[j])*M_phys*SOLARMASSINSEC,-2./3.)<LSO_max && j_max==j_max_temp){
      		j_min=max(j-j_RR,j0);
      		IEKG iekg((1.-e_fit[j]*e_fit[j])*pow(OmegaPhiTDI(v_fit[j],e_fit[j],clam,S_fit[j],M_fit[j])*M_phys*SOLARMASSINSEC,-2./3.),e_fit[j],clam,S_phys);
      		if(iekg.Stable==-1 || iekg.E>1) j_max=min(j+1,j_max_temp);
    	}
  	}

  	while(j_max-j_min>1){
    	int j_mid=(j_max+j_min)/2;
    	IEKG iekg((1.-e_fit[j_mid]*e_fit[j_mid])*pow(OmegaPhiTDI(v_fit[j_mid],e_fit[j_mid],clam,S_fit[j_mid],M_fit[j_mid])*M_phys*SOLARMASSINSEC,-2./3.),e_fit[j_mid],clam,S_phys);
    	if(iekg.Stable==-1 || iekg.E>1) j_max=j_mid;
    	else j_min=j_mid;
  	}
  	// ----------

  	// ----- evolve and fit AK from t_0 to t_start -----
  	if(backint){

  		j_end=j_max;
    	t_end=min(t_fit[j_end],vlength*timestep);
    	j_start=j0+floor((t_end-vlength*timestep)/dt_large);

    	for(int j=j0;j>j_start;j--){
      		edotp=edot;
      		vdotp=vdot;
      		edot=dedt(v_AK[j],e_AK[j],clam,mu,M0,S0);
      		vdot=dvdt(v_AK[j],e_AK[j],clam,mu,M0,S0);
      		if(j==j0){
        		edotp=edot;
        		vdotp=vdot;
      		}
      		e_AK[j-1]=e_AK[j]-(1.5*edot-.5*edotp)*dt_large;
      		v_AK[j-1]=v_AK[j]-(1.5*vdot-.5*vdotp)*dt_large;
    	}

    	for(int j=j_start;j<j0;j++){
      		double dt_undecayed=dt_large*(j-j0);
      		double dt;
      		if(j>j0-points) dt=(dt_undecayed+dt_undecayed*(j-j0+points)/points)/2.; // decaying C1 fit
      		else dt=-dt_large*points/2.;
      		double dt2=dt*dt;
      		t_fit[j]=dt_undecayed;
      		e_fit[j]=e_AK[j]+e_coeff[0]*dt+e_coeff[1]*dt2;
      		v_fit[j]=v_AK[j]+v_coeff[0]*dt+v_coeff[1]*dt2;
      		M_fit[j]=M0+(M_coeff[0]*dt+M_coeff[1]*dt2)*SOLARMASSINSEC;
      		S_fit[j]=S0+S_coeff[0]*dt+S_coeff[1]*dt2;
    	}

  	}
  	// ----------

  	// ----- interpolate AK from t_start to t_end -----
  	double *t_in,*e_in,*v_in,*M_in,*S_in;
  	t_in=(double*)malloc((j_end-j_start+1)*sizeof(double));
  	e_in=(double*)malloc((j_end-j_start+1)*sizeof(double));
  	v_in=(double*)malloc((j_end-j_start+1)*sizeof(double));
  	M_in=(double*)malloc((j_end-j_start+1)*sizeof(double));
  	S_in=(double*)malloc((j_end-j_start+1)*sizeof(double));

  	for(int j=j_start;j<=j_end;j++){
    	t_in[j-j_start]=t_fit[j];
    	e_in[j-j_start]=e_fit[j];
    	v_in[j-j_start]=v_fit[j];
    	M_in[j-j_start]=M_fit[j];
    	S_in[j-j_start]=S_fit[j];    
  	}

  	int i0=vlength-(int)(t_end/timestep);
  	tvec[0]=t_in[0];
  	for(int i=1;i<vlength;i++) tvec[i]=(i-i0)*timestep;
  	tvec[vlength]=t_in[j_end-j_start];

  	Interp(t_in,e_in,j_end-j_start+1,tvec,evec,vlength+1);
  	Interp(t_in,v_in,j_end-j_start+1,tvec,vvec,vlength+1);
  	Interp(t_in,M_in,j_end-j_start+1,tvec,Mvec,vlength+1);
  	Interp(t_in,S_in,j_end-j_start+1,tvec,Svec,vlength+1);
  	// ----------

  	// ----- evolve phases from t_0 to t_end -----
  	if(j_max==j_max_temp) *i_plunge=vlength-1;
  	else *i_plunge=i0+(int)(t_fit[j_min]/timestep);

  	gimvec[i0]=gim0;
  	Phivec[i0]=Phi0;
  	alpvec[i0]=alp0;
  	for(int i=i0;i<*i_plunge;i++){

    	gimdotm=gimdot;
    	Phidotm=Phidot;
    	alpdotm=alpdot;
    	gimdot=(dthetadm(vvec[i],evec[i],clam,Svec[i])-drdm(vvec[i],evec[i],clam,Svec[i]))/dtdm(vvec[i],evec[i],clam,Svec[i])/Mvec[i];
    	Phidot=drdm(vvec[i],evec[i],clam,Svec[i])/dtdm(vvec[i],evec[i],clam,Svec[i])/Mvec[i];
    	alpdot=(dphidm(vvec[i],evec[i],clam,Svec[i])-dthetadm(vvec[i],evec[i],clam,Svec[i]))/dtdm(vvec[i],evec[i],clam,Svec[i])/Mvec[i];
    	if(i==i0){
      		gimdotm=gimdot;
      		Phidotm=Phidot;
      		alpdotm=alpdot;
    	}
    	Ampvec[i]=pow(OmegaPhiTDI(vvec[i],evec[i],clam,Svec[i],Mvec[i])*M_phys*SOLARMASSINSEC,2./3.)*zeta;
    	nuvec[i]=Phidot/2./M_PI;
    	gimdotvec[i]=gimdot;
    	alpdotvec[i]=alpdot;
    	gimvec[i+1]=gimvec[i]+(1.5*gimdot-.5*gimdotm)*timestep;
    	Phivec[i+1]=Phivec[i]+(1.5*Phidot-.5*Phidotm)*timestep;
    	alpvec[i+1]=alpvec[i]+(1.5*alpdot-.5*alpdotm)*timestep;
    	if(i>i0){
      		gimddotvec[i-1]=(gimdotvec[i]-gimdotvec[i-1])/timestep;
      		Phiddotvec[i-1]=2.*M_PI*(nuvec[i]-nuvec[i-1])/timestep;
      		alpddotvec[i-1]=(alpdotvec[i]-alpdotvec[i-1])/timestep;
    	}
    	if(i==*i_plunge-1){
      		gimddotvec[i]=gimddotvec[i-1];
      		Phiddotvec[i]=Phiddotvec[i-1];
      		alpddotvec[i]=alpddotvec[i-1];
    	}

    	EccentricLISAMotion(0.0, 0.0, tvec[i], R, q, n);
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
                om = nn*Phidot + 2.*gimdot + mm*alpdot;
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

  	}
  	// ----------

  	// ----- evolve phases from t_0 to t_start -----
  	for(int i=i0;i>=0;i--){

    	gimdotp=gimdot;
    	Phidotp=Phidot;
    	alpdotp=alpdot;
    	gimdot=(dthetadm(vvec[i],evec[i],clam,Svec[i])-drdm(vvec[i],evec[i],clam,Svec[i]))/dtdm(vvec[i],evec[i],clam,Svec[i])/Mvec[i];
    	Phidot=drdm(vvec[i],evec[i],clam,Svec[i])/dtdm(vvec[i],evec[i],clam,Svec[i])/Mvec[i];
    	alpdot=(dphidm(vvec[i],evec[i],clam,Svec[i])-dthetadm(vvec[i],evec[i],clam,Svec[i]))/dtdm(vvec[i],evec[i],clam,Svec[i])/Mvec[i];
    	if(i==i0){
    		gimdotp=gimdot;
    		Phidotp=Phidot;
    		alpdotp=alpdot;
    	}
		Ampvec[i]=pow(OmegaPhiTDI(vvec[i],evec[i],clam,Svec[i],Mvec[i])*M_phys*SOLARMASSINSEC,2./3.)*zeta;
    	nuvec[i]=Phidot/2./M_PI;
    	gimdotvec[i]=gimdot;
    	alpdotvec[i]=alpdot;
    	if(i>0){
    		gimvec[i-1]=gimvec[i]-(1.5*gimdot-.5*gimdotp)*timestep;
    		Phivec[i-1]=Phivec[i]-(1.5*Phidot-.5*Phidotp)*timestep;
    		alpvec[i-1]=alpvec[i]-(1.5*alpdot-.5*alpdotp)*timestep;
    	}
    	if(i<i0){
      		gimddotvec[i+1]=(gimdotvec[i+1]-gimdotvec[i])/timestep;
      		Phiddotvec[i+1]=2.*M_PI*(nuvec[i+1]-nuvec[i])/timestep;
      		alpddotvec[i+1]=(alpdotvec[i+1]-alpdotvec[i])/timestep;
    	}
    	if(i==0){
    		gimddotvec[i]=gimddotvec[i+1];
      		Phiddotvec[i]=Phiddotvec[i+1];
      		alpddotvec[i]=alpddotvec[i+1];
    	}

    	EccentricLISAMotion(0.0, 0.0, tvec[i], R, q, n);
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
                om = nn*Phidot + 2.*gimdot + mm*alpdot;
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

  	}
  	// ----------

   	for (int i=0; i<3; i++){
       	free(q[i]);
       	free(n[i]);
   	}
  	free(R);
  	free(q);
  	free(n);
  	free(e_AK);
  	free(v_AK);
  	free(t_fit);
  	free(e_fit);
  	free(v_fit);
  	free(M_fit);
  	free(S_fit);
  	free(e_diff);
  	free(v_diff);
  	free(e_coeff);
  	free(v_coeff);
  	free(M_coeff);
  	free(S_coeff);
  	free(t_in);
  	free(e_in);
  	free(v_in);
  	free(M_in);
  	free(S_in);
  	return;

}


void waveformTDI(double *Xf_r, double *Xf_im, double *Yf_r, double *Yf_im, double *Zf_r, double *Zf_im, double timestep, int vlength, int nmodes, double zeta, double *par, double e_traj[], double v_map[], double M_phys, double M_map[], double S_phys, double S_map[], double dt_map[], int steps, bool backint){

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
    	Xp[i]=(complex<double>*)malloc((Nph+1)*sizeof(complex<double>));
    	Xc[i]=(complex<double>*)malloc((Nph+1)*sizeof(complex<double>));
    	Yp[i]=(complex<double>*)malloc((Nph+1)*sizeof(complex<double>));
    	Yc[i]=(complex<double>*)malloc((Nph+1)*sizeof(complex<double>));
    	Zp[i]=(complex<double>*)malloc((Nph+1)*sizeof(complex<double>));
    	Zc[i]=(complex<double>*)malloc((Nph+1)*sizeof(complex<double>));
    }

  	double *Ampvec,*tvec,*evec,*vvec,*Mvec,*Svec,*gimvec,*Phivec,*alpvec,*nuvec,*gimdotvec,*alpdotvec,*gimddotvec,*Phiddotvec,*alpddotvec;
  	Ampvec=(double *)malloc((Nph+1)*sizeof(double));
  	tvec=(double *)malloc((Nph+1)*sizeof(double));
  	evec=(double *)malloc((Nph+1)*sizeof(double));
	vvec=(double *)malloc((Nph+1)*sizeof(double));
	Mvec=(double *)malloc((Nph+1)*sizeof(double));
	Svec=(double *)malloc((Nph+1)*sizeof(double));
  	gimvec=(double *)malloc((Nph+1)*sizeof(double));
  	Phivec=(double *)malloc((Nph+1)*sizeof(double));
  	alpvec=(double *)malloc((Nph+1)*sizeof(double));
  	nuvec=(double *)malloc((Nph+1)*sizeof(double));
  	gimdotvec=(double *)malloc((Nph+1)*sizeof(double));
  	alpdotvec=(double *)malloc((Nph+1)*sizeof(double));
  	gimddotvec=(double *)malloc((Nph+1)*sizeof(double));
  	Phiddotvec=(double *)malloc((Nph+1)*sizeof(double));
  	alpddotvec=(double *)malloc((Nph+1)*sizeof(double));

	int i_plunge;
	PNevolutionTDI(Xp,Xc,Yp,Yc,Zp,Zc,Ampvec,tvec,evec,vvec,Mvec,Svec,gimvec,Phivec,alpvec,nuvec,gimdotvec,alpdotvec,gimddotvec,Phiddotvec,alpddotvec,dt_ph,Nph,par,e_traj,v_map,M_phys,M_map,S_phys,S_map,dt_map,steps,zeta,nmodes,arm,&i_plunge,backint);

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
  	free(vvec);
  	free(Mvec);
  	free(Svec);
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


void GenTDI(double *f, double *Xf_r, double *Xf_im, double *Yf_r, double *Yf_im, double *Zf_r, double *Zf_im, double timestep, int vlength, double e_traj[], double v_map[], double M_phys, double M_map[], double mu, double S_phys, double S_map[], double dist, double inc, double gam0, double Phi0, double qS, double phiS, double alp0, double qK, double phiK, double dt_map[], int steps, bool backint, bool mich, bool traj){

  double par[12];
  par[0]=mu*SOLARMASSINSEC;
  par[1]=M_map[0]*SOLARMASSINSEC;
  par[2]=S_map[0];
  par[3]=e_traj[0];
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

  waveformTDI(Xf_r,Xf_im,Yf_r,Yf_im,Zf_r,Zf_im,timestep,vlength,nmodes,zeta,par,e_traj,v_map,M_phys,M_map,S_phys,S_map,dt_map,steps,backint);
  return;

}
