#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cuComplex.h>
#include "interpolate.hh"
#include "Globals.h"

__device__
double d_dtdm(double v,double e,double Y,double q){
  double v2=v*v;
  double v3=v2*v;
  double v4=v2*v2;
  double e2=e*e;
  double e4=e2*e2;
  double Y2=Y*Y;
  double q2=q*q;
  double eq=(2.*e4*(240. + v2*(-120. + v*(42.*(-27. + 4.*q2)*v + (-6567. + 1996.*q2)*v3 +
              48.*q*(8. + 77.*v2)*Y - 4.*q2*v*(90. + 1577.*v2)*Y2))) +
     e4*e2*(560. + v2*(-360. + 960.*q*v*Y + 8816.*q*v3*Y +
           v4*(-15565. + 24.*q2*(200. - 629.*Y2)) +
           v2*(-2742. + 80.*q2*(5. - 11.*Y2)))) -
     8.*e2*(-48. + v2*(8. - 64.*q*v*Y - 688.*q*v3*Y +
           2.*v2*(99. + 16.*q2*(-1. + 2.*Y2)) + v4*(1233. + 8.*q2*(-47. + 150.*Y2))))
       + 16.*(16. + v2*(24. + v2*(27.*(2. + 5.*v2) - 48.*q*v*Y +
              4.*q2*(2. - v2 + (-2. + 3.*v2)*Y2)))))/(256.*v4);
  return eq;
}

__device__
double d_dphidm(double v,double e,double Y,double q){
  double v2=v*v;
  double v3=v2*v;
  double e2=e*e;
  double e4=e2*e2;
  double Y2=Y*Y;
  double q2=q*q;
  double eq=(16. + 8.*(3. + e2)*v2 - 16.*q*v3*(-2. + (3. + e2)*Y) -
     8.*q*v3*v2*(-6. + 15.*Y + 3.*e4*Y + 2.*e2*(-4. + 7.*Y)) +
     2.*v2*v2*(27. + 3.*e4 + 2.*q2*(-1. + Y)*(1. + 7.*Y) + 2.*e2*(9. + q2*(1. + Y2))) +
     v3*v3*(5.*e4*e2 + e4*(45. + q2*(2. + 26.*Y2)) +
        e2*(135. + 4.*q2*(-19. + 5.*Y*(-7. + 9.*Y))) + 3.*(45. + 2.*q2*(-9. + Y*(-6. + 19.*Y)))))/(16.*v);
  return eq;
}

// ----- magnitude of azimuthal angular frequency for prograde/retrograde orbits -----
__device__
double d_OmegaPhi(double v, double e, double cosiota, double s, double M){

  double omegaphi;
  if(cosiota>0) omegaphi=d_dphidm(v,e,cosiota,s)/d_dtdm(v,e,cosiota,s)/M;
  else omegaphi=d_dphidm(v,e,-cosiota,-s)/d_dtdm(v,e,-cosiota,-s)/M;

  return omegaphi;

}

__device__
void d_cross(const double *u,const double *v,double *w){
  w[0] = u[1]*v[2]-u[2]*v[1];
  w[1] = u[2]*v[0]-u[0]*v[2];
  w[2] = u[0]*v[1]-u[1]*v[0];
}

__device__
double d_dot_product(const double *u,const double *v){
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

__device__
double d_vec_norm(const double *u){
    return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
}


__device__
void d_RotCoeff(double rot[],double iota,double theta_S,double phi_S,double theta_K,double phi_K,double alpha){


  double n[3];
  double L[3];
  double S[3];
  double nxL[3];
  double nxS[3];

  n[0] = sin(theta_S)*cos(phi_S);
  n[1] = sin(theta_S)*sin(phi_S);
  n[2] = cos(theta_S);
  S[0] = sin(theta_K)*cos(phi_K);
  S[1] = sin(theta_K)*sin(phi_K);
  S[2] = cos(theta_K);
  L[0] = cos(iota)*sin(theta_K)*cos(phi_K)+sin(iota)*(sin(alpha)*sin(phi_K)-cos(alpha)*cos(theta_K)*cos(phi_K));
  L[1] = cos(iota)*sin(theta_K)*sin(phi_K)-sin(iota)*(sin(alpha)*cos(phi_K)+cos(alpha)*cos(theta_K)*sin(phi_K));
  L[2] = cos(iota)*cos(theta_K)+sin(iota)*cos(alpha)*sin(theta_K);
  d_cross(n,L,nxL);
  d_cross(n,S,nxS);

  double norm=d_vec_norm(nxL)*d_vec_norm(nxS);
  double dot,cosrot,sinrot;
  //gsl_blas_ddot(nxL,nxS,&dot);
  dot = d_dot_product(nxL,nxS);
  cosrot=dot/norm;
  //gsl_blas_ddot(L,nxS,&dot);
  dot = d_dot_product(L,nxS);
  sinrot=dot;
  //gsl_blas_ddot(S,nxL,&dot);
  dot = d_dot_product(S,nxL);
  sinrot-=dot;
  sinrot/=norm;

  rot[0]=2.*cosrot*cosrot-1.;
  rot[1]=cosrot*sinrot;
  rot[2]=-rot[1];
  rot[3]=rot[0];
}

__device__
void find_index_and_xout(int *index, double *x_out, double *x_out2, double *x_out3, double dx, double x_new, double *x_old, int length){
    double x_trans;
    *index = floor((x_new - x_old[0])/dx);
    if (*index >= length - 1) *index = length - 2;
    x_trans = (x_new - x_old[*index])/(x_old[*index+1] - x_old[*index]);

    *x_out = x_trans;
    *x_out2 = x_trans*x_trans;
    *x_out3 = x_trans*x_trans*x_trans;

    /*# if __CUDA_ARCH__>=200
    if (x_new == 1.000100e+06)
        printf("interp %d, %e %e, %e, %e, %e, %e\n", *index, dx, x_old[0], x_new, x_old[*index], x_old[*index+1], x_trans);
    #endif //*/
}

__device__
double interpolate_array(InterpArrayContainer array_container, double x, double x2, double x3, int index, double x_new){
    double coeff_0 = array_container.array[index];
    double coeff_1 = array_container.coeff_1[index];
    double coeff_2 = array_container.coeff_2[index];
    double coeff_3 = array_container.coeff_3[index];
    double return_val = coeff_0 + coeff_1*x + coeff_2*x2 + coeff_3*x3;

    /*# if __CUDA_ARCH__>=200
    if ((x_new <= 100.0))
        printf("interp2 %d, %e %e, %e, %e, %.18e, %.18e, %.18e, %.18e\n", index, x, x2, x3, x_new, coeff_0, coeff_1, coeff_2, coeff_3);
    #endif //*/

    return return_val;
}

__global__
void kernel_create_waveform(double *t, double *hI, double *hII, double *tvec, InterpArrayContainer evec, InterpArrayContainer vvec, InterpArrayContainer Mvec, InterpArrayContainer Svec, InterpArrayContainer gimvec, InterpArrayContainer Phivec, InterpArrayContainer alpvec, InterpArrayContainer nuvec, InterpArrayContainer gimdotvec, double lam, double qS, double phiS, double qK, double phiK, bool mich, int init_length, int vlength,int nmodes, int i_plunge, int i_buffer, double zeta, double M_phys, double init_dt, double timestep){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= vlength) return;

  // ----------  // TODO: calculate this first section before gpu
  double coslam=cos(lam);
  double sinlam=sin(lam);
  double cosqS=cos(qS);
  double sinqS=sin(qS);
  double cosqK=cos(qK);
  double sinqK=sin(qK);
  double cosphiK=cos(phiK);
  double sinphiK=sin(phiK);
  double halfsqrt3=sqrt(3.)/2.;
  // ----- compute waveform from t_start to t_end -----
  //for(int i=0;i<vlength;i++){
  int index;
  double x, x2, x3;
  double time = timestep*i;
  t[i]= time;
  double t_plunge=tvec[i_plunge];
  double t_zero=t_plunge+timestep*i_buffer;
  if (time <=t_zero){

    hI[i]=0.;
    hII[i]=0.;

    find_index_and_xout(&index, &x, &x2, &x3, init_dt, timestep*i, tvec, init_length);

    double e=interpolate_array(evec, x, x2, x3, index, timestep*i); //evec.array[i];
    double v=interpolate_array(vvec, x, x2, x3, index, timestep*i); //vvec.array[i];
    double M=interpolate_array(Mvec, x, x2, x3, index, timestep*i); //Mvec.array[i];
    double S=interpolate_array(Svec, x, x2, x3, index, timestep*i); //Svec.array[i];
    double gim=interpolate_array(gimvec, x, x2, x3, index, timestep*i); //gimvec.array[i];
    double Phi=interpolate_array(Phivec, x, x2, x3, index, timestep*i); //Phivec.array[i];
    double alp=interpolate_array(alpvec, x, x2, x3, index, timestep*i); //alpvec.array[i];
    double nu=interpolate_array(nuvec, x, x2, x3, index, timestep*i); //nuvec.array[i];
    double gimdot=interpolate_array(gimdotvec, x, x2, x3, index, timestep*i); //gimdotvec.array[i];

    /*# if __CUDA_ARCH__>=200
    if ((i >= 100010) && (i <= 100020))
        printf("%d %e %.18e %.18e %e %e %e %e %e\n", i, t[i], e, v, gim, Phi, alp, nu, gimdot);
    #endif //*/

    double cosalp=cos(alp);
    double sinalp=sin(alp);
    double cosqL=cosqK*coslam+sinqK*sinlam*cosalp;
    double sinqL=sqrt(1.-cosqL*cosqL);
    double phiLup=sinqK*sinphiK*coslam-cosphiK*sinlam*sinalp-cosqK*sinphiK*sinlam*cosalp;
    double phiLdown=sinqK*cosphiK*coslam+sinphiK*sinlam*sinalp-cosqK*cosphiK*sinlam*cosalp;
    double phiL=atan2(phiLup,phiLdown);
    double Ldotn=cosqL*cosqS+sinqL*sinqS*cos(phiL-phiS);
    double Ldotn2=Ldotn*Ldotn;
    double Sdotn=cosqK*cosqS+sinqK*sinqS*cos(phiK-phiS);
    double betaup=-Sdotn+coslam*Ldotn;
    double betadown=sinqS*sin(phiK-phiS)*sinlam*cosalp+(cosqK*Sdotn-cosqS)/sinqK*sinlam*sinalp;
    double beta=atan2(betaup,betadown);
    double gam=2.*(gim+beta);
    double cos2gam=cos(gam);
    double sin2gam=sin(gam);

    double orbphs,cosorbphs,sinorbphs,FplusI,FcrosI,FplusII,FcrosII;
    if(mich){
      orbphs=2.*M_PI*t[i]/year;
      cosorbphs=cos(orbphs-phiS);
      sinorbphs=sin(orbphs-phiS);
      double cosq=.5*cosqS-halfsqrt3*sinqS*cosorbphs;
      double phiw=orbphs+atan2(halfsqrt3*cosqS+.5*sinqS*cosorbphs,sinqS*sinorbphs);
      double psiup=.5*cosqK-halfsqrt3*sinqK*cos(orbphs-phiK)-cosq*(cosqK*cosqS+sinqK*sinqS*cos(phiK-phiS));
      double psidown=.5*sinqK*sinqS*sin(phiK-phiS)-halfsqrt3*cos(orbphs)*(cosqK*sinqS*sin(phiS)-cosqS*sinqK*sin(phiK))-halfsqrt3*sin(orbphs)*(cosqS*sinqK*cos(phiK)-cosqK*sinqS*cos(phiS));
      double psi=atan2(psiup,psidown);
      double cosq1=.5*(1.+cosq*cosq);
      double cos2phi=cos(2.*phiw);
      double sin2phi=sin(2.*phiw);
      double cos2psi=cos(2.*psi);
      double sin2psi=sin(2.*psi);
      FplusI=cosq1*cos2phi*cos2psi-cosq*sin2phi*sin2psi;
      FcrosI=cosq1*cos2phi*sin2psi+cosq*sin2phi*cos2psi;
      FplusII=cosq1*sin2phi*cos2psi+cosq*cos2phi*sin2psi;
      FcrosII=cosq1*sin2phi*sin2psi-cosq*cos2phi*cos2psi;
    }
    else{
      FplusI=1.;
      FcrosI=0.;
      FplusII=0.;
      FcrosII=1.;
    }

    double Amp=pow(d_OmegaPhi(v,e,coslam,S,M)*M_phys*SOLARMASSINSEC,2./3.)*zeta;
    // TODO: check making num modes to gridDim (then need to do reduction to get singular waveform)
    double fn,Doppler,nPhi;
    double ne, a, b, c, Aplus, Acros, Aplusold, Acrosold;
    double rot[4], J[5];
    for(int n=1;n<=nmodes;n++){

      if(mich){
        fn=n*nu+gimdot/M_PI;
        Doppler=2.*M_PI*fn*AUsec*sinqS*cosorbphs;
        nPhi=n*Phi+Doppler;
      }
      else nPhi=n*Phi;

      ne=n*e;
      if(n==1){
        J[0]=-1.0*j1(ne);
        J[1]=j0(ne);
        J[2]=j1(ne);
        J[3]=jn(2,ne);
        J[4]=jn(3,ne);
      }
      else{
          J[0]=jn(n-2, ne);
          J[1]=jn(n-1, ne);
          J[2]=jn(n, ne);
          J[3]=jn(n+1,ne);
          J[4]=jn(n+2,ne);
      }
      a=-n*Amp*(J[0]-2.*e*J[1]+2./n*J[2]+2.*e*J[3]-J[4])*cos(nPhi);
      b=-n*Amp*sqrt(1-e*e)*(J[0]-2.*J[2]+J[4])*sin(nPhi);
      c=2.*Amp*J[2]*cos(nPhi);
      Aplus=-(1.+Ldotn2)*(a*cos2gam-b*sin2gam)+c*(1-Ldotn2);
      Acros=2.*Ldotn*(b*cos2gam+a*sin2gam);

      // ----- rotate to NK wave frame -----
      Aplusold=Aplus;
      Acrosold=Acros;
      d_RotCoeff(rot,lam,qS,phiS,qK,phiK,alp);
      Aplus=Aplusold*rot[0]+Acrosold*rot[1];
      Acros=Aplusold*rot[2]+Acrosold*rot[3];
      // ----------

      double hnI,hnII;
      if(mich){
      	hnI=halfsqrt3*(FplusI*Aplus+FcrosI*Acros);
        hnII=halfsqrt3*(FplusII*Aplus+FcrosII*Acros);
      }
      else{
      	hnI=FplusI*Aplus+FcrosI*Acros;
        hnII=FplusII*Aplus+FcrosII*Acros;
      }

      hI[i]+=hnI;
      hII[i]+=hnII;

    }
  }

  if ((time>t_plunge) &&(i<vlength)){
    if(time<t_zero){
      hI[i]=hI[i]/(exp((t_plunge-t_zero)/(t[i]-t_plunge)+(t_plunge-t_zero)/(t[i]-t_zero))+1.);
      hII[i]=hII[i]/(exp((t_plunge-t_zero)/(t[i]-t_plunge)+(t_plunge-t_zero)/(t[i]-t_zero))+1.);
    }
    else{
      hI[i]=0.;
      hII[i]=0.;
  }
  // ----------

}

}


__global__
void likelihood_prep(cuDoubleComplex *template_channel1, cuDoubleComplex *template_channel2, double *noise_channel1_inv, double *noise_channel2_inv, int length){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= length) return;
    template_channel1[i] = cuCmul(template_channel1[i], make_cuDoubleComplex(noise_channel1_inv[i], 0.0));
    template_channel2[i] = cuCmul(template_channel2[i], make_cuDoubleComplex(noise_channel2_inv[i], 0.0));
}
