#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "KSParMap.h"
#include "AK.h"

using namespace std;

// NB M, mu in solar masses, dist in Gpc, t in seconds, nu0 in Hz;

/*   Code modified from Fortran code kindly provided by Leor Barack.

     ******************************************************************
     *** "Waveform.for calculates waveform for paper                ***
     *** of high mass-ratio inspiral.                               ***
     ***                                                            ***
     *** Precision is double .                                      ***
     *** Input parameter are dt - dutation of 1 integration,        ***
     ***                     eisco - eccentricity at isco, and      ***
     ***                     S - Spin of massive hole               ***
     *** We take t=0 at middle of integration, and parametrize      ***
     *** the system by the values of its intrinsic variables        ***
     *** at t0. We calculate e(t=0) and nu(t=0) by integrating      ***
     *** back in time from the isco.                                ***
     ******************************************************************

       *** NOTATION:
       q: theta
       gam: gamma
       lam: lambda
       alp: alpha
       S: Dimensionless! In fact, this is S/M^2 (or a/M)

       THIS VERSION HAS THE 1/2 PROBLEM FIXED
       THIS VERSION HAS THE NOISE CURVE FIXED
       THIS VERSION USES the PURELY INTRINSIC GAMMA
       THIS VERSION HAS THE BESSEL FUNCTION FIXED
       THIS VERSION HAS THE phiw PROBLEM FIXED */


/* ************************************************************
   FUNCTION ArcT
   ************************************************************
   returns ArcTan(up/down), taking into account the quadrant
   of the point (down,up) */

double ArcT(double down, double up) {
  double ArcT;
  ArcT=atan(up/down);
  if (down < 0.) ArcT=ArcT+M_PI;
  //     if ((up < 0.) && (down > 0)) ArcT=ArcT+2.*M_PI;
  return ArcT;
}

/* ************************************************************
   FUNCTION J0
   ************************************************************ */

double J0(double x)
{
  double ax,z;
  double xx,y,ans,ans1,ans2;
  
  if ((ax=fabs(x)) < 8.0) {
    y=x*x;
    ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
					    +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
    ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
					  +y*(59272.64853+y*(267.8532712+y*1.0))));
    ans=ans1/ans2;
  } else {
    z=8.0/ax;
    y=z*z;
    xx=ax-0.785398164;
    ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
				    +y*(-0.2073370639e-5+y*0.2093887211e-6)));
    ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			       +y*(-0.6911147651e-5+y*(0.7621095161e-6
						       -y*0.934935152e-7)));
    ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }
  return ans;
}

/* ************************************************************
   FUNCTION J1
   ************************************************************ */

double J1(double x)
{
  double ax,z;
  double xx,y,ans,ans1,ans2;
  
  if ((ax=fabs(x)) < 8.0) {
    y=x*x;
    ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
					      +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
    ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
					   +y*(99447.43394+y*(376.9991397+y*1.0))));
    ans=ans1/ans2;
  } else {
    z=8.0/ax;
    y=z*z;
    xx=ax-2.356194491;
    ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
			       +y*(0.2457520174e-5+y*(-0.240337019e-6))));
    ans2=0.04687499995+y*(-0.2002690873e-3
			  +y*(0.8449199096e-5+y*(-0.88228987e-6
						 +y*0.105787412e-6)));
    ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
    if (x < 0.0) ans = -ans;
  }
  return ans;
}

/* ************************************************************
   FUNCTION Jn
   ************************************************************ */

double Jn(int n, double x) {
  const int IACC=160;
  const double BIGNO=1.e10;
  const double BIGNI=1.e-10;
  int j,jsum,m;
  double ax,bj,bjm,bjp,sum,tox,Jnval;
  if(n == 0)
    return J0(x);
  if(n == 1)
    return J1(x);
  ax=fabs(x);
  if(ax == 0.)
    Jnval=0.;
  else if(ax > ((double)n)) {
    tox=2./ax;
    bjm=J0(ax);
    bj=J1(ax);
    for(j=1;j<n;j++) {
      bjp=j*tox*bj-bjm;
      bjm=bj;
      bj=bjp;
    }
    Jnval=bj;
  } else {
    tox=2./ax;
    m=2*((n+((int)(sqrt((double)(IACC*n)))))/2);
    Jnval=0.;
    jsum=0;
    sum=0.;
    bjp=0.;
    bj=1.;
    for (j=m;j>0;j--) {
      bjm=j*tox*bj-bjp;
      bjp=bj;
      bj=bjm;
      if(fabs(bj) > BIGNO) {
	bj=bj*BIGNI;
	bjp=bjp*BIGNI;
	Jnval*=BIGNI;
	sum*=BIGNI;
      }
      if(jsum != 0)
	sum+=bj;
      jsum=1-jsum;
      if(j == n)
	Jnval = bjp;
    }
    sum=2.*sum-bj;
    Jnval/=sum;
  }
  if((x < 0.) && (n & 1)) 
    Jnval*=-1.;
  /* **********************************************
     ADDED BY LEOR since NumRec does not treat n=-1
     ********************************************** */
  if(n == -1)
    Jnval*=-1.;
  return Jnval;
}


/* ************************************************************
   SUBROUTINE PNEVOLUTION
   ************************************************************

   <<< Solves coupled equations for e(t),nu(t),gam(t),Phi(t),alp(t)    >>>
   <<< from t=0 to t=tend with a rough grid [these quantities          >>>
   <<< cahnge very slowly over a few weeks - e.g., edot~-1e-6/hour]    >>>
   <<< "Initial" values are given at t=t0~dend/2, from which I evolve  >>>
   <<< to the future and to the past.                                  >>>
   <<<                                                                 >>>
   <<< This subroutine returns vectors e(i=0,1,...),nu(i=0,1,...),...  >>>
   <<< Actual time difference between grid points is
   <<<                                         timestep=tend/vlength   >>>
   <<< Actual time at step i is t=i*tend/vlength (sec)                 >>>
   <<<                                                                 >>>
   <<< The evolution equations converge quadratically                  >>>
   <<< (thought this is not necessary - check how affects running time)>>> */


void PNevolution(int vlength, double timestep, double *par, double nu0, double *gimdotvec, double *e, double *nu, double *Phi,
		 double *gim, double *alp) {
  int i,i0;
  double tend;
  double gimdot;
  double edot,nudot,Phidot,alpdot;
  double edotp,nudotp,Phidotp,gimdotp,alpdotp;
  double edotm,nudotm,Phidotm,gimdotm,alpdotm;
  double t0,mu,M,e0,Phi0,qS,phiS,gim0;
  double lam,alp0,S,qK,phiK;
  double Z,Y,hour,e2;
  double cosqS,sinqS,coslam,sinlam,cosqK,sinqK;
  double cosalp0,sinalp0;
  double SdotN,L0dotN,NcrossSdotL0,kappa0;

  hour=3600.;
  
  t0=par[1];
  mu=par[2];
  M= par[3];
  e0=par[4];
  gim0=par[5];
  Phi0=par[6];
  qS=par[7];
  phiS=par[8];
  lam=par[9];
  alp0=par[10];
  S=par[11];
  qK=par[12];
  phiK=par[13];

  cosqS=cos(qS);
  sinqS=sin(qS);
  coslam=cos(lam);
  sinlam=sin(lam);
  cosqK=cos(qK);
  sinqK=sin(qK);
  cosalp0=cos(alp0);
  sinalp0=sin(alp0);

  /*     ----------------------------------------------
	 HERE CALCULATE gam0 AS A FUNCTION OF gamtilde0
	 
	 SdotN=cosqK*cosqS+sinqK*sinqS*dcos(phiK-phiS)
	 L0dotN=SdotN*coslam+(cosqS-SdotN*cosqK)*sinlam*cosalp0/sinqK
	 p       +sinqS*sin(phiK-phiS)*sinlam*sinalp0
	 NcrossSdotL0=sinqS*dsin(phiK-phiS)*sinlam*cosalp0
	 p             -(cosqS-SdotN*cosqK)*sinlam*sinalp0/sinqK
	 kappa0=ArcT(SdotN-coslam*L0dotN,NcrossSdotL0)
	 gam0=gamtilde0+kappa0-pi/2
	 write(*,*) gam0
	 ---------------------------------------------- */

  i0=(int)(t0/timestep);

  e[i0]=e0;
  nu[i0]=nu0;
  Phi[i0]=Phi0;
  gim[i0]=gim0;
  alp[i0]=alp0;

  //     <<< EVOLVE FORWARD FROM t0 to tend >>>

  for (i=i0;i<vlength;i++) {
    Z=pow(2.*M_PI*M*nu[i],1./3.);
    e2=e[i]*e[i];
    Y=1./(1.-e2);
    //cout << i << " " << Z << " " << Y << endl;
    edotm=edot;
    nudotm=nudot;
    Phidotm=Phidot;
    gimdotm=gimdot;
    alpdotm=alpdot;
    edot=e[i]*mu/M/M*(-1./15.*pow(Y,3.5)*pow(Z,8.)*((304.+121.*e2)/Y+Z*Z*(70648.-231960.*e2-56101.*e2*e2)/56.)
		      +S*coslam*pow(Z,11.)*pow(Y,4.)*(8184.+10064.*e2+789.*e2*e2)/30.);
    nudot=96./(10.*M_PI)*mu/pow(M,3)*(pow(Z,11.)*pow(Y,4.5)*((96.+292.*e2+37.*e2*e2)/Y/96.
							   +Z*Z*(20368.-61464.*e2-163170.*e2*e2-13147.*e2*e2*e2)/5376.)
				    -pow(Z,14.)*pow(Y,5.)*S*coslam*(1168.+9688.*e2+6286.*e2*e2 +195.*e2*e2*e2)/192.);
    Phidot=2.*M_PI*nu[i];
    alpdot=8.*M_PI*M_PI*nu[i]*nu[i]*S*M*pow(Y,1.5);
    gimdot=6.*M_PI*nu[i]*Z*Z*Y*(1.+.25*Z*Z*Y*(26.-15.*e2))-3.*coslam*alpdot;
    //if (!i)
    //	    cout << "gimdot0/pi = " << gimdot/M_PI << ", alpdot0 = " << alpdot/(2.*M_PI) << endl;
    if (i == i0) {
      edotm=edot;
      nudotm=nudot;
      Phidotm=Phidot;
      gimdotm=gimdot;
      alpdotm=alpdot;
      gimdotm=gimdot;
    }
    //cout << i << " " << edot << " " << nudot << endl;
    gimdotvec[i]=gimdot;
    e[i+1]=e[i]+(1.5*edot-.5*edotm)*timestep;
    nu[i+1]=nu[i]+(1.5*nudot-.5*nudotm)*timestep;
    Phi[i+1]=Phi[i]+(1.5*Phidot-.5*Phidotm)*timestep;
    gim[i+1]=gim[i]+(1.5*gimdot-.5*gimdotm)*timestep;
    alp[i+1]=alp[i]+(1.5*alpdot-.5*alpdotm)*timestep;
  }
  gimdotvec[i]=gimdot;

  //     <<< EVOLVE BACKWARD FROM t0 to t=0 >>>

  for (i=i0;i>0;i--) {
    Z=pow(2.*M_PI*M*nu[i],1./3.);
    e2=e[i]*e[i];
    Y=1./(1.-e2);
    edotp=edot;
    nudotp=nudot;
    Phidotp=Phidot;
    gimdotp=gimdot;
    alpdotp=alpdot;

    edot=e[i]*mu/M/M*(-1./15.*pow(Y,3.5)*pow(Z,8.)*((304.+121.*e2)/Y+Z*Z*(70648.-231960.*e2-56101.*e2*e2)/56.)
		      +S*coslam*pow(Z,11.)*pow(Y,4.)*(8184.+10064.*e2+789.*e2*e2)/30.);
    nudot=96./(10.*M_PI)*mu/pow(M,3)*(pow(Z,11.)*pow(Y,4.5)*((96.+292.*e2+37.*e2*e2)/Y/96.
							   +Z*Z*(20368.-61464.*e2-163170.*e2*e2-13147.*e2*e2*e2)/5376.)
				    -pow(Z,14.)*pow(Y,5.)*S*coslam*(1168.+9688.*e2+6286.*e2*e2 +195.*e2*e2*e2)/192.);
    Phidot=2.*M_PI*nu[i];
    alpdot=8.*M_PI*M_PI*nu[i]*nu[i]*S*M*pow(Y,1.5);
    gimdot=6.*M_PI*nu[i]*Z*Z*Y*(1.+.25*Z*Z*Y*(26.-15.*e2))-3.*coslam*alpdot;
    if (!i)
	    cout << "gimdot0 = " << gimdot << ", alpdot0 = " << alpdot << endl;
    if (i == i0) {
      edotp=edot;
      nudotp=nudot;
      Phidotp=Phidot;
      gimdotp=gimdot;
      alpdotp=alpdot;
    }
    gimdotvec[i]=gimdot;
    e[i-1]=e[i]-(1.5*edot-.5*edotp)*timestep;
    nu[i-1]=nu[i]-(1.5*nudot-.5*nudotp)*timestep;
    Phi[i-1]=Phi[i]-(1.5*Phidot-.5*Phidotp)*timestep;
    gim[i-1]=gim[i]-(1.5*gimdot-.5*gimdotp)*timestep;
    alp[i-1]=alp[i]-(1.5*alpdot-.5*alpdotp)*timestep;
  }
  gimdotvec[0]=gimdot;
  //cout << "At t = 0, e = " << e[0] << " and nu = " << nu[0] << endl;
  //cout << "At t = t_end, e = " << e[vlength-1] << " and nu = " << nu[vlength-1] << endl;

  return;
}


/* ************************************************************
   SUBROUTINE WAVEFORM
   ************************************************************ */

void waveform(double tend,double *par, double nu0, int vlength, double timestep, double *hI, double *hII, int nmodes, double zeta, int WCN, bool mich, bool traj) {
  int i,i0,n;
  double t0,mu,M,qS,phiS,lam,S,qK,phiK,Sn;
  double year,fn,invsqrtS[modes+1];
  double e,nu,Phi,gim,alp,t,ne,nPhi;
  double a,b,c,An,Amp,Jm2,Jm1,Jm0,Jp1,Jp2,it1,tint;
  double cosqL,sinqL,Ldotn,phiL,BB,CC,Ldotn2;
  double gam2,Aplus,Acros,Sdotn,PhiT,PhiTdot,alpdot;
  double halfsqrt3,prefact,cosq,cosqS,sinqS,coslam,sinlam;
  double cosqK,sinqK,cosalp,sinalp,cosphiK,sinphiK;
  double orbphs,cosorbphs,sinorbphs,phiw,psi,psidown,psiup;
  double cos2phi,sin2phi,cos2psi,sin2psi,cosq1;
  double FplusI,FcrosI,FplusII,FcrosII,AUsec,Doppler;
  double hnI,hnII,cos2gam,sin2gam,h2I,h3I;
  double betaup,betadown,beta;

  double *evec,*nuvec,*alpvec,*Phivec,*gimvec,*gimdotvec;
  evec=(double *)malloc((vlength+1)*sizeof(double));
  nuvec=(double *)malloc((vlength+1)*sizeof(double));
  alpvec=(double *)malloc((vlength+1)*sizeof(double));
  Phivec=(double *)malloc((vlength+1)*sizeof(double));
  gimvec=(double *)malloc((vlength+1)*sizeof(double));
  gimdotvec=(double *)malloc((vlength+1)*sizeof(double));
  
  halfsqrt3=0.5*sqrt(3.);

  if (mich)
    prefact=halfsqrt3;
  else
    prefact=1.;

  year=31536000.;
  AUsec=499.004783702731;
  
  t0=par[1];
  mu=par[2];
  M= par[3];
  qS=par[7];
  phiS=par[8];
  lam=par[9];
  S=par[11];
  qK=par[12];
  phiK=par[13];

  coslam=cos(lam);
  sinlam=sin(lam);
  cosqS=cos(qS);
  sinqS=sin(qS);
  cosqK=cos(qK);
  sinqK=sin(qK);
  cosphiK=cos(phiK);
  sinphiK=sin(phiK);

  PNevolution(vlength,timestep,par,nu0,gimdotvec,evec,nuvec,Phivec,gimvec,alpvec);
  //cout << "Finished PNevolution. Exiting!" << endl;
  //exit(0);

  // ----- output trajectory -----
  if(traj==true){
    for(int i=0;i<vlength;i++){
      hI[i]=nuvec[i];
      hII[i]=evec[i];
    }
    free(evec);
    free(nuvec);
    free(alpvec);
    free(Phivec);
    free(gimvec);
    free(gimdotvec);
    return;
  }
  // ----------

  i0=(int)(t0/timestep);
  PhiT=0.;
  h2I=h3I=0.;

  for(i=i0;i<=vlength;i++) {
    /*       If (WCN.eq.0) then
       call NoiseInst (nuvec(i),gimdot0,nmodes,invsqrtS)
        else
        call NoiseFull (nuvec(i),gimdot0,nmodes,invsqrtS)
	endif*/

    hI[i]=0;
    hII[i]=0;

    //     t is actual time (sec) at step i
    t=timestep*(double)i;

    e=evec[i];
    nu=nuvec[i];
    Phi=Phivec[i];
    gim=gimvec[i];
    alp=alpvec[i];
    cosalp=cos(alp);
    sinalp=sin(alp);

    cosqL=cosqK*coslam+sinqK*sinlam*cosalp;
    sinqL=sqrt(1-cosqL*cosqL);
    BB=sinqK*cosphiK*coslam+sinphiK*sinlam*sinalp-cosqK*cosphiK*sinlam*cosalp;
    CC=sinqK*sinphiK*coslam-cosphiK*sinlam*sinalp-cosqK*sinphiK*sinlam*cosalp;
    phiL=ArcT(BB,CC);
    Ldotn=cosqL*cosqS+sinqL*sinqS*cos(phiL-phiS);
    Ldotn2=Ldotn*Ldotn;

    if (mich) {
      orbphs=2.*M_PI*t/year;
      cosorbphs=cos(orbphs-phiS);
      sinorbphs=sin(orbphs-phiS);
      cosq=.5*cosqS-halfsqrt3*sinqS*cosorbphs;
      phiw=orbphs+ArcT(sinqS*sinorbphs,halfsqrt3*cosqS+.5*sinqS*cosorbphs);
      psiup=.5*cosqK-halfsqrt3*sinqK*cos(orbphs-phiK)-cosq*(cosqK*cosqS+sinqK*sinqS*cos(phiK-phiS));
      psidown=.5*sinqK*sinqS*sin(phiK-phiS)-halfsqrt3*cos(orbphs)*(cosqK*sinqS*sin(phiS)-cosqS*sinqK*sin(phiK))-halfsqrt3*sin(orbphs)*(cosqS*sinqK*cos(phiK)-cosqK*sinqS*cos(phiS));
      psi=ArcT(psidown,psiup);
      cosq1=.5*(1+cosq*cosq);
      cos2phi=cos(2.*phiw);
      sin2phi=sin(2.*phiw);
      cos2psi=cos(2.*psi);
      sin2psi=sin(2.*psi);

      FplusI=cosq1*cos2phi*cos2psi-cosq*sin2phi*sin2psi;
      FcrosI=cosq1*cos2phi*sin2psi+cosq*sin2phi*cos2psi;
      FplusII=cosq1*sin2phi*cos2psi+cosq*cos2phi*sin2psi;
      FcrosII=cosq1*sin2phi*sin2psi-cosq*cos2phi*cos2psi;
    } else {
      FplusI=1.;
      FcrosI=0.;
      FplusII=0.;
      FcrosII=1.;
    }
    /*     --------------------------------------------------------
	   Calculate gamma (ext/int mixed) out of gimmel (intrinsic) 
	   -------------------------------------------------------- */

    Sdotn=cosqK*cosqS+sinqK*sinqS*cos(phiK-phiS);
    betaup=-Sdotn+coslam*Ldotn;
    betadown=sinqS*sin(phiK-phiS)*sinlam*cosalp+(cosqK*Sdotn-cosqS)/sinqK*sinlam*sinalp;
    beta=ArcT(betadown,betaup);
    gam2=2.*(gim+beta);
    cos2gam=cos(gam2);
    sin2gam=sin(gam2);
    
    Amp=pow((2.*M_PI*nu*M),2./3.)*zeta;

    for(n=1;n<nmodes+1;n++) {
      fn=n*nu+gimdotvec[i]/M_PI;
      Doppler=2.*M_PI*fn*AUsec*sinqS*cosorbphs;
      if (mich)
	      nPhi=n*Phi+Doppler;
      else
	      nPhi=n*Phi;
      ne=n*e;
      Jm2=Jn(n-2,ne);
      Jm1=Jn(n-1,ne);
      Jm0=Jn(n,ne);
      Jp1=Jn(n+1,ne);
      Jp2=Jn(n+2,ne);   
      a=-n*Amp*(Jm2-2.*e*Jm1+2./n*Jm0+2.*e*Jp1-Jp2)*cos(nPhi);
      b=-n*Amp*sqrt(1-e*e)*(Jm2-2.*Jm0+Jp2)*sin(nPhi);
      c=2.*Amp*Jm0*cos(nPhi);

      Aplus=-(1.+Ldotn2)*(a*cos2gam-b*sin2gam)+c*(1-Ldotn2);
      Acros=2.*Ldotn*(b*cos2gam+a*sin2gam);

      // ROTATION TO NK WAVE FRAME (NOT COMPUTATIONALLY OPTIMISED)
      double rot[4],Aplusold=Aplus,Acrosold=Acros;
      RotCoeff(rot,lam,qS,phiS,qK,phiK,alp);
      Aplus=Aplusold*rot[0]+Acrosold*rot[1];
      Acros=Aplusold*rot[2]+Acrosold*rot[3];

      hnI=prefact*(FplusI *Aplus+FcrosI *Acros);
      hnII=prefact*(FplusII*Aplus+FcrosII*Acros);
      if (n==2)
	h2I+=(hnI*hnI);
      if (n==3)
	h3I+=(hnI*hnI);
      hI[i]=hI[i]+hnI;
      hII[i]=hII[i]+hnII;

      //      hI(i)=hI(i)+hnI*invsqrtS(n)
      //      hII(i)=hII(i)+hnII*invsqrtS(n)
    }
  }
  //cout << "Ratio of squared integral of second and third harmonic is " << h2I/h3I << endl;
  double afact,bfact,cfact,n2fact,n3fact;
  n=2;
  ne=n*e;
  Jm2=Jn(n-2,ne);
  Jm1=Jn(n-1,ne);
  Jm0=Jn(n,ne);
  Jp1=Jn(n+1,ne);
  Jp2=Jn(n+2,ne);   
  a=-n*(Jm2-2.*e*Jm1+2./n*Jm0+2.*e*Jp1-Jp2);
  b=-n*sqrt(1-e*e)*(Jm2-2.*Jm0+Jp2);
  c=2.*Jm0;
  afact=FplusI *(-(1.+Ldotn2)*(cos2gam))+FcrosI *(2.*Ldotn*(sin2gam));
  bfact=FplusI *(-(1.+Ldotn2)*(-sin2gam))+FcrosI *(2.*Ldotn*(cos2gam));
  cfact=FplusI *((1-Ldotn2));
  n2fact=(afact*a+cfact*c)*(afact*a+cfact*c)+bfact*b*bfact*b;

  n=3;
  ne=n*e;
  Jm2=Jn(n-2,ne);
  Jm1=Jn(n-1,ne);
  Jm0=Jn(n,ne);
  Jp1=Jn(n+1,ne);
  Jp2=Jn(n+2,ne);   
  a=-n*(Jm2-2.*e*Jm1+2./n*Jm0+2.*e*Jp1-Jp2);
  b=-n*sqrt(1-e*e)*(Jm2-2.*Jm0+Jp2);
  c=2.*Jm0;
  afact=FplusI *(-(1.+Ldotn2)*(cos2gam))+FcrosI *(2.*Ldotn*(sin2gam));
  bfact=FplusI *(-(1.+Ldotn2)*(-sin2gam))+FcrosI *(2.*Ldotn*(cos2gam));
  cfact=FplusI *((1-Ldotn2));
  n3fact=(afact*a+cfact*c)*(afact*a+cfact*c)+bfact*b*bfact*b;
  //cout << "Theoretical prediction: " << n2fact/n3fact << endl;

  PhiT=0;
  
  for(i=i0-1;i>=0;i--) {

    /* If (WCN.eq.0) then
       call NoiseInst (nuvec(i),gimdot0,nmodes,invsqrtS)
       else
       call NoiseFull (nuvec(i),gimdot0,nmodes,invsqrtS)
       endif*/

    hI[i]=0.;
    hII[i]=0.;

    //     t is actual time (sec) at step i
    t=timestep*(double)i;

    e=evec[i];
    nu=nuvec[i];
    Phi=Phivec[i];
    gim=gimvec[i];
    alp=alpvec[i];
    cosalp=cos(alp);
    sinalp=sin(alp);

    cosqL=cosqK*coslam+sinqK*sinlam*cosalp;
    sinqL=sqrt(1-cosqL*cosqL);
    BB=sinqK*cosphiK*coslam+sinphiK*sinlam*sinalp-cosqK*cosphiK*sinlam*cosalp;
    CC=sinqK*sinphiK*coslam-cosphiK*sinlam*sinalp-cosqK*sinphiK*sinlam*cosalp;
    phiL=ArcT(BB,CC);
    Ldotn=cosqL*cosqS+sinqL*sinqS*cos(phiL-phiS);
    Ldotn2=Ldotn*Ldotn;

    if (mich) {
      orbphs=2.*M_PI*t/year;
      cosorbphs=cos(orbphs-phiS);
      sinorbphs=sin(orbphs-phiS);
      cosq=.5*cosqS-halfsqrt3*sinqS*cosorbphs;
      phiw=orbphs+ArcT(sinqS*sinorbphs,halfsqrt3*cosqS+.5*sinqS*cosorbphs);
      psiup=.5*cosqK-halfsqrt3*sinqK*cos(orbphs-phiK)-cosq*(cosqK*cosqS+sinqK*sinqS*cos(phiK-phiS));
      psidown=.5*sinqK*sinqS*sin(phiK-phiS)-halfsqrt3*cos(orbphs)*(cosqK*sinqS*sin(phiS)-cosqS*sinqK*sin(phiK))-halfsqrt3*sin(orbphs)*(cosqS*sinqK*cos(phiK)-cosqK*sinqS*cos(phiS));
      psi=ArcT(psidown,psiup);
      cosq1=.5*(1+cosq*cosq);
      cos2phi=cos(2.*phiw);
      sin2phi=sin(2.*phiw);
      cos2psi=cos(2.*psi);
      sin2psi=sin(2.*psi);

      FplusI=cosq1*cos2phi*cos2psi-cosq*sin2phi*sin2psi;
      FcrosI=cosq1*cos2phi*sin2psi+cosq*sin2phi*cos2psi;
      FplusII=cosq1*sin2phi*cos2psi+cosq*cos2phi*sin2psi;
      FcrosII=cosq1*sin2phi*sin2psi-cosq*cos2phi*cos2psi;
    } else {
      FplusI=1.;
      FcrosI=0.;
      FplusII=0.;
      FcrosII=1.;
    }
      
    /*     --------------------------------------------------------
	   Calculate gamma (ext/int mixed) out of gimmel (intrinsic)
	   -------------------------------------------------------- */

    Sdotn=cosqK*cosqS+sinqK*sinqS*cos(phiK-phiS);
    betaup=-Sdotn+coslam*Ldotn;
    betadown=sinqS*sin(phiK-phiS)*sinlam*cosalp+(cosqK*Sdotn-cosqS)/sinqK*sinlam*sinalp;
    beta=ArcT(betadown,betaup);
    gam2=2.*(gim+beta);
    cos2gam=cos(gam2);
    sin2gam=sin(gam2);

    Amp=pow(2.*M_PI*nu*M,2./3.)*zeta;

    for(n=1;n<nmodes+1;n++) {
      fn=n*nu+gimdotvec[i]/M_PI;
      Doppler=2.*M_PI*fn*AUsec*sinqS*cosorbphs;
      if (mich)
	      nPhi=n*Phi+Doppler;
      else
	      nPhi=n*Phi;
      ne=n*e;
      Jm2=Jn(n-2,ne);
      Jm1=Jn(n-1,ne);
      Jm0=Jn(n,ne);
      Jp1=Jn(n+1,ne);
      Jp2=Jn(n+2,ne);
      a=-n*Amp*(Jm2-2.*e*Jm1+2./n*Jm0+2.*e*Jp1-Jp2)*cos(nPhi);
      b=-n*Amp*sqrt(1-e*e)*(Jm2-2.*Jm0+Jp2)*sin(nPhi);
      c=2.*Amp*Jm0*cos(nPhi);

      Aplus=-(1.+Ldotn2)*(a*cos2gam-b*sin2gam)+c*(1-Ldotn2);
      Acros=2.*Ldotn*(b*cos2gam+a*sin2gam);

      // ROTATION TO NK WAVE FRAME (NOT COMPUTATIONALLY OPTIMISED)
      double rot[4],Aplusold=Aplus,Acrosold=Acros;
      RotCoeff(rot,lam,qS,phiS,qK,phiK,alp);
      Aplus=Aplusold*rot[0]+Acrosold*rot[1];
      Acros=Aplusold*rot[2]+Acrosold*rot[3];

      hnI=prefact*(FplusI *Aplus+FcrosI *Acros);
      hnII=prefact*(FplusII*Aplus+FcrosII*Acros);

      hI[i]=hI[i]+hnI;
      hII[i]=hII[i]+hnII;
    }
  }

  free(evec);
  free(nuvec);
  free(alpvec);
  free(Phivec);
  free(gimvec);
  free(gimdotvec);
  return;
}


void GenBCWave(double *hI, double *hII, double deltat, int vlength, double e0, double nu0, double M, double mu, double S, double dist, double inc, double gam0, double Phi0,
	       double qS, double phiS, double alp0, double qK, double phiK, bool mich, bool traj) {
  int i,n,nmodes,p,pp,m,k,iv;
  int dim,direction,comp,pointM,pointe0;
  int timestep1,nsteps1;
  double par[14],dt,t,spinres,dirres;

  double bigest,count,sqrtdet;
  double integral,StoN,eisco,nuisco,gimdot0;
  double delta,deltahalf,tend,q[14],element;
  double parplus[14],parminus[14],det;
  double normI,normII,norm,dtdays;
  double t0res,mures,Mres,e0res,gim0res,Phi0res;
  double lamres,alp0res,OmegaSres,Dres,zeta,invD,Gpsc;

  /* -------------------------------------------------------------
     SETUP OF PHYSICAL AND NUMERICAL PARAMETERS
     -------------------------------------------------------------*/

  dt=deltat*((double)vlength);
  tend=dt;

  /* Set parameters for waveform. */
  //     t0
  par[1]=0.;
  //     mu
  par[2]=mu*SOLARMASSINSEC;
  //     M
  par[3]=M*SOLARMASSINSEC;
  //     e0
  par[4]=e0;
  //     S/M^2
  par[11]=S;
  //     inc
  par[9]=inc;

  invD=1./(dist*Gpc);
  zeta=par[2]*invD;

  //     gam0
  par[5]=gam0;
  //     Phi0
  par[6]=Phi0;
  //     qS
  par[7]=qS;
  //     phiS
  par[8]=phiS;
  //     alp0
  par[10]=alp0;
  //     qK
  par[12]=qK;
  //     phiK
  par[13]=phiK;

  /* -------------------------------------------------------------
     NUMBER OF MODES TO BE CONSIDERED (as a function of e0)
     ------------------------------------------------------------- */

  /*     This is an experimental result: See Peters & Mathews fig 3
	 and Eq. (20), and see ModePower.nb */
  nmodes=(int)(30*par[4]);
  if (par[4] < 0.135) nmodes=4;


  /* --------------------------------------------------------------
     CALCULATING WAVEFORM & DERIVATIVES
     -------------------------------------------------------------- */

  //printf("Calculating waveform \n");
  waveform(tend,par,nu0,vlength,deltat,hI,hII,nmodes,zeta,0,mich,traj);
  return;
}

