#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Globals.h"
#include "NRUtil.h"
#include "NRCosmology.h"
#include "Cosmology.h"

#include "Constants.h"

const double CURVTOL = 1e-16;
const double SOLVETOL = 1e-10;
const double ZTOL = 1e-10;

Cosmology::Cosmology(Real hubble, Real OmegaM, Real OmegaL) : H0(hubble),Om(OmegaM),OL(OmegaL)
{
  Ok=1.-Om-OL;
  q0=0.5*Om-OL;
  sigma0=0.5*Om;

  if (fabs(Ok) < CURVTOL)
    curv=0;
  else {
    if (Ok > 0.)
      curv=-1;
    else
      curv=1;
  }
  if (Ok)
    a0=sqrt(1./fabs(Ok))/H0;
  else
    a0=1.;
  rhocrit=(3.*H0*H0/(8.*M_PI))/(1.e6*pc*pc*G);

  Real sgn,mag;
  Real m,nsq,a2,q,r,discrim;

  a2=Ok/Om;
  q=-a2*a2/9.;
  r=-a2*a2*a2/27.-0.5*OL/Om;
  discrim=sqrt(q*q*q+r*r);

  uc=a2/3.;
  sgn=mag=0.;
  if (r+discrim > 0.) {
    sgn=1.;
    mag=pow(r+discrim,1./3.);
  }
  else if (r+discrim < 0.) {
    sgn=-1.;
    mag=pow(-r-discrim,1./3.);
  }
  uc-=sgn*mag;
  sgn=mag=0.;
  if (r-discrim > 0.) {
    sgn=1.;
    mag=pow(r-discrim,1./3.);
  }
  else if (r-discrim < 0.) {
    sgn=-1.;
    mag=pow(discrim-r,1./3.);
  }
  uc-=sgn*mag;

  m=0.5*(uc-Ok/Om);
  nsq=OL/(uc*Om)-m*m;
  p=sqrt((m+uc)*(m+uc)+nsq);
  k=sqrt((p+m+uc)/(2.*p));

  zeropt=ellf(2.*atan(sqrt((1.+uc)/p)),k);
}

Cosmology::~Cosmology()
{

}

Real Cosmology::ProperD(Real z)
{
  Real dist,chidiff;

  chidiff=(ellf(2.*atan(sqrt((1.+z+uc)/p)),k)-zeropt)/sqrt(p*Om);

  if (curv==-1)
    dist=a0*sinh(chidiff/sqrt(Ok));
  else if (curv == 1)
    dist=a0*sin(chidiff/sqrt(-Ok));
  else
    dist=chidiff/H0;

  // Distance is expressed in Mpc. This assumes that on input, H0 is given in km/s/Mpc.

  dist*=C/1000.;

  return dist;
}

Real Cosmology::LumD(Real z)
{
  Real dist;

  dist=ProperD(z);
  dist*=(1.+z);

  return dist;
}

Real Cosmology::AngD(Real z)
{
  Real dist;

  dist=ProperD(z);
  dist/=(1.+z);

  return dist;
}

Real Cosmology::zofDp(Real d)
{
  Real arg,phi,tn,tst;
  Real *sn,*cn,*dn;
  sn=(Real *)malloc(sizeof(Real));
  cn=(Real *)malloc(sizeof(Real));
  dn=(Real *)malloc(sizeof(Real));

  if (curv==0)
    arg=1000.*d*H0*sqrt(p*Om)/C;
  else if (curv==-1)
    arg=sqrt(p*Om/Ok)*asinh(H0*1000.*d*sqrt(Ok)/C);
  else
    arg=sqrt(p*Om/Ok)*asin(H0*1000.*d*sqrt(-Ok)/C);
  arg+=zeropt;

  sncndn(arg,1.-k*k,sn,cn,dn);

  tst=arg-ellf(0.5*M_PI,k);
  if (tst > 0.)
    phi=-asin(*sn)+M_PI;
  else
    phi=asin(*sn);
  tn=tan(0.5*phi);

  return tn*tn*p-uc-1;
}

Real Cosmology::zofDL(Real d)
{
  Real x,xlow,xhigh,zhigh,zlow,z;
  bool loop=true;
  xhigh=d;
  zhigh=zofDp(d);
  xlow=d/(1.+zhigh);
  zlow=zofDp(xlow);
  while(loop) {
    x=0.5*(xlow+xhigh);
    z=zofDp(x);
    if ((1.+z)*x > d) {
      xhigh=x;
      zhigh=z;
    }
    else { 
      xlow=x;
      zlow=z;
    }
    if (fabs(zhigh-zlow)<ZTOL)
      loop=false;
  }

  return 0.5*(zhigh+zlow);
}

Real Cosmology::zofMaxDA()
{
  Real xlow,xmid,xhigh,dlow,dmid,dhigh;
  Real OFF=0.1;

  bool loop=true;

  xlow=0.;
  xmid=OFF;
  xhigh=2.*OFF;
  dlow=AngD(xlow);
  dhigh=AngD(xhigh);
  dmid=AngD(xmid);

  while(loop) {
    if ((dhigh-dmid)*(dmid-dlow) > 0.) {
      if (dhigh-dmid > 0.) {
	xlow=xmid;
	dlow=dmid;
	xmid=xhigh;
	dmid=dhigh;
	xhigh=xmid+OFF;
	dhigh=AngD(xhigh);
      }
      else {
	xhigh=xmid;
	dhigh=dmid;
	xmid=xlow;
	dmid=dlow;
	xlow=xmid-OFF;
	dlow=AngD(xhigh);
      }
    }
    else {
      OFF/=2.;
      xlow=xmid-OFF;
      dlow=AngD(xlow);
      xhigh=xmid+OFF;
      dhigh=AngD(xhigh);
    }
    if (xhigh-xlow<ZTOL)
      loop=false;
  }

  return 0.5*(xhigh+xlow);
}

Real Cosmology::zofDA(Real d, int rt)
{
  Real maxdist,zmax,zlow,dlow,zhigh,dhigh,sgn,zmid,dmid;
  bool loop;
  zmax=zofMaxDA();
  maxdist=AngD(zmax);

  if (d > maxdist) {
    cerr << "Angular diameter distance never gets that large!" << endl;
    return 0.;
  }
  else {
    if (rt == 1) {
      zlow=0.;
      dlow=0.;
      zhigh=zmax;
      dhigh=maxdist;
      sgn=-1.;
    }
    else {
      zlow=zmax;
      dlow=maxdist;
      zhigh=zmax;
      dhigh=maxdist;
      while (dhigh > d) {
	zhigh*=2.;
	dhigh=AngD(zhigh);
      }
      sgn=1.;
    }

    loop=true;
    while (loop) {
      zmid=0.5*(zlow+zhigh);
      dmid=AngD(zmid);
      if (dmid*sgn < sgn*d) {
	zhigh=zmid;
	dhigh=dmid;
      }
      else {
	zlow=zmid;
	dlow=dmid;
      }
      if (zhigh-zlow < ZTOL)
	loop=false;
    }  
    return 0.5*(zlow+zhigh);
  }
}

void Cosmology::SolveCubic(const Real p0, const Real p1, const Real p2, Real *roots, bool &cmplex)
{
  Real b0,b1,b2,discrim,tmp;
  Real *rts;
  rts=Realvector(1,3);
  int ct,cnt;
  if (p0 == 0.) {
    discrim=p2*p2-4.*p1;
    if (discrim >= 0.) {
      cmplex=false;
      roots[1]=0.;
      roots[2]=0.5*(-p2+sqrt(discrim));
      roots[3]=0.5*(-p2-sqrt(discrim));
    }
    else {
      cmplex=true;
      roots[1]=0.;
      roots[2]=-0.5*p2;
      roots[3]=0.5*sqrt(-1.*discrim);
    }
  }
  else {
    CubicSolver(p0,p1,p2,roots,cmplex);
    if (!cmplex) {
      b0=roots[1]*roots[2]*roots[3]+p0;
      b1=roots[1]*(roots[2]+roots[3])+roots[2]*roots[3]-p1;
      b2=roots[1]+roots[2]+roots[3]+p2;
    }
    else {
      b0=roots[1]*(roots[2]*roots[2]+roots[3]*roots[3])+p0;
      b1=2.*roots[1]*roots[2]+roots[2]*roots[2]+roots[3]*roots[3]-p1;
      b2=roots[1]+2.*roots[2]+p2;
    }
    if ((fabs(b0/p0) > SOLVETOL) || (fabs(b1/p1) > SOLVETOL) || (fabs(b2/p2) > SOLVETOL)) {
      b0=1./p0;
      b1=p2/p0;
      b2=p1/p0;
      CubicSolver(b0,b1,b2,rts,cmplex);
      if (!cmplex) {
	for (ct=1;ct<=3;ct++) 
	  roots[ct]=1./rts[ct];
	b0=roots[1]*roots[2]*roots[3]+p0;
	b1=roots[1]*(roots[2]+roots[3])+roots[2]*roots[3]-p1;
	b2=roots[1]+roots[2]+roots[3]+p2;
      }
      else {
	roots[1]=1./rts[1];
	roots[2]=rts[2]/(rts[2]*rts[2]+rts[3]*rts[3]);
	roots[3]=rts[3]/(rts[2]*rts[2]+rts[3]*rts[3]);
	b0=roots[1]*(roots[2]*roots[2]+roots[3]*roots[3])+p0;
	b1=2.*roots[1]*roots[2]+roots[2]*roots[2]+roots[3]*roots[3]-p1;
	b2=roots[1]+2.*roots[2]+p2;
      }
      if (fabs(b0/p0)>SOLVETOL) {
	cerr << "Possible problem in cubic equation solver, difference between -p0 and product of roots is " 
	     << b0 << endl;
	Die("Dying!");
      }
      if (fabs(b1/p1)>SOLVETOL) {
	cerr << "Possible problem in cubic equation solver, difference between p1 and sum of products of pairs of roots is " << b1 << endl;
	Die("Dying!");
      }
      if (fabs(b2/p2)>SOLVETOL) {
	cerr << "Possible problem in cubic equation solver, difference between -p2 and sum of roots is " 
	     << b2 << endl;
	Die("Dying!");
      }
    }
  } 
  if (!cmplex) {
    for (ct=1;ct<=3;ct++) {
      for (cnt=ct+1;cnt<=3;cnt++) {
	if (roots[cnt] > roots[ct]) {
	  tmp=roots[cnt];
	  roots[cnt]=roots[ct];
	  roots[ct]=tmp;
	}
      }
    }
  }
  return;
}

void Cosmology::CubicSolver(const Real p0, const Real p1, const Real p2, Real *roots, bool &cmplex)
{
  Real q=p1/3.-p2*p2/9.;
  Real rc=(p1*p2-3.*p0)/6.-p2*p2*p2/27.;
  Real tst=q*q*q+rc*rc;
  Real s1,s2,mag,arg;
  if (tst > 0.) {
    cmplex=true;
    s1=pow(rc+sqrt(tst),1./3.);
    s2=pow(rc-sqrt(tst),1./3.);
    roots[1]=s1+s2-p2/3.;
    roots[2]=-(s1+s2)/2.-p2/3.;
    roots[3]=sqrt(3.)*(s1-s2)/2.;
  }
  else {
    cmplex=false;
    mag=sqrt(-1.*q);
    if (rc > 0. )
      arg=atan(sqrt(-tst)/rc)/3.;
    else
      arg=(2.*asin(1.)-atan(-1.*sqrt(-tst)/rc))/3.;
    roots[1]=2.*mag*cos(arg)-p2/3.;
    roots[2]=-mag*cos(arg)-p2/3.+sqrt(3.)*mag*sin(arg);
    roots[3]=-mag*cos(arg)-p2/3.-sqrt(3.)*mag*sin(arg);
  }
  return;
}

void Cosmology::SolveQuartic(const Real p0, const Real p1, const Real p2, const Real p3, Real *roots, int &cmplex)
{
  roots[4]=roots[1]=roots[2]=roots[3]=0.;
  Real *cub,*poss,cub1,cub2,cub3;
  Real r1,r2,r3,r4;
  cub=Realvector(1,3);
  poss=Realvector(1,3);
  bool cubic;
  int ct,numposs=0,cnt;

  if ((p0==0.)&&(p1==0.)&&(p2==0.)) {
    cmplex=0;
    roots[1]=-p3;
    roots[2]=0.;
    roots[3]=0.;
    roots[4]=0.;
    return;
  }

  if (p0 == 0.) {
    SolveCubic(p1,p2,p3,cub,cubic);
    if (cubic) {
      roots[3]=cub[2];
      roots[4]=cub[3];
      if (cub[1]>0.) {
	roots[1]=cub[1];
	roots[2]=0.;
      }
      else {
        roots[2]=cub[1];
        roots[1]=0.;
      }
    }
    else {
      for (ct=1;ct<=3;ct++) {
       if (cub[ct] > 0.)
         roots[ct]=cub[ct];
       else
         roots[ct+1]=cub[ct];
      }
    }
  }
  else {
  SolveCubic(4.*p0*p2-p0*p3*p3-p1*p1,p1*p3-4.*p0,-p2,cub,cubic);
  if (!cubic) {
    for (ct=1;ct<=3;ct++) {
      if ((0.25*p3*p3-p2+cub[ct] >=0.)&&(0.25*cub[ct]*cub[ct]-p0 >= 0.)) {
	numposs++;
	poss[numposs]=cub[ct];
      }
    }
    for (ct=1;ct<=numposs;ct++) {
      if (fabs(0.25*p3*poss[ct]-sqrt((0.25*p3*p3+poss[ct]-p2)*
				     (0.25*poss[ct]*poss[ct]-p0))- 0.5*p1) < 1.e-10)
	cub1=poss[ct];
    }
  }
  else {
    cub1=cub[0];
  }
  Real a,b,c,d,mag,arg,discrim;
  if (p3*p3/4. +cub1-p2 > 0.) {
    a=p3/2.-sqrt(p3*p3/4. +cub1-p2);
    b=0.;
  }
  else {
    a=p3/2.;
    b=-sqrt(p2-p3*p3/4. -cub1);
  }
  if (cub1*cub1/4.-p0 > 0.) {
    c=cub1/2.-sqrt(cub1*cub1/4.-p0);
    d=0.;
  }
  else {
    c=cub1/2.;
    d=-sqrt(p0-cub1*cub1/4.);    
  }
  if (fabs(b)+fabs(d) > 0.) {
    cmplex=2;
    mag=sqrt(sqrt(raise(a*a-4.*c-b,2)+raise(2.*a*b-4.*d,2)))/2.;
    arg=0.5*atan((2.*a*b-4.*d)/(a*a-4.*c-b*b));
    r4=-0.5*b+mag*sin(arg);
    if (r4 != 0.) {
      r3=-0.5*a+mag*cos(arg);
      r1=-0.5*a-mag*cos(arg);
      r2=-0.5*b-mag*sin(arg);
      if (r2 == 0.) {
        r2=r1;
	cmplex--;
      }
    }
    else {
      cmplex--;
      r1=-0.5*a+mag*cos(arg);
      r2=r1;
      r3=-0.5*a-mag*cos(arg);
      r4=-0.5*b-mag*sin(arg);
    }
  }
  else {
    cmplex=0;
    int store=0;
    discrim=a*a-4.*c;
    if (discrim > 0.) {
      r1=-0.5*a+0.5*sqrt(discrim);
      r2=-0.5*a-0.5*sqrt(discrim);
    }
    else {
      cmplex++;
      r3=-0.5*a;
      r4=0.5*sqrt(-1.*discrim);
      store=1;
    }
    a=p3/2.+sqrt(p3*p3/4. +cub1-p2);
    c=cub1/2.+sqrt(cub1*cub1/4.-p0);
    discrim=a*a-4.*c;
    if (discrim > 0.) {
      cub2=-0.5*a+0.5*sqrt(discrim);
      cub3=-0.5*a-0.5*sqrt(discrim);
    }
    else {
      cmplex++;
      cub2=-0.5*a;
      cub3=0.5*sqrt(-1.*discrim);
    }
    if (store) {
      r1=cub2;
      r2=cub3;
    }
    else {
      r3=cub2;
      r4=cub3;
    }
  }
  Real tmp;
  
  if (cmplex == 1) {
    if (r1 < r2) {
      tmp=r1;
      r1=r2;
      r2=tmp;
    }
  }
  else if (cmplex == 0) {
    roots[1]=r1;
    roots[2]=r2;
    roots[3]=r3;
    roots[4]=r4;
    for (ct=1;ct<=4;ct++) {
      for (cnt=ct+1;cnt<=4;cnt++) {
	if (roots[cnt] > roots[ct]) {
	  tmp=roots[cnt];
	  roots[cnt]=roots[ct];
	  roots[ct]=tmp;
	}
      }
    }
  }
  Real b0,b1,b2,b3;
  if (cmplex==2) {
    b0=(roots[1]*roots[1]+roots[2]*roots[2])*(roots[3]*roots[3]+roots[4]*roots[4])-p0;
    b1=(roots[1]*roots[1]+roots[2]*roots[2])*2.*roots[3]+2.*roots[1]*(roots[3]*roots[3]+roots[4]*roots[4])+p1;
    b2=roots[1]*roots[1]+roots[2]*roots[2]+4.*roots[1]*roots[3]+roots[3]*roots[3]+roots[4]*roots[4]-p2;
    b3=2.*roots[1]+2.*roots[3]+p3;
  }
  else if (cmplex==1) {
    b0=roots[1]*roots[2]*(roots[3]*roots[3]+roots[4]*roots[4])-p0;
    b1=roots[1]*(roots[2]*2.*roots[3]+roots[3]*roots[3]+roots[4]*roots[4])+roots[2]*
      (roots[3]*roots[3]+roots[4]*roots[4])+p1;
    b2=roots[1]*(roots[2]+2.*roots[3])+2.*roots[2]*roots[3]+roots[3]*roots[3]+roots[4]*roots[4]-p2;
    b3=roots[1]+roots[2]+2.*roots[3]+p3;
  }
  else {
    b0=roots[1]*roots[2]*roots[3]*roots[4]-p0;
    b1=roots[1]*(roots[2]*(roots[3]+roots[4])+roots[3]*roots[4])+roots[2]*roots[3]*roots[4]+p1;
    b2=roots[1]*(roots[2]+roots[3]+roots[4])+roots[2]*(roots[3]+roots[4])+roots[3]*roots[4]-p2;
    b3=roots[1]+roots[2]+roots[3]+roots[4]+p3;
  }
  if (fabs(b0/p0)>SOLVETOL) {
    cerr << "Possible problem in quartic equation solver, difference between p0 and product of roots is " << b0 << endl;
    Die("Dying!");
  }
  if (fabs(b1/p1)>SOLVETOL) {
    cerr << "Possible problem in quartic equation solver, difference between -p1 and sum of products of roots in threes is " 
	 << b1 << endl;
   Die("Dying!");       ;
  }
  if (fabs(b2/p2)>SOLVETOL) {
    cerr << "Possible problem in quartic equation solver, difference between p2 and sum of products of pairs of roots is " << b2 << endl;
    Die("Dying!");
  }
  if (fabs(b3/p3)>SOLVETOL) {
    cerr << "Possible problem in quartic equation solver, difference between -p3 and sum of roots is " 
	 << b3 << endl;
    Die("Dying!");
  }
  }
  return;
}

Real Cosmology::raise(Real x, int n)
{
   Real tmp=1.;
   int i;
   for (i=0;i<fabs(n);i++)
	tmp*=x;
   if (n > 0)
	return tmp;
   else
	return 1./tmp;
}
