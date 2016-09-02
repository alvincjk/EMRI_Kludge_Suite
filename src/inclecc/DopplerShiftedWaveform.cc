#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Globals.h"
#include "CKG.h"
#include "IEKG.h"
#include "NRUtil.h"
#include "GKG.h"
#include "NRLB.h"
#include "SWSH.h"
#include "GKTraj.h"
#include "GKInsp.h"
#include "DopplerShiftedWaveform.h"

#define degrees         1.74532925199e-2           /* radians per degree */
#define SOLARMASSINSEC (G*Msun/(C*C*C)) /* 1 solar mass in units of seconds. */
#define DELTAT_TRAJ (0.001) /* sampling interval for phase-space
                             trajectory generation */

/* Note that we assume that SYSTEM_M is given in solar mass units and the 
   trajectory sampling interval dt is given in seconds. The factor 
   SOLARMASSINSEC converts from solar masses to seconds where necessary. */

DopplerShiftedWaveform::DopplerShiftedWaveform(const Real *intrinsparams, const Real *extrinsparams, double deltat, const int Nfrst, 
					       const int Nlast , const Real tzero, const Real hubb, const Real omM, const Real omL) : 
  pI(intrinsparams),pE(extrinsparams),dt(deltat),Nzero(Nfrst),Nmx(Nlast),t0(tzero),H0(hubb),OmegaM(omM),OmegaL(omL)
{
  int i,j;
  Real numfact,dtTraj,doppcs,doppsin,tSSB0;

  rtthr=sqrt(3.);
  // No factor of 0.5, since I define my Ax and A+ in a slightly different way to C+L
  numfact=rtthr*pow(G*Msun/(C*C),3)/(1.e6*C*C*pc);

  norm=numfact*exp(pI[LN_ORBIT_ETA])*pow(exp(pI[LN_SYSTEM_M]),3);

  sinthsky=sqrt(1.-pE[COS_THETA_SKY]*pE[COS_THETA_SKY]);
  csthsky=pE[COS_THETA_SKY];
  phik=pE[PHI_SPIN];
  phis=pE[PHI_SKY];
  cphik=cos(phik);
  cphis=cos(phis);
  sphik=sin(phik);
  sphis=sin(phis);
  csthk=pE[COS_THETA_SPIN];
  sinthk=sqrt(1.-csthk*csthk);

  /* Compute inclination of black-hole spin to line of sight as dot product of
     source position and spin direction in ecliptic coordinate system. */

  csinc=sinthsky*sinthk*(cphis*cphik+sphis*sphik)+csthsky*csthk;

  dtTraj=dt;

  /* NB This isn't the right thing to do, since we must only redshift the source time spacing, and not the LISA 
        time spacing. Can correct this is we want, but easier to ensure we always provide the redshifted mass and 
	never rescale time. */
  if (pE[REDSHIFT]) {
    LCDM=new Cosmology(H0,OmegaM,OmegaL);
    //if (!UseRedShiftedMass)
    //  dtTraj/=(1.+pE[REDSHIFT]);
  }

  doppcs=AU*sinthsky*cphis/C;
  doppsin=AU*sinthsky*sphis/C;

  tSSB0=t0+cos(2.*M_PI*t0/YEAR)*doppcs+sin(2.*M_PI*t0/YEAR)*doppsin;

  /* Generate inspiral, at constant intervals of detector time, which are not constant intervals of SSB time. */

  int trajlen;
  if (Nzero >= 0)
    trajlen=Nmx;
  else
    trajlen=Nmx-Nzero;
  xx=(Real *)malloc(trajlen*sizeof(Real));
  yy=(Real *)malloc(trajlen*sizeof(Real));
  zz=(Real *)malloc(trajlen*sizeof(Real));
  
  if (!xx) {
    cerr << "Memory allocation error (xx) in DopplerShiftedWaveform constructor." << endl;
    exit(-1);
  }
  if (!yy) {
    cerr << "Memory allocation error (yy) in DopplerShiftedWaveform constructor." << endl;
    exit(-1);
  }
  if (!zz) {
    cerr << "Memory allocation error (zz) in DopplerShiftedWaveform constructor." << endl;
    exit(-1);
  }
  Qdmom=(Real **)malloc(5*sizeof(Real *));
  if (!Qdmom) {
    cerr << "Memory allocation error (Qdmom) in Inspiral constructor." << endl;
    exit(-1);
  }
  
  for (i=0;i<5;i++) {
    Qdmom[i]=(Real *)malloc(trajlen*sizeof(Real));
    if (!Qdmom[i]) {
      cerr << "Memory allocation error (Qdmom[" << i << "]) in Inspiral constructor." << endl;
      exit(-1);
    }
  }

  if (Nzero >= 0) {
    x=&xx[0];
    y=&yy[0];
    z=&zz[0];
    for (i=0;i<5;i++)
      Quad[i]=&(Qdmom[i][0]);
  }
  else {
    x=&xx[-Nzero];
    y=&yy[-Nzero];
    z=&zz[-Nzero];
    for (i=0;i<5;i++)
      Quad[i]=&(Qdmom[i][-Nzero]);
  }
  
  int n = 0, nMax; /* loop index and limit */
  int Ntraj_max;   /* maximum points in phase-space trajectories */
  int Necc = 0;    /* actual number of eccentric points in above */
  int Ncirc = 0;   /* actual number of circular points in above */
  int ecc_traj, circ_traj; /* return codes from phase-space generator */
  TrajData *EccTrajBack = NULL;  /* eccentric phase-space trajectory data */
  TrajData *CircTrajBack = NULL; /* circular phase-space trajectory data */
  Real t, tLISA, dtSSB, r, ct, st;         /* integration temporary variables */
  Real Iab[6]; /* Temporary quadrupole moment tensor. */

  int iup,ival;
  bool notfound;
  Real frac;

  if (Nzero < 0) {  
    GKTraj gktrajback( cos( pI[ORBIT_I] ), pI[SYSTEM_A] );
    gktrajback.p = pI[ORBIT_P];
    gktrajback.ecc = pI[ORBIT_E];
    Ntraj_max = (int)(((1 -Nzero)*dtTraj+500.)*exp(pI[LN_ORBIT_ETA])
		       /(SOLARMASSINSEC*exp(pI[LN_SYSTEM_M]))/DELTAT_TRAJ ) + 2;
    EccTrajBack = (TrajData *)
      malloc( (size_t)( ( Ntraj_max + 1 )*sizeof(TrajData) ) );
    CircTrajBack = (TrajData *)
      malloc( (size_t)( ( Ntraj_max + 1 )*sizeof(TrajData) ) );
    if ( !EccTrajBack || !CircTrajBack ) {
      fprintf( stderr, "Error: InspiralTrajectory: Out of memory!\n" );
      if ( EccTrajBack ) free ( EccTrajBack );
      if ( CircTrajBack ) free ( CircTrajBack );
      length=0;
    }
    ecc_traj = gktrajback.Eccentric( -DELTAT_TRAJ, EccTrajBack, Ntraj_max, Necc );
    if ( ecc_traj == 1 )
      circ_traj = gktrajback.Circular( DELTAT_TRAJ, CircTrajBack, Ntraj_max, Ncirc );
    GKInsp gkinspback( EccTrajBack,Necc,CircTrajBack,Ncirc,SOLARMASSINSEC*exp(pI[LN_SYSTEM_M]),
		       exp(-1.*pI[LN_ORBIT_ETA]),pI[SYSTEM_A],-dtTraj,pI[ORBIT_P],1.0 );
    
    gkinspback.x[1] = pI[ORBIT_PSI0];

    if ((pI[ORBIT_I]*(M_PI-pI[ORBIT_I])) >= 0.)
      gkinspback.x[2] = pI[ORBIT_CHI0];
    else
      gkinspback.x[2] = M_PI+pI[ORBIT_CHI0];

    gkinspback.x[3] = 0.;
    gkinspback.x[4] = 0.;

    if ( gkinspback.Necc > 0 ) {
      gkinspback.DoingEcc = 1;
      gkinspback.ihi = 2;
      nMax = ( (int)( ( gkinspback.tecc_fin - gkinspback.tecc_init )
		      / gkinspback.deltat_est ) );

      if ( nMax > 1-Nzero )
	nMax = 1-Nzero;

      for ( n = 0; n < 1-Nzero; n++ ) {
	tLISA = t0+gkinspback.tecc_init - n*dtTraj;
	t=tLISA+cos(2.*M_PI*tLISA/YEAR)*doppcs+sin(2.*M_PI*tLISA/YEAR)*doppsin-tSSB0;
	dtSSB=dtTraj*(1.-2.*M_PI*(sin(2.*M_PI*tLISA/YEAR)*doppcs-cos(2.*M_PI*tLISA/YEAR)*doppsin)/YEAR);
	r = gkinspback.rfunc( t, gkinspback.x[1] );
	ct = gkinspback.costhetafunc( t, gkinspback.x[2] );
	st = sqrt( 1.0 - ct*ct );
	x[-n] = r*st*cos( gkinspback.x[3] ); 
	y[-n] = r*st*sin( gkinspback.x[3] );
	z[-n] = r*ct;
	gkinspback.GKInsp_QuadrupoleMoments(t, gkinspback.x, Iab);
	for (i=0;i<5;i++)
	  Quad[i][-n]=Iab[i+1];
	gkinspback.TakeAStep( gkinspback.x, t, t - dtSSB );
      }
    }
    else {
      gkinspback.DoingEcc = 0;
      gkinspback.ihi = 2;

      
      /* Continue with the integration. */
      for (n=0 ; n < 1-Nzero; n++ ) {
	tLISA = t0+gkinspback.tecc_init - n*dtTraj;
	t=tLISA+cos(2.*M_PI*tLISA/YEAR)*doppcs+sin(2.*M_PI*tLISA/YEAR)*doppsin-tSSB0;
	dtSSB=dtTraj*(1.-2.*M_PI*(sin(2.*M_PI*tLISA/YEAR)*doppcs-cos(2.*M_PI*tLISA/YEAR)*doppsin)/YEAR);
	r = gkinspback.rfunc( t, gkinspback.x[1] ); 
	ct = gkinspback.costhetafunc( t, gkinspback.x[2] );
	st = sqrt( 1.0 - ct*ct );
	x[-n] = r*st*cos( gkinspback.x[3] );
	y[-n] = r*st*sin( gkinspback.x[3] );
	z[-n] = r*ct;
	gkinspback.GKInsp_QuadrupoleMoments(t, gkinspback.x, Iab);
	for (i=0;i<5;i++)
	  Quad[i][-n]=Iab[i+1];
	gkinspback.TakeAStep( gkinspback.x, t, t - dtSSB );
      }
    }
    // gkinspback.~GKInsp();
    //gktrajback.~GKTraj();
    free( EccTrajBack );
    free( CircTrajBack );
  }

  n=0;
  Necc=0;
  Ncirc=0;
  TrajData *EccTraj=NULL;
  TrajData *CircTraj=NULL;
  
  GKTraj gktraj( cos( pI[ORBIT_I] ), pI[SYSTEM_A] );
  gktraj.p = pI[ORBIT_P];
  gktraj.ecc = pI[ORBIT_E];
  /* NB We add an extra 50000 points to this to account for the fact that the timestep may be reduced towards
        the end of the inspiral, and therefore the phase space trajectory would be too short. */
  Ntraj_max = (int)( (Nmx*dtTraj+500.)*exp(pI[LN_ORBIT_ETA])
                     /(SOLARMASSINSEC*exp(pI[LN_SYSTEM_M]))/DELTAT_TRAJ ) + 2 + 50000;
  
  EccTraj = (TrajData *)
    malloc( (size_t)( ( Ntraj_max + 1 )*sizeof(TrajData) ) );
  CircTraj = (TrajData *)
    malloc( (size_t)( ( Ntraj_max + 1 )*sizeof(TrajData) ) );
  if ( !EccTraj || !CircTraj ) {
    fprintf( stderr, "Error: InspiralTrajectory: Out of memory!\n" );
    if ( EccTraj ) free ( EccTraj );
    if ( CircTraj ) free ( CircTraj );
    length=0;
  }

  ecc_traj = gktraj.Eccentric( DELTAT_TRAJ, EccTraj, Ntraj_max, Necc );

  if ( ecc_traj == 1 )
    circ_traj = gktraj.Circular( DELTAT_TRAJ, CircTraj, Ntraj_max, Ncirc );

  /* Check for initial instability or other fatal errors. */
  if ( Necc == 0 && Ncirc == 0 ) {
    free( EccTraj );
    free( CircTraj );
    length=0;
  }

  /* Assuming that the inspiral has not reached plunge, check that the phase space trajectory is long enough
     (since the timestep might have been reduced making the trajectory end too early). */
  if (ecc_traj == 2) {
    if (EccTraj[Necc-1].t < (Nmx*dtTraj+500.)*exp(pI[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(pI[LN_SYSTEM_M]))) {
      cerr << "Error: phase space trajectory is too short!!!" << endl;
      exit(0);
    }
  } else if (circ_traj == 2) {
    if (CircTraj[Ncirc-1].t < (Nmx*dtTraj+500.)*exp(pI[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(pI[LN_SYSTEM_M]))) {
      cerr << "Error: phase space trajectory is too short!!!" << endl;
      exit(0);
    }
  }

  GKInsp gkinsp( EccTraj,Necc,CircTraj,Ncirc,SOLARMASSINSEC*exp(pI[LN_SYSTEM_M]),
                 exp(-1.*pI[LN_ORBIT_ETA]),pI[SYSTEM_A],dtTraj,pI[ORBIT_P],1.0 );

  gkinsp.x[1] = pI[ORBIT_PSI0];

  if ((pI[ORBIT_I]*(M_PI-pI[ORBIT_I])) >= 0.)
    gkinsp.x[2] = pI[ORBIT_CHI0];
  else
    gkinsp.x[2] = M_PI+pI[ORBIT_CHI0];

  gkinsp.x[3] = 0.;
  gkinsp.x[4] = 0.;


  if ( gkinsp.Necc > 0 ) {
    gkinsp.DoingEcc = 1;
    gkinsp.ihi = 2;
    nMax = ( (int)( ( gkinsp.tecc_fin - gkinsp.tecc_init )
                    / gkinsp.deltat_est ) );
    
    if ( nMax > Nmx )
      nMax = Nmx;
   
    for ( n = 0; n < Nmx; n++ ) {
      t = gkinsp.tecc_init + n*dtTraj;
      tLISA = t0+gkinsp.tecc_init + n*dtTraj;
      //cout << t << " " << tLISA << endl;
      t=tLISA+cos(2.*M_PI*tLISA/YEAR)*doppcs+sin(2.*M_PI*tLISA/YEAR)*doppsin-tSSB0;
      dtSSB=dtTraj*(1.-2.*M_PI*(sin(2.*M_PI*tLISA/YEAR)*doppcs-cos(2.*M_PI*tLISA/YEAR)*doppsin)/YEAR);
      //cout << t << " " << dtSSB << " " << tSSB0 << endl;
      r = gkinsp.rfunc( t, gkinsp.x[1] );
      ct = gkinsp.costhetafunc( t, gkinsp.x[2] );
      st = sqrt( 1.0 - ct*ct );
      x[n] = r*st*cos( gkinsp.x[3] ); 
      y[n] = r*st*sin( gkinsp.x[3] );
      z[n] = r*ct;
      
      gkinsp.GKInsp_QuadrupoleMoments(t, gkinsp.x, Iab);
      for (i=0;i<5;i++)
	Quad[i][n]=Iab[i+1];
      //cout << t << " " << dtSSB << endl;
      gkinsp.TakeAStep( gkinsp.x, t, t + dtSSB );
      if (gkinsp.errencount) {
	cerr << "Inspiral ending eccentric at t = " << t << " as LSO is approached." << endl;
	Nmx=n;
      }
    }
    if (!gkinsp.Ncirc) {
      ival=0;
      notfound=true;
      while (notfound) {
	ival++;
	if (ival == Necc-1) {
	  notfound=false;
	  iup=Necc-1;
	  frac=1.;
	} else {
	  if (EccTraj[ival].t > t*exp(pI[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(pI[LN_SYSTEM_M]))) {
	    iup=ival;
	    frac=(t*exp(pI[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(pI[LN_SYSTEM_M]))-EccTraj[ival-1].t)
	      /(EccTraj[ival].t - EccTraj[ival-1].t);
	    notfound=false;
	  }
	}
      }
    } else {
      iup=Necc-1;
      frac=1.;
    }
    p_end=(1.-frac)*EccTraj[iup-1].p+frac*EccTraj[iup].p;
    e_end=(1.-frac)*EccTraj[iup-1].ecc+frac*EccTraj[iup].ecc;
    i_end=acos((1.-frac)*EccTraj[iup-1].cosiota+frac*EccTraj[iup].cosiota);
    zmin_end=(1.-frac)*EccTraj[iup-1].zedminus+frac*EccTraj[iup].zedminus;
  }

  /* Now do the circular inspiral. */
  if ( n < Nmx && gkinsp.Ncirc > 0 ) {
    gkinsp.DoingEcc = 0;
    gkinsp.ihi = 2;

    /* Make sure we don't underrun the circular trajectory table.  If
       we do, stop integration right here. */
    if ( ( t = gkinsp.tecc_init + n*dtTraj ) < gkinsp.tcirc_init ) {
      free( EccTraj );
      free( CircTraj );
      length=n;
    }

    /* Make sure we don't overrun the eccentric trajectory table. */
    nMax = ( (int)( ( gkinsp.tcirc_fin - t ) / gkinsp.deltat_est ) );
    if ( n + nMax > Nmx )
      nMax = Nmx;
    else
      nMax += n;

    /* Continue with the integration. */
    for ( ; n < Nmx; n++ ) {
      tLISA = t0+gkinsp.tecc_init + n*dtTraj;
      t=tLISA+cos(2.*M_PI*tLISA/YEAR)*doppcs+sin(2.*M_PI*tLISA/YEAR)*doppsin-tSSB0;
      dtSSB=dtTraj*(1.-2.*M_PI*(sin(2.*M_PI*tLISA/YEAR)*doppcs-cos(2.*M_PI*tLISA/YEAR)*doppsin)/YEAR);
      r = gkinsp.rfunc( t, gkinsp.x[1] ); 
      ct = gkinsp.costhetafunc( t, gkinsp.x[2] );
      st = sqrt( 1.0 - ct*ct );
      x[n] = r*st*cos( gkinsp.x[3] );
      y[n] = r*st*sin( gkinsp.x[3] );
      z[n] = r*ct;
      gkinsp.GKInsp_QuadrupoleMoments(t, gkinsp.x, Iab);
      for (i=0;i<5;i++)
	Quad[i][n]=Iab[i+1];
      gkinsp.TakeAStep( gkinsp.x, t, t + dtSSB );
      if (gkinsp.errencount) {
	cerr << "Inspiral ending circular at t = " << t << " as LSO is approached." << endl;
	Nmx=n;
      }
    }
    ival=0;
    notfound=true;
    while (notfound) {
      ival++;
      if (ival == Ncirc-1) {
	notfound=false;
	iup=Ncirc-1;
	frac=1.;
      } else {
	if (CircTraj[ival].t > t*exp(pI[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(pI[LN_SYSTEM_M]))) {
	  iup=ival;
	  frac=(t*exp(pI[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(pI[LN_SYSTEM_M]))-CircTraj[ival-1].t)
	    /(CircTraj[ival].t - CircTraj[ival-1].t);
	  notfound=false;
	}
      }
    }
    p_end=(1.-frac)*CircTraj[iup-1].p+frac*CircTraj[iup].p;
    e_end=(1.-frac)*CircTraj[iup-1].ecc+frac*CircTraj[iup].ecc;
    i_end=acos((1.-frac)*CircTraj[iup-1].cosiota+frac*CircTraj[iup].cosiota);
    zmin_end=(1.-frac)*CircTraj[iup-1].zedminus+frac*CircTraj[iup].zedminus;
  }  
  //cout << gkinsp.x[1] << " " << gkinsp.x[2] << " " << gkinsp.x[3] << endl;
  //gkinsp.~GKInsp();
  //gktraj.~GKTraj();
  free( EccTraj );
  free( CircTraj );
  length=n;
  
}

DopplerShiftedWaveform::~DopplerShiftedWaveform()
{
  int i;
  free(xx);
  free(yy);
  free(zz);
  for (i=0;i<5;i++)
    free(Qdmom[i]);
  free(Qdmom);
}


void DopplerShiftedWaveform::ProjectQuadMomTensor(int Nmin, int Nmax, Real **QuadTensor)
{
  int i,j,k;

  Real x0,x1,x2,y0,y1,y2;
  Real cphi0,sphi0,sini;
  Real coeffs[2][5];

  cphi0=cos(pE[ORBIT_PHI0]);
  sphi0=sin(pE[ORBIT_PHI0]);

  if (fabs(csinc) < 1.) 
    sini = sqrt(1.-csinc*csinc); /* precomputation */
  else
    sini=0.;

  /* Compute projections into transverse plane and onto actual x, y axes. */
  /* The plus component is given by 0.5 * d^2 (X^2 - Y^2 )/dt^2 and the cross component by d^2(XY)/dt^2, 
     where X and Y are components relative to an X, Y basis on the sky. This basis is taken to be right handed with 
     respect to the direction of propagation of the gravitational wave (i.e. it is a left handed set with respect 
     to the direction to the source, n ).
     
     We take the X axis to be perpendicular to the spin of the black hole and to our line of sight, with positive
     x in the direction of n^S (this is also how we define phi0). The x[i] coordinates of the inspiral are computed 
     taking phi(0) = 0, so we can construct X[i] = x[i]*cos(phi0)-y[i]*sin(phi0) and 
     Y[i] = -(y[i]*cos(phi0)+x[i]*sin(phi0))*cos(i) + z[i]*sin(i), where i is the inclination angle between the line 
     of sight and the spin axis of the source. 
     
     We may use these definitions and the reduced quadrupole moment tensor in the black hole frame to generate the
     plus and cross components on the sky. We can precompute the necessary coefficients for the rotation. */


  coeffs[0][0]=0.5*(1.+cphi0*cphi0-csinc*csinc*(1.+sphi0*sphi0));
  coeffs[0][1]=0.5*(1.+sphi0*sphi0-csinc*csinc*(1.+cphi0*cphi0));
  coeffs[0][2]=-cphi0*sphi0*(1.+csinc*csinc);
  coeffs[0][3]=sini*csinc*sphi0;
  coeffs[0][4]=sini*csinc*cphi0;

  coeffs[1][0]=-cphi0*sphi0*csinc;
  coeffs[1][1]=cphi0*sphi0*csinc;
  coeffs[1][2]=(sphi0*sphi0-cphi0*cphi0)*csinc;
  coeffs[1][3]=cphi0*sini;
  coeffs[1][4]=-sphi0*sini;

  for ( i = Nmin; i < Nmax; i++ ){
    for (j=0;j<2;j++) {
      QuadTensor[i-Nmin][j]=0.;
      for (k=0;k<5;k++)
	QuadTensor[i-Nmin][j]+=coeffs[j][k]*(Quad[k][i]);
    }
  }
  
  return;
}

/* Nmin, Nmax are the start and end points of the trajectory, while t0 is the time at the 0th element of the array, 
   needed only for the modulation by the detector motion. */

int DopplerShiftedWaveform::hLISA(int Nmin, int Nmax, Real *hI, Real *hII, Real t0)
{
  int i,j,k,l,ct,num,len;
  Real dist,amp,deltat;

  len=length;

  /* NB When the Inspiral code was changed to use the geodesic equations to compute the second time derivative, it 
        became no longer necessary to divide by (dt/s)^2 when computing the corresponding waveform. Thus, we need
	to divide only by the distance here. */

  if (pE[REDSHIFT]) {
    dist=LCDM->LumD(pE[REDSHIFT]);
  } 
  else
    dist=pE[SYSTEM_DP];

  amp=norm/dist;

  num=Nmax-Nmin;
  Real t,csth,tphi,tpsi,c2psi,s2psi,tdphi,c2phi,s2phi,fonep,fonex,ftwop,ftwox,orbitalt,cost,sint;
  Real **QuadTensor;
  QuadTensor=(Real **)malloc((num)*sizeof(Real *));
  for (i=0;i<(num);i++)
    QuadTensor[i]=(Real *)malloc(2*sizeof(Real));
  
  ProjectQuadMomTensor(Nmin,Nmax,QuadTensor);

  Real hplus,hcross,frac;

  for (ct=0;ct<num;ct++) {
    t=t0+((Real)(Nmin+ct))*dt;
    orbitalt=1.*t/YEAR;
    cost=cos(2.*M_PI*orbitalt);
    sint=sin(2.*M_PI*orbitalt);
    csth=0.5*csthsky-0.5*rtthr*sinthsky*(cost*cphis+sint*sphis);
    tdphi=-0.5*(cost*cphis+sint*sphis+rtthr*csthsky/sinthsky)/(sint*cphis-cost*sphis);
    tphi=(sint-cost*tdphi)/(cost+sint*tdphi);
    c2phi=2./(1.+tphi*tphi)-1.;
    s2phi=2.*tphi/(1.+tphi*tphi);
    tpsi= ( 0.5*csthk-0.5*rtthr*sinthk*(cost*cphik+sint*sphik)
	      -csth*(csthk*csthsky+sinthk*sinthsky*cos(phik-phis) ) )/
		(0.5*sinthk*sinthsky*sin(phik-phis) -0.5*rtthr*
		 cost*(csthk*sinthsky*sphis-csthsky*sinthk*sphik)
		 -0.5*rtthr*sint*(csthsky*sinthk*cphik-csthk*sinthsky*cphis));
    c2psi=2./(1.+tpsi*tpsi)-1.;
    s2psi=2.*tpsi/(1.+tpsi*tpsi);

    if (hI) {
      fonep=0.5*(1.+csth*csth)*c2phi*c2psi-csth*s2phi*s2psi;
      fonex=0.5*(1.+csth*csth)*c2phi*s2psi +csth*s2phi*c2psi;
      hI[ct+Nmin]=amp*(fonep*QuadTensor[ct][0]+fonex*QuadTensor[ct][1]);
    }
    if (hII) {
      ftwop=0.5*(1.+csth*csth)*s2phi*c2psi+csth*c2phi*s2psi;
      ftwox=0.5*(1.+csth*csth)*s2phi*s2psi-csth*c2phi*c2psi;
      hII[ct+Nmin]=amp*(ftwop*QuadTensor[ct][0]+ftwox*QuadTensor[ct][1]);
    }
  }

  for (i=0;i<num;i++)
    free(QuadTensor[i]);
  free(QuadTensor);

  return num;
}

int DopplerShiftedWaveform::hpluscross(int Nmin, int Nmax, Real *hp, Real *hx)
{
  int i,j,k,l,ct,num;
  Real **q;
  q=(Real **)malloc(2*sizeof(Real *));
  for (i=0;i<2;i++)
    q[i]=(Real *)malloc(2*sizeof(Real));

  Real dist,amp;

  /* See comment in hLISA for a note about the factors here. */

  if (pE[REDSHIFT]) {
    dist=LCDM->LumD(pE[REDSHIFT]);
  } 
  else 
    dist=pE[SYSTEM_DP];
  
  /* NB We included a factor of sqrt(3)/2 in the definition of norm which really forms part of the detector
        response, so need to get rid of this now if we only want the h+ and hx waveforms. */

  amp=2.*norm/(dist*sqrt(3.));

  num=Nmax-Nmin;

  Real t,csth,phi,psi,fonep,ftwop,fonex,ftwox,orbitalt;
  Real **QuadTensor;
  QuadTensor=(Real **)malloc((num)*sizeof(Real *));
  for (i=0;i<(num);i++)
    QuadTensor[i]=(Real *)malloc(2*sizeof(Real));

  ProjectQuadMomTensor(Nmin,Nmax,QuadTensor);

  for (ct=0;ct<num;ct++) {
    hp[ct+Nmin]=amp*QuadTensor[ct][0];
    hx[ct+Nmin]=amp*QuadTensor[ct][1];
  }

  for (i=0;i<num;i++)
    free(QuadTensor[i]);
  free(QuadTensor);

  return num;
}
