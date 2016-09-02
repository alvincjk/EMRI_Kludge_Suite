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
#include "Inspiral.h"

#define DELTAT (0.05)     /* default sampling interval as a fraction
			     of orbital periapsis */
#define DELTAT_TRAJ (0.001) /* sampling interval for phase-space
                             trajectory generation */

#define EPSMINUS (0.9999999) /* 1 minus a small resolveable difference */
#define EPSPLUS  (1.0000001) /* 1 plus a small resolveable difference */
#define REALMIN (2.0e-38)   /* a small positive normal */
#define SOLARMASSINSEC (G*Msun/(C*C*C)) /* 1 solar mass in units of seconds. */

bool useBLCoords=true;
//bool storephase=true;
bool storephase=false;

/* Note that we assume that SYSTEM_M is given in solar mass units and the 
   trajectory sampling interval dt is given in seconds. The factor 
   SOLARMASSINSEC converts from solar masses to seconds where necessary. */

Inspiral::Inspiral(const Real *params, double deltat, const int Nfrst, const int Nlast, const bool qd, const bool datout, char *filename) : 
  p0(params),dt(deltat),Nzero(Nfrst),Nmx(Nlast),genquad(qd),outputdata(datout),filenm(filename)
{ 
  IEKG *iekg;
  Real omega[3],omold[3],psiold,chiold,phiold;
  int trajlen,i,npsi;
  if (Nzero >= 0)
    trajlen=Nmx;
  else
    trajlen=Nmx-Nzero;
  xx=(Real *)malloc(trajlen*sizeof(Real));
  yy=(Real *)malloc(trajlen*sizeof(Real));
  zz=(Real *)malloc(trajlen*sizeof(Real));
  
  if (!xx) {
    cerr << "Memory allocation error (xx) in Inspiral constructor." << endl;
    exit(-1);
  }
  if (!yy) {
    cerr << "Memory allocation error (yy) in Inspiral constructor." << endl;
    exit(-1);
  }
  if (!zz) {
    cerr << "Memory allocation error (zz) in Inspiral constructor." << endl;
    exit(-1);
  }
  if (genquad) {
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
  }

  if (Nzero >= 0) {
    x=&xx[0];
    y=&yy[0];
    z=&zz[0];
    if (genquad) {
      for (i=0;i<5;i++)
	Quad[i]=&(Qdmom[i][0]);
    }
  }
  else {
    x=&xx[-Nzero];
    y=&yy[-Nzero];
    z=&zz[-Nzero];
    if (genquad) {
      for (i=0;i<5;i++)
	Quad[i]=&(Qdmom[i][-Nzero]);
    }
  }
  //cout << "In Inspiral constructor, about to start loop " << p0[ORBIT_E] << endl;
  int n = 0, nMax; /* loop index and limit */
  int Ntraj_max;   /* maximum points in phase-space trajectories */
  int Necc = 0;    /* actual number of eccentric points in above */
  int Ncirc = 0;   /* actual number of circular points in above */
  int ecc_traj, circ_traj; /* return codes from phase-space generator */
  TrajData *EccTrajBack = NULL;  /* eccentric phase-space trajectory data */
  TrajData *CircTrajBack = NULL; /* circular phase-space trajectory data */
  Real t, r, ct, st;         /* integration temporary variables */
  Real Iab[6]; /* Temporary quadrupole moment tensor. */
  Real massfactor;
  Real psival,tKS,dtBL,speeds[3],cpsi,spsi,KSoff;

  int iup,ival;
  bool notfound;
  Real frac;

  if (Nzero < 0) {  
    //cout << Nzero << " " << cos( p0[ORBIT_I] ) << " " << p0[SYSTEM_A] << endl;
    GKTraj gktrajback( cos( p0[ORBIT_I] ), p0[SYSTEM_A] );
    gktrajback.p = p0[ORBIT_P];
    gktrajback.ecc = p0[ORBIT_E];
    Ntraj_max = (int)((1 -Nzero)*dt*exp(p0[LN_ORBIT_ETA])
		       /(SOLARMASSINSEC*exp(p0[LN_SYSTEM_M]))/DELTAT_TRAJ ) + 2;
    //cout << 1-Nzero << " " << Ntraj_max << endl;
    //cout << "Here! Ntraj_max = " << Ntraj_max << endl;
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
    //cout << Ntraj_max << endl;
    ecc_traj = gktrajback.Eccentric( -DELTAT_TRAJ, EccTrajBack, Ntraj_max, Necc );
    //cout << Necc << " " << ecc_traj << endl;
    if ( ecc_traj == 1 )
      circ_traj = gktrajback.Circular( DELTAT_TRAJ, CircTrajBack, Ntraj_max, Ncirc );

    //cout << "Creating gkinspback" << endl;
    GKInsp gkinspback( EccTrajBack,Necc,CircTrajBack,Ncirc,SOLARMASSINSEC*exp(p0[LN_SYSTEM_M]),
		       exp(-1.*p0[LN_ORBIT_ETA]),p0[SYSTEM_A],-dt,p0[ORBIT_P],1.0 );
    //cout << "Created" << endl;
    gkinspback.x[1] = p0[ORBIT_PSI0];

    if ((p0[ORBIT_I]*(M_PI-p0[ORBIT_I])) >= 0.)
      gkinspback.x[2] = p0[ORBIT_CHI0];
    else
      gkinspback.x[2] = M_PI+p0[ORBIT_CHI0];

    gkinspback.x[3] = 0.;
    gkinspback.x[4] = 0.;

    if ( gkinspback.Necc > 0 ) {
      gkinspback.DoingEcc = 1;
      gkinspback.ihi = 2;
      nMax = ( (int)( ( gkinspback.tecc_fin - gkinspback.tecc_init )
		      / gkinspback.deltat_est ) );

      if ( nMax > 1-Nzero )
	nMax = 1-Nzero;

      //KSoff=gkinspback.rfunc(gkinspback.tecc_init , gkinspback.x[1] )-gkinspback.rstarfunc(gkinspback.tecc_init,gkinspback.x[1]);
      t = gkinspback.tecc_init;
      for ( n = 0; n < 1-Nzero; n++ ) {
	//cout << n << endl;
	if (useBLCoords) {
		t = gkinspback.tecc_init - n*dt;
		r = gkinspback.rfunc( t, gkinspback.x[1] );
		ct = gkinspback.costhetafunc( t, gkinspback.x[2] );
		st = sqrt( 1.0 - ct*ct );
		x[-n] = r*st*cos( gkinspback.x[3] ); 
		y[-n] = r*st*sin( gkinspback.x[3] );
		z[-n] = r*ct;
	} else {
		tKS=gkinspback.tecc_init - n*dt;
		//t=tKS+r-rs-KSoff;
		r = gkinspback.rfunc( t, gkinspback.x[1] );
                ct = gkinspback.costhetafunc( t, gkinspback.x[2] );
                st = sqrt( 1.0 - ct*ct );
		psival=gkinspback.psifunc(t,gkinspback.x[1],gkinspback.x[3]);
		cpsi=cos(psival);
		spsi=sin(psival);
		gkinspback.getorbitalspeeds(t,speeds,gkinspback.x[1],gkinspback.x[2]);
		dtBL=dt/(1.+2.*r*speeds[0]/Kerr::Delta(r,p0[SYSTEM_A]));
		x[-n]=(r*cpsi-p0[SYSTEM_A]*spsi)*st;
		y[-n]=(r*spsi+p0[SYSTEM_A]*cpsi)*st;
		z[-n]=r*ct;
		t-=dtBL;
	}
	if (genquad) {
	  gkinspback.GKInsp_QuadrupoleMoments(t, gkinspback.x, Iab);
	  for (i=0;i<5;i++)
	    Quad[i][-n]=Iab[i+1];
	}
	//cout << "About to step, t = " << t << ", t-dt = " << t-dt << endl;
	gkinspback.TakeAStep( gkinspback.x, t, t - dt );
      }
      ival=0;
      notfound=true;
      while (notfound) {
	ival++;
	if (ival == Necc-1) {
	  notfound=false;
	  iup=Necc-1;
	  frac=1.;
	} else {
	  if (EccTrajBack[ival].t < t*exp(p0[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(p0[LN_SYSTEM_M]))) {
	    iup=ival;
	    frac=(t*exp(p0[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(p0[LN_SYSTEM_M]))-EccTrajBack[ival-1].t)
	      /(EccTrajBack[ival].t - EccTrajBack[ival-1].t);
	    notfound=false;
	  }
	}
      }
      //cout << Necc << " " << iup << " " << frac << endl;
      p_st=(1.-frac)*EccTrajBack[iup-1].p+frac*EccTrajBack[iup].p;
      e_st=(1.-frac)*EccTrajBack[iup-1].ecc+frac*EccTrajBack[iup].ecc;
      i_st=acos((1.-frac)*EccTrajBack[iup-1].cosiota+frac*EccTrajBack[iup].cosiota);
      zmin_st=(1.-frac)*EccTrajBack[iup-1].zedminus+frac*EccTrajBack[iup].zedminus;
      massfactor=(exp(p0[LN_ORBIT_ETA]-p0[LN_SYSTEM_M])/SOLARMASSINSEC);
      pdot_st=massfactor*(EccTrajBack[iup].p-EccTrajBack[iup-1].p)/(EccTrajBack[iup].t-EccTrajBack[iup-1].t);
      edot_st=massfactor*(EccTrajBack[iup].ecc-EccTrajBack[iup-1].ecc)/(EccTrajBack[iup].t-EccTrajBack[iup-1].t);
      zmindot_st=massfactor*(EccTrajBack[iup].zedminus-EccTrajBack[iup-1].zedminus)/(EccTrajBack[iup].t-EccTrajBack[iup-1].t);
    }
    else {
      gkinspback.DoingEcc = 0;
      gkinspback.ihi = 2;

      
      /* Continue with the integration. */
      for (n=0 ; n < 1-Nzero; n++ ) {
	t = gkinspback.tecc_init - n*dt;
	r = gkinspback.rfunc( t, gkinspback.x[1] ); 
	ct = gkinspback.costhetafunc( t, gkinspback.x[2] );
	st = sqrt( 1.0 - ct*ct );
	x[-n] = r*st*cos( gkinspback.x[3] );
	y[-n] = r*st*sin( gkinspback.x[3] );
	z[-n] = r*ct;
	if (genquad) {
	  gkinspback.GKInsp_QuadrupoleMoments(t, gkinspback.x, Iab);
	  for (i=0;i<5;i++)
	    Quad[i][-n]=Iab[i+1];
	}
	gkinspback.TakeAStep( gkinspback.x, t, t - dt );
      }
      ival=0;
      notfound=true;
      while (notfound) {
	ival++;
	if (ival == Necc-1) {
	  notfound=false;
	  iup=Necc-1;
	  frac=1.;
	} else {
	  if (CircTrajBack[ival].t < t*exp(p0[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(p0[LN_SYSTEM_M]))) {
	    iup=ival;
	    frac=(t*exp(p0[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(p0[LN_SYSTEM_M]))-CircTrajBack[ival-1].t)
	      /(CircTrajBack[ival].t - CircTrajBack[ival-1].t);
	    notfound=false;
	  }
	}
      }
      p_st=(1.-frac)*CircTrajBack[iup-1].p+frac*CircTrajBack[iup].p;
      e_st=0.;
      i_st=acos((1.-frac)*CircTrajBack[iup-1].cosiota+frac*CircTrajBack[iup].cosiota);
      zmin_st=(1.-frac)*CircTrajBack[iup-1].zedminus+frac*CircTrajBack[iup].zedminus;
      massfactor=(exp(p0[LN_ORBIT_ETA]-p0[LN_SYSTEM_M])/SOLARMASSINSEC);
      pdot_st=massfactor*(CircTrajBack[iup].p-CircTrajBack[iup-1].p)/(CircTrajBack[iup].t-CircTrajBack[iup-1].t);
      edot_st=0.;
      idot_st=massfactor*(acos(CircTrajBack[iup].cosiota)-acos(CircTrajBack[iup-1].cosiota))/(CircTrajBack[iup].t-CircTrajBack[iup-1].t);
      zmindot_st=massfactor*(CircTrajBack[iup].zedminus-CircTrajBack[iup-1].zedminus)/(CircTrajBack[iup].t-CircTrajBack[iup-1].t);
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
  
  GKTraj gktraj( cos( p0[ORBIT_I] ), p0[SYSTEM_A] );
  gktraj.p = p0[ORBIT_P];
  gktraj.ecc = p0[ORBIT_E];
  /* NB We add an extra 50000 points to this to account for the fact that the timestep may be reduced towards
     the end of the inspiral, and therefore the phase space trajectory would be too short. */
  Ntraj_max = (int)( Nmx*dt*exp(p0[LN_ORBIT_ETA])
                     /(SOLARMASSINSEC*exp(p0[LN_SYSTEM_M]))/DELTAT_TRAJ ) + 2 + 50000;
  //cout << Nmx << " " << Ntraj_max << endl;
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
    if (EccTraj[Necc-1].t < (Nmx*dt)*exp(p0[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(p0[LN_SYSTEM_M]))) {
      cerr << "Error: phase space trajectory is too short!!!" << endl;
      exit(0);
    }
  } else if (circ_traj == 2) {
    if (CircTraj[Ncirc-1].t < (Nmx*dt)*exp(p0[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(p0[LN_SYSTEM_M]))) {
      cerr << "Error: phase space trajectory is too short!!!" << endl;
      exit(0);
    }
  }

  //cout << ecc_traj << " " << EccTraj[Necc-1].t << " " << 
  //  EccTraj[Necc-1].t*(SOLARMASSINSEC*exp(p0[LN_SYSTEM_M]))/exp(p0[LN_ORBIT_ETA]) << endl;
 
  GKInsp gkinsp( EccTraj,Necc,CircTraj,Ncirc,SOLARMASSINSEC*exp(p0[LN_SYSTEM_M]),
                 exp(-1.*p0[LN_ORBIT_ETA]),p0[SYSTEM_A],dt,p0[ORBIT_P],1.0 );

  gkinsp.x[1] = p0[ORBIT_PSI0];

  if ((p0[ORBIT_I]*(M_PI-p0[ORBIT_I])) >= 0.)
    gkinsp.x[2] = p0[ORBIT_CHI0];
  else
    gkinsp.x[2] = M_PI+p0[ORBIT_CHI0];

  gkinsp.x[3] = 0.;
  gkinsp.x[4] = 0.;
  int Nrad=0;
  Real orbelts[3],orbspeeds[3],orbconsts[3];
  FILE *outfile;
  bool writeout=false;
  if (filenm) {
	writeout=true;
  	outfile=fopen(filenm,"w");
	printf("Writing file!!\n");
  }
  if ( gkinsp.Necc > 0 ) {
    gkinsp.DoingEcc = 1;
    gkinsp.ihi = 2;
    nMax = ( (int)( ( gkinsp.tecc_fin - gkinsp.tecc_init )
                    / gkinsp.deltat_est ) );
    
    if ( nMax > Nmx )
      nMax = Nmx;
   
    for ( n = 0; n < Nmx; n++ ) {
      t = gkinsp.tecc_init + n*dt;
      r = gkinsp.rfunc( t, gkinsp.x[1] );
      ct = gkinsp.costhetafunc( t, gkinsp.x[2] );
      st = sqrt( 1.0 - ct*ct );
      if (storephase) {
	x[n]=gkinsp.x[1];
	y[n]=gkinsp.x[2];
	z[n]=gkinsp.x[3];
      } else {
      	x[n] = r*st*cos( gkinsp.x[3] ); 
      	y[n] = r*st*sin( gkinsp.x[3] );
      	z[n] = r*ct;
      }
      gkinsp.getorbitalelements(t,orbelts);
      gkinsp.getorbitalspeeds(t,orbspeeds,gkinsp.x[1],gkinsp.x[2]);
      gkinsp.getorbitalconsts(t,orbconsts);
      iekg = new IEKG(orbelts[0],orbelts[1],orbelts[2],p0[SYSTEM_A]);
      iekg->Frequencies(omega);
      delete(iekg);
      if (writeout) {
	if (!n) {
		psiold=gkinsp.x[1];
		chiold=gkinsp.x[2];
		phiold=gkinsp.x[3];
		omold[0]=omega[0];
		omold[1]=omega[1];
		omold[2]=omega[2];
		psiold=0.;
		chiold=0.;
		phiold=0.;
		npsi=1;
        }
	//fprintf(outfile,"%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",t,r,acos(ct),gkinsp.x[3],gkinsp.x[1]-2.*M_PI*omega[0]*t,gkinsp.x[2]-2.*M_PI*omega[1]*t,gkinsp.x[3]-2.*M_PI*omega[2]*t,(gkinsp.x[1]-2.*M_PI*omega[0]*t-psiold)/dt,(gkinsp.x[2]-2.*M_PI*omega[1]*t-chiold)/dt,(gkinsp.x[3]-2.*M_PI*omega[2]*t-phiold)/dt,2.*M_PI*(omega[0]-omold[0])*t/dt,2.*M_PI*(omega[1]-omold[1])*t/dt,2.*M_PI*(omega[2]-omold[2])*t/dt);
	if (gkinsp.x[1] > ((double)npsi)*M_PI) {
		//fprintf(outfile,"%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",t,r,acos(ct),gkinsp.x[3],gkinsp.x[1]-2.*M_PI*omega[0]*t,gkinsp.x[2]-2.*M_PI*omega[1]*t,gkinsp.x[3]-2.*M_PI*omega[2]*t,psiold,chiold,phiold);
		//fprintf(outfile,"%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",t,r,acos(ct),gkinsp.x[3],gkinsp.x[1]-2.*M_PI*omega[0]*t,gkinsp.x[2]-2.*M_PI*omega[1]*t,gkinsp.x[3]-2.*M_PI*omega[2]*t,(omega[0]-omold[0])/dt,(omega[1]-omold[1])/dt,(omega[2]-omold[2])/dt);
		//fprintf(outfile,"%12.8e %12.8e %12.8e %12.8e\n",t,r,acos(ct),gkinsp.x[3]);
		//fprintf(outfile,"%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",t,orbconsts[0],orbconsts[1],orbconsts[2],gkinsp.x[1]-2.*M_PI*omega[0]*t,gkinsp.x[2]-2.*M_PI*omega[1]*t,gkinsp.x[3]-2.*M_PI*omega[2]*t,psiold,chiold,phiold);
		npsi+=2;
	}
	//psiold=gkinsp.x[1]-2.*M_PI*omega[0]*t;
	//chiold=gkinsp.x[2]-2.*M_PI*omega[1]*t;
	//phiold=gkinsp.x[3]-2.*M_PI*omega[2]*t;
	psiold+=-2.*M_PI*(omega[0]-omold[0])*t;
	chiold+=-2.*M_PI*(omega[1]-omold[1])*t;
	phiold+=-2.*M_PI*(omega[2]-omold[2])*t;
	omold[0]=omega[0];
	omold[1]=omega[1];
	omold[2]=omega[2];
	//fprintf(outfile,"%10.8e %10.8e %10.8e %10.8e\n",t,x[n],y[n],z[n]);
        fprintf(outfile,"%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",t,gkinsp.x[1],gkinsp.x[2],gkinsp.x[3],r,ct,x[n],y[n],z[n]);
        //fprintf(outfile,"%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",t,r,acos(ct),orbspeeds[0],orbspeeds[1],orbspeeds[2],gkinsp.x[3],orbconsts[0],orbconsts[1],orbconsts[2]);
      	//fprintf(outfile,"%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",t,gkinsp.x[1],gkinsp.x[2],orbspeeds[0],orbspeeds[1],orbspeeds[2],gkinsp.x[3],orbconsts[0],orbconsts[1],orbconsts[2]);
      }
      //printf("%12.8e %12.8e %12.8e %12.8e \n",t/(SOLARMASSINSEC*exp(p0[LN_SYSTEM_M])),r,acos(ct),gkinsp.x[3]);
      //cout << x[n] << " " << y[n] << endl;
      //ra=gkinsp.rfunc(t,M_PI);
      //rp=gkinsp.rfunc(t,0.);
      //cout << ((ra-rp)/(ra+rp)-0.3)/t << " " << (2.*ra*rp/(ra+rp)-25.083)/t << endl;
      if (genquad) {
	gkinsp.GKInsp_QuadrupoleMoments(t, gkinsp.x, Iab);
	for (i=0;i<5;i++)
	  Quad[i][n]=Iab[i+1];
      }
      gkinsp.TakeAStep( gkinsp.x, t, t + dt );
      if (gkinsp.x[1] > ((double)Nrad)*2.*M_PI) {
	//cout << t << endl;
	Nrad++;
      }
      if (gkinsp.errencount) {
	cerr << "Inspiral ending eccentric at t = " << t << " as LSO is approached." << endl;
	Nmx=n;
      }
    }
    if (writeout)
    	fclose(outfile);

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
	  if (EccTraj[ival].t > t*exp(p0[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(p0[LN_SYSTEM_M]))) {
	    iup=ival;
	    frac=(t*exp(p0[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(p0[LN_SYSTEM_M]))-EccTraj[ival-1].t)
	      /(EccTraj[ival].t - EccTraj[ival-1].t);
	    notfound=false;
	  }
	}
      }
    } else {
      iup=Necc-1;
      frac=1.;
    }
    //cout << Necc << " " << iup << " " << frac << endl;
    p_end=(1.-frac)*EccTraj[iup-1].p+frac*EccTraj[iup].p;
    e_end=(1.-frac)*EccTraj[iup-1].ecc+frac*EccTraj[iup].ecc;
    i_end=acos((1.-frac)*EccTraj[iup-1].cosiota+frac*EccTraj[iup].cosiota);
    zmin_end=(1.-frac)*EccTraj[iup-1].zedminus+frac*EccTraj[iup].zedminus;
    massfactor=(exp(p0[LN_ORBIT_ETA]-p0[LN_SYSTEM_M])/SOLARMASSINSEC);
    pdot_end=massfactor*(EccTraj[iup].p-EccTraj[iup-1].p)/(EccTraj[iup].t-EccTraj[iup-1].t);
    edot_end=massfactor*(EccTraj[iup].ecc-EccTraj[iup-1].ecc)/(EccTraj[iup].t-EccTraj[iup-1].t);
    idot_end=massfactor*(acos(EccTraj[iup].cosiota)-acos(EccTraj[iup-1].cosiota))/(EccTraj[iup].t-EccTraj[iup-1].t);
    zmindot_end=massfactor*(EccTraj[iup].zedminus-EccTraj[iup-1].zedminus)/(EccTraj[iup].t-EccTraj[iup-1].t);
  }

  /* Now do the circular inspiral. */
  if ( n < Nmx && gkinsp.Ncirc > 0 ) {
    gkinsp.DoingEcc = 0;
    gkinsp.ihi = 2;

    /* Make sure we don't underrun the circular trajectory table.  If
       we do, stop integration right here. */
    if ( ( t = gkinsp.tecc_init + n*dt ) < gkinsp.tcirc_init ) {
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
      t = gkinsp.tecc_init + n*dt;
      r = gkinsp.rfunc( t, gkinsp.x[1] ); 
      ct = gkinsp.costhetafunc( t, gkinsp.x[2] );
      st = sqrt( 1.0 - ct*ct );
      x[n] = r*st*cos( gkinsp.x[3] );
      y[n] = r*st*sin( gkinsp.x[3] );
      z[n] = r*ct;
      if (genquad) {
	gkinsp.GKInsp_QuadrupoleMoments(t, gkinsp.x, Iab);
	for (i=0;i<5;i++)
	  Quad[i][n]=Iab[i+1];
      }
      //std::cout << x[n] << " " << y[n] << " " << z[n] << " " << Quad[0][n] << " " << Quad[1][n] << " "<< Quad[2][n] << " "<< Quad[3][n] << " "<< Quad[4][n] << std::endl;
      gkinsp.TakeAStep( gkinsp.x, t, t + dt );
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
	if (CircTraj[ival].t > t*exp(p0[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(p0[LN_SYSTEM_M]))) {
	  iup=ival;
	  frac=(t*exp(p0[LN_ORBIT_ETA])/(SOLARMASSINSEC*exp(p0[LN_SYSTEM_M]))-CircTraj[ival-1].t)
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
  //gkinsp.~GKInsp();
  //gktraj.~GKTraj();
  free( EccTraj );
  free( CircTraj );
  length=n;
}

Inspiral::~Inspiral()
{
  int i;
  free(xx);
  free(yy);
  free(zz);
  if (genquad) {
    for (i=0;i<5;i++)
      free(Qdmom[i]);
    free(Qdmom);
  }
}

