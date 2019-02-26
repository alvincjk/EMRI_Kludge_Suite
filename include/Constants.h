#ifndef _CONSTANTS_H
#define _CONSTANTS_H

const double G = 6.67259e-11;
const double C = 299792458.0;

const double Msun = 1.98892e30;
const double MSUN_CM = G*Msun/C/C*100.;
const double MSUN_SEC = G*Msun/(C*C*C);
const double SOLARMASSINSEC = G*Msun/(C*C*C); // duplicate definition (for historical reasons)

const double pc = 3.0856775807e16;
const double GPCINSEC = pc*1.e9;
const double Gpc = GPCINSEC/C; // names are messed up (for historical reasons)

const double year = 31558149.8; // sidereal year
const double YEAR = 31558149.8;

const double AU = 1.4959787066e11;
const double AUsec = AU/C;

#endif
