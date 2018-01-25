#ifndef _CONSTANTS_H
#define _CONSTANTS_H

const double G = 6.6726e-11;
const double C = 299792458.0;

const double Msun = 1.9889e30;
const double MSUN_CM = 147661.2814476609;
const double MSUN_SEC = MSUN_CM/(100.*C);
const double SOLARMASSINSEC = G*Msun/(C*C*C);

const double pc = 3.08567818585e16;
const double Gpc = 1.02938e17;
const double GPCINSEC = 3.08567818585e25;

const double year = 31536000.0;
const double YEAR = 31536000.0;

const double AU = 1.4959787066e11;
const double AUsec = 499.004783702731;

#endif
