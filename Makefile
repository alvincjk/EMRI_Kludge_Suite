# EMRI Kludge Suite: Makefile

#############################################################################

# Paths and variables

TOP = .
BIN = $(TOP)/bin
INC = $(TOP)/include
LIB = $(TOP)/lib

SRC = $(TOP)/src
CIRCSRC = $(SRC)/circ
EXECSRC = $(SRC)/exec
IEKGSRC = $(SRC)/inclecc
NRSRC = $(SRC)/numrec
UTILSRC = $(SRC)/utility
KSSRC = $(SRC)/suite
ALLSRCS = $(CIRCSRC):$(EXECSRC):$(IEKGSRC):$(NRSRC):$(UTILSRC):$(KSSRC)

VPATH = $(BIN):$(INC):$(LIB):$(ALLSRCS)

CC = g++

AR = ar rv

SYSLIBS = -lm -lgsl -lgslcblas -lfftw3

CFLAGS = -O3 -Wall -Wno-unused -Wno-uninitialized -Wno-deprecated 

#############################################################################

# Executable files

all : AAK_Waveform AK_Waveform NK_Waveform

AAK_Waveform : AAK_Waveform.cc -lKS -lIEKG -lGKG -lCirc -lLB -lRRGW -lNR Globals.h GKTrajFast.h KSParMap.h KSTools.h AAK.h AAKpy.h
	$(CC) $(EXECSRC)/AAK_Waveform.cc -o $(BIN)/AAK_Waveform $(CFLAGS) -I$(INC) -L$(LIB) -lKS -lIEKG -lCirc -lGKG -lLB -lRRGW -lNR $(SYSLIBS)

AK_Waveform : AK_Waveform.cc -lKS -lIEKG -lGKG -lCirc -lLB -lRRGW -lNR Globals.h IEKG.h KSParMap.h KSTools.h AK.h
	$(CC) $(EXECSRC)/AK_Waveform.cc -o $(BIN)/AK_Waveform $(CFLAGS) -I$(INC) -L$(LIB) -lKS -lIEKG -lCirc -lGKG -lLB -lRRGW -lNR $(SYSLIBS)

NK_Waveform : NK_Waveform.cc -lKS -lIEKG -lGKG -lCirc -lLB -lRRGW -lNR Globals.h GKInsp.h GKR.h IEKG.h CKG.h Inspiral.h Cosmology.h KSParMap.h KSTools.h DopplerShiftedWaveform.h
	$(CC) $(EXECSRC)/NK_Waveform.cc -o $(BIN)/NK_Waveform $(CFLAGS) -I$(INC) -L$(LIB) -lKS -lIEKG -lCirc -lGKG -lLB -lRRGW -lNR $(SYSLIBS)

#############################################################################

# Object files

CIRCOBJS = CKG.o CKR.o

GKGOBJS = GKG.o

IEKGOBJS = GKR.o GKInsp.o GKTraj.o IEKG.o IEKGMNewt.o Cosmology.o Inspiral.o Waveform.o DopplerShiftedWaveform.o WDBinary.o

LBOBJS = LBFJ.o LBMNewt.o

NROBJS = NRElle.o NREllf.o NREllpi.o NRFactrl.o NRGammln.o NRIndexx.o \
	NRLaguer.o NRLUbksb.o NRLUdcmp.o NRPythag.o NRRunge.o NRUtil.o \
	NRchder.o NRchebev.o NRchebft.o NRpzextr.o NRran2.o NRrc.o NRrd.o \
	NRrf.o NRrj.o NRrzextr.o NRspline.o NRsplint.o NRtqli.o NRtred2.o \
	NRzroot.o NRsncndn.o NRsvdcmp.o NRGaussj.o NRdfour1.o

RRGWOBJS = RRGW.o

KSOBJS = AAK.o AK.o GKTrajFast.o KSParMap.o KSTools.o AAKpy.o

.INTERMEDIATE : $(CIRCOBJS) $(GKGOBJS) $(IEKGOBJS) $(LBOBJS) $(NROBJS) $(RRGWOBJS) $(KSOBJS)

%.o : %.cc
	$(CC) $(CFLAGS) -I$(INC) -c $< -o $@

#############################################################################

# Library files

-lCirc : $(CIRCOBJS) Globals.h CKG.h CKR.h Integral.h NRCKG.h NRUtil.h SWSH.h SNT.h
	$(AR) $(LIB)/libCirc.a $(CIRCOBJS)

-lGKG : $(GKGOBJS) Globals.h GKG.h NRUtil.h
	$(AR) $(LIB)/libGKG.a $(GKGOBJS)

-lLB : $(LBOBJS) Globals.h GKG.h LB.h NRLB.h NRUtil.h
	$(AR) $(LIB)/libLB.a $(LBOBJS)

-lIEKG : $(IEKGOBJS) Globals.h IEKG.h GKR.h GKInsp.h GKTraj.h NRRoot.h NRUtil.h NRIEKG.h \
	NRMM.h Cosmology.h Inspiral.h Waveform.h WDBinary.h
	$(AR) $(LIB)/libIEKG.a $(IEKGOBJS)

-lNR : $(NROBJS) Globals.h NRUtil.h
	$(AR) $(LIB)/libNR.a $(NROBJS)

-lRRGW : $(RRGWOBJS) Globals.h RRGW.h SWSH.h
	$(AR) $(LIB)/libRRGW.a $(RRGWOBJS)

-lKS : $(KSOBJS) Globals.h AAK.h AK.h GKTrajFast.h KSParMap.h KSTools.h AAKpy.h
	$(AR) $(LIB)/libKS.a $(KSOBJS)

#############################################################################

# Make clean

.PHONY: dummy

clean : dummy
	$(RM) $(BIN)/*
	$(RM) $(LIB)/*

