# Makefile

tsProg = vbHmmTs
tsObjs = vbHmmTs_Main.o vbHmm_Common.o vbHmmTs.o vbHmmTsDataHandler.o rand.o

tsFretProg = vbHmmTsFret
tsFretObjs = vbHmmTsFret_Main.o vbHmm_Common.o vbHmmTsFret.o vbHmmTsFretDataHandler.o rand.o

pcProg = vbHmmPc
pcObjs = vbHmmPc_Main.o vbHmm_Common.o vbHmmPc.o vbHmmPcDataHandler.o rand.o

pcFretProg = vbHmmPcFret
pcFretObjs = vbHmmPcFret_Main.o vbHmm_Common.o vbHmmPcFret.o vbHmmPcFretDataHandler.o rand.o

gaussProg = vbHmmGauss
gaussObjs = vbHmmGauss_Main.o vbHmm_Common.o vbHmmGauss.o vbHmmGaussDataHandler.o rand.o

tempProg = vbHmm_template
tempObjs = vbHmm_template_Main.o vbHmm_Common.o vbHmm_template.o vbHmm_Template_DataHandler.o rand.o

## Compiler dependent variables
## Uncomment either 3 lines below
# for Clang
CC = clang
CCFLAGS = -O3 -I/usr/local/include
LDFLAGS = -lm -lgsl -lgslcblas -L/usr/local/lib
# for general gcc
#CC = gcc
#CCFLAGS = -O3 -fopenmp
#LDFLAGS = -lm -lgsl -lgslcblas -fopenmp
# for intel compiler
#CC = icc
#CCFLAGS = -O3 -openmp
#LDFLAGS = -lm -lgsl -lgslcblas -openmp -static-intel

.SUFFIXES: .c .o

.PHONY: all
all: $(tsProg) $(tsFretProg) $(pcProg) $(pcFretProg) $(tempProg)

#$(program) : $(objs)
#	$(CC) $^ $(LDFLAGS) -o $(program)

ts : $(tsProg)
$(tsProg) : $(tsObjs)
	$(CC) $^ $(LDFLAGS) -o $(tsProg)

fret : $(tsFretProg)
$(tsFretProg) : $(tsFretObjs)
	$(CC) $^ $(LDFLAGS) -o $(tsFretProg)

pc : $(pcProg)
$(pcProg) : $(pcObjs)
	$(CC) $^ $(LDFLAGS) -o $(pcProg)

pcFret : $(pcFretProg)
$(pcFretProg) : $(pcFretObjs)
	$(CC) $^ $(LDFLAGS) -o $(pcFretProg)

gauss : $(gaussProg)
$(gaussProg) : $(gaussObjs)
	$(CC) $^ $(LDFLAGS) -o $(gaussProg)

template : $(tempProg)
$(tempProg) : $(tempObjs)
	$(CC) $^ $(LDFLAGS) -o $(tempProg)

.c.o:
	$(CC) $(CCFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) $(tsProg) $(tsObjs) $(tsFretProg) $(tsFretObjs) $(pcProg) $(pcObjs) $(pcFretProg) $(pcFretObjs) $(gaussProg) $(gaussObjs) $(tempProg) $(tempObjs)

vbHmm_Common.o : vbHmm_Common.h
rand.o : rand.h

vbHmmTs.o            : vbHmm_Common.h vbHmmTs.h rand.h
vbHmmTs_Main.o       : vbHmm_Common.h vbHmmTs.h vbHmmTsDataHandler.h vbHmmTs_Main.h
vbHmmTsDataHandler.o : vbHmm_Common.h vbHmmTs.h vbHmmTsDataHandler.h

vbHmmTsFret.o            : vbHmm_Common.h vbHmmTsFret.h rand.h
vbHmmTsFret_Main.o       : vbHmm_Common.h vbHmmTsFret.h vbHmmTsFretDataHandler.h vbHmmTsFret_Main.h
vbHmmTsFretDataHandler.o : vbHmm_Common.h vbHmmTsFret.h vbHmmTsFretDataHandler.h

vbHmmPc.o            : vbHmm_Common.h vbHmmPc.h rand.h
vbHmmPc_Main.o       : vbHmm_Common.h vbHmmPc.h vbHmmPcDataHandler.h vbHmmPc_Main.h
vbHmmPcDataHandler.o : vbHmm_Common.h vbHmmPc.h vbHmmPcDataHandler.h

vbHmmPcFret.o            : vbHmm_Common.h vbHmmPcFret.h rand.h
vbHmmPcFret_Main.o       : vbHmm_Common.h vbHmmPcFret.h vbHmmPcFretDataHandler.h vbHmmPcFret_Main.h
vbHmmPcFretDataHandler.o : vbHmm_Common.h vbHmmPcFret.h vbHmmPcFretDataHandler.h

vbHmmGauss.o            : vbHmm_Common.h vbHmmGauss.h rand.h
vbHmmGauss_Main.o       : vbHmm_Common.h vbHmmGauss.h vbHmmGaussDataHandler.h vbHmmGauss_Main.h
vbHmmGaussDataHandler.o : vbHmm_Common.h vbHmmGauss.h vbHmmGaussDataHandler.h

vbHmm_template.o             : vbHmm_Common.h vbHmm_template.h rand.h
vbHmm_template_Main.o        : vbHmm_Common.h vbHmm_template.h vbHmm_Template_DataHandler.h vbHmm_template_Main.h
vbHmm_Template_DataHandler.o : vbHmm_Common.h vbHmm_template.h vbHmm_Template_DataHandler.h

#
