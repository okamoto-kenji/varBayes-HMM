# Makefile

gaussDiffProg = vbHmmGaussDiff
gaussDiffObjs = vbHmm_Class.o vbHmmCore.o vbHmmModel_Common.o vbHmmData_GaussDiffusion.o vbHmmModel_GaussDiffusion.o vbHmmGaussDiffusion_Main.o mathUtils.o rand.o

templateProg = vbHmm_template
templateObjs = vbHmm_Class.o vbHmmCore.o vbHmmModel_Common.o vbHmmData_template.o vbHmmModel_template.o vbHmmTemplate_Main.o mathUtils.o rand.o

INCLUDES = /usr/local/include
LIBS = /usr/local/lib

# for gcc
#CC = g++
#CCFLAGS = -O3 -fopenmp
#CXXFLAGS = -O3 -fopenmp
#LDFLAGS = -lm -lgsl -lgslcblas -fopenmp

# for clang
CC = clang++
#CCFLAGS = -O3
CXXFLAGS = -O3 -I$(INCLUDES)
LDFLAGS = -lgsl -lgslcblas -L$(LIBS)

# for intel compiler
#CC = icc
#CCFLAGS = -O3 -openmp
#CXXFLAGS = -O3 -openmp
#LDFLAGS = -lm -lgsl -lgslcblas -static-intel

.SUFFIXES: .c .cpp .o

$(program) : $(objs)
	$(CC) $^ $(LDFLAGS) -o $(program)

gaussDiff : $(gaussDiffProg)
$(gaussDiffProg) : $(gaussDiffObjs)
	$(CC) $^ $(LDFLAGS) -o $(gaussDiffProg)

temp : $(templateProg)
$(templateProg) : $(templateObjs)
	$(CC) $^ $(LDFLAGS) -o $(templateProg)

.c.cpp.o:
	$(CC) $(CCFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) $(gaussDiffProg) $(gaussDiffObjs) $(templateProg) $(templateObjs)

vbHmm_Class.o : vbHmm_Class.h
vbHmmModel_Common.o : vbHmmModel_Common.h vbHmm_Class.h mathUtils.h rand.h
vbHmmCore.o : vbHmmCore.h vbHmmModel_Common.h vbHmm_Class.h mathUtils.h rand.h

vbHmmData_GaussDiffusion.o : vbHmm_GaussDiffusion.h vbHmmModel_Common.h vbHmmCore.h vbHmm_Class.h mathUtils.h rand.h
vbHmmModel_GaussDiffusion.o : vbHmm_GaussDiffusion.h vbHmmModel_Common.h vbHmmCore.h vbHmm_Class.h mathUtils.h rand.h
vbHmmGaussDiffusion_Main.o : vbHmmGaussDiffusion_Main.h vbHmm_GaussDiffusion.h vbHmmModel_Common.h vbHmmCore.h vbHmm_Class.h mathUtils.h rand.h

vbHmmData_template.o : vbHmm_template.h vbHmmModel_Common.h vbHmmCore.h vbHmm_Class.h mathUtils.h rand.h
vbHmmModel_template.o : vbHmm_template.h vbHmmModel_Common.h vbHmmCore.h vbHmm_Class.h mathUtils.h rand.h
vbHmmTemplate_Main.o : vbHmmTemplate_Main.h vbHmm_template.h vbHmmModel_Common.h vbHmmCore.h vbHmm_Class.h mathUtils.h rand.h


mathUtils.o : mathUtils.h
rand.o : rand.h

#
