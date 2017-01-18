include $(WORK)/cantera-2.2.1/cantera-2.2.1/include/cantera/Cantera.mak

CXX   = mpicxx -std=c++11
F90   = mpif90
FLAG  = -g -O0 -w #-O2 -w
LIBS  = -L$(WORK)/sundials-2.6.1/sundials-2.6.1/lib -lsundials_cvode -lsundials_nvecserial
LIBS += $(CANTERA_LIBS)
LIBS += -L$(WORK)/ceq -lceq
LIBS += -L$(WORK)/lapack-3.6.1/ -llapack -lrefblas
INCS  = -I.
INCS += -I$(WORK)/sundials-2.6.1/sundials-2.6.1/include/
INCS += $(CANTERA_INCLUDES)
INCS += -I$(WORK)/ceq/
INCS += -I$(WORK)/eigen-3.3.1/eigen-3.3.1/include/eigen3/
LINK  = -lifcore

all: preproc postproc redkin

preproc: preprocessor.o
	$(CXX) $(FLAG) -o preproc preprocessor.o

postproc: postprocessor.o
	$(CXX) $(FLAG) -o postproc postprocessor.o

redkin: ceqct.o driver.o
	$(CXX) $(FLAG) -o redkin ceqct.o driver.o $(LIBS) $(LINK)

preprocessor.o: preprocessor.cpp
	$(CXX) $(FLAG) -c $<

postprocessor.o: postprocessor.cpp
	$(CXX) $(FLAG) -c $<

ceqct.o: ceqct.f90
	$(F90) $(FLAG) -c $< $(INCS)

driver.o: driver.cpp ceqct.o
	$(CXX) $(FLAG) -c $< $(INCS)

clean:
	rm *.o

real-clean:
	rm *~
