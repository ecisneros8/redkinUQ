include $(WORK)/cantera-2.2.1/cantera-2.2.1/include/cantera/Cantera.mak

CXX   = mpicxx -std=c++11
F90   = mpif90
FLAG  = -g -O0 -w #-O2 -w
LIBS  = $(CANTERA_LIBS)
LIBS += -L$(WORK)/ceq -lceq
LIBS += -L$(WORK)/lapack-3.6.0/ -llapack -lrefblas
INCS  = -I.
INCS += $(CANTERA_INCLUDES)
INCS += -I$(WORK)/ceq/
INCS += -I$(WORK)/eigen-3.3.1/eigen-3.3.1/include/eigen3/
LINK  = -lifcore

all: redkin

redkin: ceqct.o driver.o
	$(CXX) $(FLAG) -o redkin ceqct.o driver.o $(LIBS) $(LINK)

ceqct.o: ceqct.f90
	$(F90) $(FLAG) -c $< $(INCS)

driver.o: driver.cpp ceqct.o
	$(CXX) $(FLAG) -c $< $(INCS)

clean:
	rm *.o

real-clean:
	rm *~
