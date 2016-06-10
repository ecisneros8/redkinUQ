CXX   = mpicxx -std=c++0x
F90   = mpif90
FLAG  = -g -ftrapuv -fp-stack-check -fstack-protector-all -traceback -C -O0 #-O2 -w
LIBS  = -L$(WORK)/CEQ -lceq
LIBS += -L$(WORK)/lapack-3.6.0/ -llapack -lrefblas
LIBS += -L$(WORK)/sundials-2.6.1/sundials-2.6.1/lib -lsundials_cvode -lsundials_nvecserial
INCS  = -I$(WORK)/CEQ
INCS += -I$(WORK)/sundials-2.6.1/sundials-2.6.1/include/
INCS += -I$(WORK)/cppad-2.0.0/cppad-2.0.0/include/
INCS += -I$(WORK)/eigen-3.2.8/eigen-3.2.8/include/eigen3/
LINK  = -lifcore

all: simex

simex: localCEQ.o driver.o
	$(CXX) $(FLAG) -o simex localCEQ.o driver.o $(LIBS) $(LINK)

localCEQ.o: localCEQ.f90
	$(F90) $(FLAG) -c $< $(INCS)

driver.o: driver.cpp localCEQ.o
	$(CXX) $(FLAG) -c $< $(INCS)

clean:
	rm *.o

dist-clean:
	rm *~
