#include "uqSolver.h"

using namespace std;
using namespace mech;

/* run rcce */
int main(int argc, char **argv) {

  /* declarations */
  cout.precision(6);
  cout.setf(ios::scientific);

  int            vio  = 1;
  int            ver  = 0;
  double         t0   = 1.0e-12;
  double         tf   = 1.0e-02;
  double         dt   = 1.0e-08;
  double         Ti   = 1000.0;
  double         p    = OneAtm;
  vector<double> xi(m_kk,0.0);

  /* MPI INIT */
  MPI_Init(&argc, &argv);
  mpiEnvironment myEnv(MPI_COMM_WORLD,argv[1]);

  /* set initial conditions */
  if(mech_name == "MAASPOPE") {
    xi[0] = 0.30e0; // H2
    xi[2] = 0.15e0; // O2
    xi[8] = 0.55e0; // N2
  } else if(mech_name == "GRI12") {
    xi[3]  = 0.190e0; // O2
    xi[13] = 0.095e0; // CH4
    xi[30] = 0.715e0; // N2
  }

  IdealGasMix myGas(Ti, p, xi);

  uqSolver solver(myEnv,myGas);
  solver.run(vio, ver, t0, tf, dt);

  MPI_Finalize();
  return(0);

  //mpiEnvironment myEnv(MPI_COMM_WORLD,argv[1]);
  //IdealGasMix myGas(Ti, p, xi);
  //mpiEnvironment* myEnv = new mpiEnvironment(MPI_COMM_WORLD,argv[1]);
  //IdealGasMix* myGas = new IdealGasMix(Ti, p, xi);

  // int    ncs = 2;
  // double rtm;
  // vector<int> csi(2,0);
  // vector<double> Qk(2,0.0);
  // csi[0] = 1;
  // csi[1] = 4;
  // Integrator* integ = newIntegrator(myGas,"RCCE");
  // integ->initialize(vio, ver, ncs, t0, tf, dt);
  // integ->setInitialConditions(Ti, csi);
  // integ->integrate(rtm, Qk);

}
