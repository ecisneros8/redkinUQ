#include "uqSolver.h"

using namespace std;

/* run rcce */
int main(int argc, char **argv) {

  /* declarations */
  cout.precision(6);
  cout.setf(ios::scientific);

  std::string    fuel    = "CH4";
  int            mod;
  int            vio  = 1;
  int            ver  = 0;
  double         oar  = 0.21;
  double         phi  = 0.50;
  double         nux  = 0.00;
  double         t0   = 1.0e-12;
  double         tf   = 1.0e-03;
  double         dt   = 1.0e-09;
  double         Ti   = 1644.0;
  double         p    = 3.76 * OneAtm;
  vector<double> xi(m_kk,0.0);

  /* MPI INIT */
  MPI_Init(&argc, &argv);
  mpiEnvironment myEnv(MPI_COMM_WORLD,argv[1]);

  /* ideal gas mixture */
  std::string fuel = "H2O";
  std::string mech = "sanDiego";
  std::string ctif = "ctis/"+mech+".cti";
  Cantera::IdealGasMix myGas(ctif,"gas");

  /* uq solver */
  uqSolver solver(myEnv,myGas);
  solver.run(fuel, vio, ver, t0, tf, dt);

  MPI_Finalize();
  return(0);

}
