#include "uqSolver.h"

using namespace std;
using namespace mech;

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

  /* set initial conditions */
  if(fuel == "H2O") {
    
    /* hydrogen-air combustion */
    /* nux */
    nux   = 0.50e0;
    /* H2 */
    xi[0] = oar * phi / (nux + oar * phi);
    /* O2 */
    xi[2] = nux * xi[0] / phi;
    /* N2 */
    xi[8] = 1.0 - (1.0 + nux/phi) * xi[0];

  } else if(fuel == "CH4") {

    /* methane-air combustion */
    /* nux */
    nux    = 2.00e0;
    /* CH4 */
    xi[13] = oar * phi / (nux + oar * phi);
    /* O2  */
    xi[3]  = nux * xi[13] / phi;
    /* N2  */
    xi[47] = 1.0 - (1.0 + nux/phi) * xi[13];

  }

  IdealGasMix myGas(Ti, p, xi);

  /* uq solver */
  uqSolver solver(myEnv,myGas);
  solver.run(fuel, vio, ver, t0, tf, dt);

  MPI_Finalize();
  return(0);

}
