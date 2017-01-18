#include "uqSolver.h"

using namespace std;

/* run rcce */
int main(int argc, char **argv) {

  /* declarations */
  cout.precision(6);
  cout.setf(ios::scientific);

  std::string    fuel = "CH4";
  std::string    mech = "gri30";
  std::string    ctif = "ctis/"+mech+".cti";
  int            mod;
  int            vio  = 0;
  int            ver  = 0;
  double         t0   = 1.0e-14;
  double         tf   = 1.0e-03;
  double         dt   = 1.0e-09;  

  /* MPI INIT */
  MPI_Init(&argc, &argv);
  mpiEnvironment myEnv(MPI_COMM_WORLD,argv[1],argv[2]);

  /* ideal gas mixture */
  if(myEnv.worldRank() == 0) { std::cout << "Reading cti file..." << std::endl; }
  Cantera::IdealGasMix myGas(ctif,"gas");

  /* uq solver */
  if(myEnv.worldRank() == 0) { std::cout << "Initializing uqSolver..." << std::endl; }
  uqSolver solver(myEnv,myGas);

  if(myEnv.worldRank() == 0) { std::cout << "Run uqSolver..." << std::endl; }
  solver.run(fuel, vio, ver, t0, tf, dt);

  MPI_Finalize();
  return(0);

  /* run single realization at single scenario */
  //Cantera::IdealGasMix myGas(ctif,"gas");
  //vector<int> csi(7);
  //csi[0] = myGas.speciesIndex("H2O2");
  //csi[1] = myGas.speciesIndex("CH2");
  //csi[2] = myGas.speciesIndex("CH3");
  //csi[3] = myGas.speciesIndex("CH4");
  //csi[4] = myGas.speciesIndex("CH2O");
  //csi[5] = myGas.speciesIndex("CH2OH");
  //csi[6] = myGas.speciesIndex("CH3O");
  //int    flag;
  //double p   = 3.76 * Cantera::OneAtm;
  //double T0  = 1644;
  //double phi = 0.5;
  //double mout;
  //rcceInteg integ;
  //integ.setIdealGasMix(myGas);
  //integ.initialize(vio, ver, t0, tf, dt);
  //integ.setInitialConditions(fuel, p, T0, phi, csi);
  //cout << "Integrating:" << endl;
  //integ.integrate(flag, mout);
  //cout << "Done!" << endl;

}
