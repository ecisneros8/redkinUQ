#ifndef MECH_INTEGRATOR
#define MECH_INTEGRATOR

#include <iostream>

using namespace std;

namespace mech
{

  class Integrator
  {
  public:

    Integrator() {};
    
    virtual int integratorType() {
      cout << "integratorType called from base class" << endl;
      return(-1);
    };

    virtual IdealGasMix& gas() {
      cout << "gas called from base class" << endl;
    }

    virtual void initialize(int& vio, int& ver, double& t0, double& dt, double& rtol,
			    vector<double>& atol, vector<double>& z0,
			    vector<string>& trackSpecies) {
      cout << "initialize called from base class" << endl;
    };

    virtual void initialize(int& vio, int& ver, int& ncs, 
			    double& t0, double& tf, double& dt) {
      cout << "initialize called from base class" << endl;
    };

    virtual void setInitialConditions(double& To, vector<int>& csi) {
      cout << "setInitialConditions called from base class" << endl;
    };
    
    virtual void integrate(double& RTM, vector<double>& QOI) {
      cout << "integrate called from base class" << endl;
    };
    
  };

  IdealGasMix* m_gas;
  
}

#endif
