#ifndef MECH_IDEALGASMIX
#define MECH_IDEALGASMIX

#include "GasThermo.h"
#include "GasKinetics.h"
#include "mech_defs.h"

namespace mech
{

  class IdealGasMix : public GasThermo, public GasKinetics
  {
  public:

    IdealGasMix(double& T, double& p, vector<double>& x) {
      addThermo(this);
      setState_TPX(T, p, x);
    }
    
  };

}

#endif
