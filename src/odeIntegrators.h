#include "CVodeInt.h"
#include "rcceInt.h"

namespace mech
{

  Integrator* newIntegrator(IdealGasMix& gas, const string& iType) {
    if (iType == "CVODE") {
      return new CVodeInt(gas);
    } else if (iType == "RCCE") {
      return new rcceInt(gas);
    }
    return(0);
  };
  
}
