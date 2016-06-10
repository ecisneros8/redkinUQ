#ifndef MECH_GASKINETICS
#define MECH_GASKINETICS

#include <cmath>
#include <vector>
#include <iostream>
#include "GasThermo.h"
#include "mech_defs.h"
#include <cppad/cppad.hpp>

using namespace std;

using std::vector;
using CppAD::AD;
using CppAD::Value;
using CppAD::Var2Par;

namespace mech
{

  class GasKinetics
  {
  public:

    // Constructor
  GasKinetics() :
      m_nrev(18),
      m_nirrev(1),
      m_ntbrxn(5)
	{
	  a_z.resize(m_kk);
	  a_fz.resize(m_kk);
	}

      /*
       * Non-templated routines:
       * These routines are inteded to be either for computational 
       * (e.g. getJacobian) or for informational purposes.
       */

      // Add thermo pointer
      void addThermo(GasThermo* thermo);

      // Source
      void getSource(vector<double>& z, vector<double>& fz);

      // Jacobian
      void getJacobian(vector<double>& z, vector<double>& dfz);

      GasThermo& thermo();

      // Informational routines
      int nReactions();
      int nReversibleReactions();
      int nIrreversibleReactions();
      int nThreeBodyReactions();

      /*
       * Templated routines:
       * Species rate laws have been hardcoded to
       * aid in the computation of the Jacobian through
       * Automatic Differentiation package CppAD.
       */
      
      // Rates of Progress
      template <class Type>
	void updateRateConstants(Type& T, vector<Type>& c,
				 vector<Type>& kfwd, vector<Type> & krev);
      
      // Net production rates
      template <class Type>
	void getNetProductionRates(Type& T, vector<Type>& z, vector<Type>& fz);
      
      /*
       * Member variables
       */
      int m_nrev;
      int m_nirrev;
      int m_ntbrxn;
      
  protected:
      GasThermo*           m_thermo;
      vector< AD<double> > a_z;
      vector< AD<double> > a_fz;
  };

  /*
   * Non-templated routines:
   * These routines are inteded to be either for computational 
   * (e.g. getJacobian) or for informational purposes.
   */

  void GasKinetics::addThermo(GasThermo* thermo) {

    m_thermo = thermo;
    
  };

  void GasKinetics::getSource(vector<double>& z, vector<double>& fz) {

    double T;
    thermo().getTemperature(z, T);
    getNetProductionRates(T, z, fz);
    
  }
  
  void GasKinetics::getJacobian(vector<double>& z, vector<double>& dfz) {

    copy(z.begin(), z.end(), a_z.begin());

    /* AD Stuff */
    CppAD::Independent(a_z);
    AD<double> a_T;
    thermo().getTemperature(a_z, a_T);
    getNetProductionRates(a_T, a_z, a_fz);
    CppAD::ADFun<double> tapef(a_z, a_fz);

    dfz = tapef.Jacobian(z);
    
  };

  GasThermo& GasKinetics::thermo() {
    return *m_thermo;
  };

  int GasKinetics::nReactions() {
    return m_ii;
  };

  int GasKinetics::nReversibleReactions() {
    return m_nrev;
  };

  int GasKinetics::nIrreversibleReactions() {
    return m_nirrev;
  };

  int GasKinetics::nThreeBodyReactions() {
    return m_ntbrxn;
  };

  /*
   * Templated routines:
   * Species rate laws have been hardcoded to
   * aid in the computation of the Jacobian through
   * Automatic Differentiation package CppAD.
   */
  
  template <class Type>
    void GasKinetics::updateRateConstants(Type& T, vector<Type>& C,
			     vector<Type>& kfwd, vector<Type>& krev) {
    
    Type         tlog = log(T);
    Type         rt   = 1.0 / T;
    vector<Type> keqs(m_ii, 0.0);
    
    thermo().getEquilibriumConstants(T, keqs);
    
    kfwd[0]  =  exp(2.602158320349e+01 - 8.455147086424267e+03 * rt);
    kfwd[1]  =  exp(3.923951576293e+00 + 2.6700e+00 * tlog - 3.163163134750473e+03 * rt);
    kfwd[2]  =  exp(1.151292546497e+01 + 1.6000e+00 * tlog - 1.659758603024963e+03 * rt);
    kfwd[3]  =  exp(1.422097566607e+01 + 1.1400e+00 * tlog - 5.051439226597713e+01 * rt);
    kfwd[4]  =  exp(2.821880778083e+01 + -1.0000e+00 * tlog);
    kfwd[5]  =  exp(3.762981884827e+01 + -2.0000e+00 * tlog);
    kfwd[6]  =  exp(2.639314675993e+01 + -1.0000e+00 * tlog);
    kfwd[7]  =  exp(2.846393023886e+01 + -8.0000e-01 * tlog);
    kfwd[8]  =  exp(2.573390113104e+01 - 5.051439226597713e+02 * rt);
    kfwd[9]  =  exp(2.394214166181e+01 - 3.487898513603183e+02 * rt);
    kfwd[10] =  exp(2.412446321861e+01 - 8.659610102738936e+02 * rt);
    kfwd[11] =  exp(2.361363759484e+01 - -2.044630163146693e+02 * rt);
    kfwd[12] =  exp(2.481761039917e+01);
    kfwd[13] =  exp(1.933697147583e+01 - -6.254162851978120e+02 * rt);
    kfwd[14] =  exp(3.802001648425e+01 + -2.0000e+00 * tlog);
    kfwd[15] =  exp(2.125389408801e+01 - 1.888276091847240e+03 * rt);
    kfwd[16] =  exp(2.302585092994e+01 - 1.804085438070612e+03 * rt);
    kfwd[17] =  exp(2.405547034712e+01 - 3.223299316019493e+03 * rt);
    kfwd[18] =  exp(2.240966479052e+01 - 5.051439226597713e+02 * rt);
    
    kfwd[4] = (1.0e+00 * C[0] + 1.0e+00 * C[1] + 3.5e-01 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 6.5e+00 * C[7] + 5.0e-01 * C[8]) * kfwd[4];
    kfwd[5] = (1.0e+00 * C[0] + 1.0e+00 * C[1] + 3.5e-01 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 6.5e+00 * C[7] + 5.0e-01 * C[8]) * kfwd[5];
    kfwd[6] = (1.0e+00 * C[0] + 1.0e+00 * C[1] + 3.5e-01 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 6.5e+00 * C[7] + 5.0e-01 * C[8]) * kfwd[6];
    kfwd[7] = (1.0e+00 * C[0] + 1.0e+00 * C[1] + 3.5e-01 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 6.5e+00 * C[7] + 5.0e-01 * C[8]) * kfwd[7];
    kfwd[14] = (1.0e+00 * C[0] + 1.0e+00 * C[1] + 3.5e-01 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 6.5e+00 * C[7] + 5.0e-01 * C[8]) * kfwd[14];
    
    for(int i = 0; i < m_ii; ++i) { krev[i] = kfwd[i] * keqs[i]; }
    
  };
  
  template <class Type>
    void GasKinetics::getNetProductionRates(Type& T, vector<Type>& z, vector<Type>& fz) {

    Type           W;
    Type           rho;
    vector<double> mw = thermo().molecularWeights();
    vector<Type>   C(m_kk,    0.0);
    vector<Type>   kfwd(m_ii, 0.0);
    vector<Type>   krev(m_ii, 0.0);
    vector<Type>   Rfwd(m_ii, 0.0);
    vector<Type>   Rrev(m_ii, 0.0);
    vector<Type>   Rnet(m_ii, 0.0);
    vector<Type>   omega(m_kk,0.0);

    thermo().getMeanMolecularWeightZ(z, W);
    thermo().getDensity(T, W, rho);
    
    for(int k = 0; k < m_kk; ++k) { C[k] = rho * z[k]; }
    updateRateConstants(T, C, kfwd, krev);

    Rfwd[0]  = kfwd[0] * C[1] * C[2];
    Rfwd[1]  = kfwd[1] * C[0] * C[3];
    Rfwd[2]  = kfwd[2] * C[0] * C[4];
    Rfwd[3]  = kfwd[3] * C[4] * C[4];
    Rfwd[4]  = kfwd[4] * C[1] * C[1];
    Rfwd[5]  = kfwd[5] * C[1] * C[4];
    Rfwd[6]  = kfwd[6] * C[3] * C[3];
    Rfwd[7]  = kfwd[7] * C[1] * C[2];
    Rfwd[8]  = kfwd[8] * C[1] * C[5];
    Rfwd[9]  = kfwd[9] * C[1] * C[5];
    Rfwd[10] = kfwd[10] * C[1] * C[5];
    Rfwd[11] = kfwd[11] * C[5] * C[3];
    Rfwd[12] = kfwd[12] * C[5] * C[4];
    Rfwd[13] = kfwd[13] * C[5] * C[5];
    Rfwd[14] = kfwd[14] * C[4] * C[4];
    Rfwd[15] = kfwd[15] * C[1] * C[6];
    Rfwd[16] = kfwd[16] * C[1] * C[6];
    Rfwd[17] = kfwd[17] * C[6] * C[3];
    Rfwd[18] = kfwd[18] * C[6] * C[4];

    Rrev[0]  = krev[0] * C[3] * C[4];
    Rrev[1]  = krev[1] * C[1] * C[4];
    Rrev[2]  = krev[2] * C[1] * C[7];
    Rrev[3]  = krev[3] * C[7] * C[3];
    Rrev[4]  = krev[4] * C[0];
    Rrev[5]  = krev[5] * C[7];
    Rrev[6]  = krev[6] * C[2];
    Rrev[7]  = krev[7] * C[5];
    Rrev[8]  = krev[8] * C[4] * C[4];
    Rrev[9]  = krev[9] * C[0] * C[2];
    Rrev[10] = krev[10] * C[7] * C[3];
    Rrev[11] = krev[11] * C[2] * C[4];
    Rrev[12] = krev[12] * C[7] * C[2];
    Rrev[14] = krev[14] * C[6];
    Rrev[15] = krev[15] * C[0] * C[5];
    Rrev[16] = krev[16] * C[7] * C[4];
    Rrev[17] = krev[17] * C[5] * C[4];
    Rrev[18] = krev[18] * C[7] * C[5];
    
    for(int i = 0; i < m_ii; ++i) { Rnet[i] = Rfwd[i] - Rrev[i]; }

    omega[0] = (  + Rnet[4] + Rnet[9] + Rnet[15] - Rnet[1] - Rnet[2] );
    omega[1] = (  + Rnet[1] + Rnet[2] - Rnet[0] - Rnet[4] - Rnet[4] - Rnet[5] - Rnet[7]
		  - Rnet[8] - Rnet[9] - Rnet[10] - Rnet[15] - Rnet[16] );
    omega[2] = (  + Rnet[6] + Rnet[9] + Rnet[11] + Rnet[12] + Rnet[13] - Rnet[0] - Rnet[7] );
    omega[3] = (  + Rnet[0] + Rnet[3] + Rnet[10] - Rnet[1] - Rnet[6] - Rnet[6] - Rnet[11]
      - Rnet[17] );
    omega[4] = (  + Rnet[0] + Rnet[1] + Rnet[8] + Rnet[8] + Rnet[11] + Rnet[16] + Rnet[17]
		  - Rnet[2] - Rnet[3] - Rnet[3] - Rnet[5] - Rnet[12] - Rnet[14] - Rnet[14]
		  - Rnet[18] );
    omega[5] = (  + Rnet[7] + Rnet[15] + Rnet[17] + Rnet[18] - Rnet[8] - Rnet[9] - Rnet[10]
		  - Rnet[11] - Rnet[12] - Rnet[13] - Rnet[13] );
    omega[6] = (  + Rnet[14] + Rnet[13] - Rnet[15] - Rnet[16] - Rnet[17] - Rnet[18] );
    omega[7] = (  + Rnet[5] + Rnet[2] + Rnet[3] + Rnet[10] + Rnet[12] + Rnet[16] + Rnet[18] );

    for(int k = 0; k < m_kk; ++k) { fz[k] = omega[k] / rho; }

  };

}
    
#endif
