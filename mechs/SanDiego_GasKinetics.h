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
        m_nrev(24),
	m_nirrev(0),
	m_nfall(2),
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
    
    kfwd[0] =  exp(3.119206719853e+01 + -7.0000e-01 * tlog - 8.589852132466873e+03 * rt);
    kfwd[1] =  exp(3.923951576293e+00 + 2.6700e+00 * tlog - 3.165568582001234e+03 * rt);
    kfwd[2] =  exp(1.397251430677e+01 + 1.3000e+00 * tlog - 1.829342634203600e+03 * rt);
    kfwd[3] =  exp(6.551080335043e+00 + 2.3300e+00 * tlog - 7.320978707690542e+03 * rt);
    kfwd[4] =  exp(2.789338538040e+01 + -1.0000e+00 * tlog);
    kfwd[5] =  exp(3.822765584902e+01 + -2.0000e+00 * tlog);
    kfwd[6] =  exp(2.254296467486e+01 + -5.0000e-01 * tlog);
    kfwd[7] =  exp(2.918070902396e+01 + -1.0000e+00 * tlog);
    kfwd[8] =  exp(2.280270737863e+01);
    kfwd[10] =  exp(2.498312483765e+01 - 1.484160953719423e+02 * rt);
    kfwd[11] =  exp(2.353266853231e+01 - 4.140977442184744e+02 * rt);
    kfwd[12] =  exp(2.415725304143e+01 - 8.659610102738936e+02 * rt);
    kfwd[13] =  exp(2.371899811050e+01);
    kfwd[14] =  exp(2.683251341971e+01 - 5.500055138864605e+03 * rt);
    kfwd[15] =  exp(2.273816885749e+01 - -5.508474204242268e+02 * rt);
    kfwd[17] =  exp(1.908336871703e+01 - -7.090055771617505e+02 * rt);
    kfwd[18] =  exp(2.535799482518e+01 - 5.556583149257484e+03 * rt);
    kfwd[19] =  exp(2.385876005288e+01 - 4.000619595102850e+03 * rt);
    kfwd[20] =  exp(2.302585092994e+01 - 1.804085438070612e+03 * rt);
    kfwd[21] =  exp(2.505268252135e+01 - 3.659887992032581e+03 * rt);
    kfwd[22] =  exp(2.127715095017e+01 - 1.599622421755942e+02 * rt);
    kfwd[23] =  exp(9.172638504792e+00 + 2.0000e+00 * tlog - 2.008548454385281e+03 * rt);

    getFalloffRates(T, C, kfwd);
    
    kfwd[4] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 1.2e+01 * C[7] + 1.0e+00 * C[8]) * kfwd[4];
    kfwd[5] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 1.2e+01 * C[7] + 1.0e+00 * C[8]) * kfwd[5];
    kfwd[6] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 1.2e+01 * C[7] + 1.0e+00 * C[8]) * kfwd[6];
    kfwd[7] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 1.2e+01 * C[7] + 1.0e+00 * C[8]) * kfwd[7];
    kfwd[8] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 1.2e+01 * C[7] + 1.0e+00 * C[8]) * kfwd[8];
    
    for(int i = 0; i < m_ii; ++i) { krev[i] = kfwd[i] * keqs[i]; }
    
  };

  template<class Type>
    void GasKinetics::getFalloffRates(Type& T, vector<Type>& C, vector<Type>& kfwd) {

    Type         tlog = log(T);
    Type         rt   = 1.0 / T;
    Type         lpr;
    Type         cc;
    Type         nn;
    Type         f1;
    Type         lgf;
    vector<Type> pr(m_nfall,   0.0);
    vector<Type> khi(m_nfall,  0.0);
    vector<Type> klo(m_nfall,  0.0);
    vector<Type> work(m_nfall, 0.0);
    
    khi[0] =  exp(2.226013305655e+01 + 4.4000e-01 * tlog);
    khi[1] =  exp(2.528239208443e+01 + -2.7000e-01 * tlog);

    klo[0] =  exp(3.168280606373e+01 + -1.4000e+00 * tlog);
    klo[1] =  exp(4.476434744662e+01 + -3.2000e+00 * tlog);

    pr[0] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 1.6e+01 * C[7] + 1.0e+00 * C[8]) * klo[0] / khi[0];
    pr[1] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 6.0e+00 * C[7] + 1.0e+00 * C[8]) * klo[1] / khi[1];

    work[0] = 0.5;
    work[1] = 0.5;

    for(int r = 0; r < m_nfall; ++r) {
      lpr     =  log10(pr[r]);
      cc      = -0.40 - 0.67 * log10(work[r]);
      nn      =  0.75 - 1.27 * log10(work[r]);
      f1      =  (lpr + cc)/(nn - 0.14 * (lpr + cc));
      work[r] =  log10(work[r])/(1 + f1 * f1);
      work[r] =  pow(10.0, work[r]);
      work[r] =  (pr[r] * work[r])/(1 + pr[r]);
    }

    kfwd[9] = khi[0] * work[0];
    kfwd[16] = khi[1] * work[1];

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

    Rfwd[0] = kfwd[0] * C[1] * C[2];
    Rfwd[1] = kfwd[1] * C[0] * C[3];
    Rfwd[2] = kfwd[2] * C[0] * C[4];
    Rfwd[3] = kfwd[3] * C[7] * C[3];
    Rfwd[4] = kfwd[4] * C[1] * C[1];
    Rfwd[5] = kfwd[5] * C[1] * C[4];
    Rfwd[6] = kfwd[6] * C[3] * C[3];
    Rfwd[7] = kfwd[7] * C[1] * C[3];
    Rfwd[8] = kfwd[8] * C[3] * C[4];
    Rfwd[9] = kfwd[9] * C[1] * C[2];
    Rfwd[10] = kfwd[10] * C[1] * C[5];
    Rfwd[11] = kfwd[11] * C[1] * C[5];
    Rfwd[12] = kfwd[12] * C[1] * C[5];
    Rfwd[13] = kfwd[13] * C[5] * C[3];
    Rfwd[14] = kfwd[14] * C[5] * C[4];
    Rfwd[15] = kfwd[15] * C[5] * C[4];
    Rfwd[16] = kfwd[16] * C[4] * C[4];
    Rfwd[17] = kfwd[17] * C[5] * C[5];
    Rfwd[18] = kfwd[18] * C[5] * C[5];
    Rfwd[19] = kfwd[19] * C[1] * C[6];
    Rfwd[20] = kfwd[20] * C[1] * C[6];
    Rfwd[21] = kfwd[21] * C[6] * C[4];
    Rfwd[22] = kfwd[22] * C[6] * C[4];
    Rfwd[23] = kfwd[23] * C[6] * C[3];

    Rrev[0] = krev[0] * C[3] * C[4];
    Rrev[1] = krev[1] * C[1] * C[4];
    Rrev[2] = krev[2] * C[1] * C[7];
    Rrev[3] = krev[3] * C[4] * C[4];
    Rrev[4] = krev[4] * C[0];
    Rrev[5] = krev[5] * C[7];
    Rrev[6] = krev[6] * C[2];
    Rrev[7] = krev[7] * C[4];
    Rrev[8] = krev[8] * C[5];
    Rrev[9] = krev[9] * C[5];
    Rrev[10] = krev[10] * C[4] * C[4];
    Rrev[11] = krev[11] * C[0] * C[2];
    Rrev[12] = krev[12] * C[7] * C[3];
    Rrev[13] = krev[13] * C[2] * C[4];
    Rrev[14] = krev[14] * C[7] * C[2];
    Rrev[15] = krev[15] * C[7] * C[2];
    Rrev[16] = krev[16] * C[6];
    Rrev[17] = krev[17] * C[6] * C[2];
    Rrev[18] = krev[18] * C[6] * C[2];
    Rrev[19] = krev[19] * C[0] * C[5];
    Rrev[20] = krev[20] * C[7] * C[4];
    Rrev[21] = krev[21] * C[7] * C[5];
    Rrev[22] = krev[22] * C[7] * C[5];
    Rrev[23] = krev[23] * C[5] * C[4];
    
    for(int i = 0; i < m_ii; ++i) { Rnet[i] = Rfwd[i] - Rrev[i]; }

    omega[0] =  + Rnet[4] + Rnet[11] + Rnet[19] - Rnet[1] - Rnet[2] ;
    omega[1] =  + Rnet[1] + Rnet[2] - Rnet[0] - Rnet[4] - Rnet[4] - Rnet[5] - Rnet[7]
      - Rnet[9] - Rnet[10] - Rnet[11] - Rnet[12] - Rnet[19] - Rnet[20] ;
    omega[2] =  + Rnet[6] + Rnet[11] + Rnet[13] + Rnet[14] + Rnet[15] + Rnet[17] + Rnet[18]
      - Rnet[0] - Rnet[9] ;
    omega[3] =  + Rnet[0] + Rnet[12] - Rnet[1] - Rnet[3] - Rnet[6] - Rnet[6] - Rnet[7]
      - Rnet[8] - Rnet[13] - Rnet[23] ;
    omega[4] =  + Rnet[7] + Rnet[0] + Rnet[1] + Rnet[3] + Rnet[3] + Rnet[10] + Rnet[10]
      + Rnet[13] + Rnet[20] + Rnet[23] - Rnet[2] - Rnet[5] - Rnet[8] - Rnet[14]
      - Rnet[15] - Rnet[16] - Rnet[16] - Rnet[21] - Rnet[22] ;
    omega[5] =  + Rnet[8] + Rnet[9] + Rnet[19] + Rnet[21] + Rnet[22] + Rnet[23] - Rnet[10]
      - Rnet[11] - Rnet[12] - Rnet[13] - Rnet[14] - Rnet[15] - Rnet[17] - Rnet[17]
      - Rnet[18] - Rnet[18] ;
    omega[6] =  + Rnet[16] + Rnet[17] + Rnet[18] - Rnet[19] - Rnet[20] - Rnet[21] - Rnet[22]
      - Rnet[23] ;
    omega[7] =  + Rnet[5] + Rnet[2] + Rnet[12] + Rnet[14] + Rnet[15] + Rnet[20] + Rnet[21]
      + Rnet[22] - Rnet[3] ;

    for(int k = 0; k < m_kk; ++k) { fz[k] = omega[k] / rho; }

  };

}
    
#endif
