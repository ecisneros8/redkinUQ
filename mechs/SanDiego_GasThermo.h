#ifndef MECH_GASTHERMO
#define MECH_GASTHERMO

#include <cmath>
#include <vector>
#include <iostream>
#include <numeric>
#include <map>
#include "mech_defs.h"
#include <cppad/cppad.hpp>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

using std::vector;
using CppAD::AD;
using CppAD::Value;
using CppAD::Var2Par;

namespace mech
{

  class GasThermo
  {
  public:
    
    // Constructor
  GasThermo() {
	  // resize vectors
    m_speciesNames.resize(m_kk),
      	  m_mw.resize(m_kk);
	  m_cp0_R.resize(m_kk);
	  m_h0_RT.resize(m_kk);
	  m_s0_R.resize(m_kk);
	  m_g0_RT.resize(m_kk);
	  m_z.resize(m_kk);
	  m_Emat.resize(m_kk,m_mm);

	  // species names
	  m_speciesNames[0] = "H2";
	  m_speciesNames[1] = "H";
	  m_speciesNames[2] = "O2";
	  m_speciesNames[3] = "O";
	  m_speciesNames[4] = "OH";
	  m_speciesNames[5] = "HO2";
	  m_speciesNames[6] = "H2O2";
	  m_speciesNames[7] = "H2O";
	  m_speciesNames[8] = "N2";
	  //m_speciesNames[8]  = "TM";
	  //m_speciesNames[9]  = "FO";
	  //m_speciesNames[10] = "AV";

	  // species indices 
	  m_speciesIndex["H2"]   = 0;
	  m_speciesIndex["H"]    = 1;
	  m_speciesIndex["O2"]   = 2;
	  m_speciesIndex["O"]    = 3;
	  m_speciesIndex["OH"]   = 4;
	  m_speciesIndex["HO2"]  = 5;
	  m_speciesIndex["H2O2"] = 6;
	  m_speciesIndex["H2O"]  = 7;
	  m_speciesIndex["N2"]   = 8;
	  //m_speciesNames[8]  = "TM";
	  //m_speciesNames[9]  = "FO";
	  //m_speciesNames[10] = "AV";

	  // set molecular weights
	  m_mw[0] = 2.015880e+00;
	  m_mw[1] = 1.007940e+00;
	  m_mw[2] = 3.199880e+01;
	  m_mw[3] = 1.599940e+01;
	  m_mw[4] = 1.700734e+01;
	  m_mw[5] = 3.300674e+01;
	  m_mw[6] = 3.401468e+01;
	  m_mw[7] = 1.801528e+01;
	  m_mw[8] = 2.801348e+01;

	  // set element matrix
	  m_Emat(0,0) = 2.0e+00;
	  m_Emat(0,1) = 0.0e+00;
	  m_Emat(0,2) = 0.0e+00;

	  m_Emat(1,0) = 1.0e+00;
	  m_Emat(1,1) = 0.0e+00;
	  m_Emat(1,2) = 0.0e+00;

	  m_Emat(2,0) = 0.0e+00;
	  m_Emat(2,1) = 2.0e+00;
	  m_Emat(2,2) = 0.0e+00;

	  m_Emat(3,0) = 0.0e+00;
	  m_Emat(3,1) = 1.0e+00;
	  m_Emat(3,2) = 0.0e+00;

	  m_Emat(4,0) = 1.0e+00;
	  m_Emat(4,1) = 1.0e+00;
	  m_Emat(4,2) = 0.0e+00;

	  m_Emat(5,0) = 1.0e+00;
	  m_Emat(5,1) = 2.0e+00;
	  m_Emat(5,2) = 0.0e+00;

	  m_Emat(6,0) = 2.0e+00;
	  m_Emat(6,1) = 2.0e+00;
	  m_Emat(6,2) = 0.0e+00;

	  m_Emat(7,0) = 2.0e+00;
	  m_Emat(7,1) = 1.0e+00;
	  m_Emat(7,2) = 0.0e+00;

	  m_Emat(8,0) = 0.0e+00;
	  m_Emat(8,1) = 0.0e+00;
	  m_Emat(8,2) = 2.0e+00;
	  
	}

      /*
       * Non-templated routines:
       * These routines are inteded to be either for computational 
       * (e.g. set-routines) or for informational purposes.
       */

      // State Routines
      void setState_TPX(double& T, double& p, vector<double>& x);
      void setState_TPZ(double& T, double& p, vector<double>& z);
      void setState_TRZ(double& T, double& rho, vector<double>& z);
      void setTemperature(double& T);
      void setPressure(double& p);
      void setDensity(double& rho);
      void setEnthalpyMass(double& h);
      void setMoleFractions(vector<double>& x);
      void setSpecificMoles(vector<double>& z);
      void getCombinations(int offset, int k, vector< vector<int> >& v);
      
      // Informational routines
      int            nSpecies();
      double         temperature();
      double         pressure();
      double         density();
      double         enthalpy_mass();
      double         cp_mass();
      vector<double> elementMatrix();
      vector<double> specificMoles();
      vector<double> moleFractions();
      vector<double> massFractions();
      vector<double> molecularWeights();

      /*
       * Templated routines:
       * Species thermo (NASA7) have been hardcoded to
       * aid in the computation of the Jacobian through
       * Automatic Differentiation package CppAD.
       */

      // Obtain temperature from composition and enthalpy
      template <class Type>
	void getTemperature(vector<Type>& z, Type& T);

      // Obtain density from temperature and composition (Ideal Gas Eq. of State)
      template <class Type>
	void getDensity(Type& T, vector<Type>& z, Type& rho);

      // Obtain mass-based enthalpy from temperature and composition
      template <class Type>
	void getEnthalpyMass(Type& T, vector<Type>& z, Type& h);

      // Obtain mass-based specific heat from temperature and composition
      template <class Type>
	void getSpecificHeatMass(Type& T, vector<Type>& z, Type& cpmass);

      // Obtain mean molecular weight from temperature & density
      template <class Type>
	void getMeanMolecularWeightTR(Type& T, Type& rho, Type& W);
      
      // Obtain mean molecular weight from mass fractions
      template <class Type>
	void getMeanMolecularWeightZ(vector<Type>& z, Type& W);

      // Obtain density using Ideal Gas Eq. of State
      template <class Type>
	void getDensity(Type& T, Type& W, Type& rho);

      // Obtain equilibrium constants 
      template <class Type>
	void getEquilibriumConstants(Type& T, vector<Type>& keqs);

      // Species specific heats at constant pressure
      template <class Type>
	void getSpecificHeats_R(Type& T, vector<Type>& cp0_R);

      // Species normalized enthalpies
      template <class Type>
	void getEnthalpies_RT(Type& T, vector<Type>& h0_RT);

      // Derivatives of NASA7 polynomials
      template <class Type>
	void getEnthalpiesDerivatives(Type& T, vector<Type>& dh0dT);

      // Species normalized entropies
      template <class Type>
	void getEntropies_R(Type& T, vector<Type>& s0_R);

      // Species normalized Gibbs free energy
      template <class Type>
	void getGibbsFunctions_RT(Type& T, vector<Type>& g0_RT);
      
      /*
       * Member variables
       */
      double          m_told;
      double          m_temp;
      double          m_pres;
      double          m_dens;
      double          m_enthalpy;
      double          m_cpmass;
      double          m_mmw;
      vector<string>  m_speciesNames;
      vector<int>     m_combination;
      vector<double>  m_mw;
      vector<double>  m_cp0_R;
      vector<double>  m_h0_RT;
      vector<double>  m_s0_R;
      vector<double>  m_g0_RT;
      vector<double>  m_z;
      MatrixXd        m_Emat;
      map<string,int> m_speciesIndex;
      
  };
     
  /*
   * Non-templated routines:
   * These routines are inteded to be either for computational 
   * (e.g. set-routines) or for informational purposes.
   */

  // State routines
  
  void GasThermo::setState_TPX(double& T, double& p, vector<double>& x) {
    // Set TPX
    setTemperature(T);
    setPressure(p);
    setMoleFractions(x);
    // Get rho, bulk mass-based enthalpy and specific heat
    getDensity(m_temp, m_mmw, m_dens);
    getEnthalpyMass(m_temp, m_z, m_enthalpy);
    getSpecificHeatMass(m_temp, m_z, m_cpmass);
  };

  void GasThermo::setState_TPZ(double& T, double& p, vector<double>& z) {
    // Set TPZ
    setTemperature(T);
    setPressure(p);
    setSpecificMoles(z);
    // Get rho, bulk mass-based enthalpy and specific heat
    getDensity(m_temp, m_mmw, m_dens);
    getEnthalpyMass(m_temp, m_z, m_enthalpy);
    getSpecificHeatMass(m_temp, m_z, m_cpmass);
  };

  void GasThermo::setState_TRZ(double& T, double& rho, vector<double>& z) {
    // Set TPZ
    setTemperature(T);
    setDensity(rho);
    setSpecificMoles(z);
    // Get rho, bulk mass-based enthalpy and specific heat
    getEnthalpyMass(m_temp, m_z, m_enthalpy);
    getSpecificHeatMass(m_temp, m_z, m_cpmass);
  };

  void GasThermo::setTemperature(double& T) {
    m_temp = T;
    m_told = T;
  };

  void GasThermo::setPressure(double& p) {
    m_pres = p;
  };

  void GasThermo::setDensity(double& rho) {
    m_dens = rho;
  };

  void GasThermo::setEnthalpyMass(double& h) {
    m_enthalpy = h;
  };

  void GasThermo::setMoleFractions(vector<double>& x) {
    
    double         norm   = 0.0;
    double         sum    = 0.0;
    double         invSum = 0.0;
    vector<double> zm(m_kk, 0.0);

    // Ignore negative mole fractions
    for (int k = 0; k < m_kk; k++) {
        double xk = max(x[k], 0.0);
        m_z[k]  = xk;
        norm   += xk;
        sum    += m_mw[k] * xk;
    }
    
    invSum = 1.0/sum;
    m_mmw  = sum/norm;
    for (int k=0; k < m_kk; k++) { m_z[k] = m_z[k] * invSum; }
    
  };

  void GasThermo::setSpecificMoles(vector<double>& z) {

    double         norm;
    double         invNorm;
    vector<double> zm(m_kk, 0.0);

    // Ignore negative mass fractions
    for (int k = 0; k < m_kk; k++) { m_z[k] = max(z[k], 0.0); }

    // Set mean molecular weight
    m_mmw = 1.0 / accumulate(m_z.begin(), m_z.end(), 0.0);
    
  };

  void GasThermo::getCombinations(int offset, int k, vector< vector<int> >& v) {
    if (k == 0) {
      v.push_back(m_combination);
      return;
    }
    for (int i = offset; i <= m_kk-2; ++i) {  
      string speciesName  = m_speciesNames[i];
      int    speciesIndex = m_speciesIndex[speciesName];
      m_combination.push_back(speciesIndex);
      getCombinations(i+1, k-1, v);
      m_combination.pop_back();
    }
  }

  // Informational routines

  int GasThermo::nSpecies() {
    return m_kk;
  };

  double GasThermo::temperature() {
    if(m_temp == m_told) {
      return m_temp;
    } else {
      getTemperature(m_z, m_temp);
      return m_temp;
    }
  };

  double GasThermo::pressure() {
    return m_pres;
  };

  double GasThermo::density() {
    return m_dens;
  };

  double GasThermo::enthalpy_mass() {
    return m_enthalpy;
  };

  double GasThermo::cp_mass() {
    return m_cpmass;
  };

  vector<double> GasThermo::elementMatrix() {
    vector<double> Emat(m_kk*m_mm,0);
    for(int k = 0; k < m_kk; ++k) {
      for(int m = 0; m < m_mm; ++m) {
	Emat[m*m_kk+k] = m_Emat(k,m);
      }
    }
    return Emat;
  }

  vector<double> GasThermo::specificMoles() {
    return m_z;
  };

  vector<double> GasThermo::moleFractions() {
    vector<double> x(m_kk, 0.0);
    for(int k = 0; k < m_kk; ++k) { x[k] = m_mmw * m_z[k]; }
    return x;
  };

  vector<double> GasThermo::massFractions() {
    vector<double> y(m_kk, 0.0);
    for(int k = 0; k < m_kk; ++k) { y[k] = m_mw[k] * m_z[k]; }
    return y;
  }
  
  vector<double> GasThermo::molecularWeights() {
    return m_mw;
  };

  /*
   * Templated routines:
   * Species thermo (NASA7) have been hardcoded to
   * aid in the computation of the Jacobian through
   * Automatic Differentiation package CppAD.
   */
  
  template <class Type>
    void GasThermo::getTemperature(vector<Type>& z, Type& T) {

    double       tol   = 1.0e-06;
    double       Told  = m_told;
    int          niter = 50;
    Type         RT;
    Type         To;
    Type         Tp;
    Type         dT = 1.0;
    Type         FZ = 0.0;
    Type         JZ = 0.0;
    vector<Type> hi(m_kk,    0.0);
    vector<Type> dhidT(m_kk, 0.0);
    
    To = Told;
    Tp = Told;
    
    for(int i = 0; i < niter; ++i) { 

      RT = GasConstant * To;
      getEnthalpies_RT(To, hi);
      getEnthalpiesDerivatives(To, dhidT);
      
      for(int k = 0; k < m_kk; ++k) { hi[k]    = RT * hi[k]; }
      for(int k = 0; k < m_kk; ++k) { dhidT[k] = GasConstant * dhidT[k]; }
      
      for(int k = 0; k < m_kk; ++k) { FZ -= hi[k] * z[k]; }
      for(int k = 0; k < m_kk; ++k) { JZ -= dhidT[k] * z[k]; }

      FZ =  m_enthalpy + FZ;
      JZ =  1.0 / JZ;
      dT = -FZ * JZ;
      Tp =  To + dT;
      To =  Tp;
      
      if( (fabs(dT) < tol)) {
	T      = Tp;
	return;
      }

      FZ = 0.0;
      JZ = 0.0;

    }

    T      = Tp;

  };

  template <class Type>
    void GasThermo::getDensity(Type& T, Type& W, Type& rho) {

    rho = (m_pres * W)/(GasConstant * T);
    
  }

  template <class Type>
    void GasThermo::getEnthalpyMass(Type& T, vector<Type>& z, Type& h) {

    Type         RT;
    vector<Type> hi(m_kk, 0.0);

    RT = GasConstant * T;
    h  = 0.0;
    getEnthalpies_RT(T, hi);

    for(int k = 0; k < m_kk; ++k) { hi[k] = RT * hi[k]; }
    for(int k = 0; k < m_kk; ++k) { h += hi[k] * z[k]; }
    
  };

  template <class Type>
    void GasThermo::getSpecificHeatMass(Type& T, vector<Type>& z, Type& cpmass) {

    vector<Type> cpi(m_kk, 0.0);

    cpmass = 0.0;
    getSpecificHeats_R(T, cpi);
    for(int k = 0; k < m_kk; ++k) { cpi[k]  = GasConstant * cpi[k]; }
    for(int k = 0; k < m_kk; ++k) { cpmass += cpi[k] * z[k]; }

  };

  template <class Type>
    void GasThermo::getMeanMolecularWeightTR(Type& T, Type& rho, Type& W) {

    double p  = m_pres;
    Type   RT = GasConstant * T;

    W = rho * RT / p;
    
  };

  template <class Type>
    void GasThermo::getMeanMolecularWeightZ(vector<Type>& z, Type& W) {

    W   = 0.0;
    for(int k = 0; k < m_kk; ++k) { W += z[k]; }
    W   = 1.0/W;
    
  };

  template <class Type>
    void GasThermo::getSpecificHeats_R(Type& T, vector<Type>& cp0_R) {

    Type tt0 = T;
    Type tt1 = T * tt0;
    Type tt2 = T * tt1;
    Type tt3 = T * tt2;

    if(tt0 < 1.0000e+03) {
      cp0_R[0] = 3.298007330271616e+00 + 8.258468556414628e-04 * tt0 + -8.167383830175186e-07 * tt1 + -9.201733687992775e-11 * tt2 + 4.124019474420322e-13 * tt3;
    } else {
      cp0_R[0] = 2.991423098020802e+00 + 7.000645004017671e-04 * tt0 + -5.633838501075994e-08 * tt1 + -9.231549368983405e-12 * tt2 + 1.582749414838935e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[1] = 2.500000000000000e+00 + 0.000000000000000e+00 * tt0 + 0.000000000000000e+00 * tt1 + 0.000000000000000e+00 * tt2 + 0.000000000000000e+00 * tt3;
    } else {
      cp0_R[1] = 2.500000000000000e+00 + 0.000000000000000e+00 * tt0 + 0.000000000000000e+00 * tt1 + 0.000000000000000e+00 * tt2 + 0.000000000000000e+00 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[2] = 3.212715989395561e+00 + 1.129189930372770e-03 * tt0 + -5.802181361464514e-07 * tt1 + 1.319054799445732e-09 * tt2 + -8.789124843560576e-13 * tt3;
    } else {
      cp0_R[2] = 3.697578198388313e+00 + 6.135198631579162e-04 * tt0 + -1.258843839277266e-07 * tt1 + 1.775286094678163e-11 * tt2 + -1.136439853728577e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[3] = 2.946268214996042e+00 + -1.636929514520428e-03 * tt0 + 2.417689047111522e-06 * tt1 + -1.599085904220948e-09 * tt2 + 3.875781637558715e-13 * tt3;
    } else {
      cp0_R[3] = 2.542059346810753e+00 + -2.755084819709930e-05 * tt0 + -3.102750735284210e-09 * tt1 + 4.551064604029875e-12 * tt2 + -4.368053603400201e-16 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[4] = 3.637469954817199e+00 + 1.835107588557445e-04 * tt0 + -1.671890821934399e-06 * tt1 + 2.382387984290979e-09 * tt2 + -8.412272879086942e-13 * tt3;
    } else {
      cp0_R[4] = 2.882730010358268e+00 + 1.013973808582920e-03 * tt0 + -2.276873525826253e-07 * tt1 + 2.174674530723264e-11 * tt2 + -5.126235449665062e-16 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[5] = 2.980113629643359e+00 + 4.995528635916930e-03 * tt0 + -3.787831766740593e-06 * tt1 + 2.350618553607452e-09 * tt2 + -8.074751784591024e-13 * tt3;
    } else {
      cp0_R[5] = 4.072190882104516e+00 + 2.131295860812287e-03 * tt0 + -5.308143605520907e-07 * tt1 + 6.112265234996528e-11 * tt2 + -2.841160746631655e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[6] = 3.388694479335974e+00 + 6.569676580383018e-03 * tt0 + -1.497067806981738e-07 * tt1 + -4.624463890229065e-09 * tt2 + 2.470987660509708e-12 * tt3;
    } else {
      cp0_R[6] = 4.573166050599029e+00 + 4.336137813959980e-03 * tt0 + -1.474689772850899e-06 * tt1 + 2.348905095864344e-10 * tt2 + -1.431655199308401e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[7] = 3.386728712311508e+00 + 3.475860768086567e-03 * tt0 + -6.357077086401517e-06 * tt1 + 6.971269292764603e-09 * tt2 + -2.507661553195187e-12 * tt3;
    } else {
      cp0_R[7] = 2.672145205625812e+00 + 3.056292898549634e-03 * tt0 + -8.730259974096839e-07 * tt1 + 1.200996455260830e-10 * tt2 + -6.391618725872017e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[8] = 3.298616281368291e+00 + 1.408708904146039e-03 * tt0 + -3.964481129109908e-06 * tt1 + 5.642920880408571e-09 * tt2 + -2.445407041148433e-12 * tt3;
    } else {
      cp0_R[8] = 2.926639911210682e+00 + 1.487977101178227e-03 * tt0 + -5.684761849244810e-07 * tt1 + 1.009704225872734e-10 * tt2 + -6.753354387142974e-15 * tt3;
    };


  };

  template <class Type>
    void GasThermo::getEnthalpies_RT(Type& T, vector<Type>& h0_RT) {

    Type tt0 = T;
    Type tt1 = T * tt0;
    Type tt2 = T * tt1;
    Type tt3 = T * tt2;
    Type tt4 = 1.0 / T;

    if(tt0 < 1.0000e+03) {
      h0_RT[0] = 3.298007330271616e+00 + 8.258468556414628e-04 * tt0 * 0.50 + -8.167383830175186e-07 * tt1 * OneThird + -9.201733687992775e-11 * tt2 * 0.25 + 4.124019474420322e-13 * tt3 * 0.20 + -1.012509746630069e+03 * tt4;
    } else {
      h0_RT[0] = 2.991423098020802e+00 + 7.000645004017671e-04 * tt0 * 0.50 + -5.633838501075994e-08 * tt1 * OneThird + -9.231549368983405e-12 * tt2 * 0.25 + 1.582749414838935e-15 * tt3 * 0.20 + -8.350336100339579e+02 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[1] = 2.500000000000000e+00 + 0.000000000000000e+00 * tt0 * 0.50 + 0.000000000000000e+00 * tt1 * OneThird + 0.000000000000000e+00 * tt2 * 0.25 + 0.000000000000000e+00 * tt3 * 0.20 + 2.547162000000000e+04 * tt4;
    } else {
      h0_RT[1] = 2.500000000000000e+00 + 0.000000000000000e+00 * tt0 * 0.50 + 0.000000000000000e+00 * tt1 * OneThird + 0.000000000000000e+00 * tt2 * 0.25 + 0.000000000000000e+00 * tt3 * 0.20 + 2.547162000000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[2] = 3.212715989395561e+00 + 1.129189930372770e-03 * tt0 * 0.50 + -5.802181361464514e-07 * tt1 * OneThird + 1.319054799445732e-09 * tt2 * 0.25 + -8.789124843560576e-13 * tt3 * 0.20 + -1.005227964894372e+03 * tt4;
    } else {
      h0_RT[2] = 3.697578198388313e+00 + 6.135198631579162e-04 * tt0 * 0.50 + -1.258843839277266e-07 * tt1 * OneThird + 1.775286094678163e-11 * tt2 * 0.25 + -1.136439853728577e-15 * tt3 * 0.20 + -1.233929448628332e+03 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[3] = 2.946268214996042e+00 + -1.636929514520428e-03 * tt0 * 0.50 + 2.417689047111522e-06 * tt1 * OneThird + -1.599085904220948e-09 * tt2 * 0.25 + 3.875781637558715e-13 * tt3 * 0.20 + 2.914765544214541e+04 * tt4;
    } else {
      h0_RT[3] = 2.542059346810753e+00 + -2.755084819709930e-05 * tt0 * 0.50 + -3.102750735284210e-09 * tt1 * OneThird + 4.551064604029875e-12 * tt2 * 0.25 + -4.368053603400201e-16 * tt3 * 0.20 + 2.923079932806830e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[4] = 3.637469954817199e+00 + 1.835107588557445e-04 * tt0 * 0.50 + -1.671890821934399e-06 * tt1 * OneThird + 2.382387984290979e-09 * tt2 * 0.25 + -8.412272879086942e-13 * tt3 * 0.20 + 3.606761664188268e+03 * tt4;
    } else {
      h0_RT[4] = 2.882730010358268e+00 + 1.013973808582920e-03 * tt0 * 0.50 + -2.276873525826253e-07 * tt1 * OneThird + 2.174674530723264e-11 * tt2 * 0.25 + -5.126235449665062e-16 * tt3 * 0.20 + 3.886886304206210e+03 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[5] = 2.980113629643359e+00 + 4.995528635916930e-03 * tt0 * 0.50 + -3.787831766740593e-06 * tt1 * OneThird + 2.350618553607452e-09 * tt2 * 0.25 + -8.074751784591024e-13 * tt3 * 0.20 + 1.762129276620481e+02 * tt4;
    } else {
      h0_RT[5] = 4.072190882104516e+00 + 2.131295860812287e-03 * tt0 * 0.50 + -5.308143605520907e-07 * tt1 * OneThird + 6.112265234996528e-11 * tt2 * 0.25 + -2.841160746631655e-15 * tt3 * 0.20 + -1.579732342044106e+02 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[6] = 3.388694479335974e+00 + 6.569676580383018e-03 * tt0 * 0.50 + -1.497067806981738e-07 * tt1 * OneThird + -4.624463890229065e-09 * tt2 * 0.25 + 2.470987660509708e-12 * tt3 * 0.20 + -1.766313519508993e+04 * tt4;
    } else {
      h0_RT[6] = 4.573166050599029e+00 + 4.336137813959980e-03 * tt0 * 0.50 + -1.474689772850899e-06 * tt1 * OneThird + 2.348905095864344e-10 * tt2 * 0.25 + -1.431655199308401e-14 * tt3 * 0.20 + -1.800695414321054e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[7] = 3.386728712311508e+00 + 3.475860768086567e-03 * tt0 * 0.50 + -6.357077086401517e-06 * tt1 * OneThird + 6.971269292764603e-09 * tt2 * 0.25 + -2.507661553195187e-12 * tt3 * 0.20 + -3.020809909662633e+04 * tt4;
    } else {
      h0_RT[7] = 2.672145205625812e+00 + 3.056292898549634e-03 * tt0 * 0.50 + -8.730259974096839e-07 * tt1 * OneThird + 1.200996455260830e-10 * tt2 * 0.25 + -6.391618725872017e-15 * tt3 * 0.20 + -2.989921025992035e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[8] = 3.298616281368291e+00 + 1.408708904146039e-03 * tt0 * 0.50 + -3.964481129109908e-06 * tt1 * OneThird + 5.642920880408571e-09 * tt2 * 0.25 + -2.445407041148433e-12 * tt3 * 0.20 + -1.020894198687962e+03 * tt4;
    } else {
      h0_RT[8] = 2.926639911210682e+00 + 1.487977101178227e-03 * tt0 * 0.50 + -5.684761849244810e-07 * tt1 * OneThird + 1.009704225872734e-10 * tt2 * 0.25 + -6.753354387142974e-15 * tt3 * 0.20 + -9.227966980051905e+02 * tt4;
    };

  };

  template <class Type>
    void GasThermo::getEnthalpiesDerivatives(Type& T, vector<Type>& dh0dT) {

    Type tt0 = T;
    Type tt1 = T * tt0;
    Type tt2 = T * tt1;
    Type tt3 = T * tt2;

    if(tt0 < 1.0000e+03) {
      dh0dT[0] = 3.298007330271616e+00 + 8.258468556414628e-04 * tt0 + -8.167383830175186e-07 * tt1 + -9.201733687992775e-11 * tt2 + 4.124019474420322e-13 * tt3;
    } else {
      dh0dT[0] = 2.991423098020802e+00 + 7.000645004017671e-04 * tt0 + -5.633838501075994e-08 * tt1 + -9.231549368983405e-12 * tt2 + 1.582749414838935e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[1] = 2.500000000000000e+00 + 0.000000000000000e+00 * tt0 + 0.000000000000000e+00 * tt1 + 0.000000000000000e+00 * tt2 + 0.000000000000000e+00 * tt3;
    } else {
      dh0dT[1] = 2.500000000000000e+00 + 0.000000000000000e+00 * tt0 + 0.000000000000000e+00 * tt1 + 0.000000000000000e+00 * tt2 + 0.000000000000000e+00 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[2] = 3.212715989395561e+00 + 1.129189930372770e-03 * tt0 + -5.802181361464514e-07 * tt1 + 1.319054799445732e-09 * tt2 + -8.789124843560576e-13 * tt3;
    } else {
      dh0dT[2] = 3.697578198388313e+00 + 6.135198631579162e-04 * tt0 + -1.258843839277266e-07 * tt1 + 1.775286094678163e-11 * tt2 + -1.136439853728577e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[3] = 2.946268214996042e+00 + -1.636929514520428e-03 * tt0 + 2.417689047111522e-06 * tt1 + -1.599085904220948e-09 * tt2 + 3.875781637558715e-13 * tt3;
    } else {
      dh0dT[3] = 2.542059346810753e+00 + -2.755084819709930e-05 * tt0 + -3.102750735284210e-09 * tt1 + 4.551064604029875e-12 * tt2 + -4.368053603400201e-16 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[4] = 3.637469954817199e+00 + 1.835107588557445e-04 * tt0 + -1.671890821934399e-06 * tt1 + 2.382387984290979e-09 * tt2 + -8.412272879086942e-13 * tt3;
    } else {
      dh0dT[4] = 2.882730010358268e+00 + 1.013973808582920e-03 * tt0 + -2.276873525826253e-07 * tt1 + 2.174674530723264e-11 * tt2 + -5.126235449665062e-16 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[5] = 2.980113629643359e+00 + 4.995528635916930e-03 * tt0 + -3.787831766740593e-06 * tt1 + 2.350618553607452e-09 * tt2 + -8.074751784591024e-13 * tt3;
    } else {
      dh0dT[5] = 4.072190882104516e+00 + 2.131295860812287e-03 * tt0 + -5.308143605520907e-07 * tt1 + 6.112265234996528e-11 * tt2 + -2.841160746631655e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[6] = 3.388694479335974e+00 + 6.569676580383018e-03 * tt0 + -1.497067806981738e-07 * tt1 + -4.624463890229065e-09 * tt2 + 2.470987660509708e-12 * tt3;
    } else {
      dh0dT[6] = 4.573166050599029e+00 + 4.336137813959980e-03 * tt0 + -1.474689772850899e-06 * tt1 + 2.348905095864344e-10 * tt2 + -1.431655199308401e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[7] = 3.386728712311508e+00 + 3.475860768086567e-03 * tt0 + -6.357077086401517e-06 * tt1 + 6.971269292764603e-09 * tt2 + -2.507661553195187e-12 * tt3;
    } else {
      dh0dT[7] = 2.672145205625812e+00 + 3.056292898549634e-03 * tt0 + -8.730259974096839e-07 * tt1 + 1.200996455260830e-10 * tt2 + -6.391618725872017e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[8] = 3.298616281368291e+00 + 1.408708904146039e-03 * tt0 + -3.964481129109908e-06 * tt1 + 5.642920880408571e-09 * tt2 + -2.445407041148433e-12 * tt3;
    } else {
      dh0dT[8] = 2.926639911210682e+00 + 1.487977101178227e-03 * tt0 + -5.684761849244810e-07 * tt1 + 1.009704225872734e-10 * tt2 + -6.753354387142974e-15 * tt3;
    };

  };

  template <class Type>
    void GasThermo::getEntropies_R(Type& T, vector<Type>& s0_R) {

    Type tt0 = T;
    Type tt1 = T * tt0;
    Type tt2 = T * tt1;
    Type tt3 = T * tt2;
    Type tt4 = 1.0 / T;
    Type tt5 = log(T);

    if(tt0 < 1.0000e+03) {
      s0_R[0] = 3.298007330271616e+00 * tt5 + 8.258468556414628e-04 * tt0 + -8.167383830175186e-07 * tt1 * 0.50 + -9.201733687992775e-11 * tt2 * OneThird + 4.124019474420322e-13 * tt3 * 0.25 + -3.293611841861683e+00;
    } else {
      s0_R[0] = 2.991423098020802e+00 * tt5 +  7.000645004017671e-04 * tt0 + -5.633838501075994e-08 * tt1 * 0.50 + -9.231549368983405e-12 * tt2 * OneThird + 1.582749414838935e-15 * tt3 * 0.25 + -1.355111099838972e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[1] = 2.500000000000000e+00 * tt5 + 0.000000000000000e+00 * tt0 + 0.000000000000000e+00 * tt1 * 0.50 + 0.000000000000000e+00 * tt2 * OneThird + 0.000000000000000e+00 * tt3 * 0.25 + -4.601176000000000e-01;
    } else {
      s0_R[1] = 2.500000000000000e+00 * tt5 +  0.000000000000000e+00 * tt0 + 0.000000000000000e+00 * tt1 * 0.50 + 0.000000000000000e+00 * tt2 * OneThird + 0.000000000000000e+00 * tt3 * 0.25 + -4.601176000000000e-01;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[2] = 3.212715989395561e+00 * tt5 + 1.129189930372770e-03 * tt0 + -5.802181361464514e-07 * tt1 * 0.50 + 1.319054799445732e-09 * tt2 * OneThird + -8.789124843560576e-13 * tt3 * 0.25 + 6.035646054074739e+00;
    } else {
      s0_R[2] = 3.697578198388313e+00 * tt5 +  6.135198631579162e-04 * tt0 + -1.258843839277266e-07 * tt1 * 0.50 + 1.775286094678163e-11 * tt2 * OneThird + -1.136439853728577e-15 * tt3 * 0.25 + 3.189163063139013e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[3] = 2.946268214996042e+00 * tt5 + -1.636929514520428e-03 * tt0 + 2.417689047111522e-06 * tt1 * 0.50 + -1.599085904220948e-09 * tt2 * OneThird + 3.875781637558715e-13 * tt3 * 0.25 + 2.964654488320935e+00;
    } else {
      s0_R[3] = 2.542059346810753e+00 * tt5 +  -2.755084819709930e-05 * tt0 + -3.102750735284210e-09 * tt1 * 0.50 + 4.551064604029875e-12 * tt2 * OneThird + -4.368053603400201e-16 * tt3 * 0.25 + 4.920305749943392e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[4] = 3.637469954817199e+00 * tt5 + 1.835107588557445e-04 * tt0 + -1.671890821934399e-06 * tt1 * 0.50 + 2.382387984290979e-09 * tt2 * OneThird + -8.412272879086942e-13 * tt3 * 0.25 + 1.358017185765906e+00;
    } else {
      s0_R[4] = 2.882730010358268e+00 * tt5 +  1.013973808582920e-03 * tt0 + -2.276873525826253e-07 * tt1 * 0.50 + 2.174674530723264e-11 * tt2 * OneThird + -5.126235449665062e-16 * tt3 * 0.25 + 5.595712983861354e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[5] = 2.980113629643359e+00 * tt5 + 4.995528635916930e-03 * tt0 + -3.787831766740593e-06 * tt1 * 0.50 + 2.350618553607452e-09 * tt2 * OneThird + -8.074751784591024e-13 * tt3 * 0.25 + 9.222101780643662e+00;
    } else {
      s0_R[5] = 4.072190882104516e+00 * tt5 +  2.131295860812287e-03 * tt0 + -5.308143605520907e-07 * tt1 * 0.50 + 6.112265234996528e-11 * tt2 * OneThird + -2.841160746631655e-15 * tt3 * 0.25 + 3.476030242900233e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[6] = 3.388694479335974e+00 * tt5 + 6.569676580383018e-03 * tt0 + -1.497067806981738e-07 * tt1 * 0.50 + -4.624463890229065e-09 * tt2 * OneThird + 2.470987660509708e-12 * tt3 * 0.25 + 6.785608772166807e+00;
    } else {
      s0_R[6] = 4.573166050599029e+00 * tt5 +  4.336137813959980e-03 * tt0 + -1.474689772850899e-06 * tt1 * 0.50 + 2.348905095864344e-10 * tt2 * OneThird + -1.431655199308401e-14 * tt3 * 0.25 + 5.011405386567667e-01;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[7] = 3.386728712311508e+00 * tt5 + 3.475860768086567e-03 * tt0 + -6.357077086401517e-06 * tt1 * 0.50 + 6.971269292764603e-09 * tt2 * OneThird + -2.507661553195187e-12 * tt3 * 0.25 + 2.590699531648500e+00;
    } else {
      s0_R[7] = 2.672145205625812e+00 * tt5 +  3.056292898549634e-03 * tt0 + -8.730259974096839e-07 * tt1 * 0.50 + 1.200996455260830e-10 * tt2 * OneThird + -6.391618725872017e-15 * tt3 * 0.25 + 6.862815579400046e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[8] = 3.298616281368291e+00 * tt5 + 1.408708904146039e-03 * tt0 + -3.964481129109908e-06 * tt1 * 0.50 + 5.642920880408571e-09 * tt2 * OneThird + -2.445407041148433e-12 * tt3 * 0.25 + 3.950623591964732e+00;
    } else {
      s0_R[8] = 2.926639911210682e+00 * tt5 +  1.487977101178227e-03 * tt0 + -5.684761849244810e-07 * tt1 * 0.50 + 1.009704225872734e-10 * tt2 * OneThird + -6.753354387142974e-15 * tt3 * 0.25 + 5.980528055036107e+00;
    };

  };

  template <class Type>
    void GasThermo::getGibbsFunctions_RT(Type& T,vector<Type>& g0_RT) {

    vector<Type> h0_RT(m_kk, 0.0);
    vector<Type> s0_R(m_kk,  0.0);

    getEnthalpies_RT(T, h0_RT);
    getEntropies_R(T, s0_R);
    for(int k = 0; k < m_kk; ++k) { g0_RT[k] = h0_RT[k] - s0_R[k]; }

  };

  template <class Type>
    void GasThermo::getEquilibriumConstants(Type& T, vector<Type>& keqs) {

    double       p0 = OneAtm;
    Type         RT = GasConstant * T;
    Type         C0 = p0 / RT;
    vector<Type> g0_RT(m_kk, 0.0);

    getGibbsFunctions_RT(T, g0_RT);
    for(int k = 0; k < m_kk; ++k) { g0_RT[k] = exp(g0_RT[k]); }

    keqs[0] = (g0_RT[3] * g0_RT[4]);
    keqs[1] = (g0_RT[1] * g0_RT[4]);
    keqs[2] = (g0_RT[1] * g0_RT[7]);
    keqs[3] = (g0_RT[4] * g0_RT[4]);
    keqs[4] = (C0 * g0_RT[0]);
    keqs[5] = (C0 * g0_RT[7]);
    keqs[6] = (C0 * g0_RT[2]);
    keqs[7] = (C0 * g0_RT[4]);
    keqs[8] = (C0 * g0_RT[5]);
    keqs[9] = (C0 * g0_RT[5]);
    keqs[10] = (g0_RT[4] * g0_RT[4]);
    keqs[11] = (g0_RT[0] * g0_RT[2]);
    keqs[12] = (g0_RT[7] * g0_RT[3]);
    keqs[13] = (g0_RT[2] * g0_RT[4]);
    keqs[14] = (g0_RT[7] * g0_RT[2]);
    keqs[15] = (g0_RT[7] * g0_RT[2]);
    keqs[16] = (C0 * g0_RT[6]);
    keqs[17] = (g0_RT[6] * g0_RT[2]);
    keqs[18] = (g0_RT[6] * g0_RT[2]);
    keqs[19] = (g0_RT[0] * g0_RT[5]);
    keqs[20] = (g0_RT[7] * g0_RT[4]);
    keqs[21] = (g0_RT[7] * g0_RT[5]);
    keqs[22] = (g0_RT[7] * g0_RT[5]);
    keqs[23] = (g0_RT[5] * g0_RT[4]);

    keqs[0] /= (g0_RT[1] * g0_RT[2]);
    keqs[1] /= (g0_RT[0] * g0_RT[3]);
    keqs[2] /= (g0_RT[0] * g0_RT[4]);
    keqs[3] /= (g0_RT[7] * g0_RT[3]);
    keqs[4] /= (g0_RT[1] * g0_RT[1]);
    keqs[5] /= (g0_RT[1] * g0_RT[4]);
    keqs[6] /= (g0_RT[3] * g0_RT[3]);
    keqs[7] /= (g0_RT[1] * g0_RT[3]);
    keqs[8] /= (g0_RT[3] * g0_RT[4]);
    keqs[9] /= (g0_RT[1] * g0_RT[2]);
    keqs[10] /= (g0_RT[1] * g0_RT[5]);
    keqs[11] /= (g0_RT[1] * g0_RT[5]);
    keqs[12] /= (g0_RT[1] * g0_RT[5]);
    keqs[13] /= (g0_RT[5] * g0_RT[3]);
    keqs[14] /= (g0_RT[5] * g0_RT[4]);
    keqs[15] /= (g0_RT[5] * g0_RT[4]);
    keqs[16] /= (g0_RT[4] * g0_RT[4]);
    keqs[17] /= (g0_RT[5] * g0_RT[5]);
    keqs[18] /= (g0_RT[5] * g0_RT[5]);
    keqs[19] /= (g0_RT[1] * g0_RT[6]);
    keqs[20] /= (g0_RT[1] * g0_RT[6]);
    keqs[21] /= (g0_RT[6] * g0_RT[4]);
    keqs[22] /= (g0_RT[6] * g0_RT[4]);
    keqs[23] /= (g0_RT[6] * g0_RT[3]);

  };

}

#endif
