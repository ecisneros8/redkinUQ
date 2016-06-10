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
          m_speciesNames.resize(m_kk);
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
	  m_speciesNames[2] = "O";
	  m_speciesNames[3] = "O2";
	  m_speciesNames[4] = "OH";
	  m_speciesNames[5] = "H2O";
	  m_speciesNames[6] = "HO2";
	  m_speciesNames[7] = "H2O2";
	  m_speciesNames[8] = "C";
	  m_speciesNames[9] = "CH";
	  m_speciesNames[10] = "CH2";
	  m_speciesNames[11] = "CH2(S)";
	  m_speciesNames[12] = "CH3";
	  m_speciesNames[13] = "CH4";
	  m_speciesNames[14] = "CO";
	  m_speciesNames[15] = "CO2";
	  m_speciesNames[16] = "HCO";
	  m_speciesNames[17] = "CH2O";
	  m_speciesNames[18] = "CH2OH";
	  m_speciesNames[19] = "CH3O";
	  m_speciesNames[20] = "CH3OH";
	  m_speciesNames[21] = "C2H";
	  m_speciesNames[22] = "C2H2";
	  m_speciesNames[23] = "C2H3";
	  m_speciesNames[24] = "C2H4";
	  m_speciesNames[25] = "C2H5";
	  m_speciesNames[26] = "C2H6";
	  m_speciesNames[27] = "HCCO";
	  m_speciesNames[28] = "CH2CO";
	  m_speciesNames[29] = "HCCOH";
	  m_speciesNames[30] = "N2";
	  m_speciesNames[31] = "AR";

	  // species indices
	  m_speciesIndex["H2"] = 0;
	  m_speciesIndex["H"] = 1;
	  m_speciesIndex["O"] = 2;
	  m_speciesIndex["O2"] = 3;
	  m_speciesIndex["OH"] = 4;
	  m_speciesIndex["H2O"] = 5;
	  m_speciesIndex["HO2"] = 6;
	  m_speciesIndex["H2O2"] = 7;
	  m_speciesIndex["C"] = 8;
	  m_speciesIndex["CH"] = 9;
	  m_speciesIndex["CH2"] = 10;
	  m_speciesIndex["CH2(S)"] = 11;
	  m_speciesIndex["CH3"] = 12;
	  m_speciesIndex["CH4"] = 13;
	  m_speciesIndex["CO"] = 14;
	  m_speciesIndex["CO2"] = 15;
	  m_speciesIndex["HCO"] = 16;
	  m_speciesIndex["CH2O"] = 17;
	  m_speciesIndex["CH2OH"] = 18;
	  m_speciesIndex["CH3O"] = 19;
	  m_speciesIndex["CH3OH"] = 20;
	  m_speciesIndex["C2H"] = 21;
	  m_speciesIndex["C2H2"] = 22;
	  m_speciesIndex["C2H3"] = 23;
	  m_speciesIndex["C2H4"] = 24;
	  m_speciesIndex["C2H5"] = 25;
	  m_speciesIndex["C2H6"] = 26;
	  m_speciesIndex["HCCO"] = 27;
	  m_speciesIndex["CH2CO"] = 28;
	  m_speciesIndex["HCCOH"] = 29;
	  m_speciesIndex["N2"] = 30;
	  m_speciesIndex["AR"] = 31;

	  // set molecular weights
	  m_mw[0] = 2.015880e+00;
	  m_mw[1] = 1.007940e+00;
	  m_mw[2] = 1.599940e+01;
	  m_mw[3] = 3.199880e+01;
	  m_mw[4] = 1.700734e+01;
	  m_mw[5] = 1.801528e+01;
	  m_mw[6] = 3.300674e+01;
	  m_mw[7] = 3.401468e+01;
	  m_mw[8] = 1.201100e+01;
	  m_mw[9] = 1.301894e+01;
	  m_mw[10] = 1.402688e+01;
	  m_mw[11] = 1.402688e+01;
	  m_mw[12] = 1.503482e+01;
	  m_mw[13] = 1.604276e+01;
	  m_mw[14] = 2.801040e+01;
	  m_mw[15] = 4.400980e+01;
	  m_mw[16] = 2.901834e+01;
	  m_mw[17] = 3.002628e+01;
	  m_mw[18] = 3.103422e+01;
	  m_mw[19] = 3.103422e+01;
	  m_mw[20] = 3.204216e+01;
	  m_mw[21] = 2.502994e+01;
	  m_mw[22] = 2.603788e+01;
	  m_mw[23] = 2.704582e+01;
	  m_mw[24] = 2.805376e+01;
	  m_mw[25] = 2.906170e+01;
	  m_mw[26] = 3.006964e+01;
	  m_mw[27] = 4.102934e+01;
	  m_mw[28] = 4.203728e+01;
	  m_mw[29] = 4.203728e+01;
	  m_mw[30] = 2.801348e+01;
	  m_mw[31] = 3.994800e+01;

	  // set element matrix
	  m_Emat(0,0) = 0.0e+00;
	  m_Emat(0,1) = 2.0e+00;
	  m_Emat(0,2) = 0.0e+00;
	  m_Emat(0,3) = 0.0e+00;
	  m_Emat(0,4) = 0.0e+00;

	  m_Emat(1,0) = 0.0e+00;
	  m_Emat(1,1) = 1.0e+00;
	  m_Emat(1,2) = 0.0e+00;
	  m_Emat(1,3) = 0.0e+00;
	  m_Emat(1,4) = 0.0e+00;

	  m_Emat(2,0) = 1.0e+00;
	  m_Emat(2,1) = 0.0e+00;
	  m_Emat(2,2) = 0.0e+00;
	  m_Emat(2,3) = 0.0e+00;
	  m_Emat(2,4) = 0.0e+00;

	  m_Emat(3,0) = 2.0e+00;
	  m_Emat(3,1) = 0.0e+00;
	  m_Emat(3,2) = 0.0e+00;
	  m_Emat(3,3) = 0.0e+00;
	  m_Emat(3,4) = 0.0e+00;

	  m_Emat(4,0) = 1.0e+00;
	  m_Emat(4,1) = 1.0e+00;
	  m_Emat(4,2) = 0.0e+00;
	  m_Emat(4,3) = 0.0e+00;
	  m_Emat(4,4) = 0.0e+00;

	  m_Emat(5,0) = 1.0e+00;
	  m_Emat(5,1) = 2.0e+00;
	  m_Emat(5,2) = 0.0e+00;
	  m_Emat(5,3) = 0.0e+00;
	  m_Emat(5,4) = 0.0e+00;

	  m_Emat(6,0) = 2.0e+00;
	  m_Emat(6,1) = 1.0e+00;
	  m_Emat(6,2) = 0.0e+00;
	  m_Emat(6,3) = 0.0e+00;
	  m_Emat(6,4) = 0.0e+00;

	  m_Emat(7,0) = 2.0e+00;
	  m_Emat(7,1) = 2.0e+00;
	  m_Emat(7,2) = 0.0e+00;
	  m_Emat(7,3) = 0.0e+00;
	  m_Emat(7,4) = 0.0e+00;

	  m_Emat(8,0) = 0.0e+00;
	  m_Emat(8,1) = 0.0e+00;
	  m_Emat(8,2) = 1.0e+00;
	  m_Emat(8,3) = 0.0e+00;
	  m_Emat(8,4) = 0.0e+00;

	  m_Emat(9,0) = 0.0e+00;
	  m_Emat(9,1) = 1.0e+00;
	  m_Emat(9,2) = 1.0e+00;
	  m_Emat(9,3) = 0.0e+00;
	  m_Emat(9,4) = 0.0e+00;

	  m_Emat(10,0) = 0.0e+00;
	  m_Emat(10,1) = 2.0e+00;
	  m_Emat(10,2) = 1.0e+00;
	  m_Emat(10,3) = 0.0e+00;
	  m_Emat(10,4) = 0.0e+00;

	  m_Emat(11,0) = 0.0e+00;
	  m_Emat(11,1) = 2.0e+00;
	  m_Emat(11,2) = 1.0e+00;
	  m_Emat(11,3) = 0.0e+00;
	  m_Emat(11,4) = 0.0e+00;

	  m_Emat(12,0) = 0.0e+00;
	  m_Emat(12,1) = 3.0e+00;
	  m_Emat(12,2) = 1.0e+00;
	  m_Emat(12,3) = 0.0e+00;
	  m_Emat(12,4) = 0.0e+00;

	  m_Emat(13,0) = 0.0e+00;
	  m_Emat(13,1) = 4.0e+00;
	  m_Emat(13,2) = 1.0e+00;
	  m_Emat(13,3) = 0.0e+00;
	  m_Emat(13,4) = 0.0e+00;

	  m_Emat(14,0) = 1.0e+00;
	  m_Emat(14,1) = 0.0e+00;
	  m_Emat(14,2) = 1.0e+00;
	  m_Emat(14,3) = 0.0e+00;
	  m_Emat(14,4) = 0.0e+00;

	  m_Emat(15,0) = 2.0e+00;
	  m_Emat(15,1) = 0.0e+00;
	  m_Emat(15,2) = 1.0e+00;
	  m_Emat(15,3) = 0.0e+00;
	  m_Emat(15,4) = 0.0e+00;

	  m_Emat(16,0) = 1.0e+00;
	  m_Emat(16,1) = 1.0e+00;
	  m_Emat(16,2) = 1.0e+00;
	  m_Emat(16,3) = 0.0e+00;
	  m_Emat(16,4) = 0.0e+00;

	  m_Emat(17,0) = 1.0e+00;
	  m_Emat(17,1) = 2.0e+00;
	  m_Emat(17,2) = 1.0e+00;
	  m_Emat(17,3) = 0.0e+00;
	  m_Emat(17,4) = 0.0e+00;

	  m_Emat(18,0) = 1.0e+00;
	  m_Emat(18,1) = 3.0e+00;
	  m_Emat(18,2) = 1.0e+00;
	  m_Emat(18,3) = 0.0e+00;
	  m_Emat(18,4) = 0.0e+00;

	  m_Emat(19,0) = 1.0e+00;
	  m_Emat(19,1) = 3.0e+00;
	  m_Emat(19,2) = 1.0e+00;
	  m_Emat(19,3) = 0.0e+00;
	  m_Emat(19,4) = 0.0e+00;

	  m_Emat(20,0) = 1.0e+00;
	  m_Emat(20,1) = 4.0e+00;
	  m_Emat(20,2) = 1.0e+00;
	  m_Emat(20,3) = 0.0e+00;
	  m_Emat(20,4) = 0.0e+00;

	  m_Emat(21,0) = 0.0e+00;
	  m_Emat(21,1) = 1.0e+00;
	  m_Emat(21,2) = 2.0e+00;
	  m_Emat(21,3) = 0.0e+00;
	  m_Emat(21,4) = 0.0e+00;

	  m_Emat(22,0) = 0.0e+00;
	  m_Emat(22,1) = 2.0e+00;
	  m_Emat(22,2) = 2.0e+00;
	  m_Emat(22,3) = 0.0e+00;
	  m_Emat(22,4) = 0.0e+00;

	  m_Emat(23,0) = 0.0e+00;
	  m_Emat(23,1) = 3.0e+00;
	  m_Emat(23,2) = 2.0e+00;
	  m_Emat(23,3) = 0.0e+00;
	  m_Emat(23,4) = 0.0e+00;

	  m_Emat(24,0) = 0.0e+00;
	  m_Emat(24,1) = 4.0e+00;
	  m_Emat(24,2) = 2.0e+00;
	  m_Emat(24,3) = 0.0e+00;
	  m_Emat(24,4) = 0.0e+00;

	  m_Emat(25,0) = 0.0e+00;
	  m_Emat(25,1) = 5.0e+00;
	  m_Emat(25,2) = 2.0e+00;
	  m_Emat(25,3) = 0.0e+00;
	  m_Emat(25,4) = 0.0e+00;

	  m_Emat(26,0) = 0.0e+00;
	  m_Emat(26,1) = 6.0e+00;
	  m_Emat(26,2) = 2.0e+00;
	  m_Emat(26,3) = 0.0e+00;
	  m_Emat(26,4) = 0.0e+00;

	  m_Emat(27,0) = 1.0e+00;
	  m_Emat(27,1) = 1.0e+00;
	  m_Emat(27,2) = 2.0e+00;
	  m_Emat(27,3) = 0.0e+00;
	  m_Emat(27,4) = 0.0e+00;

	  m_Emat(28,0) = 1.0e+00;
	  m_Emat(28,1) = 2.0e+00;
	  m_Emat(28,2) = 2.0e+00;
	  m_Emat(28,3) = 0.0e+00;
	  m_Emat(28,4) = 0.0e+00;

	  m_Emat(29,0) = 1.0e+00;
	  m_Emat(29,1) = 2.0e+00;
	  m_Emat(29,2) = 2.0e+00;
	  m_Emat(29,3) = 0.0e+00;
	  m_Emat(29,4) = 0.0e+00;

	  m_Emat(30,0) = 0.0e+00;
	  m_Emat(30,1) = 0.0e+00;
	  m_Emat(30,2) = 0.0e+00;
	  m_Emat(30,3) = 2.0e+00;
	  m_Emat(30,4) = 0.0e+00;

	  m_Emat(31,0) = 0.0e+00;
	  m_Emat(31,1) = 0.0e+00;
	  m_Emat(31,2) = 0.0e+00;
	  m_Emat(31,3) = 0.0e+00;
	  m_Emat(31,4) = 1.0e+00;
	  
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
    for (int i = offset; i <= m_kk-3; ++i) {  
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

    T = Tp;

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
      cp0_R[0] = 2.344331120000000e+00 + 7.980520749999999e-03 * tt0 + -1.947815100000000e-05 * tt1 + 2.015720940000000e-08 * tt2 + -7.376117610000001e-12 * tt3;
    } else {
      cp0_R[0] = 3.337279200000000e+00 + -4.940247310000000e-05 * tt0 + 4.994567780000000e-07 * tt1 + -1.795663940000000e-10 * tt2 + 2.002553760000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[1] = 2.500000000000000e+00 + 7.053328190000000e-13 * tt0 + -1.995919640000000e-15 * tt1 + 2.300816320000000e-18 * tt2 + -9.277323320000001e-22 * tt3;
    } else {
      cp0_R[1] = 2.500000010000000e+00 + -2.308429730000000e-11 * tt0 + 1.615619480000000e-14 * tt1 + -4.735152350000000e-18 * tt2 + 4.981973570000000e-22 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[2] = 3.168267100000000e+00 + -3.279318840000000e-03 * tt0 + 6.643063960000000e-06 * tt1 + -6.128066240000000e-09 * tt2 + 2.112659710000000e-12 * tt3;
    } else {
      cp0_R[2] = 2.569420780000000e+00 + -8.597411370000000e-05 * tt0 + 4.194845890000000e-08 * tt1 + -1.001777990000000e-11 * tt2 + 1.228336910000000e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[3] = 3.782456360000000e+00 + -2.996734160000000e-03 * tt0 + 9.847302010000000e-06 * tt1 + -9.681295090000001e-09 * tt2 + 3.243728370000000e-12 * tt3;
    } else {
      cp0_R[3] = 3.282537840000000e+00 + 1.483087540000000e-03 * tt0 + -7.579666690000000e-07 * tt1 + 2.094705550000000e-10 * tt2 + -2.167177940000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[4] = 3.992015430000000e+00 + -2.401317520000000e-03 * tt0 + 4.617938410000000e-06 * tt1 + -3.881133330000000e-09 * tt2 + 1.364114700000000e-12 * tt3;
    } else {
      cp0_R[4] = 3.092887670000000e+00 + 5.484297160000000e-04 * tt0 + 1.265052280000000e-07 * tt1 + -8.794615559999999e-11 * tt2 + 1.174123760000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[5] = 4.198640560000000e+00 + -2.036434100000000e-03 * tt0 + 6.520402110000000e-06 * tt1 + -5.487970620000000e-09 * tt2 + 1.771978170000000e-12 * tt3;
    } else {
      cp0_R[5] = 3.033992490000000e+00 + 2.176918040000000e-03 * tt0 + -1.640725180000000e-07 * tt1 + -9.704198700000000e-11 * tt2 + 1.682009920000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[6] = 4.301798010000000e+00 + -4.749120510000000e-03 * tt0 + 2.115828910000000e-05 * tt1 + -2.427638940000000e-08 * tt2 + 9.292251240000000e-12 * tt3;
    } else {
      cp0_R[6] = 4.017210900000000e+00 + 2.239820130000000e-03 * tt0 + -6.336581500000000e-07 * tt1 + 1.142463700000000e-10 * tt2 + -1.079085350000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[7] = 4.276112690000000e+00 + -5.428224169999999e-04 * tt0 + 1.673357010000000e-05 * tt1 + -2.157708130000000e-08 * tt2 + 8.624543630000000e-12 * tt3;
    } else {
      cp0_R[7] = 4.165002850000000e+00 + 4.908316940000000e-03 * tt0 + -1.901392250000000e-06 * tt1 + 3.711859860000000e-10 * tt2 + -2.879083050000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[8] = 2.554239550000000e+00 + -3.215377240000000e-04 * tt0 + 7.337922450000000e-07 * tt1 + -7.322348890000000e-10 * tt2 + 2.665214460000000e-13 * tt3;
    } else {
      cp0_R[8] = 2.492668880000000e+00 + 4.798892840000000e-05 * tt0 + -7.243350200000001e-08 * tt1 + 3.742910290000000e-11 * tt2 + -4.872778930000000e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[9] = 3.489816650000000e+00 + 3.238355410000000e-04 * tt0 + -1.688990650000000e-06 * tt1 + 3.162173270000000e-09 * tt2 + -1.406090670000000e-12 * tt3;
    } else {
      cp0_R[9] = 2.878464730000000e+00 + 9.709136810000000e-04 * tt0 + 1.444456550000000e-07 * tt1 + -1.306878490000000e-10 * tt2 + 1.760793830000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[10] = 3.762678670000000e+00 + 9.688721430000000e-04 * tt0 + 2.794898410000000e-06 * tt1 + -3.850911530000000e-09 * tt2 + 1.687417190000000e-12 * tt3;
    } else {
      cp0_R[10] = 2.874101130000000e+00 + 3.656392920000000e-03 * tt0 + -1.408945970000000e-06 * tt1 + 2.601795490000000e-10 * tt2 + -1.877275670000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[11] = 4.198604110000000e+00 + -2.366614190000000e-03 * tt0 + 8.232962200000000e-06 * tt1 + -6.688159810000000e-09 * tt2 + 1.943147370000000e-12 * tt3;
    } else {
      cp0_R[11] = 2.292038420000000e+00 + 4.655886370000000e-03 * tt0 + -2.011919470000000e-06 * tt1 + 4.179060000000000e-10 * tt2 + -3.397163650000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[12] = 3.673590400000000e+00 + 2.010951750000000e-03 * tt0 + 5.730218560000000e-06 * tt1 + -6.871174250000000e-09 * tt2 + 2.543857340000000e-12 * tt3;
    } else {
      cp0_R[12] = 2.285717720000000e+00 + 7.239900370000000e-03 * tt0 + -2.987143480000000e-06 * tt1 + 5.956846440000000e-10 * tt2 + -4.671543940000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[13] = 5.149876130000000e+00 + -1.367097880000000e-02 * tt0 + 4.918005990000000e-05 * tt1 + -4.847430260000000e-08 * tt2 + 1.666939560000000e-11 * tt3;
    } else {
      cp0_R[13] = 7.485149500000000e-02 + 1.339094670000000e-02 * tt0 + -5.732858090000000e-06 * tt1 + 1.222925350000000e-09 * tt2 + -1.018152300000000e-13 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[14] = 3.579533470000000e+00 + -6.103536800000000e-04 * tt0 + 1.016814330000000e-06 * tt1 + 9.070058840000000e-10 * tt2 + -9.044244990000000e-13 * tt3;
    } else {
      cp0_R[14] = 2.715185610000000e+00 + 2.062527430000000e-03 * tt0 + -9.988257710000001e-07 * tt1 + 2.300530080000000e-10 * tt2 + -2.036477160000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[15] = 2.356773520000000e+00 + 8.984596770000000e-03 * tt0 + -7.123562690000000e-06 * tt1 + 2.459190220000000e-09 * tt2 + -1.436995480000000e-13 * tt3;
    } else {
      cp0_R[15] = 3.857460290000000e+00 + 4.414370260000000e-03 * tt0 + -2.214814040000000e-06 * tt1 + 5.234901880000001e-10 * tt2 + -4.720841640000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[16] = 4.221185840000000e+00 + -3.243925320000000e-03 * tt0 + 1.377994460000000e-05 * tt1 + -1.331440930000000e-08 * tt2 + 4.337688650000000e-12 * tt3;
    } else {
      cp0_R[16] = 2.772174380000000e+00 + 4.956955260000000e-03 * tt0 + -2.484456130000000e-06 * tt1 + 5.891617780000000e-10 * tt2 + -5.335087110000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[17] = 4.793723150000000e+00 + -9.908333690000000e-03 * tt0 + 3.732200080000000e-05 * tt1 + -3.792852610000000e-08 * tt2 + 1.317726520000000e-11 * tt3;
    } else {
      cp0_R[17] = 1.760690080000000e+00 + 9.200000820000000e-03 * tt0 + -4.422588130000000e-06 * tt1 + 1.006412120000000e-09 * tt2 + -8.838556400000001e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[18] = 3.863889180000000e+00 + 5.596723040000000e-03 * tt0 + 5.932717910000000e-06 * tt1 + -1.045320120000000e-08 * tt2 + 4.369672780000000e-12 * tt3;
    } else {
      cp0_R[18] = 3.692665690000000e+00 + 8.645767970000001e-03 * tt0 + -3.751011200000000e-06 * tt1 + 7.872346360000000e-10 * tt2 + -6.485542010000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[19] = 2.106083740868644e+00 + 7.217545611840916e-03 * tt0 + 5.335825777358757e-06 * tt1 + -7.374562574086794e-09 * tt2 + 2.074342158630007e-12 * tt3;
    } else {
      cp0_R[19] = 3.771000099870984e+00 + 7.871073375136668e-03 * tt0 + -2.656064859783170e-06 * tt1 + 3.943402155606197e-10 * tt2 + -2.111411617357040e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[20] = 5.715395820000000e+00 + -1.523091290000000e-02 * tt0 + 6.524411550000000e-05 * tt1 + -7.108068890000000e-08 * tt2 + 2.613526980000000e-11 * tt3;
    } else {
      cp0_R[20] = 1.789707910000000e+00 + 1.409382920000000e-02 * tt0 + -6.365008350000000e-06 * tt1 + 1.381710850000000e-09 * tt2 + -1.170602200000000e-13 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[21] = 2.889657330000000e+00 + 1.340996110000000e-02 * tt0 + -2.847695010000000e-05 * tt1 + 2.947910450000000e-08 * tt2 + -1.093315110000000e-11 * tt3;
    } else {
      cp0_R[21] = 3.167806520000000e+00 + 4.752219020000000e-03 * tt0 + -1.837870770000000e-06 * tt1 + 3.041902520000000e-10 * tt2 + -1.772327700000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[22] = 8.086810940000000e-01 + 2.336156290000000e-02 * tt0 + -3.551718150000000e-05 * tt1 + 2.801524370000000e-08 * tt2 + -8.500729740000000e-12 * tt3;
    } else {
      cp0_R[22] = 4.147569640000000e+00 + 5.961666640000000e-03 * tt0 + -2.372948520000000e-06 * tt1 + 4.674121710000000e-10 * tt2 + -3.612352130000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[23] = 3.212466450000000e+00 + 1.514791620000000e-03 * tt0 + 2.592094120000000e-05 * tt1 + -3.576578470000000e-08 * tt2 + 1.471508730000000e-11 * tt3;
    } else {
      cp0_R[23] = 3.016724000000000e+00 + 1.033022920000000e-02 * tt0 + -4.680823490000000e-06 * tt1 + 1.017632880000000e-09 * tt2 + -8.626070410000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[24] = 3.959201480000000e+00 + -7.570522470000000e-03 * tt0 + 5.709902920000000e-05 * tt1 + -6.915887530000000e-08 * tt2 + 2.698843730000000e-11 * tt3;
    } else {
      cp0_R[24] = 2.036111160000000e+00 + 1.464541510000000e-02 * tt0 + -6.710779150000000e-06 * tt1 + 1.472229230000000e-09 * tt2 + -1.257060610000000e-13 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[25] = 4.306465680000000e+00 + -4.186588920000000e-03 * tt0 + 4.971428070000000e-05 * tt1 + -5.991266060000000e-08 * tt2 + 2.305090040000000e-11 * tt3;
    } else {
      cp0_R[25] = 1.954656420000000e+00 + 1.739727220000000e-02 * tt0 + -7.982066680000000e-06 * tt1 + 1.752176890000000e-09 * tt2 + -1.496415760000000e-13 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[26] = 4.291424920000000e+00 + -5.501542700000000e-03 * tt0 + 5.994382880000000e-05 * tt1 + -7.084662850000000e-08 * tt2 + 2.686857710000000e-11 * tt3;
    } else {
      cp0_R[26] = 1.071881500000000e+00 + 2.168526770000000e-02 * tt0 + -1.002560670000000e-05 * tt1 + 2.214120010000000e-09 * tt2 + -1.900028900000000e-13 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[27] = 2.251845286753941e+00 + 1.765401529121510e-02 * tt0 + -2.372620998746318e-05 * tt1 + 1.727224682516638e-08 * tt2 + -5.064957679024219e-12 * tt3;
    } else {
      cp0_R[27] = 5.628172209693201e+00 + 4.085394191788521e-03 * tt0 + -1.593486747480799e-06 * tt1 + 2.862686614342911e-10 * tt2 + -1.940857878719715e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[28] = 2.135836300000000e+00 + 1.811887210000000e-02 * tt0 + -1.739474740000000e-05 * tt1 + 9.343975680000000e-09 * tt2 + -2.014576150000000e-12 * tt3;
    } else {
      cp0_R[28] = 4.511297320000000e+00 + 9.003597449999999e-03 * tt0 + -4.169396350000000e-06 * tt1 + 9.233458820000000e-10 * tt2 + -7.948382010000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[29] = 1.242373300000000e+00 + 3.107220100000000e-02 * tt0 + -5.086686400000000e-05 * tt1 + 4.313713100000000e-08 * tt2 + -1.401459400000000e-11 * tt3;
    } else {
      cp0_R[29] = 5.923829100000000e+00 + 6.792360000000000e-03 * tt0 + -2.565856400000000e-06 * tt1 + 4.498784100000000e-10 * tt2 + -2.994010100000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[30] = 3.298616281368291e+00 + 1.408708904146039e-03 * tt0 + -3.964481129109908e-06 * tt1 + 5.642920880408571e-09 * tt2 + -2.445407041148433e-12 * tt3;
    } else {
      cp0_R[30] = 2.926639911210682e+00 + 1.487977101178227e-03 * tt0 + -5.684761849244810e-07 * tt1 + 1.009704225872734e-10 * tt2 + -6.753354387142974e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[31] = 2.500000000000000e+00 + 0.000000000000000e+00 * tt0 + 0.000000000000000e+00 * tt1 + 0.000000000000000e+00 * tt2 + 0.000000000000000e+00 * tt3;
    } else {
      cp0_R[31] = 2.500000000000000e+00 + 0.000000000000000e+00 * tt0 + 0.000000000000000e+00 * tt1 + 0.000000000000000e+00 * tt2 + 0.000000000000000e+00 * tt3;
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
      h0_RT[0] = 2.344331120000000e+00 + 7.980520749999999e-03 * tt0 * 0.50 + -1.947815100000000e-05 * tt1 * OneThird + 2.015720940000000e-08 * tt2 * 0.25 + -7.376117610000001e-12 * tt3 * 0.20 + -9.179351730000000e+02 * tt4;
    } else {
      h0_RT[0] = 3.337279200000000e+00 + -4.940247310000000e-05 * tt0 * 0.50 + 4.994567780000000e-07 * tt1 * OneThird + -1.795663940000000e-10 * tt2 * 0.25 + 2.002553760000000e-14 * tt3 * 0.20 + -9.501589220000000e+02 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[1] = 2.500000000000000e+00 + 7.053328190000000e-13 * tt0 * 0.50 + -1.995919640000000e-15 * tt1 * OneThird + 2.300816320000000e-18 * tt2 * 0.25 + -9.277323320000001e-22 * tt3 * 0.20 + 2.547365990000000e+04 * tt4;
    } else {
      h0_RT[1] = 2.500000010000000e+00 + -2.308429730000000e-11 * tt0 * 0.50 + 1.615619480000000e-14 * tt1 * OneThird + -4.735152350000000e-18 * tt2 * 0.25 + 4.981973570000000e-22 * tt3 * 0.20 + 2.547365990000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[2] = 3.168267100000000e+00 + -3.279318840000000e-03 * tt0 * 0.50 + 6.643063960000000e-06 * tt1 * OneThird + -6.128066240000000e-09 * tt2 * 0.25 + 2.112659710000000e-12 * tt3 * 0.20 + 2.912225920000000e+04 * tt4;
    } else {
      h0_RT[2] = 2.569420780000000e+00 + -8.597411370000000e-05 * tt0 * 0.50 + 4.194845890000000e-08 * tt1 * OneThird + -1.001777990000000e-11 * tt2 * 0.25 + 1.228336910000000e-15 * tt3 * 0.20 + 2.921757910000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[3] = 3.782456360000000e+00 + -2.996734160000000e-03 * tt0 * 0.50 + 9.847302010000000e-06 * tt1 * OneThird + -9.681295090000001e-09 * tt2 * 0.25 + 3.243728370000000e-12 * tt3 * 0.20 + -1.063943560000000e+03 * tt4;
    } else {
      h0_RT[3] = 3.282537840000000e+00 + 1.483087540000000e-03 * tt0 * 0.50 + -7.579666690000000e-07 * tt1 * OneThird + 2.094705550000000e-10 * tt2 * 0.25 + -2.167177940000000e-14 * tt3 * 0.20 + -1.088457720000000e+03 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[4] = 3.992015430000000e+00 + -2.401317520000000e-03 * tt0 * 0.50 + 4.617938410000000e-06 * tt1 * OneThird + -3.881133330000000e-09 * tt2 * 0.25 + 1.364114700000000e-12 * tt3 * 0.20 + 3.615080560000000e+03 * tt4;
    } else {
      h0_RT[4] = 3.092887670000000e+00 + 5.484297160000000e-04 * tt0 * 0.50 + 1.265052280000000e-07 * tt1 * OneThird + -8.794615559999999e-11 * tt2 * 0.25 + 1.174123760000000e-14 * tt3 * 0.20 + 3.858657000000000e+03 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[5] = 4.198640560000000e+00 + -2.036434100000000e-03 * tt0 * 0.50 + 6.520402110000000e-06 * tt1 * OneThird + -5.487970620000000e-09 * tt2 * 0.25 + 1.771978170000000e-12 * tt3 * 0.20 + -3.029372670000000e+04 * tt4;
    } else {
      h0_RT[5] = 3.033992490000000e+00 + 2.176918040000000e-03 * tt0 * 0.50 + -1.640725180000000e-07 * tt1 * OneThird + -9.704198700000000e-11 * tt2 * 0.25 + 1.682009920000000e-14 * tt3 * 0.20 + -3.000429710000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[6] = 4.301798010000000e+00 + -4.749120510000000e-03 * tt0 * 0.50 + 2.115828910000000e-05 * tt1 * OneThird + -2.427638940000000e-08 * tt2 * 0.25 + 9.292251240000000e-12 * tt3 * 0.20 + 2.948080400000000e+02 * tt4;
    } else {
      h0_RT[6] = 4.017210900000000e+00 + 2.239820130000000e-03 * tt0 * 0.50 + -6.336581500000000e-07 * tt1 * OneThird + 1.142463700000000e-10 * tt2 * 0.25 + -1.079085350000000e-14 * tt3 * 0.20 + 1.118567130000000e+02 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[7] = 4.276112690000000e+00 + -5.428224169999999e-04 * tt0 * 0.50 + 1.673357010000000e-05 * tt1 * OneThird + -2.157708130000000e-08 * tt2 * 0.25 + 8.624543630000000e-12 * tt3 * 0.20 + -1.770258210000000e+04 * tt4;
    } else {
      h0_RT[7] = 4.165002850000000e+00 + 4.908316940000000e-03 * tt0 * 0.50 + -1.901392250000000e-06 * tt1 * OneThird + 3.711859860000000e-10 * tt2 * 0.25 + -2.879083050000000e-14 * tt3 * 0.20 + -1.786178770000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[8] = 2.554239550000000e+00 + -3.215377240000000e-04 * tt0 * 0.50 + 7.337922450000000e-07 * tt1 * OneThird + -7.322348890000000e-10 * tt2 * 0.25 + 2.665214460000000e-13 * tt3 * 0.20 + 8.544388320000000e+04 * tt4;
    } else {
      h0_RT[8] = 2.492668880000000e+00 + 4.798892840000000e-05 * tt0 * 0.50 + -7.243350200000001e-08 * tt1 * OneThird + 3.742910290000000e-11 * tt2 * 0.25 + -4.872778930000000e-15 * tt3 * 0.20 + 8.545129530000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[9] = 3.489816650000000e+00 + 3.238355410000000e-04 * tt0 * 0.50 + -1.688990650000000e-06 * tt1 * OneThird + 3.162173270000000e-09 * tt2 * 0.25 + -1.406090670000000e-12 * tt3 * 0.20 + 7.079729340000000e+04 * tt4;
    } else {
      h0_RT[9] = 2.878464730000000e+00 + 9.709136810000000e-04 * tt0 * 0.50 + 1.444456550000000e-07 * tt1 * OneThird + -1.306878490000000e-10 * tt2 * 0.25 + 1.760793830000000e-14 * tt3 * 0.20 + 7.101243640000001e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[10] = 3.762678670000000e+00 + 9.688721430000000e-04 * tt0 * 0.50 + 2.794898410000000e-06 * tt1 * OneThird + -3.850911530000000e-09 * tt2 * 0.25 + 1.687417190000000e-12 * tt3 * 0.20 + 4.600404010000000e+04 * tt4;
    } else {
      h0_RT[10] = 2.874101130000000e+00 + 3.656392920000000e-03 * tt0 * 0.50 + -1.408945970000000e-06 * tt1 * OneThird + 2.601795490000000e-10 * tt2 * 0.25 + -1.877275670000000e-14 * tt3 * 0.20 + 4.626360400000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[11] = 4.198604110000000e+00 + -2.366614190000000e-03 * tt0 * 0.50 + 8.232962200000000e-06 * tt1 * OneThird + -6.688159810000000e-09 * tt2 * 0.25 + 1.943147370000000e-12 * tt3 * 0.20 + 5.049681630000000e+04 * tt4;
    } else {
      h0_RT[11] = 2.292038420000000e+00 + 4.655886370000000e-03 * tt0 * 0.50 + -2.011919470000000e-06 * tt1 * OneThird + 4.179060000000000e-10 * tt2 * 0.25 + -3.397163650000000e-14 * tt3 * 0.20 + 5.092599970000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[12] = 3.673590400000000e+00 + 2.010951750000000e-03 * tt0 * 0.50 + 5.730218560000000e-06 * tt1 * OneThird + -6.871174250000000e-09 * tt2 * 0.25 + 2.543857340000000e-12 * tt3 * 0.20 + 1.644499880000000e+04 * tt4;
    } else {
      h0_RT[12] = 2.285717720000000e+00 + 7.239900370000000e-03 * tt0 * 0.50 + -2.987143480000000e-06 * tt1 * OneThird + 5.956846440000000e-10 * tt2 * 0.25 + -4.671543940000000e-14 * tt3 * 0.20 + 1.677558430000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[13] = 5.149876130000000e+00 + -1.367097880000000e-02 * tt0 * 0.50 + 4.918005990000000e-05 * tt1 * OneThird + -4.847430260000000e-08 * tt2 * 0.25 + 1.666939560000000e-11 * tt3 * 0.20 + -1.024664760000000e+04 * tt4;
    } else {
      h0_RT[13] = 7.485149500000000e-02 + 1.339094670000000e-02 * tt0 * 0.50 + -5.732858090000000e-06 * tt1 * OneThird + 1.222925350000000e-09 * tt2 * 0.25 + -1.018152300000000e-13 * tt3 * 0.20 + -9.468344590000001e+03 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[14] = 3.579533470000000e+00 + -6.103536800000000e-04 * tt0 * 0.50 + 1.016814330000000e-06 * tt1 * OneThird + 9.070058840000000e-10 * tt2 * 0.25 + -9.044244990000000e-13 * tt3 * 0.20 + -1.434408600000000e+04 * tt4;
    } else {
      h0_RT[14] = 2.715185610000000e+00 + 2.062527430000000e-03 * tt0 * 0.50 + -9.988257710000001e-07 * tt1 * OneThird + 2.300530080000000e-10 * tt2 * 0.25 + -2.036477160000000e-14 * tt3 * 0.20 + -1.415187240000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[15] = 2.356773520000000e+00 + 8.984596770000000e-03 * tt0 * 0.50 + -7.123562690000000e-06 * tt1 * OneThird + 2.459190220000000e-09 * tt2 * 0.25 + -1.436995480000000e-13 * tt3 * 0.20 + -4.837196970000000e+04 * tt4;
    } else {
      h0_RT[15] = 3.857460290000000e+00 + 4.414370260000000e-03 * tt0 * 0.50 + -2.214814040000000e-06 * tt1 * OneThird + 5.234901880000001e-10 * tt2 * 0.25 + -4.720841640000000e-14 * tt3 * 0.20 + -4.875916600000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[16] = 4.221185840000000e+00 + -3.243925320000000e-03 * tt0 * 0.50 + 1.377994460000000e-05 * tt1 * OneThird + -1.331440930000000e-08 * tt2 * 0.25 + 4.337688650000000e-12 * tt3 * 0.20 + 3.839564960000000e+03 * tt4;
    } else {
      h0_RT[16] = 2.772174380000000e+00 + 4.956955260000000e-03 * tt0 * 0.50 + -2.484456130000000e-06 * tt1 * OneThird + 5.891617780000000e-10 * tt2 * 0.25 + -5.335087110000000e-14 * tt3 * 0.20 + 4.011918150000000e+03 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[17] = 4.793723150000000e+00 + -9.908333690000000e-03 * tt0 * 0.50 + 3.732200080000000e-05 * tt1 * OneThird + -3.792852610000000e-08 * tt2 * 0.25 + 1.317726520000000e-11 * tt3 * 0.20 + -1.430895670000000e+04 * tt4;
    } else {
      h0_RT[17] = 1.760690080000000e+00 + 9.200000820000000e-03 * tt0 * 0.50 + -4.422588130000000e-06 * tt1 * OneThird + 1.006412120000000e-09 * tt2 * 0.25 + -8.838556400000001e-14 * tt3 * 0.20 + -1.399583230000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[18] = 3.863889180000000e+00 + 5.596723040000000e-03 * tt0 * 0.50 + 5.932717910000000e-06 * tt1 * OneThird + -1.045320120000000e-08 * tt2 * 0.25 + 4.369672780000000e-12 * tt3 * 0.20 + -3.193913670000000e+03 * tt4;
    } else {
      h0_RT[18] = 3.692665690000000e+00 + 8.645767970000001e-03 * tt0 * 0.50 + -3.751011200000000e-06 * tt1 * OneThird + 7.872346360000000e-10 * tt2 * 0.25 + -6.485542010000000e-14 * tt3 * 0.20 + -3.242506270000000e+03 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[19] = 2.106083740868644e+00 + 7.217545611840916e-03 * tt0 * 0.50 + 5.335825777358757e-06 * tt1 * OneThird + -7.374562574086794e-09 * tt2 * 0.25 + 2.074342158630007e-12 * tt3 * 0.20 + 9.786125620403312e+02 * tt4;
    } else {
      h0_RT[19] = 3.771000099870984e+00 + 7.871073375136668e-03 * tt0 * 0.50 + -2.656064859783170e-06 * tt1 * OneThird + 3.943402155606197e-10 * tt2 * 0.25 + -2.111411617357040e-14 * tt3 * 0.20 + 1.277614246529541e+02 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[20] = 5.715395820000000e+00 + -1.523091290000000e-02 * tt0 * 0.50 + 6.524411550000000e-05 * tt1 * OneThird + -7.108068890000000e-08 * tt2 * 0.25 + 2.613526980000000e-11 * tt3 * 0.20 + -2.564276560000000e+04 * tt4;
    } else {
      h0_RT[20] = 1.789707910000000e+00 + 1.409382920000000e-02 * tt0 * 0.50 + -6.365008350000000e-06 * tt1 * OneThird + 1.381710850000000e-09 * tt2 * 0.25 + -1.170602200000000e-13 * tt3 * 0.20 + -2.537487470000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[21] = 2.889657330000000e+00 + 1.340996110000000e-02 * tt0 * 0.50 + -2.847695010000000e-05 * tt1 * OneThird + 2.947910450000000e-08 * tt2 * 0.25 + -1.093315110000000e-11 * tt3 * 0.20 + 6.683939320000001e+04 * tt4;
    } else {
      h0_RT[21] = 3.167806520000000e+00 + 4.752219020000000e-03 * tt0 * 0.50 + -1.837870770000000e-06 * tt1 * OneThird + 3.041902520000000e-10 * tt2 * 0.25 + -1.772327700000000e-14 * tt3 * 0.20 + 6.712106500000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[22] = 8.086810940000000e-01 + 2.336156290000000e-02 * tt0 * 0.50 + -3.551718150000000e-05 * tt1 * OneThird + 2.801524370000000e-08 * tt2 * 0.25 + -8.500729740000000e-12 * tt3 * 0.20 + 2.642898070000000e+04 * tt4;
    } else {
      h0_RT[22] = 4.147569640000000e+00 + 5.961666640000000e-03 * tt0 * 0.50 + -2.372948520000000e-06 * tt1 * OneThird + 4.674121710000000e-10 * tt2 * 0.25 + -3.612352130000000e-14 * tt3 * 0.20 + 2.593599920000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[23] = 3.212466450000000e+00 + 1.514791620000000e-03 * tt0 * 0.50 + 2.592094120000000e-05 * tt1 * OneThird + -3.576578470000000e-08 * tt2 * 0.25 + 1.471508730000000e-11 * tt3 * 0.20 + 3.485984680000000e+04 * tt4;
    } else {
      h0_RT[23] = 3.016724000000000e+00 + 1.033022920000000e-02 * tt0 * 0.50 + -4.680823490000000e-06 * tt1 * OneThird + 1.017632880000000e-09 * tt2 * 0.25 + -8.626070410000000e-14 * tt3 * 0.20 + 3.461287390000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[24] = 3.959201480000000e+00 + -7.570522470000000e-03 * tt0 * 0.50 + 5.709902920000000e-05 * tt1 * OneThird + -6.915887530000000e-08 * tt2 * 0.25 + 2.698843730000000e-11 * tt3 * 0.20 + 5.089775930000000e+03 * tt4;
    } else {
      h0_RT[24] = 2.036111160000000e+00 + 1.464541510000000e-02 * tt0 * 0.50 + -6.710779150000000e-06 * tt1 * OneThird + 1.472229230000000e-09 * tt2 * 0.25 + -1.257060610000000e-13 * tt3 * 0.20 + 4.939886140000000e+03 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[25] = 4.306465680000000e+00 + -4.186588920000000e-03 * tt0 * 0.50 + 4.971428070000000e-05 * tt1 * OneThird + -5.991266060000000e-08 * tt2 * 0.25 + 2.305090040000000e-11 * tt3 * 0.20 + 1.284162650000000e+04 * tt4;
    } else {
      h0_RT[25] = 1.954656420000000e+00 + 1.739727220000000e-02 * tt0 * 0.50 + -7.982066680000000e-06 * tt1 * OneThird + 1.752176890000000e-09 * tt2 * 0.25 + -1.496415760000000e-13 * tt3 * 0.20 + 1.285752000000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[26] = 4.291424920000000e+00 + -5.501542700000000e-03 * tt0 * 0.50 + 5.994382880000000e-05 * tt1 * OneThird + -7.084662850000000e-08 * tt2 * 0.25 + 2.686857710000000e-11 * tt3 * 0.20 + -1.152220550000000e+04 * tt4;
    } else {
      h0_RT[26] = 1.071881500000000e+00 + 2.168526770000000e-02 * tt0 * 0.50 + -1.002560670000000e-05 * tt1 * OneThird + 2.214120010000000e-09 * tt2 * 0.25 + -1.900028900000000e-13 * tt3 * 0.20 + -1.142639320000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[27] = 2.251845286753941e+00 + 1.765401529121510e-02 * tt0 * 0.50 + -2.372620998746318e-05 * tt1 * OneThird + 1.727224682516638e-08 * tt2 * 0.25 + -5.064957679024219e-12 * tt3 * 0.20 + 2.005943757764270e+04 * tt4;
    } else {
      h0_RT[27] = 5.628172209693201e+00 + 4.085394191788521e-03 * tt0 * 0.50 + -1.593486747480799e-06 * tt1 * OneThird + 2.862686614342911e-10 * tt2 * 0.25 + -1.940857878719715e-14 * tt3 * 0.20 + 1.932723151197489e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[28] = 2.135836300000000e+00 + 1.811887210000000e-02 * tt0 * 0.50 + -1.739474740000000e-05 * tt1 * OneThird + 9.343975680000000e-09 * tt2 * 0.25 + -2.014576150000000e-12 * tt3 * 0.20 + -7.042918040000000e+03 * tt4;
    } else {
      h0_RT[28] = 4.511297320000000e+00 + 9.003597449999999e-03 * tt0 * 0.50 + -4.169396350000000e-06 * tt1 * OneThird + 9.233458820000000e-10 * tt2 * 0.25 + -7.948382010000000e-14 * tt3 * 0.20 + -7.551053110000000e+03 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[29] = 1.242373300000000e+00 + 3.107220100000000e-02 * tt0 * 0.50 + -5.086686400000000e-05 * tt1 * OneThird + 4.313713100000000e-08 * tt2 * 0.25 + -1.401459400000000e-11 * tt3 * 0.20 + 8.031614300000000e+03 * tt4;
    } else {
      h0_RT[29] = 5.923829100000000e+00 + 6.792360000000000e-03 * tt0 * 0.50 + -2.565856400000000e-06 * tt1 * OneThird + 4.498784100000000e-10 * tt2 * 0.25 + -2.994010100000000e-14 * tt3 * 0.20 + 7.264626000000000e+03 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[30] = 3.298616281368291e+00 + 1.408708904146039e-03 * tt0 * 0.50 + -3.964481129109908e-06 * tt1 * OneThird + 5.642920880408571e-09 * tt2 * 0.25 + -2.445407041148433e-12 * tt3 * 0.20 + -1.020894198687962e+03 * tt4;
    } else {
      h0_RT[30] = 2.926639911210682e+00 + 1.487977101178227e-03 * tt0 * 0.50 + -5.684761849244810e-07 * tt1 * OneThird + 1.009704225872734e-10 * tt2 * 0.25 + -6.753354387142974e-15 * tt3 * 0.20 + -9.227966980051905e+02 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[31] = 2.500000000000000e+00 + 0.000000000000000e+00 * tt0 * 0.50 + 0.000000000000000e+00 * tt1 * OneThird + 0.000000000000000e+00 * tt2 * 0.25 + 0.000000000000000e+00 * tt3 * 0.20 + -7.453750000000000e+02 * tt4;
    } else {
      h0_RT[31] = 2.500000000000000e+00 + 0.000000000000000e+00 * tt0 * 0.50 + 0.000000000000000e+00 * tt1 * OneThird + 0.000000000000000e+00 * tt2 * 0.25 + 0.000000000000000e+00 * tt3 * 0.20 + -7.453750000000000e+02 * tt4;
    };

  };

  template <class Type>
    void GasThermo::getEnthalpiesDerivatives(Type& T, vector<Type>& dh0dT) {

    Type tt0 = T;
    Type tt1 = T * tt0;
    Type tt2 = T * tt1;
    Type tt3 = T * tt2;

    if(tt0 < 1.0000e+03) {
      dh0dT[0] = 2.344331120000000e+00 + 7.980520749999999e-03 * tt0 + -1.947815100000000e-05 * tt1 + 2.015720940000000e-08 * tt2 + -7.376117610000001e-12 * tt3;
    } else {
      dh0dT[0] = 3.337279200000000e+00 + -4.940247310000000e-05 * tt0 + 4.994567780000000e-07 * tt1 + -1.795663940000000e-10 * tt2 + 2.002553760000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[1] = 2.500000000000000e+00 + 7.053328190000000e-13 * tt0 + -1.995919640000000e-15 * tt1 + 2.300816320000000e-18 * tt2 + -9.277323320000001e-22 * tt3;
    } else {
      dh0dT[1] = 2.500000010000000e+00 + -2.308429730000000e-11 * tt0 + 1.615619480000000e-14 * tt1 + -4.735152350000000e-18 * tt2 + 4.981973570000000e-22 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[2] = 3.168267100000000e+00 + -3.279318840000000e-03 * tt0 + 6.643063960000000e-06 * tt1 + -6.128066240000000e-09 * tt2 + 2.112659710000000e-12 * tt3;
    } else {
      dh0dT[2] = 2.569420780000000e+00 + -8.597411370000000e-05 * tt0 + 4.194845890000000e-08 * tt1 + -1.001777990000000e-11 * tt2 + 1.228336910000000e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[3] = 3.782456360000000e+00 + -2.996734160000000e-03 * tt0 + 9.847302010000000e-06 * tt1 + -9.681295090000001e-09 * tt2 + 3.243728370000000e-12 * tt3;
    } else {
      dh0dT[3] = 3.282537840000000e+00 + 1.483087540000000e-03 * tt0 + -7.579666690000000e-07 * tt1 + 2.094705550000000e-10 * tt2 + -2.167177940000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[4] = 3.992015430000000e+00 + -2.401317520000000e-03 * tt0 + 4.617938410000000e-06 * tt1 + -3.881133330000000e-09 * tt2 + 1.364114700000000e-12 * tt3;
    } else {
      dh0dT[4] = 3.092887670000000e+00 + 5.484297160000000e-04 * tt0 + 1.265052280000000e-07 * tt1 + -8.794615559999999e-11 * tt2 + 1.174123760000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[5] = 4.198640560000000e+00 + -2.036434100000000e-03 * tt0 + 6.520402110000000e-06 * tt1 + -5.487970620000000e-09 * tt2 + 1.771978170000000e-12 * tt3;
    } else {
      dh0dT[5] = 3.033992490000000e+00 + 2.176918040000000e-03 * tt0 + -1.640725180000000e-07 * tt1 + -9.704198700000000e-11 * tt2 + 1.682009920000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[6] = 4.301798010000000e+00 + -4.749120510000000e-03 * tt0 + 2.115828910000000e-05 * tt1 + -2.427638940000000e-08 * tt2 + 9.292251240000000e-12 * tt3;
    } else {
      dh0dT[6] = 4.017210900000000e+00 + 2.239820130000000e-03 * tt0 + -6.336581500000000e-07 * tt1 + 1.142463700000000e-10 * tt2 + -1.079085350000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[7] = 4.276112690000000e+00 + -5.428224169999999e-04 * tt0 + 1.673357010000000e-05 * tt1 + -2.157708130000000e-08 * tt2 + 8.624543630000000e-12 * tt3;
    } else {
      dh0dT[7] = 4.165002850000000e+00 + 4.908316940000000e-03 * tt0 + -1.901392250000000e-06 * tt1 + 3.711859860000000e-10 * tt2 + -2.879083050000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[8] = 2.554239550000000e+00 + -3.215377240000000e-04 * tt0 + 7.337922450000000e-07 * tt1 + -7.322348890000000e-10 * tt2 + 2.665214460000000e-13 * tt3;
    } else {
      dh0dT[8] = 2.492668880000000e+00 + 4.798892840000000e-05 * tt0 + -7.243350200000001e-08 * tt1 + 3.742910290000000e-11 * tt2 + -4.872778930000000e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[9] = 3.489816650000000e+00 + 3.238355410000000e-04 * tt0 + -1.688990650000000e-06 * tt1 + 3.162173270000000e-09 * tt2 + -1.406090670000000e-12 * tt3;
    } else {
      dh0dT[9] = 2.878464730000000e+00 + 9.709136810000000e-04 * tt0 + 1.444456550000000e-07 * tt1 + -1.306878490000000e-10 * tt2 + 1.760793830000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[10] = 3.762678670000000e+00 + 9.688721430000000e-04 * tt0 + 2.794898410000000e-06 * tt1 + -3.850911530000000e-09 * tt2 + 1.687417190000000e-12 * tt3;
    } else {
      dh0dT[10] = 2.874101130000000e+00 + 3.656392920000000e-03 * tt0 + -1.408945970000000e-06 * tt1 + 2.601795490000000e-10 * tt2 + -1.877275670000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[11] = 4.198604110000000e+00 + -2.366614190000000e-03 * tt0 + 8.232962200000000e-06 * tt1 + -6.688159810000000e-09 * tt2 + 1.943147370000000e-12 * tt3;
    } else {
      dh0dT[11] = 2.292038420000000e+00 + 4.655886370000000e-03 * tt0 + -2.011919470000000e-06 * tt1 + 4.179060000000000e-10 * tt2 + -3.397163650000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[12] = 3.673590400000000e+00 + 2.010951750000000e-03 * tt0 + 5.730218560000000e-06 * tt1 + -6.871174250000000e-09 * tt2 + 2.543857340000000e-12 * tt3;
    } else {
      dh0dT[12] = 2.285717720000000e+00 + 7.239900370000000e-03 * tt0 + -2.987143480000000e-06 * tt1 + 5.956846440000000e-10 * tt2 + -4.671543940000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[13] = 5.149876130000000e+00 + -1.367097880000000e-02 * tt0 + 4.918005990000000e-05 * tt1 + -4.847430260000000e-08 * tt2 + 1.666939560000000e-11 * tt3;
    } else {
      dh0dT[13] = 7.485149500000000e-02 + 1.339094670000000e-02 * tt0 + -5.732858090000000e-06 * tt1 + 1.222925350000000e-09 * tt2 + -1.018152300000000e-13 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[14] = 3.579533470000000e+00 + -6.103536800000000e-04 * tt0 + 1.016814330000000e-06 * tt1 + 9.070058840000000e-10 * tt2 + -9.044244990000000e-13 * tt3;
    } else {
      dh0dT[14] = 2.715185610000000e+00 + 2.062527430000000e-03 * tt0 + -9.988257710000001e-07 * tt1 + 2.300530080000000e-10 * tt2 + -2.036477160000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[15] = 2.356773520000000e+00 + 8.984596770000000e-03 * tt0 + -7.123562690000000e-06 * tt1 + 2.459190220000000e-09 * tt2 + -1.436995480000000e-13 * tt3;
    } else {
      dh0dT[15] = 3.857460290000000e+00 + 4.414370260000000e-03 * tt0 + -2.214814040000000e-06 * tt1 + 5.234901880000001e-10 * tt2 + -4.720841640000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[16] = 4.221185840000000e+00 + -3.243925320000000e-03 * tt0 + 1.377994460000000e-05 * tt1 + -1.331440930000000e-08 * tt2 + 4.337688650000000e-12 * tt3;
    } else {
      dh0dT[16] = 2.772174380000000e+00 + 4.956955260000000e-03 * tt0 + -2.484456130000000e-06 * tt1 + 5.891617780000000e-10 * tt2 + -5.335087110000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[17] = 4.793723150000000e+00 + -9.908333690000000e-03 * tt0 + 3.732200080000000e-05 * tt1 + -3.792852610000000e-08 * tt2 + 1.317726520000000e-11 * tt3;
    } else {
      dh0dT[17] = 1.760690080000000e+00 + 9.200000820000000e-03 * tt0 + -4.422588130000000e-06 * tt1 + 1.006412120000000e-09 * tt2 + -8.838556400000001e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[18] = 3.863889180000000e+00 + 5.596723040000000e-03 * tt0 + 5.932717910000000e-06 * tt1 + -1.045320120000000e-08 * tt2 + 4.369672780000000e-12 * tt3;
    } else {
      dh0dT[18] = 3.692665690000000e+00 + 8.645767970000001e-03 * tt0 + -3.751011200000000e-06 * tt1 + 7.872346360000000e-10 * tt2 + -6.485542010000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[19] = 2.106083740868644e+00 + 7.217545611840916e-03 * tt0 + 5.335825777358757e-06 * tt1 + -7.374562574086794e-09 * tt2 + 2.074342158630007e-12 * tt3;
    } else {
      dh0dT[19] = 3.771000099870984e+00 + 7.871073375136668e-03 * tt0 + -2.656064859783170e-06 * tt1 + 3.943402155606197e-10 * tt2 + -2.111411617357040e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[20] = 5.715395820000000e+00 + -1.523091290000000e-02 * tt0 + 6.524411550000000e-05 * tt1 + -7.108068890000000e-08 * tt2 + 2.613526980000000e-11 * tt3;
    } else {
      dh0dT[20] = 1.789707910000000e+00 + 1.409382920000000e-02 * tt0 + -6.365008350000000e-06 * tt1 + 1.381710850000000e-09 * tt2 + -1.170602200000000e-13 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[21] = 2.889657330000000e+00 + 1.340996110000000e-02 * tt0 + -2.847695010000000e-05 * tt1 + 2.947910450000000e-08 * tt2 + -1.093315110000000e-11 * tt3;
    } else {
      dh0dT[21] = 3.167806520000000e+00 + 4.752219020000000e-03 * tt0 + -1.837870770000000e-06 * tt1 + 3.041902520000000e-10 * tt2 + -1.772327700000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[22] = 8.086810940000000e-01 + 2.336156290000000e-02 * tt0 + -3.551718150000000e-05 * tt1 + 2.801524370000000e-08 * tt2 + -8.500729740000000e-12 * tt3;
    } else {
      dh0dT[22] = 4.147569640000000e+00 + 5.961666640000000e-03 * tt0 + -2.372948520000000e-06 * tt1 + 4.674121710000000e-10 * tt2 + -3.612352130000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[23] = 3.212466450000000e+00 + 1.514791620000000e-03 * tt0 + 2.592094120000000e-05 * tt1 + -3.576578470000000e-08 * tt2 + 1.471508730000000e-11 * tt3;
    } else {
      dh0dT[23] = 3.016724000000000e+00 + 1.033022920000000e-02 * tt0 + -4.680823490000000e-06 * tt1 + 1.017632880000000e-09 * tt2 + -8.626070410000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[24] = 3.959201480000000e+00 + -7.570522470000000e-03 * tt0 + 5.709902920000000e-05 * tt1 + -6.915887530000000e-08 * tt2 + 2.698843730000000e-11 * tt3;
    } else {
      dh0dT[24] = 2.036111160000000e+00 + 1.464541510000000e-02 * tt0 + -6.710779150000000e-06 * tt1 + 1.472229230000000e-09 * tt2 + -1.257060610000000e-13 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[25] = 4.306465680000000e+00 + -4.186588920000000e-03 * tt0 + 4.971428070000000e-05 * tt1 + -5.991266060000000e-08 * tt2 + 2.305090040000000e-11 * tt3;
    } else {
      dh0dT[25] = 1.954656420000000e+00 + 1.739727220000000e-02 * tt0 + -7.982066680000000e-06 * tt1 + 1.752176890000000e-09 * tt2 + -1.496415760000000e-13 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[26] = 4.291424920000000e+00 + -5.501542700000000e-03 * tt0 + 5.994382880000000e-05 * tt1 + -7.084662850000000e-08 * tt2 + 2.686857710000000e-11 * tt3;
    } else {
      dh0dT[26] = 1.071881500000000e+00 + 2.168526770000000e-02 * tt0 + -1.002560670000000e-05 * tt1 + 2.214120010000000e-09 * tt2 + -1.900028900000000e-13 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[27] = 2.251845286753941e+00 + 1.765401529121510e-02 * tt0 + -2.372620998746318e-05 * tt1 + 1.727224682516638e-08 * tt2 + -5.064957679024219e-12 * tt3;
    } else {
      dh0dT[27] = 5.628172209693201e+00 + 4.085394191788521e-03 * tt0 + -1.593486747480799e-06 * tt1 + 2.862686614342911e-10 * tt2 + -1.940857878719715e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[28] = 2.135836300000000e+00 + 1.811887210000000e-02 * tt0 + -1.739474740000000e-05 * tt1 + 9.343975680000000e-09 * tt2 + -2.014576150000000e-12 * tt3;
    } else {
      dh0dT[28] = 4.511297320000000e+00 + 9.003597449999999e-03 * tt0 + -4.169396350000000e-06 * tt1 + 9.233458820000000e-10 * tt2 + -7.948382010000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[29] = 1.242373300000000e+00 + 3.107220100000000e-02 * tt0 + -5.086686400000000e-05 * tt1 + 4.313713100000000e-08 * tt2 + -1.401459400000000e-11 * tt3;
    } else {
      dh0dT[29] = 5.923829100000000e+00 + 6.792360000000000e-03 * tt0 + -2.565856400000000e-06 * tt1 + 4.498784100000000e-10 * tt2 + -2.994010100000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[30] = 3.298616281368291e+00 + 1.408708904146039e-03 * tt0 + -3.964481129109908e-06 * tt1 + 5.642920880408571e-09 * tt2 + -2.445407041148433e-12 * tt3;
    } else {
      dh0dT[30] = 2.926639911210682e+00 + 1.487977101178227e-03 * tt0 + -5.684761849244810e-07 * tt1 + 1.009704225872734e-10 * tt2 + -6.753354387142974e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[31] = 2.500000000000000e+00 + 0.000000000000000e+00 * tt0 + 0.000000000000000e+00 * tt1 + 0.000000000000000e+00 * tt2 + 0.000000000000000e+00 * tt3;
    } else {
      dh0dT[31] = 2.500000000000000e+00 + 0.000000000000000e+00 * tt0 + 0.000000000000000e+00 * tt1 + 0.000000000000000e+00 * tt2 + 0.000000000000000e+00 * tt3;
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
      s0_R[0] = 2.344331120000000e+00 * tt5 + 7.980520749999999e-03 * tt0 + -1.947815100000000e-05 * tt1 * 0.50 + 2.015720940000000e-08 * tt2 * OneThird + -7.376117610000001e-12 * tt3 * 0.25 + 6.830102380000000e-01;
    } else {
      s0_R[0] = 3.337279200000000e+00 * tt5 +  -4.940247310000000e-05 * tt0 + 4.994567780000000e-07 * tt1 * 0.50 + -1.795663940000000e-10 * tt2 * OneThird + 2.002553760000000e-14 * tt3 * 0.25 + -3.205023310000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[1] = 2.500000000000000e+00 * tt5 + 7.053328190000000e-13 * tt0 + -1.995919640000000e-15 * tt1 * 0.50 + 2.300816320000000e-18 * tt2 * OneThird + -9.277323320000001e-22 * tt3 * 0.25 + -4.466828530000000e-01;
    } else {
      s0_R[1] = 2.500000010000000e+00 * tt5 +  -2.308429730000000e-11 * tt0 + 1.615619480000000e-14 * tt1 * 0.50 + -4.735152350000000e-18 * tt2 * OneThird + 4.981973570000000e-22 * tt3 * 0.25 + -4.466829140000000e-01;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[2] = 3.168267100000000e+00 * tt5 + -3.279318840000000e-03 * tt0 + 6.643063960000000e-06 * tt1 * 0.50 + -6.128066240000000e-09 * tt2 * OneThird + 2.112659710000000e-12 * tt3 * 0.25 + 2.051933460000000e+00;
    } else {
      s0_R[2] = 2.569420780000000e+00 * tt5 +  -8.597411370000000e-05 * tt0 + 4.194845890000000e-08 * tt1 * 0.50 + -1.001777990000000e-11 * tt2 * OneThird + 1.228336910000000e-15 * tt3 * 0.25 + 4.784338640000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[3] = 3.782456360000000e+00 * tt5 + -2.996734160000000e-03 * tt0 + 9.847302010000000e-06 * tt1 * 0.50 + -9.681295090000001e-09 * tt2 * OneThird + 3.243728370000000e-12 * tt3 * 0.25 + 3.657675730000000e+00;
    } else {
      s0_R[3] = 3.282537840000000e+00 * tt5 +  1.483087540000000e-03 * tt0 + -7.579666690000000e-07 * tt1 * 0.50 + 2.094705550000000e-10 * tt2 * OneThird + -2.167177940000000e-14 * tt3 * 0.25 + 5.453231290000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[4] = 3.992015430000000e+00 * tt5 + -2.401317520000000e-03 * tt0 + 4.617938410000000e-06 * tt1 * 0.50 + -3.881133330000000e-09 * tt2 * OneThird + 1.364114700000000e-12 * tt3 * 0.25 + -1.039254580000000e-01;
    } else {
      s0_R[4] = 3.092887670000000e+00 * tt5 +  5.484297160000000e-04 * tt0 + 1.265052280000000e-07 * tt1 * 0.50 + -8.794615559999999e-11 * tt2 * OneThird + 1.174123760000000e-14 * tt3 * 0.25 + 4.476696100000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[5] = 4.198640560000000e+00 * tt5 + -2.036434100000000e-03 * tt0 + 6.520402110000000e-06 * tt1 * 0.50 + -5.487970620000000e-09 * tt2 * OneThird + 1.771978170000000e-12 * tt3 * 0.25 + -8.490322080000000e-01;
    } else {
      s0_R[5] = 3.033992490000000e+00 * tt5 +  2.176918040000000e-03 * tt0 + -1.640725180000000e-07 * tt1 * 0.50 + -9.704198700000000e-11 * tt2 * OneThird + 1.682009920000000e-14 * tt3 * 0.25 + 4.966770100000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[6] = 4.301798010000000e+00 * tt5 + -4.749120510000000e-03 * tt0 + 2.115828910000000e-05 * tt1 * 0.50 + -2.427638940000000e-08 * tt2 * OneThird + 9.292251240000000e-12 * tt3 * 0.25 + 3.716662450000000e+00;
    } else {
      s0_R[6] = 4.017210900000000e+00 * tt5 +  2.239820130000000e-03 * tt0 + -6.336581500000000e-07 * tt1 * 0.50 + 1.142463700000000e-10 * tt2 * OneThird + -1.079085350000000e-14 * tt3 * 0.25 + 3.785102150000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[7] = 4.276112690000000e+00 * tt5 + -5.428224169999999e-04 * tt0 + 1.673357010000000e-05 * tt1 * 0.50 + -2.157708130000000e-08 * tt2 * OneThird + 8.624543630000000e-12 * tt3 * 0.25 + 3.435050740000000e+00;
    } else {
      s0_R[7] = 4.165002850000000e+00 * tt5 +  4.908316940000000e-03 * tt0 + -1.901392250000000e-06 * tt1 * 0.50 + 3.711859860000000e-10 * tt2 * OneThird + -2.879083050000000e-14 * tt3 * 0.25 + 2.916156620000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[8] = 2.554239550000000e+00 * tt5 + -3.215377240000000e-04 * tt0 + 7.337922450000000e-07 * tt1 * 0.50 + -7.322348890000000e-10 * tt2 * OneThird + 2.665214460000000e-13 * tt3 * 0.25 + 4.531308480000000e+00;
    } else {
      s0_R[8] = 2.492668880000000e+00 * tt5 +  4.798892840000000e-05 * tt0 + -7.243350200000001e-08 * tt1 * 0.50 + 3.742910290000000e-11 * tt2 * OneThird + -4.872778930000000e-15 * tt3 * 0.25 + 4.801503730000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[9] = 3.489816650000000e+00 * tt5 + 3.238355410000000e-04 * tt0 + -1.688990650000000e-06 * tt1 * 0.50 + 3.162173270000000e-09 * tt2 * OneThird + -1.406090670000000e-12 * tt3 * 0.25 + 2.084011080000000e+00;
    } else {
      s0_R[9] = 2.878464730000000e+00 * tt5 +  9.709136810000000e-04 * tt0 + 1.444456550000000e-07 * tt1 * 0.50 + -1.306878490000000e-10 * tt2 * OneThird + 1.760793830000000e-14 * tt3 * 0.25 + 5.484979990000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[10] = 3.762678670000000e+00 * tt5 + 9.688721430000000e-04 * tt0 + 2.794898410000000e-06 * tt1 * 0.50 + -3.850911530000000e-09 * tt2 * OneThird + 1.687417190000000e-12 * tt3 * 0.25 + 1.562531850000000e+00;
    } else {
      s0_R[10] = 2.874101130000000e+00 * tt5 +  3.656392920000000e-03 * tt0 + -1.408945970000000e-06 * tt1 * 0.50 + 2.601795490000000e-10 * tt2 * OneThird + -1.877275670000000e-14 * tt3 * 0.25 + 6.171193240000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[11] = 4.198604110000000e+00 * tt5 + -2.366614190000000e-03 * tt0 + 8.232962200000000e-06 * tt1 * 0.50 + -6.688159810000000e-09 * tt2 * OneThird + 1.943147370000000e-12 * tt3 * 0.25 + -7.691189670000000e-01;
    } else {
      s0_R[11] = 2.292038420000000e+00 * tt5 +  4.655886370000000e-03 * tt0 + -2.011919470000000e-06 * tt1 * 0.50 + 4.179060000000000e-10 * tt2 * OneThird + -3.397163650000000e-14 * tt3 * 0.25 + 8.626501690000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[12] = 3.673590400000000e+00 * tt5 + 2.010951750000000e-03 * tt0 + 5.730218560000000e-06 * tt1 * 0.50 + -6.871174250000000e-09 * tt2 * OneThird + 2.543857340000000e-12 * tt3 * 0.25 + 1.604564330000000e+00;
    } else {
      s0_R[12] = 2.285717720000000e+00 * tt5 +  7.239900370000000e-03 * tt0 + -2.987143480000000e-06 * tt1 * 0.50 + 5.956846440000000e-10 * tt2 * OneThird + -4.671543940000000e-14 * tt3 * 0.25 + 8.480071790000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[13] = 5.149876130000000e+00 * tt5 + -1.367097880000000e-02 * tt0 + 4.918005990000000e-05 * tt1 * 0.50 + -4.847430260000000e-08 * tt2 * OneThird + 1.666939560000000e-11 * tt3 * 0.25 + -4.641303760000000e+00;
    } else {
      s0_R[13] = 7.485149500000000e-02 * tt5 +  1.339094670000000e-02 * tt0 + -5.732858090000000e-06 * tt1 * 0.50 + 1.222925350000000e-09 * tt2 * OneThird + -1.018152300000000e-13 * tt3 * 0.25 + 1.843731800000000e+01;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[14] = 3.579533470000000e+00 * tt5 + -6.103536800000000e-04 * tt0 + 1.016814330000000e-06 * tt1 * 0.50 + 9.070058840000000e-10 * tt2 * OneThird + -9.044244990000000e-13 * tt3 * 0.25 + 3.508409280000000e+00;
    } else {
      s0_R[14] = 2.715185610000000e+00 * tt5 +  2.062527430000000e-03 * tt0 + -9.988257710000001e-07 * tt1 * 0.50 + 2.300530080000000e-10 * tt2 * OneThird + -2.036477160000000e-14 * tt3 * 0.25 + 7.818687720000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[15] = 2.356773520000000e+00 * tt5 + 8.984596770000000e-03 * tt0 + -7.123562690000000e-06 * tt1 * 0.50 + 2.459190220000000e-09 * tt2 * OneThird + -1.436995480000000e-13 * tt3 * 0.25 + 9.901052220000000e+00;
    } else {
      s0_R[15] = 3.857460290000000e+00 * tt5 +  4.414370260000000e-03 * tt0 + -2.214814040000000e-06 * tt1 * 0.50 + 5.234901880000001e-10 * tt2 * OneThird + -4.720841640000000e-14 * tt3 * 0.25 + 2.271638060000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[16] = 4.221185840000000e+00 * tt5 + -3.243925320000000e-03 * tt0 + 1.377994460000000e-05 * tt1 * 0.50 + -1.331440930000000e-08 * tt2 * OneThird + 4.337688650000000e-12 * tt3 * 0.25 + 3.394372430000000e+00;
    } else {
      s0_R[16] = 2.772174380000000e+00 * tt5 +  4.956955260000000e-03 * tt0 + -2.484456130000000e-06 * tt1 * 0.50 + 5.891617780000000e-10 * tt2 * OneThird + -5.335087110000000e-14 * tt3 * 0.25 + 9.798344920000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[17] = 4.793723150000000e+00 * tt5 + -9.908333690000000e-03 * tt0 + 3.732200080000000e-05 * tt1 * 0.50 + -3.792852610000000e-08 * tt2 * OneThird + 1.317726520000000e-11 * tt3 * 0.25 + 6.028129000000000e-01;
    } else {
      s0_R[17] = 1.760690080000000e+00 * tt5 +  9.200000820000000e-03 * tt0 + -4.422588130000000e-06 * tt1 * 0.50 + 1.006412120000000e-09 * tt2 * OneThird + -8.838556400000001e-14 * tt3 * 0.25 + 1.365632300000000e+01;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[18] = 3.863889180000000e+00 * tt5 + 5.596723040000000e-03 * tt0 + 5.932717910000000e-06 * tt1 * 0.50 + -1.045320120000000e-08 * tt2 * OneThird + 4.369672780000000e-12 * tt3 * 0.25 + 5.473022430000000e+00;
    } else {
      s0_R[18] = 3.692665690000000e+00 * tt5 +  8.645767970000001e-03 * tt0 + -3.751011200000000e-06 * tt1 * 0.50 + 7.872346360000000e-10 * tt2 * OneThird + -6.485542010000000e-14 * tt3 * 0.25 + 5.810432150000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[19] = 2.106083740868644e+00 * tt5 + 7.217545611840916e-03 * tt0 + 5.335825777358757e-06 * tt1 * 0.50 + -7.374562574086794e-09 * tt2 * OneThird + 2.074342158630007e-12 * tt3 * 0.25 + 1.315267015071675e+01;
    } else {
      s0_R[19] = 3.771000099870984e+00 * tt5 +  7.871073375136668e-03 * tt0 + -2.656064859783170e-06 * tt1 * 0.50 + 3.943402155606197e-10 * tt2 * OneThird + -2.111411617357040e-14 * tt3 * 0.25 + 2.928482743514920e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[20] = 5.715395820000000e+00 * tt5 + -1.523091290000000e-02 * tt0 + 6.524411550000000e-05 * tt1 * 0.50 + -7.108068890000000e-08 * tt2 * OneThird + 2.613526980000000e-11 * tt3 * 0.25 + -1.504098230000000e+00;
    } else {
      s0_R[20] = 1.789707910000000e+00 * tt5 +  1.409382920000000e-02 * tt0 + -6.365008350000000e-06 * tt1 * 0.50 + 1.381710850000000e-09 * tt2 * OneThird + -1.170602200000000e-13 * tt3 * 0.25 + 1.450236230000000e+01;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[21] = 2.889657330000000e+00 * tt5 + 1.340996110000000e-02 * tt0 + -2.847695010000000e-05 * tt1 * 0.50 + 2.947910450000000e-08 * tt2 * OneThird + -1.093315110000000e-11 * tt3 * 0.25 + 6.222964380000000e+00;
    } else {
      s0_R[21] = 3.167806520000000e+00 * tt5 +  4.752219020000000e-03 * tt0 + -1.837870770000000e-06 * tt1 * 0.50 + 3.041902520000000e-10 * tt2 * OneThird + -1.772327700000000e-14 * tt3 * 0.25 + 6.635894750000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[22] = 8.086810940000000e-01 * tt5 + 2.336156290000000e-02 * tt0 + -3.551718150000000e-05 * tt1 * 0.50 + 2.801524370000000e-08 * tt2 * OneThird + -8.500729740000000e-12 * tt3 * 0.25 + 1.393970510000000e+01;
    } else {
      s0_R[22] = 4.147569640000000e+00 * tt5 +  5.961666640000000e-03 * tt0 + -2.372948520000000e-06 * tt1 * 0.50 + 4.674121710000000e-10 * tt2 * OneThird + -3.612352130000000e-14 * tt3 * 0.25 + -1.230281210000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[23] = 3.212466450000000e+00 * tt5 + 1.514791620000000e-03 * tt0 + 2.592094120000000e-05 * tt1 * 0.50 + -3.576578470000000e-08 * tt2 * OneThird + 1.471508730000000e-11 * tt3 * 0.25 + 8.510540250000000e+00;
    } else {
      s0_R[23] = 3.016724000000000e+00 * tt5 +  1.033022920000000e-02 * tt0 + -4.680823490000000e-06 * tt1 * 0.50 + 1.017632880000000e-09 * tt2 * OneThird + -8.626070410000000e-14 * tt3 * 0.25 + 7.787323780000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[24] = 3.959201480000000e+00 * tt5 + -7.570522470000000e-03 * tt0 + 5.709902920000000e-05 * tt1 * 0.50 + -6.915887530000000e-08 * tt2 * OneThird + 2.698843730000000e-11 * tt3 * 0.25 + 4.097330960000000e+00;
    } else {
      s0_R[24] = 2.036111160000000e+00 * tt5 +  1.464541510000000e-02 * tt0 + -6.710779150000000e-06 * tt1 * 0.50 + 1.472229230000000e-09 * tt2 * OneThird + -1.257060610000000e-13 * tt3 * 0.25 + 1.030536930000000e+01;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[25] = 4.306465680000000e+00 * tt5 + -4.186588920000000e-03 * tt0 + 4.971428070000000e-05 * tt1 * 0.50 + -5.991266060000000e-08 * tt2 * OneThird + 2.305090040000000e-11 * tt3 * 0.25 + 4.707209240000000e+00;
    } else {
      s0_R[25] = 1.954656420000000e+00 * tt5 +  1.739727220000000e-02 * tt0 + -7.982066680000000e-06 * tt1 * 0.50 + 1.752176890000000e-09 * tt2 * OneThird + -1.496415760000000e-13 * tt3 * 0.25 + 1.346243430000000e+01;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[26] = 4.291424920000000e+00 * tt5 + -5.501542700000000e-03 * tt0 + 5.994382880000000e-05 * tt1 * 0.50 + -7.084662850000000e-08 * tt2 * OneThird + 2.686857710000000e-11 * tt3 * 0.25 + 2.666823160000000e+00;
    } else {
      s0_R[26] = 1.071881500000000e+00 * tt5 +  2.168526770000000e-02 * tt0 + -1.002560670000000e-05 * tt1 * 0.50 + 2.214120010000000e-09 * tt2 * OneThird + -1.900028900000000e-13 * tt3 * 0.25 + 1.511561070000000e+01;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[27] = 2.251845286753941e+00 * tt5 + 1.765401529121510e-02 * tt0 + -2.372620998746318e-05 * tt1 * 0.50 + 1.727224682516638e-08 * tt2 * OneThird + -5.064957679024219e-12 * tt3 * 0.25 + 1.248990921083339e+01;
    } else {
      s0_R[27] = 5.628172209693201e+00 * tt5 +  4.085394191788521e-03 * tt0 + -1.593486747480799e-06 * tt1 * 0.50 + 2.862686614342911e-10 * tt2 * OneThird + -1.940857878719715e-14 * tt3 * 0.25 + -3.930065989049631e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[28] = 2.135836300000000e+00 * tt5 + 1.811887210000000e-02 * tt0 + -1.739474740000000e-05 * tt1 * 0.50 + 9.343975680000000e-09 * tt2 * OneThird + -2.014576150000000e-12 * tt3 * 0.25 + 1.221564800000000e+01;
    } else {
      s0_R[28] = 4.511297320000000e+00 * tt5 +  9.003597449999999e-03 * tt0 + -4.169396350000000e-06 * tt1 * 0.50 + 9.233458820000000e-10 * tt2 * OneThird + -7.948382010000000e-14 * tt3 * 0.25 + 6.322472050000000e-01;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[29] = 1.242373300000000e+00 * tt5 + 3.107220100000000e-02 * tt0 + -5.086686400000000e-05 * tt1 * 0.50 + 4.313713100000000e-08 * tt2 * OneThird + -1.401459400000000e-11 * tt3 * 0.25 + 1.387431900000000e+01;
    } else {
      s0_R[29] = 5.923829100000000e+00 * tt5 +  6.792360000000000e-03 * tt0 + -2.565856400000000e-06 * tt1 * 0.50 + 4.498784100000000e-10 * tt2 * OneThird + -2.994010100000000e-14 * tt3 * 0.25 + -7.601774200000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[30] = 3.298616281368291e+00 * tt5 + 1.408708904146039e-03 * tt0 + -3.964481129109908e-06 * tt1 * 0.50 + 5.642920880408571e-09 * tt2 * OneThird + -2.445407041148433e-12 * tt3 * 0.25 + 3.950623591964732e+00;
    } else {
      s0_R[30] = 2.926639911210682e+00 * tt5 +  1.487977101178227e-03 * tt0 + -5.684761849244810e-07 * tt1 * 0.50 + 1.009704225872734e-10 * tt2 * OneThird + -6.753354387142974e-15 * tt3 * 0.25 + 5.980528055036107e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[31] = 2.500000000000000e+00 * tt5 + 0.000000000000000e+00 * tt0 + 0.000000000000000e+00 * tt1 * 0.50 + 0.000000000000000e+00 * tt2 * OneThird + 0.000000000000000e+00 * tt3 * 0.25 + 4.366000000000000e+00;
    } else {
      s0_R[31] = 2.500000000000000e+00 * tt5 +  0.000000000000000e+00 * tt0 + 0.000000000000000e+00 * tt1 * 0.50 + 0.000000000000000e+00 * tt2 * OneThird + 0.000000000000000e+00 * tt3 * 0.25 + 4.366000000000000e+00;
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

    keqs[0] = (C0 * g0_RT[3]);
    keqs[1] = (C0 * g0_RT[4]);
    keqs[2] = (g0_RT[1] * g0_RT[4]);
    keqs[3] = (g0_RT[3] * g0_RT[4]);
    keqs[4] = (g0_RT[6] * g0_RT[4]);
    keqs[5] = (g0_RT[1] * g0_RT[14]);
    keqs[6] = (g0_RT[1] * g0_RT[16]);
    keqs[7] = (g0_RT[0] * g0_RT[14]);
    keqs[8] = (g0_RT[1] * g0_RT[16]);
    keqs[9] = (g0_RT[17] * g0_RT[1]);
    keqs[10] = (g0_RT[12] * g0_RT[4]);
    keqs[11] = (C0 * g0_RT[15]);
    keqs[12] = (g0_RT[14] * g0_RT[4]);
    keqs[13] = (g0_RT[1] * g0_RT[15]);
    keqs[14] = (g0_RT[16] * g0_RT[4]);
    keqs[15] = (g0_RT[17] * g0_RT[4]);
    keqs[16] = (g0_RT[17] * g0_RT[4]);
    keqs[17] = (g0_RT[18] * g0_RT[4]);
    keqs[18] = (g0_RT[19] * g0_RT[4]);
    keqs[19] = (g0_RT[9] * g0_RT[14]);
    keqs[20] = (g0_RT[1] * g0_RT[27]);
    keqs[21] = (g0_RT[21] * g0_RT[4]);
    keqs[22] = (g0_RT[10] * g0_RT[14]);
    keqs[23] = (g0_RT[1] * g0_RT[28]);
    keqs[24] = (g0_RT[12] * g0_RT[16]);
    keqs[25] = (g0_RT[17] * g0_RT[12]);
    keqs[26] = (g0_RT[25] * g0_RT[4]);
    keqs[27] = (g0_RT[1] * g0_RT[14] * g0_RT[14]);
    keqs[28] = (g0_RT[27] * g0_RT[4]);
    keqs[29] = (g0_RT[10] * g0_RT[15]);
    keqs[30] = (g0_RT[15] * g0_RT[2]);
    keqs[31] = (g0_RT[6] * g0_RT[16]);
    keqs[32] = (C0 * g0_RT[6]);
    keqs[33] = (C0 * g0_RT[6] * g0_RT[3]);
    keqs[34] = (C0 * g0_RT[5] * g0_RT[6]);
    keqs[35] = (C0 * g0_RT[30] * g0_RT[6]);
    keqs[36] = (C0 * g0_RT[31] * g0_RT[6]);
    keqs[37] = (g0_RT[2] * g0_RT[4]);
    keqs[38] = (C0 * g0_RT[0]);
    keqs[39] = (C0 * g0_RT[0] * g0_RT[0]);
    keqs[40] = (C0 * g0_RT[0] * g0_RT[5]);
    keqs[41] = (C0 * g0_RT[0] * g0_RT[15]);
    keqs[42] = (C0 * g0_RT[5]);
    keqs[43] = (g0_RT[5] * g0_RT[2]);
    keqs[44] = (g0_RT[0] * g0_RT[3]);
    keqs[45] = (g0_RT[4] * g0_RT[4]);
    keqs[46] = (g0_RT[0] * g0_RT[6]);
    keqs[47] = (g0_RT[5] * g0_RT[4]);
    keqs[48] = (g0_RT[0] * g0_RT[8]);
    keqs[49] = (C0 * g0_RT[12]);
    keqs[50] = (g0_RT[0] * g0_RT[9]);
    keqs[51] = (C0 * g0_RT[13]);
    keqs[52] = (g0_RT[0] * g0_RT[12]);
    keqs[53] = (C0 * g0_RT[17]);
    keqs[54] = (g0_RT[0] * g0_RT[14]);
    keqs[55] = (C0 * g0_RT[18]);
    keqs[56] = (C0 * g0_RT[19]);
    keqs[57] = (g0_RT[0] * g0_RT[16]);
    keqs[58] = (C0 * g0_RT[20]);
    keqs[59] = (g0_RT[0] * g0_RT[17]);
    keqs[60] = (g0_RT[12] * g0_RT[4]);
    keqs[61] = (g0_RT[11] * g0_RT[5]);
    keqs[62] = (C0 * g0_RT[20]);
    keqs[63] = (g0_RT[1] * g0_RT[18]);
    keqs[64] = (g0_RT[0] * g0_RT[17]);
    keqs[65] = (g0_RT[12] * g0_RT[4]);
    keqs[66] = (g0_RT[11] * g0_RT[5]);
    keqs[67] = (g0_RT[0] * g0_RT[18]);
    keqs[68] = (g0_RT[0] * g0_RT[19]);
    keqs[69] = (C0 * g0_RT[22]);
    keqs[70] = (C0 * g0_RT[23]);
    keqs[71] = (C0 * g0_RT[24]);
    keqs[72] = (g0_RT[0] * g0_RT[22]);
    keqs[73] = (C0 * g0_RT[25]);
    keqs[74] = (g0_RT[0] * g0_RT[23]);
    keqs[75] = (C0 * g0_RT[26]);
    keqs[76] = (g0_RT[0] * g0_RT[24]);
    keqs[77] = (g0_RT[0] * g0_RT[25]);
    keqs[78] = (g0_RT[11] * g0_RT[14]);
    keqs[79] = (g0_RT[0] * g0_RT[27]);
    keqs[80] = (g0_RT[12] * g0_RT[14]);
    keqs[81] = (g0_RT[1] * g0_RT[28]);
    keqs[82] = (C0 * g0_RT[17]);
    keqs[83] = (g0_RT[1] * g0_RT[5]);
    keqs[84] = (C0 * g0_RT[7]);
    keqs[85] = (g0_RT[5] * g0_RT[2]);
    keqs[86] = (g0_RT[5] * g0_RT[3]);
    keqs[87] = (g0_RT[5] * g0_RT[6]);
    keqs[88] = (g0_RT[5] * g0_RT[6]);
    keqs[89] = (g0_RT[1] * g0_RT[14]);
    keqs[90] = (g0_RT[1] * g0_RT[16]);
    keqs[91] = (g0_RT[17] * g0_RT[1]);
    keqs[92] = (g0_RT[5] * g0_RT[9]);
    keqs[93] = (g0_RT[17] * g0_RT[1]);
    keqs[94] = (C0 * g0_RT[20]);
    keqs[95] = (g0_RT[10] * g0_RT[5]);
    keqs[96] = (g0_RT[11] * g0_RT[5]);
    keqs[97] = (g0_RT[5] * g0_RT[12]);
    keqs[98] = (g0_RT[1] * g0_RT[15]);
    keqs[99] = (g0_RT[5] * g0_RT[14]);
    keqs[100] = (g0_RT[5] * g0_RT[16]);
    keqs[101] = (g0_RT[17] * g0_RT[5]);
    keqs[102] = (g0_RT[17] * g0_RT[5]);
    keqs[103] = (g0_RT[18] * g0_RT[5]);
    keqs[104] = (g0_RT[5] * g0_RT[19]);
    keqs[105] = (g0_RT[1] * g0_RT[27]);
    keqs[106] = (g0_RT[1] * g0_RT[28]);
    keqs[107] = (g0_RT[1] * g0_RT[29]);
    keqs[108] = (g0_RT[21] * g0_RT[5]);
    keqs[109] = (g0_RT[12] * g0_RT[14]);
    keqs[110] = (g0_RT[5] * g0_RT[22]);
    keqs[111] = (g0_RT[5] * g0_RT[23]);
    keqs[112] = (g0_RT[25] * g0_RT[5]);
    keqs[113] = (g0_RT[5] * g0_RT[27]);
    keqs[114] = (g0_RT[3] * g0_RT[7]);
    keqs[115] = (g0_RT[3] * g0_RT[7]);
    keqs[116] = (g0_RT[17] * g0_RT[4]);
    keqs[117] = (g0_RT[13] * g0_RT[3]);
    keqs[118] = (g0_RT[19] * g0_RT[4]);
    keqs[119] = (g0_RT[15] * g0_RT[4]);
    keqs[120] = (g0_RT[16] * g0_RT[7]);
    keqs[121] = (g0_RT[14] * g0_RT[2]);
    keqs[122] = (g0_RT[1] * g0_RT[21]);
    keqs[123] = (g0_RT[1] * g0_RT[22]);
    keqs[124] = (g0_RT[16] * g0_RT[2]);
    keqs[125] = (g0_RT[1] * g0_RT[10]);
    keqs[126] = (g0_RT[17] * g0_RT[1]);
    keqs[127] = (g0_RT[1] * g0_RT[22]);
    keqs[128] = (g0_RT[1] * g0_RT[23]);
    keqs[129] = (g0_RT[1] * g0_RT[24]);
    keqs[130] = (C0 * g0_RT[27]);
    keqs[131] = (g0_RT[14] * g0_RT[16]);
    keqs[132] = (g0_RT[1] * g0_RT[28]);
    keqs[133] = (g0_RT[14] * g0_RT[22]);
    keqs[134] = (g0_RT[16] * g0_RT[4]);
    keqs[135] = (g0_RT[1] * g0_RT[12]);
    keqs[136] = (g0_RT[0] * g0_RT[22]);
    keqs[137] = (g0_RT[1] * g0_RT[24]);
    keqs[138] = (g0_RT[12] * g0_RT[12]);
    keqs[139] = (C0 * g0_RT[28]);
    keqs[140] = (g0_RT[14] * g0_RT[23]);
    keqs[141] = (g0_RT[10] * g0_RT[30]);
    keqs[142] = (g0_RT[10] * g0_RT[31]);
    keqs[143] = (g0_RT[1] * g0_RT[14] * g0_RT[4]);
    keqs[144] = (g0_RT[5] * g0_RT[14]);
    keqs[145] = (g0_RT[1] * g0_RT[12]);
    keqs[146] = (C0 * g0_RT[20]);
    keqs[147] = (g0_RT[10] * g0_RT[5]);
    keqs[148] = (g0_RT[1] * g0_RT[24]);
    keqs[149] = (g0_RT[12] * g0_RT[12]);
    keqs[150] = (g0_RT[10] * g0_RT[14]);
    keqs[151] = (g0_RT[10] * g0_RT[15]);
    keqs[152] = (g0_RT[17] * g0_RT[14]);
    keqs[153] = (g0_RT[25] * g0_RT[12]);
    keqs[154] = (g0_RT[19] * g0_RT[2]);
    keqs[155] = (g0_RT[17] * g0_RT[4]);
    keqs[156] = (g0_RT[13] * g0_RT[6]);
    keqs[157] = (C0 * g0_RT[26]);
    keqs[158] = (g0_RT[1] * g0_RT[25]);
    keqs[159] = (g0_RT[14] * g0_RT[13]);
    keqs[160] = (g0_RT[13] * g0_RT[16]);
    keqs[161] = (g0_RT[18] * g0_RT[13]);
    keqs[162] = (g0_RT[19] * g0_RT[13]);
    keqs[163] = (g0_RT[13] * g0_RT[23]);
    keqs[164] = (g0_RT[25] * g0_RT[13]);
    keqs[165] = (g0_RT[1] * g0_RT[5] * g0_RT[14]);
    keqs[166] = (g0_RT[1] * g0_RT[14]);
    keqs[167] = (g0_RT[14] * g0_RT[6]);
    keqs[168] = (g0_RT[17] * g0_RT[6]);
    keqs[169] = (g0_RT[17] * g0_RT[6]);
    keqs[170] = (g0_RT[14] * g0_RT[16]);
    keqs[171] = (g0_RT[1] * g0_RT[22]);
    keqs[172] = (g0_RT[17] * g0_RT[16]);
    keqs[173] = (g0_RT[0] * g0_RT[22]);
    keqs[174] = (g0_RT[24] * g0_RT[6]);
    keqs[175] = (g0_RT[14] * g0_RT[14] * g0_RT[4]);
    keqs[176] = (g0_RT[14] * g0_RT[14] * g0_RT[22]);

    keqs[0] /= (g0_RT[2] * g0_RT[2]);
    keqs[1] /= (g0_RT[1] * g0_RT[2]);
    keqs[2] /= (g0_RT[0] * g0_RT[2]);
    keqs[3] /= (g0_RT[6] * g0_RT[2]);
    keqs[4] /= (g0_RT[7] * g0_RT[2]);
    keqs[5] /= (g0_RT[9] * g0_RT[2]);
    keqs[6] /= (g0_RT[10] * g0_RT[2]);
    keqs[7] /= (g0_RT[11] * g0_RT[2]);
    keqs[8] /= (g0_RT[11] * g0_RT[2]);
    keqs[9] /= (g0_RT[12] * g0_RT[2]);
    keqs[10] /= (g0_RT[13] * g0_RT[2]);
    keqs[11] /= (g0_RT[14] * g0_RT[2]);
    keqs[12] /= (g0_RT[16] * g0_RT[2]);
    keqs[13] /= (g0_RT[16] * g0_RT[2]);
    keqs[14] /= (g0_RT[17] * g0_RT[2]);
    keqs[15] /= (g0_RT[18] * g0_RT[2]);
    keqs[16] /= (g0_RT[19] * g0_RT[2]);
    keqs[17] /= (g0_RT[20] * g0_RT[2]);
    keqs[18] /= (g0_RT[20] * g0_RT[2]);
    keqs[19] /= (g0_RT[21] * g0_RT[2]);
    keqs[20] /= (g0_RT[22] * g0_RT[2]);
    keqs[21] /= (g0_RT[22] * g0_RT[2]);
    keqs[22] /= (g0_RT[22] * g0_RT[2]);
    keqs[23] /= (g0_RT[23] * g0_RT[2]);
    keqs[24] /= (g0_RT[24] * g0_RT[2]);
    keqs[25] /= (g0_RT[25] * g0_RT[2]);
    keqs[26] /= (g0_RT[26] * g0_RT[2]);
    keqs[27] /= (C0 * g0_RT[27] * g0_RT[2]);
    keqs[28] /= (g0_RT[28] * g0_RT[2]);
    keqs[29] /= (g0_RT[28] * g0_RT[2]);
    keqs[30] /= (g0_RT[14] * g0_RT[3]);
    keqs[31] /= (g0_RT[17] * g0_RT[3]);
    keqs[32] /= (g0_RT[1] * g0_RT[3]);
    keqs[33] /= (g0_RT[1] * g0_RT[3] * g0_RT[3]);
    keqs[34] /= (g0_RT[1] * g0_RT[5] * g0_RT[3]);
    keqs[35] /= (g0_RT[1] * g0_RT[30] * g0_RT[3]);
    keqs[36] /= (g0_RT[1] * g0_RT[31] * g0_RT[3]);
    keqs[37] /= (g0_RT[1] * g0_RT[3]);
    keqs[38] /= (g0_RT[1] * g0_RT[1]);
    keqs[39] /= (g0_RT[0] * g0_RT[1] * g0_RT[1]);
    keqs[40] /= (g0_RT[1] * g0_RT[1] * g0_RT[5]);
    keqs[41] /= (g0_RT[1] * g0_RT[1] * g0_RT[15]);
    keqs[42] /= (g0_RT[1] * g0_RT[4]);
    keqs[43] /= (g0_RT[1] * g0_RT[6]);
    keqs[44] /= (g0_RT[1] * g0_RT[6]);
    keqs[45] /= (g0_RT[1] * g0_RT[6]);
    keqs[46] /= (g0_RT[1] * g0_RT[7]);
    keqs[47] /= (g0_RT[1] * g0_RT[7]);
    keqs[48] /= (g0_RT[1] * g0_RT[9]);
    keqs[49] /= (g0_RT[1] * g0_RT[10]);
    keqs[50] /= (g0_RT[1] * g0_RT[11]);
    keqs[51] /= (g0_RT[1] * g0_RT[12]);
    keqs[52] /= (g0_RT[1] * g0_RT[13]);
    keqs[53] /= (g0_RT[1] * g0_RT[16]);
    keqs[54] /= (g0_RT[1] * g0_RT[16]);
    keqs[55] /= (g0_RT[17] * g0_RT[1]);
    keqs[56] /= (g0_RT[17] * g0_RT[1]);
    keqs[57] /= (g0_RT[17] * g0_RT[1]);
    keqs[58] /= (g0_RT[1] * g0_RT[18]);
    keqs[59] /= (g0_RT[1] * g0_RT[18]);
    keqs[60] /= (g0_RT[1] * g0_RT[18]);
    keqs[61] /= (g0_RT[1] * g0_RT[18]);
    keqs[62] /= (g0_RT[1] * g0_RT[19]);
    keqs[63] /= (g0_RT[1] * g0_RT[19]);
    keqs[64] /= (g0_RT[1] * g0_RT[19]);
    keqs[65] /= (g0_RT[1] * g0_RT[19]);
    keqs[66] /= (g0_RT[1] * g0_RT[19]);
    keqs[67] /= (g0_RT[20] * g0_RT[1]);
    keqs[68] /= (g0_RT[20] * g0_RT[1]);
    keqs[69] /= (g0_RT[1] * g0_RT[21]);
    keqs[70] /= (g0_RT[1] * g0_RT[22]);
    keqs[71] /= (g0_RT[1] * g0_RT[23]);
    keqs[72] /= (g0_RT[1] * g0_RT[23]);
    keqs[73] /= (g0_RT[1] * g0_RT[24]);
    keqs[74] /= (g0_RT[1] * g0_RT[24]);
    keqs[75] /= (g0_RT[1] * g0_RT[25]);
    keqs[76] /= (g0_RT[1] * g0_RT[25]);
    keqs[77] /= (g0_RT[1] * g0_RT[26]);
    keqs[78] /= (g0_RT[1] * g0_RT[27]);
    keqs[79] /= (g0_RT[1] * g0_RT[28]);
    keqs[80] /= (g0_RT[1] * g0_RT[28]);
    keqs[81] /= (g0_RT[1] * g0_RT[29]);
    keqs[82] /= (g0_RT[0] * g0_RT[14]);
    keqs[83] /= (g0_RT[0] * g0_RT[4]);
    keqs[84] /= (g0_RT[4] * g0_RT[4]);
    keqs[85] /= (g0_RT[4] * g0_RT[4]);
    keqs[86] /= (g0_RT[6] * g0_RT[4]);
    keqs[87] /= (g0_RT[7] * g0_RT[4]);
    keqs[88] /= (g0_RT[7] * g0_RT[4]);
    keqs[89] /= (g0_RT[8] * g0_RT[4]);
    keqs[90] /= (g0_RT[9] * g0_RT[4]);
    keqs[91] /= (g0_RT[10] * g0_RT[4]);
    keqs[92] /= (g0_RT[10] * g0_RT[4]);
    keqs[93] /= (g0_RT[11] * g0_RT[4]);
    keqs[94] /= (g0_RT[12] * g0_RT[4]);
    keqs[95] /= (g0_RT[12] * g0_RT[4]);
    keqs[96] /= (g0_RT[12] * g0_RT[4]);
    keqs[97] /= (g0_RT[13] * g0_RT[4]);
    keqs[98] /= (g0_RT[14] * g0_RT[4]);
    keqs[99] /= (g0_RT[16] * g0_RT[4]);
    keqs[100] /= (g0_RT[17] * g0_RT[4]);
    keqs[101] /= (g0_RT[18] * g0_RT[4]);
    keqs[102] /= (g0_RT[19] * g0_RT[4]);
    keqs[103] /= (g0_RT[20] * g0_RT[4]);
    keqs[104] /= (g0_RT[20] * g0_RT[4]);
    keqs[105] /= (g0_RT[21] * g0_RT[4]);
    keqs[106] /= (g0_RT[22] * g0_RT[4]);
    keqs[107] /= (g0_RT[22] * g0_RT[4]);
    keqs[108] /= (g0_RT[22] * g0_RT[4]);
    keqs[109] /= (g0_RT[22] * g0_RT[4]);
    keqs[110] /= (g0_RT[23] * g0_RT[4]);
    keqs[111] /= (g0_RT[24] * g0_RT[4]);
    keqs[112] /= (g0_RT[26] * g0_RT[4]);
    keqs[113] /= (g0_RT[28] * g0_RT[4]);
    keqs[114] /= (g0_RT[6] * g0_RT[6]);
    keqs[115] /= (g0_RT[6] * g0_RT[6]);
    keqs[116] /= (g0_RT[10] * g0_RT[6]);
    keqs[117] /= (g0_RT[12] * g0_RT[6]);
    keqs[118] /= (g0_RT[12] * g0_RT[6]);
    keqs[119] /= (g0_RT[14] * g0_RT[6]);
    keqs[120] /= (g0_RT[17] * g0_RT[6]);
    keqs[121] /= (g0_RT[8] * g0_RT[3]);
    keqs[122] /= (g0_RT[8] * g0_RT[10]);
    keqs[123] /= (g0_RT[8] * g0_RT[12]);
    keqs[124] /= (g0_RT[9] * g0_RT[3]);
    keqs[125] /= (g0_RT[0] * g0_RT[9]);
    keqs[126] /= (g0_RT[5] * g0_RT[9]);
    keqs[127] /= (g0_RT[10] * g0_RT[9]);
    keqs[128] /= (g0_RT[12] * g0_RT[9]);
    keqs[129] /= (g0_RT[9] * g0_RT[13]);
    keqs[130] /= (g0_RT[9] * g0_RT[14]);
    keqs[131] /= (g0_RT[9] * g0_RT[15]);
    keqs[132] /= (g0_RT[17] * g0_RT[9]);
    keqs[133] /= (g0_RT[9] * g0_RT[27]);
    keqs[134] /= (g0_RT[10] * g0_RT[3]);
    keqs[135] /= (g0_RT[0] * g0_RT[10]);
    keqs[136] /= (g0_RT[10] * g0_RT[10]);
    keqs[137] /= (g0_RT[10] * g0_RT[12]);
    keqs[138] /= (g0_RT[10] * g0_RT[13]);
    keqs[139] /= (g0_RT[10] * g0_RT[14]);
    keqs[140] /= (g0_RT[10] * g0_RT[27]);
    keqs[141] /= (g0_RT[11] * g0_RT[30]);
    keqs[142] /= (g0_RT[11] * g0_RT[31]);
    keqs[143] /= (C0 * g0_RT[11] * g0_RT[3]);
    keqs[144] /= (g0_RT[11] * g0_RT[3]);
    keqs[145] /= (g0_RT[0] * g0_RT[11]);
    keqs[146] /= (g0_RT[11] * g0_RT[5]);
    keqs[147] /= (g0_RT[11] * g0_RT[5]);
    keqs[148] /= (g0_RT[11] * g0_RT[12]);
    keqs[149] /= (g0_RT[11] * g0_RT[13]);
    keqs[150] /= (g0_RT[11] * g0_RT[14]);
    keqs[151] /= (g0_RT[11] * g0_RT[15]);
    keqs[152] /= (g0_RT[11] * g0_RT[15]);
    keqs[153] /= (g0_RT[11] * g0_RT[26]);
    keqs[154] /= (g0_RT[12] * g0_RT[3]);
    keqs[155] /= (g0_RT[12] * g0_RT[3]);
    keqs[156] /= (g0_RT[12] * g0_RT[7]);
    keqs[157] /= (g0_RT[12] * g0_RT[12]);
    keqs[158] /= (g0_RT[12] * g0_RT[12]);
    keqs[159] /= (g0_RT[12] * g0_RT[16]);
    keqs[160] /= (g0_RT[17] * g0_RT[12]);
    keqs[161] /= (g0_RT[20] * g0_RT[12]);
    keqs[162] /= (g0_RT[20] * g0_RT[12]);
    keqs[163] /= (g0_RT[12] * g0_RT[24]);
    keqs[164] /= (g0_RT[26] * g0_RT[12]);
    keqs[165] /= (C0 * g0_RT[5] * g0_RT[16]);
    keqs[166] /= (C0 * g0_RT[16]);
    keqs[167] /= (g0_RT[16] * g0_RT[3]);
    keqs[168] /= (g0_RT[18] * g0_RT[3]);
    keqs[169] /= (g0_RT[19] * g0_RT[3]);
    keqs[170] /= (g0_RT[21] * g0_RT[3]);
    keqs[171] /= (g0_RT[0] * g0_RT[21]);
    keqs[172] /= (g0_RT[23] * g0_RT[3]);
    keqs[173] /= (C0 * g0_RT[24]);
    keqs[174] /= (g0_RT[25] * g0_RT[3]);
    keqs[175] /= (C0 * g0_RT[27] * g0_RT[3]);
    keqs[176] /= (C0 * g0_RT[27] * g0_RT[27]);

  };

}

#endif
