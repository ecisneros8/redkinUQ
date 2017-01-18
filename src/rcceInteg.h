#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include <Eigen/Dense>
#include "cantera/IdealGasMix.h"
#include "ceqInterface.h"

#define Ith(v,i)    NV_Ith_S(v,i-1)       
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

using namespace Eigen;

class rcceInteg
{
 public:

  rcceInteg() { };
    
  Cantera::IdealGasMix& gas() { return *m_gas; }
    
  void setIdealGasMix(Cantera::IdealGasMix& gas);
  
  void initialize(int& vio, int& ver, double& t0, double& tf, double& dt);
  
  void setInitialConditions(std::string& fuel, double& p, double& T0,
			    double& phi, std::vector<int>& csi);
  
  void integrate(int& flag, double& mout);
  
  double modelOutput(std::vector<double>& all_T);
  
  static int rcce_rhs(realtype t, N_Vector z, N_Vector zdot, void* f_data);
  
  void rcce_jac();

  double                       m_Teq;      // equilibrium temperature
  
 private:
  std::string                  m_filename;    
  void*                        m_cvode_mem;
  void*                        m_cvode_data;
  int                          m_cvode_flag;
  int                          m_vio;
  int                          m_ver;    
  int                          m_method;
  int                          m_flag;     // realizability of constraints
  int                          m_ncs;      // # constraint species
  int                          m_ng;       // # linear constraint
  int                          m_mm;       // # elements
  static int                   m_nc;       // # constraints
  static int                   m_kk;       // # species 
  int                          m_iter;
  int                          m_neq;
  double                       m_t0;       // init time
  double                       m_tf;       // final time
  double                       m_dt;       // time step
  static double                m_T;        // temperature
  double                       m_T0;       // initial temperature
  double                       m_rtol;
  N_Vector                     m_atol;
  N_Vector                     m_statevec;
  static double                m_p0;
  static double                m_HR;       // molar enthalpy over gas constant
  VectorXd                     m_z;
  VectorXd                     m_r;
  static VectorXd              m_fz;       // chemistry source
  static VectorXd              m_R;        // Source term
  MatrixXd                     m_Emat;     // elemental constraints
  MatrixXd                     m_CS;       // species constraints
  MatrixXd                     m_BG;       // general constraints
  static MatrixXd              m_C;        // constraints matrix
  static Cantera::IdealGasMix* m_gas;
  
};
  
int                   rcceInteg::m_nc;
int                   rcceInteg::m_kk;
double                rcceInteg::m_T;
double                rcceInteg::m_p0;
double                rcceInteg::m_HR;
VectorXd              rcceInteg::m_fz;
VectorXd              rcceInteg::m_R;
MatrixXd              rcceInteg::m_C;
Cantera::IdealGasMix* rcceInteg::m_gas;

void rcceInteg::setIdealGasMix(Cantera::IdealGasMix& gas) {

  m_gas = &gas;
  m_kk  = m_gas->nSpecies();
  m_mm  = m_gas->nElements();

}

void rcceInteg::initialize(int& vio, int& ver, double& t0, double& tf, double& dt) {
      
  // Integrator options
  m_vio = vio;
  m_ver = ver;
  m_t0  = t0;
  m_tf  = tf;
  m_dt  = dt;
  m_method = CV_BDF;
  m_iter   = CV_NEWTON;
    
};

void rcceInteg::setInitialConditions(std::string& fuel, double& p, double& T0,
				     double& phi, std::vector<int>& csi) {

  /********************************************/
  /* set thermo state of the mixture          */
  /********************************************/
  int     XF;
  int     XO;
  int     XN;
  double  nux;
  double  oar = 0.21;
  double* x0  = new double[m_kk];

  if(fuel == "H2O") {
    /* nux */
    nux = 0.5;
    /* H2 */
    XF = 0;
    /* O2 */
    XO = 2;
    /* N2 */
    XN = 8;
  } else if(fuel == "CH4") {
    /* nux */
    nux = 2.0;
    /* CH4 */
    XF = 13;
    /* O2 */
    XO = 3;
    /* N2 */
    XN = 47;
  }

  for(int k = 0; k < m_kk; ++k) { x0[k] = 0.0; }
  
  x0[XF] = oar * phi / (nux + oar * phi);
  x0[XO] = nux * x0[XF] / phi;
  x0[XN] = 1.0 - (1.0 + nux/phi) * x0[XF];

  m_gas->setState_TPX(T0, p, x0);
  m_gas->equilibrate("HP", "auto");
  m_Teq = m_gas->temperature();
  m_gas->setState_TPX(T0, p, x0);

  m_p0 = p;

  /********************************************/
  /* Fill constraints matrix                  */
  /********************************************/
  m_ncs = csi.size();
  m_ng  = 0;
  m_nc  = m_mm + m_ncs + m_ng;
  m_neq = m_nc;

  m_r.setZero(m_nc);
  m_z.setZero(m_kk);
  m_fz.setZero(m_kk);
  m_R.setZero(m_nc);

  m_Emat.setZero(m_kk,m_mm);
  m_C.setZero(m_kk,m_nc);
  m_CS.setZero(m_kk,m_ncs);
  if(m_ng == 0) {
    m_BG.setZero(m_kk,1);
  } else {
    m_BG.setZero(m_kk,m_ng);
  }

  for(int k = 0; k < m_kk; ++k) {
    for(int m = 0; m < m_mm; ++m) {
      m_Emat(k,m) = m_gas->nAtoms(k,m);
    }
  }
    
  for(int i = 0; i < m_ncs; ++i) {
    int k     = csi[i];
    m_CS(k,i) = 1.0;
  }
    
  if(m_ng == 0) {
    m_C << m_CS, m_Emat;
  } else {
    for(int k = 0; k < m_kk; ++k) { m_BG(k,0) = 1.0; }
    m_C << m_CS, m_Emat, m_BG;
  }
  
  /********************************************/
  /* Initialize ceq                           */
  /********************************************/
  double  mmw  = 0.0;
  double  patm = m_gas->pressure() / Cantera::OneAtm;
  double* ther = new double[m_kk * 15];
  double* Bg   = new double[m_kk * m_ng];
  double* Em   = new double[m_kk * m_mm];
  double* hRT  = new double[m_kk];
  double* wts  = new double[m_kk];
  double* z0   = new double[m_kk];
  double* zc   = new double[m_kk];
  double* coef = new double[15];
    
  m_gas->getMolecularWeights(wts);
  m_gas->getMoleFractions(x0);
  m_gas->getEnthalpy_RT(hRT);
  for(int k = 0; k < m_kk; ++k) { mmw += x0[k] * wts[k]; }
    
  mmw  = 1.0 / mmw;
  m_HR = 0.0;
  for(int k = 0; k < m_kk; ++k) { z0[k]  = x0[k] * mmw; }
  for(int k = 0; k < m_kk; ++k) { m_HR  += z0[k] * hRT[k] * T0; }
  
  // these loops put together the element matrix
  for(int k = 0; k < m_kk; ++k) {
    for(int m = 0; m < m_mm; ++m) {
      Em[m*m_kk+k] = m_Emat(k,m);
    }
  }
    
  // these loops put together the general
  // constraint matrix
  for(int k = 0; k < m_kk; ++k) {
    for(int i = 0; i < m_ng; ++i) {
      Bg[i*m_kk+k] = m_BG(k,i);
    }
  }

  // get the termo parameterization from cantera
  for(int k = 0; k < m_kk; ++k) {
    int    type;
    double minTemp;
    double maxTemp;
    double refPressure;
    m_gas->thermo().speciesThermo().reportParams(k,type,coef,minTemp,maxTemp,refPressure);
    for(int i = 0; i < 15; ++i) {
      ther[i*m_kk+k] = coef[i];
    }
  }

  // now that everything ready, initialize
  ceqInitialize(m_ncs,m_mm,m_ng,m_kk,&csi[0],T0,patm,m_HR,wts,x0,ther,Em,Bg,m_T,zc);

  // set temperature and reduced representation
  m_T0  = m_T;
  m_z   = VectorXd::Map(zc, m_kk);
  m_r   = m_C.transpose() * m_z;

  /********************************************/
  /* check whether constraints are realizable */
  /********************************************/
  if(m_T < 0.0) {
    m_flag = 1; // not realizable
  } else {
    m_flag = 0; // realizable
  }

  /********************************************/
  /* Initialize CVODE                         */
  /********************************************/
  double atol;
  if(fuel == "H2O") {
    atol = 1.0e-10;
  } else if(fuel == "CH4") {
    atol = 1.0e-08;
  }

  m_cvode_mem  = NULL;
  m_cvode_data = NULL;
  m_atol       = NULL;
  m_statevec   = NULL;

  m_statevec = N_VNew_Serial(m_neq);
  m_atol     = N_VNew_Serial(m_neq);

  m_rtol = 1.0e-08;
  for(int i = 0; i < m_neq; ++i) { Ith(m_statevec,i+1) = m_r(i); }
  for(int i = 0; i < m_neq; ++i) { Ith(m_atol,i+1)     = atol; }

  m_cvode_mem  = CVodeCreate(m_method, m_iter);
  m_cvode_flag = CVodeInit(m_cvode_mem, rcce_rhs, m_t0, m_statevec);
  m_cvode_flag = CVodeSVtolerances(m_cvode_mem, m_rtol, m_atol);
  m_cvode_flag = CVDense(m_cvode_mem, m_neq);
  m_cvode_flag = CVDlsSetDenseJacFn(m_cvode_mem, NULL);
  m_cvode_flag = CVodeSetMinStep(m_cvode_mem, 1.0e-10);

  /********************************************/
  /* output files                             */
  /********************************************/
  std::string        temps;
  std::string        ndims;
  std::string        specs;
  int                iT0;
  std::vector<int>   ici;
  std::ostringstream ostrt;
  std::ostringstream ostrd;

  iT0 = int(T0);
  ostrt << iT0;
  ostrd << m_ncs+m_ng;
  temps = ostrt.str();
  ndims = ostrd.str();

  for(int i = 0; i < csi.size(); ++i) { 
    ici.push_back(csi[i]);
    std::ostringstream ostrc;
    ostrc << ici[i] + 1;
    specs = specs+ostrc.str()+"-";
  }
    
  m_filename = "outs/temp/"+ndims+"D/";
  m_filename = m_filename+"rcce";
  m_filename = m_filename+"."+temps+"K";
  m_filename = m_filename+"."+specs+"S";
  m_filename = m_filename+".dat";

};

void rcceInteg::integrate(int& flag, double& mout) {

  std::ofstream  out;
  int            i      = 0;
  int            iprint = (int) (1.0e-6/m_dt);
  realtype       t      = m_t0;
  double         tout   = m_t0 + m_dt;
  std::vector<double> all_T;

  if(m_vio) {
    out.precision(8);
    out.setf(std::ios::scientific);
    out.open(m_filename.c_str());
    out << t << "\t" << m_T << "\t" << m_Teq << std::endl;
  }

  if(m_ver == 1 && i%iprint == 0) {
    std::cout << std::endl;
    std::cout << t << "\t" << m_T << "\t" << m_Teq << std::endl;
  }

  all_T.push_back(m_T);
 
  if(m_flag == 0) {

    while(1) {

      /* integrate */
      m_cvode_flag = CVode(m_cvode_mem, tout, m_statevec, &t, CV_NORMAL);
	
      if(m_cvode_flag == 0) {
	tout += m_dt;
	i    += 1;

	/* save temperature */
	all_T.push_back(m_T);
	  
	/* output */
	if(m_ver == 1 && i%iprint == 0) {
	  std::cout << std::endl;
	  std::cout << t << "\t" << m_T << "\t" << m_Teq << std::endl;
	}
	if(m_vio == 1 && i%iprint == 0) {
	  out << t << "\t" << m_T << "\t" << m_Teq << std::endl;
	}

      } else {

	//std::cout << "Dump CVode! \n" << std::endl;
	mout = -1.0e10;;
	break;

      }
      
      if(tout < m_tf && fabs(m_T-m_Teq) < 0.05) break;

      if(tout >= m_tf) break;

    }

    if(tout >= m_tf) {
      mout = modelOutput(all_T);
    } else {
      mout = -1.0e10;
    }

  } else {

    mout = -1.0e10;
      
  }

  out.close();
  N_VDestroy_Serial(m_statevec);
  N_VDestroy_Serial(m_atol);
  CVodeFree(&m_cvode_mem);

};

double rcceInteg::modelOutput(std::vector<double>& all_T) {

  int                 nt = all_T.size();
  double              w; //  = m_dt / 3.0;
  double              t  = m_t0;
  double              OneSixth = 1.0/6.0;
  double              mout;
  std::vector<double> fnum(nt,0.0);
  std::vector<double> fden(nt,0.0);
  std::string         mode = "trap";

  /* put the functions together */
  for(int i = 0; i < nt; ++i) {
    double T = all_T[i];
    fden[i]  = pow(T - m_T0, 8.0) * pow(m_Teq - T, 8.0);
    fnum[i]  = fden[i] * log(t);
    t += m_dt;
  }

  /* integrate */
  double qnum = 0.0;
  double qden = 0.0;

  if(mode == "trap") {

    /* integrate with trapezoidal rule */
    t = m_t0;
    for(int i = 0; i < nt-1; ++i) {
      w     = log(t + m_dt) - log(t);
      qnum += w * fnum[i+1] + fnum[i];
      qden += w * fden[i+1] + fden[i];
      t    += m_dt;
    }

  } else if(mode == "simp") {

    /* integrate with simpsons rule */
    t = m_t0;
    for(int i = 0; i < nt-1; i+=2) {
      w     = log(t+2*m_dt) - log(t);
      qnum += OneSixth * w * (fnum[i+2] + 4.0 * fnum[i+1] + fnum[i]);
      qden += OneSixth * w * (fden[i+2] + 4.0 * fden[i+1] + fden[i]);
      t    += m_dt;
    }

  }

  /* the ratio between integrals is the ignition time */
  mout = qnum/qden;
  mout = exp(mout);
  mout = mout * 1.0e6;

  return(mout);
    
}

int rcceInteg::rcce_rhs(realtype t, N_Vector r, N_Vector rdot, void* f_data) {

  int            flag;
  double*        local_r    = NV_DATA_S(r);
  double*        local_rdot = NV_DATA_S(rdot);
  double*        wts        = new double[m_kk];
  double*        x          = new double[m_kk];
  double*        omega      = new double[m_kk];
  double*        fz         = new double[m_kk];

  /***********************************/
  /* perform species reconstruction  */
  /***********************************/
  m_gas->getMolecularWeights(wts);
  speciesReconstruction(local_r, m_T, x, flag);

  /***********************************/
  /* set state and get full source   */
  /***********************************/  
  m_gas->setState_TPX(m_T, m_p0, x);
  m_gas->getNetProductionRates(omega);

  double rho  = m_gas->density();
  double irho = 1.0 / rho;
  for(int k = 0; k < m_kk; ++k) { fz[k] = omega[k] * irho; }

  /***********************************/
  /* transfer data to eigen objects  */
  /* to perform MV product           */
  /***********************************/
  m_fz = VectorXd::Map(fz, m_kk);
  m_R  = m_C.transpose() * m_fz;
  
  /***********************************/
  /* dump data in pointer to rdot    */
  /***********************************/
  for(int i = 0; i < m_nc; ++i) { local_rdot[i] = m_R(i); }

  /**/
  if(flag == 1) {
    return(1);
  } else {
    return(0);
  }
  /**/
      
};

void rcceInteg::rcce_jac() {

};

