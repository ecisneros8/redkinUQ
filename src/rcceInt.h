#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>

using namespace std;
using namespace std::chrono;

extern "C" {
  void localceq_mp_ceqinit_(int& ncs, const int& ne, int& ng, const int& ns,
			    int& nc, int* csi, double* Ti, double* pi, 
			    double* hi, double* mw, double* xi, 
			    double* Ein, double* Bg, double* Tad, 
			    double* Teq, double* zeq);

  void localceq_mp_ceqrecon_(int& nc, const int& ns,
			      double* p, double* hr, double* c, double* T, 
			      double* z, int& flag);
}

namespace mech {

  class rcceInt : public Integrator
  {
  public:

  rcceInt(IdealGasMix& gas) :
    m_type(2),
      m_method(CV_BDF),
      m_iter(CV_NEWTON) 
      { m_gas = &gas; };
    
    int integratorType() { return m_type; }

    IdealGasMix& gas() { return *m_gas; }
    
    void initialize(int& vio, int& ver, double& t0, double& tf, double& dt);

    void setInitialConditions(double& p, double& To, vector<int>& csi);

    void integrate(int& flag, double& QoI);

    static int rcce_rhs(realtype t, N_Vector z, N_Vector zdot, void* f_data);

    void rcce_jac();

  private:
    static bool     m_stop;
    string          m_filename;    
    void*           m_cvode_mem;
    void*           m_cvode_data;
    int             m_cvode_flag;
    int             m_vio;
    int             m_ver;
    int             m_type;    
    int             m_method;
    int             m_flag;     // realizability of constraints
    int             m_ncs;      // # constraint species
    int             m_ng;       // # linear constraint
    static int      m_nc;       // # constraints
    int             m_iter;
    int             m_neq;
    double          m_t0;       // init time
    double          m_tf;       // final time
    double          m_dt;       // time step
    static double   m_T;        // temperature
    double          m_Tin;      // initial temperature
    double          m_Tor;      // original (scenario) temperature
    double          m_Teq;      // equilibrium temperature
    double          m_rtol;
    N_Vector        m_atol;
    N_Vector        m_statevec;
    static double   m_HR;       // molar enthalpy over gas constant
    static VectorXd m_z;        // specific moles
    VectorXd        m_c;
    static VectorXd m_fz;       // chemistry source
    static VectorXd m_R;        // Source term
    MatrixXd        m_CS;       // species constraints
    MatrixXd        m_BG;       // general constraints
    static MatrixXd m_C;        // constraints matrix
     
  };

  bool     rcceInt::m_stop;
  int      rcceInt::m_nc;
  double   rcceInt::m_T;
  double   rcceInt::m_HR;
  VectorXd rcceInt::m_z;
  VectorXd rcceInt::m_fz;
  VectorXd rcceInt::m_R;
  MatrixXd rcceInt::m_C;

  void rcceInt::initialize(int& vio, int& ver, double& t0, double& tf, double& dt) {
    
    // Integrator options
    m_vio = vio;
    m_ver = ver;
    m_t0  = t0;
    m_tf  = tf;
    m_dt  = dt;
    
  };

  void rcceInt::setInitialConditions(double& p, double& To, vector<int>& csi) {

    // set pressure
    m_gas->setPressure(p);
    
    // number of constraints
    m_ncs = csi.size();
    m_ng  = 0;
    m_nc  = m_mm + m_ncs + m_ng;
    m_neq = m_nc;

    m_c.setZero(m_nc);
    m_fz.setZero(m_kk);
    m_R.setZero(m_nc);

    // Fill constraints matrix
    m_C.setZero(m_kk,m_nc);
    m_CS.setZero(m_kk,m_ncs);
    if(m_ng == 0) {
      m_BG.setZero(m_kk,1);
    } else {
      m_BG.setZero(m_kk,m_ng);
    }
    
    for(int i = 0; i < m_ncs; ++i) {
      int k     = csi[i];
      m_CS(k,i) = 1;
    }
    
    if(m_ng == 0) {
      m_C << m_CS, m_gas->m_Emat;
    } else {
      for(int k = 0; k < m_kk; ++k) { m_BG(k,0) = 1.0; }
      m_C << m_CS, m_gas->m_Emat, m_BG;
    }
    
    // Initialize CEQ
    double         patm;
    vector<double> Bg(m_kk*m_ng, 0.0);
    vector<double> hi(m_kk, 0.0);
    vector<double> zc(m_kk, 0.0);
    vector<double> mw = m_gas->molecularWeights();
    vector<double> xi = m_gas->moleFractions();
    vector<double> zi = m_gas->specificMoles();
    vector<double> Em = m_gas->elementMatrix();

    patm  = m_gas->m_pres / OneAtm;
    m_Tor = To;
    m_Tin = To;
    m_HR  = 0.0;
    m_gas->getEnthalpies_RT(To, hi);
    for(int k = 0; k < m_kk; ++k) { m_HR += zi[k] * hi[k] * To; }

    for(int k = 0; k < m_kk; ++k) {
      for(int i = 0; i < m_ng; ++i) {
	Bg[i*m_kk+k] = m_BG(k,i);
      }
    }

    localceq_mp_ceqinit_(m_ncs, m_mm, m_ng, m_kk, m_nc, &csi[0],
			 &To, &patm, &m_HR, &mw[0], &xi[0], &Em[0], &Bg[0], 
			 &m_Teq, &m_T, &zc[0]);
 
    m_Tin  = m_T;
    m_z    = VectorXd::Map(&zc[0], zc.size());
    m_c    = m_C.transpose() * m_z;
    
    if(m_T < 0.0) {
      m_flag = 1; // the constraints are not realizable
    } else {
      m_flag = 0; // the constraints are realizable
    }

    // Initialize CVODE
    m_cvode_mem  = NULL;
    m_cvode_data = NULL;
    m_atol       = NULL;
    m_statevec   = NULL;

    m_statevec = N_VNew_Serial(m_neq);
    m_atol     = N_VNew_Serial(m_neq);

    m_rtol = 1.0e-06;
    for(int i = 0; i < m_neq; ++i) { Ith(m_statevec,i+1) = m_c(i); }
    for(int i = 0; i < m_neq; ++i) { Ith(m_atol,i+1)     = 1.0e-08; }

    m_cvode_mem  = CVodeCreate(m_method, m_iter);
    m_cvode_flag = CVodeInit(m_cvode_mem, rcce_rhs, m_t0, m_statevec);
    m_cvode_flag = CVodeSVtolerances(m_cvode_mem, m_rtol, m_atol);
    m_cvode_flag = CVDense(m_cvode_mem, m_neq);
    m_cvode_flag = CVDlsSetDenseJacFn(m_cvode_mem, NULL);
    m_cvode_flag = CVodeSetMinStep(m_cvode_mem, 1.0e-14);

    // Output file name
    string        temps;
    string        ndims;
    string        specs;
    int           iTo;
    vector<int>   ici;
    ostringstream ostrt;
    ostringstream ostrd;
    //ostringstream ostrc; 

    iTo = int(To);
    ostrt << iTo;
    ostrd << m_ncs+m_ng;
    temps = ostrt.str();
    ndims = ostrd.str();

    for(int i = 0; i < csi.size(); ++i) { 
      ici.push_back(csi[i]);
      ostringstream ostrc;
      ostrc << ici[i] + 1;
      specs = specs+ostrc.str()+"-";
    }
    //specs = ostrc.str(); 
    
    m_filename = "outs/"+ndims+"D/";
    m_filename = m_filename+"rcce";
    m_filename = m_filename+"."+temps+"K";
    m_filename = m_filename+"."+specs+"S";
    m_filename = m_filename+".dat";
    
  };

  void rcceInt::integrate(int& flag, double& QoI) {

    ofstream out;
    int      i  = 0;
    int      iprint = (int) (1.0e-6/m_dt);
    realtype t = m_t0;
    double   tout = m_t0 + m_dt;
    double   funt = 0.0;
    double   qnum = 0.0;
    double   qden = 0.0;
    double   qoin = 0.0;
    double   qoid = 0.0;
    double   old_time;
    double   old_temp;
    double   old_funt;
    VectorXd statevec(m_nc);

    if(m_vio) {
      for(int k = 0; k < m_nc; ++k) { statevec(k) = Ith(m_statevec,k+1); }
      out.precision(8);
      out.setf(ios::scientific);
      out.open(m_filename.c_str());
      out << t
	  << "\t" << m_T
	  << "\t" << m_Teq << endl;
    }

    if(m_ver == 1 && i%iprint == 0) {
      cout << endl;
      cout << t
	   << "\t" << m_T
	   << "\t" << m_Teq << endl;
    }
 
    if(m_flag == 0) {

      while(1) {
	
	/* previous soln */
	old_time = t;
	old_temp = m_T;
	old_funt = pow(m_Teq - old_temp, 8.0) * pow(old_temp - m_Tin, 8.0);

	/* integrate */
	m_cvode_flag = CVode(m_cvode_mem, tout, m_statevec, &t, CV_NORMAL);
	
	if(m_cvode_flag == 0) {
	  tout += m_dt;
	  i    += 1;

	  /* qoi */
	  funt  = pow(m_Teq - m_T, 8.0) * pow(m_T - m_Tin, 8.0);
	  qnum  = funt * log(t) + old_funt * log(old_time);
	  qden  = funt + old_funt;
	  qoin += (log(t) - log(old_time)) * qnum;
	  qoid += (log(t) - log(old_time)) * qden;

	  for(int k = 0; k < m_nc; ++k) { statevec(k) = Ith(m_statevec,k+1); }
	  
	  /* output */
	  if(m_ver == 1 && i%iprint == 0) {
	    cout << endl;
	    cout << t
		 << "\t" << m_T
		 << "\t" << m_Teq
		 << "\t" << qoin/qoid << endl;
	  }
	  if(m_vio == 1 && i%iprint == 0) {
	    out << t
		<< "\t" << m_T
		<< "\t" << m_Teq << endl;
	  }

	} else {

	  std::cout << "Dump CVode! \n" << std::endl;
	  qoin = -1.0e10;
	  qoid =  1.0e00;
	  break;

	}

	if(tout >= m_tf) break;

      }

      QoI = qoin/qoid;

    } else {

      QoI = -1.0e10;
      
    }

    out.close();
    N_VDestroy_Serial(m_statevec);
    N_VDestroy_Serial(m_atol);
    CVodeFree(&m_cvode_mem);

  };

  int rcceInt::rcce_rhs(realtype t, N_Vector c, N_Vector cdot, void* f_data) {

    int            flag;
    double         patm       = m_gas->m_pres / OneAtm;
    double*        local_c    = NV_DATA_S(c);
    double*        local_cdot = NV_DATA_S(cdot);
    vector<double> ci(m_nc,0.0);
    vector<double> zi(m_kk,0.0);
    vector<double> fz(m_kk,0.0);

    for(int i = 0; i < m_nc; ++i) { ci[i] = local_c[i]; }
    localceq_mp_ceqrecon_(m_nc, m_kk, &patm, &m_HR, &ci[0], &m_T, &zi[0], flag);

    m_gas->getNetProductionRates(m_T, zi, fz);
    m_fz = VectorXd::Map(&fz[0], fz.size());
    m_R  = m_C.transpose() * m_fz;
    for(int i = 0; i < m_nc; ++i) { local_cdot[i] = m_R(i); }

    /*std::cout << "------------------------------------" << std::endl;
    std::cout << t << "\t" << m_R.transpose() << std::endl;
    std::cout << "------------------------------------" << std::endl;*/

    if(flag == 1) {
      return(1);
    } else {
      return(0);
    }
      
  };

  void rcceInt::rcce_jac() {

  };

  
}
