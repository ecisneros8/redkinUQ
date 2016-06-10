#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include "Integrator.h"

#define Ith(v,i)    NV_Ith_S(v,i-1)       
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

namespace mech
{

  class CVodeInt : public Integrator
  {
  public:
    
  CVodeInt(IdealGasMix& gas) :
      m_type(0),
      m_method(CV_BDF),
      m_iter(CV_NEWTON)
      {
	m_gas = &gas;
	m_neq = m_gas->nSpecies();
      };

    int integratorType() { return m_type; }

    void initialize(int& vio, int& ver, double& t0, double& tf, double& dt, 
		    double& rtol, vector<double>& atol, vector<double>& z0,
		    vector<string>& trackSpecies);
    
    void integrate(double& RTM, vector<double>& QOI);

    static int cvode_rhs(realtype t, N_Vector z, N_Vector zdot, void* f_data);

    static int cvode_jac(long int N, realtype t,
			 N_Vector z, N_Vector fz, DlsMat J, void *user_data,
			 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    
  private:
    void*       m_cvode_mem;
    void*       m_cvode_data;
    int         m_cvode_flag;
    int         m_vio;
    int         m_ver;
    int         m_type;
    int         m_method;
    int         m_iter;
    int         m_neq;
    string      m_filename;
    double      m_t0;
    double      m_tf;
    double      m_dt;
    double      m_rtol;
    N_Vector    m_atol;
    N_Vector    m_statevec;
    vector<int> m_trackSpecies;
    
  };

  void CVodeInt::initialize(int& vio, int& ver, double& t0, double& tf, double& dt, 
			    double& rtol, vector<double>& atol, vector<double>& z0,
			    vector<string>& trackSpecies) {

    // filename
    int Ti = int(m_gas->m_temp);
    m_filename = "outs/full/gas.1500.dat";
    
    m_cvode_mem  = NULL;
    m_cvode_data = NULL;
    m_atol       = NULL;
    m_statevec   = NULL;

    m_statevec = N_VNew_Serial(m_neq);
    m_atol     = N_VNew_Serial(m_neq);

    m_vio      = vio;
    m_ver      = ver;
    m_t0       = t0;
    m_tf       = tf;
    m_dt       = dt;
    m_rtol     = rtol;
    for(int i = 0; i < m_neq; ++i) { Ith(m_statevec, i+1) = z0[i];   }
    for(int i = 0; i < m_neq; ++i) { Ith(m_atol, i+1)     = atol[i]; }
    for(int k = 0; k < trackSpecies.size(); ++k) {
      m_trackSpecies.push_back(m_gas->m_speciesIndex[trackSpecies[k]]);
    }

    m_cvode_mem  = CVodeCreate(m_method, m_iter);
    m_cvode_flag = CVodeInit(m_cvode_mem, cvode_rhs, t0, m_statevec);
    m_cvode_flag = CVodeSVtolerances(m_cvode_mem, m_rtol, m_atol);
    m_cvode_flag = CVDense(m_cvode_mem, m_neq);
    m_cvode_flag = CVDlsSetDenseJacFn(m_cvode_mem, cvode_jac);
    
  };

  void CVodeInt::integrate(double& RTM, vector<double>& QoI) {

    int            i = 0;
    realtype       t;
    ofstream       out;
    double         tout = m_t0 + m_dt;
    double*        local_z;
    double         Temp;
    double         allz;
    vector<double> z(m_kk,0.0);

    if(m_vio) {

      local_z = NV_DATA_S(m_statevec);
      for(int k = 0; k < m_kk; ++k) { z[k]  = local_z[k]; }
      
      out.open(m_filename.c_str());
      out.precision(8);
      out.setf(ios::scientific);
      out << m_t0
	  << "\t" << m_gas->m_temp;      
      for(int k = 0; k < m_trackSpecies.size(); ++k) {
	out << "\t" << z[m_trackSpecies[k]];
      }      
      out << endl;
      
    }

    while(1) {
      m_cvode_flag = CVode(m_cvode_mem, tout, m_statevec, &t, CV_NORMAL);
      if(m_cvode_flag == 0) {
	tout += m_dt;
	i    += 1;

	// output temperature
	local_z = NV_DATA_S(m_statevec);
	allz    = 0.0;
	for(int k = 0; k < m_kk; ++k) { z[k]  = local_z[k]; } 
	m_gas->getTemperature(z, Temp);
	
	if(m_ver == 1) {	  
	  cout << t
	       << "\t" << Temp;
	  for(int k = 0; k < m_trackSpecies.size(); ++k) {
	    cout << "\t" << z[m_trackSpecies[k]];
	  }
	  cout << endl;
	}
	
	if(m_vio == 1 && i%1000 == 0) {
	  out << t
	      << "\t" << Temp;
	  for(int k = 0; k < m_trackSpecies.size(); ++k) {
	    out << "\t" << z[m_trackSpecies[k]];
	  }
	  out << endl;
	}
	
      };
      if(tout >= m_tf) break;
    }

    N_VDestroy_Serial(m_statevec);
    N_VDestroy_Serial(m_atol);
    CVodeFree(&m_cvode_mem);

  };

  int CVodeInt::cvode_rhs(realtype t, N_Vector z, N_Vector zdot, void* f_data) {

    double*        local_z    = NV_DATA_S(z);
    double*        local_zdot = NV_DATA_S(zdot);
    vector<double> zz(m_kk, 0.0);
    vector<double> fz(m_kk, 0.0);

    for(int k = 0; k < m_kk; ++k) { zz[k] = local_z[k]; }    
    m_gas->getSource(zz, fz);
    copy(fz.begin(), fz.end(), local_zdot);

    return(0);
    
  };

  int CVodeInt::cvode_jac(long int N, realtype t,
			  N_Vector z, N_Vector fz, DlsMat J, void *user_data,
			  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

    double*        local_z = NV_DATA_S(z);
    vector<double> zz(m_kk, 0.0);
    vector<double> dfz(m_kk*m_kk,0.0);

    for(int k = 0; k < m_kk; ++k) { zz[k] = local_z[k]; };
    m_gas->getJacobian(zz, dfz);

    for(int i = 0; i < m_kk; ++i) {
      for(int j = 0; j < m_kk; ++j) {
	IJth(J,i+1,j+1) = dfz[m_kk * i + j];
      }
    }

    return(0);

  };
  
}
