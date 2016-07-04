#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;

extern "C" {
  void localceq_mp_ceqinit_(int& ncs, const int& ne, int& ng, const int& ns,
			    int& nc, int* csi, double* Ti, double* hi, double* mw,
			    double* xi, double* Ein, double* Tad, double* Teq,
			    double* zeq);

  void localceq_mp_ceqcalc_(int& nc, const int& ns,
			    double* hr, double* c, double* T, double* z);
}

namespace mech {

  class rcceInt : public Integrator
  {
  public:

  rcceInt(IdealGasMix& gas) :
    m_type(2) { m_gas = &gas; };
    
    int integratorType() { return m_type; }

    IdealGasMix& gas() { return *m_gas; }
    
    void initialize(int& vio, int& ver, int& ncs, double& t0, double& tf, double& dt);

    void setInitialConditions(double& To, vector<int>& csi);

    void integrate(double& RTM, vector<double>& QOI);

    void rcce_rhs(double& t, VectorXd& statevec, VectorXd& fz);

    void rcce_jac();

  private:
    string      m_filename;
    int         m_vio;
    int         m_ver;
    int         m_type;     
    int         m_ncs;      // # constraint species
    int         m_ng;       // # linear constraint
    int         m_nc;       // # constraints
    double      m_t0;       // init time
    double      m_tf;       // final time
    double      m_dt;       // time step
    double      m_T;        // temperature
    double      m_Tin;      // initial temperature
    double      m_Tor;      // original (scenario) temperature
    double      m_Teq;      // equilibrium temperature
    double      m_HR;       // molar enthalpy over gas constant
    VectorXd    m_z;        // specific moles
    VectorXd    m_statevec; // state vector
    VectorXd    m_fz;       // chemistry source
    VectorXd    m_R;        // Source term
    MatrixXd    m_CS;       // species constraints
    MatrixXd    m_BG;       // general constraints
    MatrixXd    m_C;        // constraints matrix
    vector<int> m_csi;
    
  };

  void rcceInt::initialize(int& vio, int& ver, double& t0, double& tf, double& dt) {
    
    // Integrator options
    m_vio = vio;
    m_ver = ver;
    m_t0  = t0;
    m_tf  = tf;
    m_dt  = dt;
    
  };

  void rcceInt::setInitialConditions(double& To, vector<int>& csi) {

    m_ncs = 0;
    m_ng  = 0;
    for(int i = 0; i < csi.size(); ++i) {
      if(csi[i] <  8) m_ncs += 1;
      if(csi[i] >= 8) m_ng  += 1;
    }
    m_nc  = m_mm + m_ncs + m_ng;

    m_statevec.setZero(m_nc);
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
      for(int i = m_ncs; i < m_nc; ++i) {
	if(csi[i] == TM) {
	  for(int k = 0; k < m_kk; ++k) { m_BG(k,i-m_ncs) = 1.0; }
	} else if(csi[i] == FO) {
	  m_BG(3,i-m_ncs) = 1.0;
	  m_BG(4,i-m_ncs) = 1.0;
	  m_BG(7,i-m_ncs) = 1.0;
	} else if(csi[i] == AV) {
	  m_BG(1,i-m_ncs) = 1.0;
	  m_BG(3,i-m_ncs) = 2.0;
	  m_BG(4,i-m_ncs) = 1.0;
	}
      }
      m_C << m_CS, m_gas->m_Emat, m_BG;
    }
    
    // Initialize CEQ
    vector<double> Bg(m_kk*m_ng, 0.0);
    vector<double> hi(m_kk, 0.0);
    vector<double> zc(m_kk, 0.0);
    vector<double> mw = m_gas->molecularWeights();
    vector<double> xi = m_gas->moleFractions();
    vector<double> zi = m_gas->specificMoles();
    vector<double> Em = m_gas->elementMatrix();

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
			 &To, &m_HR, &mw[0], &xi[0], &Em[0], &Bg[0], 
			 &m_Teq, &m_T, &zc[0]);
 
    m_Tin      = m_T;
    m_z        = VectorXd::Map(&zc[0], zc.size());
    m_statevec = m_C.transpose() * m_z;

    if(m_T < 0.0) {
      m_flag = 1; // the constraints are not realizable
    } else {
      m_flag = 0; // the constraints are realizable
    }
    
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

  void rcceInt::integrate(double& RTM, vector<double>& QOI) {

    ofstream out;
    int      i  = 0;
    double   tn = m_t0;
    double   funt = 0.0;
    double   qnum = 0.0;
    double   qden = 0.0;
    double   qoin = 0.0;
    double   qoid = 0.0;
    double   qoit = 0.0;
    double   old_time;
    double   old_temp;
    double   old_funt;
    VectorXd  u(m_nc);
    VectorXd k1(m_nc);
    VectorXd k2(m_nc);
    VectorXd k3(m_nc);
    VectorXd k4(m_nc);
    
    if(m_vio) {
      out.precision(8);
      out.setf(ios::scientific);
      out.open(m_filename.c_str());
      out << tn    << "\t"
	  << m_T   << "\t" 
	  << m_Teq << "\t"
	  << m_Tin << "\t"
	  << funt  << "\t"
	  << qoin  << "\t"
	  << qoid << endl;
    }

    high_resolution_clock::time_point start_clock = high_resolution_clock::now();
    
    while(1) {

      /* previous soln */
      old_time = tn;
      old_temp = m_T;
      old_funt = pow(m_Teq - old_temp, 3.0) * pow(old_temp - m_Tin, 2.0);

      /* 4th Order Runge-Kutta */
      u  = m_statevec;
      rcce_rhs(tn, u, k1);
      k1 = m_dt * k1;
      
      u  = m_statevec + 0.5 * k1;
      rcce_rhs(tn, u, k2);
      k2 = m_dt * k2;

      u  = m_statevec + 0.5 * k2;
      rcce_rhs(tn, u, k3);
      k3 = m_dt * k3;

      u  = m_statevec + k3;
      rcce_rhs(tn, u, k4);
      k4 = m_dt * k4;

      m_statevec = m_statevec + OneSixth * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
      
      tn = tn + m_dt;
      i  = i  + 1;      
   
      /* qoi */
      funt  = pow(m_Teq - m_T, 3.0) * pow(m_T - m_Tin, 2.0);
      qnum  = funt * log(tn) + old_funt * log(old_time);
      qden  = funt + old_funt;
      qoin += (log(tn) - log(old_time)) * qnum;
      qoid += (log(tn) - log(old_time)) * qden;
      qoit += 0.5 * (tn - old_time) * (m_T + old_temp); 

      /* output */
      if(m_ver == 1 /* && i%100 == 0 */) {
	cout << endl;
	cout << tn
	     << "\t" << m_Tor
	     << "\t" << m_T
	     << "\t" << m_Teq
	     << "\t" << m_Tin
	     << "\t" << qoit << endl;
	cout << endl;
      }
      if(m_vio == 1) {
	out << tn    << "\t"
	    << m_T   << "\t" 
	    << m_Teq << "\t"
	    << m_Tin << "\t"
	    << funt  << "\t" 
	    << qoin  << "\t" 
	    << qoid  << endl;
      }
      
      if(tn >= m_tf || abs(m_T - m_Teq) < 0.05) {
	if(tn < m_tf) {
	  qoit += m_Teq *(m_tf-tn);
	}
	if(qoin == 0.0 && qoid == 0 && abs(m_T - m_Teq < 0.05)) {
	  qoid = 1.0;
	}
	out << m_tf << "\t" 
	    << m_T  << "\t" 
	    << m_Teq << "\t"
	    << m_Tin << "\t"
	    << funt << "\t" 
	    << qoin << "\t" 
	    << qoid << endl;
	out.close();
	break;
      }

    }

    high_resolution_clock::time_point end_clock = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end_clock - start_clock).count();
    RTM = duration;

    QOI[0] = qoit;
    QOI[1] = qoin/qoid;

  };

  void rcceInt::rcce_rhs(double& t, VectorXd& statevec, VectorXd& R) {

    vector<double> ci(m_nc,0.0);
    vector<double> zi(m_kk,0.0);
    vector<double> fz(m_kk,0.0);

    VectorXd::Map(&ci[0], statevec.size()) = statevec;
    localceq_mp_ceqcalc_(m_nc, m_kk, &m_HR, &ci[0], &m_T, &zi[0]);

    m_gas->getNetProductionRates(m_T, zi, fz);
    m_fz = VectorXd::Map(&fz[0], fz.size());
    m_R  = m_C.transpose() * m_fz;
    R    = m_R;
     
  };

  void rcceInt::rcce_jac() {

  };

  
}
