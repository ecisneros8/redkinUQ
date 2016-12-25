#include "IdealGasMix.h"
#include "rcceInteg.h"
#include "mpiEnvironment.h"
#include <cfloat>
#include <vector>

using namespace std;

class uqSolver {

 public:
    
  uqSolver(mpiEnvironment& env, IdealGasMix& gas) {

    rcceInteg* integ = new rcceInteg(gas);
    m_env = &env;

  };

  Integrator&     integ();
  mpiEnvironment& env();
  void            loadData(string& fuel);
  vector<double>  realization(vector<int>& csi);
  vector<double>  realizationFromFiles(vector<int>& csi);
  double          likelihood(double& sigma);
  double          gaussian(double& x);
  double          adaptiveQuadrature(double& a, double& b, double& tol);
  double          adaptiveSimpsons(int& level, int& maxLevel, 
				   double& a, double& b, double& tol,
				   double& IW);
  double          simpsonsRule(double& a, double& b);
  void            posteriorRepresentation(int& rank, int& jump, int& ncsi);
  void            getCombinations(int offset, int k, vector< vector<int> >& v);
  void            run(string& fuel, int& vio, int& ver,
		      double& t0, double& tf, double& dt);    

 protected:
  string                   m_fileroot; 
  string                   m_filename;
  string                   m_ndims;
  string                   m_fuel;
  int                      m_ncs;      // model class
  int                      m_Nk;    
  double                   m_ddm;
  vector<double>           m_pk;
  vector<double>           m_Sk;       // scenarios
  vector<double>           m_Dk;       // data
  vector< vector<double> > m_yk;
  vector<int>              m_combination;
  vector< vector<int> >    m_allcsi;
  Integrator*              m_integ;
  mpiEnvironment*          m_env;

};

Integrator& uqSolver::integ() {
  return *m_integ;
}; 

mpiEnvironment& uqSolver::env() {
  return *m_env;
};

void uqSolver::loadData(string& fuel) {

  /* open data files */
  if(fuel == "H2O") {
      
    ifstream inp;
    string   filename = "inps/tign.Bhaskaran.dat";
    inp.open(filename.c_str());

    /* read in data */
    while(!inp.eof()) {
      double scen;
      double data;
      inp >> scen >> data;
      m_Sk.push_back(scen);
      m_Dk.push_back(data);
    }
    m_Sk.pop_back();
    m_Dk.pop_back();
    m_Nk = m_Dk.size();
 
  } else if(fuel == "CH4") {

    ifstream inp;
    string   filename = "inps/tign.Zhukov.dat";
    inp.open(filename.c_str());

    /* read in data */
    while(!inp.eof()) {
      double pres;
      double Temp;
      double data;
      inp >> pres >> Temp >> data;
      m_pk.push_back(pres);
      m_Sk.push_back(Temp);
      m_Dk.push_back(data);
    }
    m_pk.pop_back();
    m_Sk.pop_back();
    m_Dk.pop_back();
    m_Nk = m_Sk.size();

  }

  /* resize model output vector */
  m_yk.resize(m_Nk);

  /* save fuel */
  m_fuel = fuel;

};

vector<double> uqSolver::realization(vector<int>& csi) {

  /* declarations */
  ofstream       out;
  string         filename;
  string         ndims;
  string         specs;
  int            flag;
  double         p;
  double         T0;
  double         phi;
  double         yk;
  vector<double> mout;

  /* output file */
  specs = "";
  for(int i = 0; i < csi.size(); ++i) { 
    ostringstream ostrs;
    ostrs << csi[i] + 1;
    specs = specs+ostrs.str()+"-";
  }

  filename = "outs/mout/"+m_ndims+"D/";
  filename = filename+"rcce";
  filename = filename+"."+specs+"S";
  filename = filename+".dat";

  cout.precision(4);
  out.precision(8);
  out.setf(ios::scientific);
  out.open(filename.c_str());

  /* run realization at every data point */
  for(int k = 0; k < m_Nk; ++k) {

    /* (1) set scenario */
    if(m_fuel == "H2O") {
      p   = 2.50 * Cantera::OneAtm;
      T0  = m_Sk[k];
      phi = 1.0;
    } else if(m_fuel == "CH4") {
      p   = m_pk[k] * Cantera::OneAtm;
      T0  = m_Sk[k];
      phi = 0.5;
    }
      
    /* (2) integrate */
    integ().setInitialConditions(fuel, p, T0, phi, csi);
    integ().integrate(flag, yk);
    mout.push_back(yk);

    /* (3) write output to file */
    out << Ti << "\t" << yk << std::endl;

  }

  /* close output file */
  out.close();

  /* return */
  return(mout);
    
}

vector<double> uqSolver::realizationFromFiles(vector<int>& csi) {
    
  ifstream       inp;
  string         filename;
  string         specs;
  ostringstream  ostrn;
  vector<double> mout;

  /* input file */
  specs = "";
  for(int i = 0; i < csi.size(); ++i) {
    ostringstream ostrs;
    ostrs << csi[i] + 1;
    specs = specs + ostrs.str() + "-";
  }

  filename = "outs/mout/"+m_ndims+"D/";
  filename = filename+"rcce";
  filename = filename+"."+specs+"S";
  filename = filename+".dat";

  /* retrieve model outputs */
  inp.open(filename.c_str());
  while(!inp.eof()) {
    double temp;
    double yk;
    inp >> temp >> yk;
    yk = exp(yk) * 1.0e6;
    mout.push_back(yk);
  }
  inp.close();
  mout.pop_back();

  return(mout);
     
}

double uqSolver::likelihood(double& sigma) {

  /* declaration */
  double g;
  double lh;
  double llh;

  g = 1.0 / sigma;
  for(int k = 0; k < m_Nk-1; ++k) { g *= 1.0 / sigma; }

  llh = -0.5 * m_ddm / (sigma * sigma);
  llh =  llh; // - m_Nk * log(sigma);
  lh  =  g * exp(llh);

  /* return */
  return(lh);

}

double uqSolver::adaptiveQuadrature(double& a, double& b, double& tol) {

  int maxLevel = 10;
  int myLevel  = 0;
  double IW = simpsonsRule(a, b);
  double I  = adaptiveSimpsons(myLevel, maxLevel, a, b, tol, IW);
  return(I);

}

double uqSolver::adaptiveSimpsons(int& level, int& maxLevel, 
				  double& a, double& b, double& tol, 
				  double& IW) {

  double c  = 0.5 * (a + b);
  double IL = simpsonsRule(a, c);
  double IR = simpsonsRule(c, b);
  double II = IL + IR;
  double I;

  if(level == maxLevel || fabs(II-IW) <= 15.0 * tol) {
    I = II + (II-IW)/15.0;
    return(I);
  }
    
  int    newlev = level + 1;
  double newtol = 0.5 * tol;
  I = adaptiveSimpsons(newlev,maxLevel,a,c,newtol,IL) + 
    adaptiveSimpsons(newlev,maxLevel,c,b,newtol,IR);
  return(I);

}

double uqSolver::simpsonsRule(double& a, double& b) {

  double c = 0.5 * (a + b);
  double fa = likelihood(a);
  double fb = likelihood(b);
  double fc = likelihood(c);
  double I  = OneSixth * (b-a) * (fa + 4.0 * fc + fb);
    
  return(I);

}
  
void uqSolver::posteriorRepresentation(int& rank, int& jump, int& ncsi) {

  /* declarations */
  ofstream                 out;
  string                   filename;
  double                   IL;
  double                   IR;
  double                   sigopt;
  double                   sigmin = 1.0e-03;
  double                   sigmax = 1.0;
  double                   tol    = 1.0e-06;
  vector<double>           prep;
  ostringstream            orank;

  /* output file */
  orank    << rank;
  filename = "outs/prep/"+m_ndims+"D/";
  filename = filename+"gri30.prep."+orank.str()+".dat";
  out.precision(6);
  out.setf(ios::scientific);
  out.open(filename.c_str());

  /* compute statistics */
  for(int i = 0; i < ncsi; ++i) {

    /* data mismatch */
    m_ddm = 0.0;
    for(int k = 0; k < m_Nk; ++k) {
      double d = (m_yk[i][k] - m_Dk[k]) / m_Dk[k];
      m_ddm += d*d;
    }
      
    sigopt = sqrt(2.0 * m_ddm / m_Nk);
      
    IL = adaptiveQuadrature(sigmin, sigopt, tol);
    IR = adaptiveQuadrature(sigopt, sigmax, tol);
    prep.push_back(IL+IR);

    /* output */
    cout << rank << "\t" << i*jump+rank << "\t";
    for(int s = 0; s < m_ncs; ++s) {
      cout << integ().gas().speciesName(m_allcsi[i*jump+rank-1][s]) << "\t";
      out  << integ().gas().speciesName(m_allcsi[i*jump+rank-1][s]) << "\t";
    }
    cout << prep[i] << endl;
    out  << prep[i] << endl;
      
  }

  cout << "Done with rank " << rank << endl;
  cout << endl;
    
  /* close output file */
  out.close();

}

void uqSolver::getCombinations(int offset, int k, vector< vector<int> >& v) {

  if(k == 0) {
    v.push_back(m_combination);
    return;
  }

  int lim;
  if(m_fuel == "H2O") {
    lim = m_kk - 2;
  } else if(m_fuel == "CH4") {
    lim = 20;
  }
  
  for (int i = offset; i <= lim; ++i) {
    string speciesName  = gas().speciesName(i);
    int    speciesIndex = gas().speciesIndex(speciesName);
    m_combinations.push_back(speciesIndex);
    getCombinations(i+1, k-1, v);
    m_combination.pop_back();
  }
  
}

void uqSolver::run(string& fuel, int& vio, int& ver,
		   double& t0, double& tf, double& dt) {

  ofstream              out;
  string                filename;
  int                   ncsi   = 0;
  int                   tag    = 0;
  int                   rank   = env().worldRank();
  int                   root   = env().worldRoot();
  int                   size   = env().worldSize();
  int                   jump   = env().worldSize()-1;
  vector<int>           csi;
  vector<int>           local_size;
  vector< vector<int> > local_csi;
  ostringstream         ostrncs;
  ostringstream         orank;

  /* model class */
  m_ncs   = env().worldOpts();
  ostrncs << m_ncs;
  m_ndims = ostrncs.str();

  /* load data */
  loadData(fuel);

  /* initialize integrator */
  integ().initialize(vio, ver, t0, tf, dt);

  /* access parameter space */
  getCombinations(0, m_ncs, m_allcsi);
  csi.resize(m_ncs);
    
  if(rank != root) {

    /* distribute combinations across processes */
    for(vector< vector<int> >::iterator it = m_allcsi.begin();
	it != m_allcsi.end(); it+=env().worldSize()-1) {
      if(it+rank-1 >= m_allcsi.end()) break;
      local_csi.push_back(*(it+rank-1));
    }

    /* number of local combinations */
    ncsi = local_csi.size();

    /* resize vector of model output vectors */
    m_yk.resize(ncsi);

    /* talk to root */
    msg[0] = rank;
    msg[1] = ncsi;
    env().sendMessage(&msg[0], 2 * sizeof(int), MPI_INT, root, tag);

  } else {

    filename = "outs/logs/"+m_ndims+"D/";
    filename = filename+"gri30.run.root.log";
    out.open(filename.c_str());

    /* greetings */
    out << endl;
    out << "Greetings from root!" << endl;
    out << "World size is " << env().worldSize() << endl;
    out << "Working with " << m_ncs << "D RCCE" << endl;
    out << "Integrators initialized with:" << endl;
    out << "\t vio = " << vio << endl;
    out << "\t ver = " << ver << endl;
    out << "\t t0  = " << t0 << endl;
    out << "\t tf  = " << tf << endl;
    out << "\t dt  = " << dt << endl;
    out << "Size of parameter space is " << m_allcsi.size() << endl;
    out << "Working with " << m_Nk << " data points" << endl;
    out << endl;
      
    for(int source = 1; source < env().worldSize(); ++source) {
      env().receiveMessage(&msg[0], 2 * sizeof(int), MPI_INT, source, tag);
      out << "Rank " << msg[0] << " is alive!" << endl;
      out << "Working with " << msg[1] << " combinations" << endl;
      out << endl;
      local_size.push_back(msg[1]);
    }
    out << endl;
    out.close();

  }

  /* run */
  if(rank != root) {

    /* open log file */
    orank    << rank;
    filename = "outs/logs/"+m_ndims+"D/";
    filename = filename+"gri30.run."+orank.str()+".log";
      
    out.open(filename.c_str());
      
    /* model outputs */
    for(int i = 0; i < ncsi; ++i) {
	
      /* print m.id. to log */
      out  << "Working with the " << i
	   << "th local combination, " << i*jump+rank
	   << "th overall..." << endl;
      out  << "Represented Species: " << endl;
      cout << "Represented Species: " << endl;
      for(int s = 0; s < m_ncs; ++s) {
	cout << integ().gas().speciesName( m_allcsi[i*jump+rank-1][s] ) << "\t";
	out  << integ().gas().speciesName( m_allcsi[i*jump+rank-1][s] ) << "\t";
      }
      cout << endl;
      out << endl;

      /* compute ignition times */
      csi  = local_csi[i];
      //vector<double> yi = realization(csi);
      vector<double> yi = realizationFromFiles(csi);
      cout << "Done!" << endl;
      cout << endl;
      m_yk[i] = yi;

    }

    out.close();

    /* compute and write posteriors */
    posteriorRepresentation(rank, jump, ncsi);

  }

}


