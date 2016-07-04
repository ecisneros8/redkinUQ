#include "IdealGasMix.h"
#include "odeIntegrators.h"
#include "mpiEnvironment.h"
#include <vector>

using namespace std;

namespace mech 
{

  class uqSolver {

  public:
    
    uqSolver(mpiEnvironment& env, IdealGasMix& gas) {

      m_integ = newIntegrator(gas,"RCCE");
      m_env   = &env;

    };

    Integrator&     integ();
    mpiEnvironment& env();
    void loadData();
    void likelihood(vector<int>& csi);
    void likelihoodFromFiles(vector<int>& csi);
    void run(int& vio, int& ver, double& t0, double& tf, double& dt);    

  protected:
    string          m_fileroot; 
    string          m_filename;
    string          m_ndims;
    int             m_modid;    // model class
    int             m_Nk;       // number of data points
    double          m_sig;      // sigma
    double          m_isig;     // inverse sigma
    double          m_llh;      // log likelihood
    double          m_lev;      // log evidence      
    double          m_lpr;      // log prior
    double          m_lps;      // log posterior
    vector<double>  m_Dk;       // data
    vector<double>  m_Qk;       // model class prediction
    vector<double>  m_vllh;     // all log likelihoods
    Integrator*     m_integ;    // pointer to integrator
    mpiEnvironment* m_env;      // point to mpi environment

  };

  Integrator& uqSolver::integ() {
    return *m_integ;
  }; 

  mpiEnvironment& uqSolver::env() {
    return *m_env;
  };

  void uqSolver::loadData() {

    ifstream inp;
    string   filename = "inps/rcce.logt.Dk.dat";
    inp.open(filename.c_str());

    while(!inp.eof()) {
      double data;
      inp >> data;
      m_Dk.push_back(data);
    }
    m_Dk.pop_back();
    m_Nk = m_Dk.size();

  };

  void uqSolver::likelihood(vector<int>& csi) {

    ofstream       out;
    string         filename;
    string         ndims;
    string         specs;
    double         Ti = 1000;
    double         rtm;
    vector<double> Qk(2,0.0);
    ostringstream  ostrn;
    //ostringstream  ostrs;

    /* output file */
    ostrn << m_modid;
    specs = "";
    for(int i = 0; i < csi.size(); ++i) { 
      ostringstream ostrs;
      ostrs << csi[i] + 1;
      specs = specs+ostrs.str()+"-";
    }
    ndims    = ostrn.str();
    //specs    = ostrs.str();

    filename = "outs/qois/"+ndims+"D/";
    filename = filename+"rcce";
    filename = filename+"."+specs+"S";
    filename = filename+".dat";

    out.precision(8);
    out.setf(ios::scientific);
    out.open(filename.c_str());

    /* set sigma */
    m_sig  = 0.2;
    m_isig = 1.0 / m_sig;

    /* run realization at every data point */
    m_llh = 0.0;
    for(int k = 0; k < m_Nk; ++k) {
      
      /* (1) integrate */
      integ().setInitialConditions(Ti, csi);
      integ().integrate(rtm, Qk);
      
      /* (2) compute log likelihood */
      m_llh -= pow( (m_Dk[k] - Qk[1]) * m_isig, 2.0 );
      m_Qk.push_back(Qk[1]);

      /* (3) write prediction to file */
      out << Ti << "\t" << Qk[i] << endl;

      /* (4) move on to next scenario */
      Ti = Ti + 100;

    }
    
    m_llh = 0.5 * m_llh;

    /* close output file */
    out.close();
    
  };

  void uqSolver::likelihoodFromFiles(vector<int>& csi) {

    ifstream       inp;
    string         filename;
    string         ndims;
    string         specs;
    int            flag = 0;
    vector<double> Qk;
    ostringstream  ostrn;

    /* output file */
    ostrn << m_modid;
    specs = "";
    for(int i = 0; i < csi.size(); ++i) { 
      ostringstream ostrs;
      ostrs << csi[i] + 1;
      specs = specs+ostrs.str()+"-";
    }
    ndims    = ostrn.str();
    m_ndims  = ndims;

    filename = "outs/qois/"+ndims+"D/";
    filename = filename+"rcce";
    filename = filename+"."+specs+"S";
    filename = filename+".dat";

    /* read in simulation data */
    inp.open(filename.c_str());
    while(!inp.eof()) {
      double scen;
      double data;
      inp >> scen >> data;
      Qk.push_back(data);
    }

    /* close output file */
    inp.close();

    /* set sigma */
    m_sig  = 0.02;
    m_isig = 1.0 / m_sig;

    /* run realization at every data point */
    m_llh = 0.0;
    for(int k = 0; k < m_Nk; ++k) {
      
      /* (1) compute log likelihood */
      m_llh -= pow( (m_Dk[k] - Qk[k]) * m_isig, 2.0 );

    }
    
    m_llh = 0.5 * m_llh;
    
  };

  void uqSolver::run(int& vio, int& ver, double& t0, double& tf, double& dt) {

    int    tag = 0;
    int    rank = env().worldRank();
    int    root = env().worldRoot();
    int    size = env().worldSize();
    double llh;
    vector<int>              msg(2,0);
    vector<int>              csi;
    vector<int>              local_size;
    vector< vector<int> >    local_csi;
    vector< vector<int> >    all_csi;
    vector<double>::iterator mllh;

    /* model class */
    m_modid = env().worldOpts();

    /* load data */
    loadData();
 
    /* initialize integrator */
    integ().initialize(vio, ver, m_modid, t0, tf, dt);

    /* access parameter space */
    integ().gas().getCombinations(0, m_modid, all_csi);
    csi.resize(m_modid);
    
    /* greetings */
    if(rank != root) {

      for(vector< vector<int> >::iterator it = all_csi.begin();
	  it != all_csi.end(); it+=env().worldSize()-1) {
	if(it+rank-1 >= all_csi.end()) break;
	local_csi.push_back(*(it+rank-1));
      }

      msg[0] = rank;
      msg[1] = local_csi.size();
      env().sendMessage(&msg[0], 2 * sizeof(int), MPI_INT, root, tag);

    } else {

      cout << endl;
      cout << "Greetings from root!" << endl;
      cout << "World size is " << env().worldSize() << endl;
      cout << "Working with " << m_modid << "D RCCE" << endl;
      cout << "Integrators initialized with:" << endl;
      cout << "\t vio = " << vio << endl;
      cout << "\t ver = " << ver << endl; 
      cout << "\t t0  = " << t0 << endl;
      cout << "\t tf  = " << tf << endl;
      cout << "\t dt  = " << dt << endl;
      cout << "Size of parameter space is " << all_csi.size() << endl;
      cout << "Working with " << m_Nk << " data points" << endl;
      cout << endl;
      
      for(int source = 1; source < env().worldSize(); ++source) {
	env().receiveMessage(&msg[0], 2 * sizeof(int), MPI_INT, source, tag);
	cout << "Rank " << msg[0] << " is alive!" << endl;
	cout << "Working with " << msg[1] << " combinations" << endl;
	cout << endl;
	local_size.push_back(msg[1]);
      }
      cout << endl;

    }

    if(rank != root) {

      for(int i = 0; i < local_csi.size(); ++i) {
	
	/* compute likelihood */
	csi = local_csi[i];
	likelihood(csi);
	
	/* send likelihood to root */
	env().sendMessage(&m_llh, sizeof(double), MPI_DOUBLE, root, tag);

      }

    } else {

      for(int source = 1; source < env().worldSize(); ++source) {	

	for(int i = 0; i < local_size[source-1]; ++i) {

	  /* receive llh from source */
	  cout << "Receiving likelihood from Rank " << source << endl;
	  env().receiveMessage(&llh, sizeof(double), MPI_DOUBLE, source, tag);
	  m_vllh.push_back(llh);
	
	  /* print to screen */
	  int jump = env().worldSize()-1;
	  cout << source << "\t" << i*jump+source << "\t";
	  for(int s = 0; s < m_modid; ++s) {
	    cout << integ().gas().m_speciesNames[all_csi[i*jump+source-1][s]] 
		 << "\t";
	  }
	  cout << llh << endl;
	  cout << endl;

	}

      };

      /* compute evidence */
      mllh  = max_element(m_vllh.begin(), m_vllh.end());
      m_lev = 0.0;
      for(int i = 0; i < m_vllh.size(); ++i) { m_lev += exp(m_vllh[i] - *mllh); }
      m_lev = log(m_lev) + *mllh;

      /* impose occam penalty */
      m_lpr = -log(all_csi.size());
      m_lps =  m_lev + m_lpr;

      /* final output */
      cout << endl;
      cout << m_modid << "\t"
    	   << m_lev   << "\t"
    	   << m_lpr   << "\t"
    	   << m_lps   << endl;
      cout << endl;

      ostringstream ostrn;
      ostrn << m_modid;
      m_ndims = ostrn.str();
      m_filename = "outs/pdfs/rcce.llh.";
      m_filename = m_filename+m_ndims+"D.bin";
      ofstream outs;
      outs.open(m_filename.c_str(), ios::out | ios::binary);
      outs.write((char*) &m_vllh[0], m_vllh.size() * sizeof(double));
      outs.close();
      
      m_filename = "outs/pdfs/rcce.pdf.";
      m_filename = m_filename+m_ndims+"D.bin";
      outs.open(m_filename.c_str());
      outs << m_modid << "\t"
    	   << m_lev   << "\t"
    	   << m_lpr   << endl;
      outs.close();
    }

  };

}

