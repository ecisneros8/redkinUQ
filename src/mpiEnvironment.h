#include <sstream>
#include <mpi.h>

using namespace std;

class mpiEnvironment 
{
 public:
    
 mpiEnvironment(MPI_Comm inputComm, const char* OPT0, const char* OPT1) :
  m_worldRoot(0),
    m_worldRank(-1),
    m_worldSize(-1) {

    m_worldComm = inputComm;
    MPI_Comm_size(inputComm,&m_worldSize);
    MPI_Comm_rank(inputComm,&m_worldRank);
   
    m_worldOpts.resize(2);

    /* OPT0: number of constrained species */
    istringstream istr0(OPT0);
    istr0 >> m_worldOpts[0];

    /* OPT1: combination file id */
    istringstream istr1(OPT1);
    istr1 >> m_worldOpts[1];

  };

  int  worldRank();
  int  worldSize();
  int  worldRoot();
  std::vector<int> worldOpts();
  void sendMessage(void* buf, int count, MPI_Datatype dtype, 
		   int dest, int tag);
  void receiveMessage(void* buf, int count, MPI_Datatype dtype, 
		      int source, int tag);
    
 protected:
  MPI_Comm   m_worldComm;
  MPI_Status m_status;
  int        m_worldSize;
  int        m_worldRank;
  int        m_worldRoot;
  std::vector<int> m_worldOpts;

};

int mpiEnvironment::worldRank() {
  return m_worldRank;
};

int mpiEnvironment::worldRoot() {
  return m_worldRoot;
};

int mpiEnvironment::worldSize() {
  return m_worldSize;
};

std::vector<int> mpiEnvironment::worldOpts() {
  return m_worldOpts;
};

void mpiEnvironment::sendMessage(void* buf, int count, MPI_Datatype dtype, 
				 int dest, int tag) {
  
  MPI_Send(buf, count, dtype, dest, tag, m_worldComm); 
  
};

void mpiEnvironment::receiveMessage(void* buf, int count, MPI_Datatype dtype, 
				    int source, int tag) {

  MPI_Recv(buf, count, dtype, source, tag, m_worldComm, &m_status);

};



