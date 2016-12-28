#include <sstream>
#include <mpi.h>

using namespace std;

class mpiEnvironment 
{
 public:
    
 mpiEnvironment(MPI_Comm inputComm, const char* options) :
  m_worldRoot(0),
    m_worldRank(-1),
    m_worldSize(-1) {

    m_worldComm = inputComm;
    MPI_Comm_size(inputComm,&m_worldSize);
    MPI_Comm_rank(inputComm,&m_worldRank);

    istringstream istr(options);
    istr >> m_worldOpts;

  };

  int  worldRank();
  int  worldSize();
  int  worldRoot();
  int  worldOpts();
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
  int        m_worldOpts;

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

int mpiEnvironment::worldOpts() {
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



