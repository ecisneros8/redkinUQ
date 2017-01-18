#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

vector<int> agents;
vector<int> combination;

void choose(int offset, int k, vector< vector<int> >& v) {

  if (k == 0) {
    v.push_back(combination);
    return;
  }

  for(int i = offset; i <= agents.size() - k; ++i) {
    combination.push_back(agents[i]);
    choose(i+1, k-1, v);
    combination.pop_back();
  }

}

int main(int argc, char **argv) {

  int n  = 21;            // number of eligible species
  int k  = atoi(argv[1]); // number of constrained species
  int np = atoi(argv[2]); // number of processors in node
  int nc;                // total number of combinations
  int nf;                // total number of output files
  vector< vector<int> > v;

  cout << endl;
  cout << "Greetings!" << endl;
  cout << "Creating combinations..." << endl;
  for(int i = 0; i < n; ++i) { agents.push_back(i); }
  choose(0, k, v);

  nc = v.size();
  nf = (int) ceil(nc/(double)np);
  
  cout << "\t Number of eligible species:   " << n  << endl;
  cout << "\t Number of constraint species: " << k  << endl;
  cout << "\t Total number of combinations: " << nc << endl;
  cout << "\t Number of combinations/file : " << np << endl;
  cout << "\t Total number of files:        " << nf << endl;
  cout << endl;

  /* write combination files */
  ofstream out;
  int      count = 0;
  string   root  = "inps/comb/"+std::to_string(k)+"D/"+"csi.";
  string   filename;

  cout << "Creating input files..." << endl;
  for(int f = 0; f < nf; ++f) {
    filename = root+std::to_string(static_cast<long long>(f))+".dat";
    cout << "\t Current file: " << filename << endl;
    out.open(filename.c_str());
    int top = min((f+1)*np, nc);
    for(int i = f*np; i < top; ++i) {
      for(int j = 0; j < k; ++j) {
	out << v[i][j] << "\t";
      }
      out << endl;
    }
    out.close();
    ++count;
  }

  cout << endl;
  cout << "Done writing combination files..." << endl;
  cout << "\t Files created: " << count << endl;
  cout << endl;
  cout << endl;

  /* save number of files for use in job submission */
  out.open("jobs.sh");
  out << count << endl;
  out.close();

  return(0);

}
