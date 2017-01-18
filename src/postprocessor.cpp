#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>

using namespace std;

template <typename T>
vector<int> sort_indices(const vector<T>& v) {

  vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  sort(idx.begin(), idx.end(), [&v](int i1, int i2) { return v[i1] > v[i2]; }); 
  return(idx);

}

int main(int argc, char **argv) {

  int            M   = atoi(argv[1]); // number of constrained species
  int            ND  = atoi(argv[2]); // number of CSI directories
  int            NF  = atoi(argv[3]); // number of folders in each CSI directory
  int            NCK = atoi(argv[4]); // total number of combinations
  int            fexists;
  string         root;
  string         fileroot;
  string         filename;
  string         line;
  ifstream       inp;
  vector<string> allnames;
  vector<double> allp;

  root = "outs/prep/"+to_string(M)+"D/";
  
  /***********************************/
  /* assemble posterior distribution */
  /***********************************/
  cout << "Assembling posterior distribution..." << endl;
  cout << "Looping over CSI directories and P-rep files..." << endl;
  
  for(int i = 0; i <= ND; ++i) { 

    /* CSI folder */
    fileroot = root+"CSI"+to_string(i)+"/";
    
    for(int f = 1; f <= NF; ++f) {
      
      /* open file */
      fexists  = 0;
      filename = fileroot+"gri30.prep."+to_string(f)+".dat";
      inp.open(filename.c_str());
  
      /* read in data */
      while(inp) {

        fexists = 1;

        while(getline(inp,line)) {
	  
	  int     count = 0;
	  string  temp;
	  string::reverse_iterator rit = line.rbegin();
      
          while(*rit != '	') {
	    temp += *rit;
	    ++count;
	    ++rit;
	  }
          
	  string prep = string(temp.rbegin(), temp.rend());
	  string name = line.substr(0,line.size()-count);
	  double p    = atof(prep.c_str());
	  
	  allnames.push_back(name);
	  allp.push_back(p); 
          
	}

	inp.close();

      }

      if(fexists == 0) {
        cout << endl;	
	cout << filename << " does not exist..." << endl;
      }

    }

  }

  cout << endl;

  /******************************/
  /* sort and store in new file */
  /******************************/
  ofstream    outone;
  ofstream    outtwo;
  ofstream    outevd;
  string      fileone = root+"gri30."+to_string(M)+"D.preps.dat";
  string      filetwo = root+"gri30."+to_string(M)+"D.names.dat";
  string      fileevd = "outs/gri30.evd.dat";
  vector<int> idx = sort_indices(allp);

  outone.precision(6);
  outtwo.precision(6);
  outone.setf(ios::scientific);
  outtwo.setf(ios::scientific);
  outone.open(fileone);
  outtwo.open(filetwo);

  /* sort and save posterior & names */
  outone << "c pc" << endl;
  int lim = min(30,NCK);
  for(int i = 0; i < lim; ++i) {
    outone << i+1 << "\t" << allp[idx[i]] << endl;
    outtwo << allnames[idx[i]] << "\t" << allp[idx[i]] << endl;
    cout << allnames[idx[i]] << "\t" << allp[idx[i]] << endl;
  }
  cout << endl;
  outone.close();
  outtwo.close();

 /* compute evidence and save */
 double evd = 0.0;
 for(int i = 0; i < allp.size(); ++i) { evd += allp[i]; }
 evd = evd / NCK;

 outevd.precision(6);
 outevd.setf(ios::scientific);
 outevd.open(fileevd.c_str(), ofstream::out | ofstream::app);
 outevd << M << "\t" << evd << endl;
 cout << M << "\t" << endl;
 outevd.close();
  
}
