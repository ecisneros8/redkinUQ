extern "C" {

  void ceqct_mp_ceqinit_(int& ncs, int& ne, int& ng, int& ns, int* csi,
			 double* T0, double* p0, double* h0,
			 double* wts, double* x0, double* ther, double* Ein, double* Bg,
			 double* Tceq, double* zceq);

  void ceqct_mp_ceqrecon_(double* ri, double* T, double* x, int& flag);

}

void ceqInitialize(int& ncs, int& ne, int& ng, int& ns, int* csi,
		   double& T0, double& p0, double& h0,
		   double* wts, double* x0, double* ther, double* Ein, double* Bg,
		   double& Tceq, double* zceq) {

  ceqct_mp_ceqinit_(ncs, ne, ng, ns, csi, &T0, &p0, &h0, wts, x0, ther, Ein, Bg, &Tceq, zceq);
  
}

void speciesReconstruction(double* ri, double& T, double* x, int& flag) {

  ceqct_mp_ceqrecon_(ri, &T, x, flag);
  
}
