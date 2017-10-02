// Tell emacs that this is -*-c++-*- mode
//================================================
// Class containing all the dumping stuff
//================================================
//
#ifndef DUMPER_H
#define DUMPER_H

//#include <omp.h>
#include "nr3.h"
#include "gridfunction.h"

class dumper {
private:
  int dump_step;
  char * filestem;
public:
  //================================================
  // Constructor
  //================================================
  dumper(int step, char * filename_i): 
    dump_step(step), filestem(filename_i) {

  };
  //================================================
  // Destructor
  //================================================
  ~dumper() {};
  //================================================
  // "Conditional dump"
  //================================================
  void cond_dump(double time, int timestep, 
  		 gf3d * psi, gf3d * psiana) {
    if (timestep % dump_step == 0)  // time to dump?
      dump(time,timestep,psi,psiana);
  };
  //================================================
  // forced dump
  //================================================
  void dump(double time, int timestep, gf3d * psi, gf3d * psiana) 
  {
    //================================================
    // first right out rays
    //================================================
    ofstream outfile;
    ostringstream rayfilename;
    rayfilename << filestem << "_rays_" << setfill('0') << setw(4)
	     << timestep << ends;
    outfile.open(rayfilename.str().c_str());
    //    cout << " ... dumping file " << filename.str() << endl;
    outfile.setf(ios::left);
    outfile << "# Data at time " << time << endl;
    outfile << "# " << setw(14) << "r" << setw(16) << "psi"  << setw(16) 
	    << "psi_ana"  << endl;
    outfile << "#==================================================" << endl;
    outfile.setf(ios::right);
    int N_r = psi->dim1();
    int N_theta = psi->dim2();
    int N_phi = psi->dim3();
    //    cout << " Dimensions of gridfuction in dumper : (" << N_r << "," 
    //	 << N_theta << "," << N_phi << ") " << endl;
    int nt = N_theta/2;
    int np = N_phi/2;
    for (int i = 1; i < N_r; i++) {
      outfile << setprecision(8) << setw(16) << psi->r(i) 
	      << setw(16) << (*psi)[i][nt][np] 
	      << setw(16) << (*psiana)[i][nt][np] << endl; 
    }
    outfile.close();
    //================================================
    // Now surfaces
    //================================================
    ostringstream surffilename;
    surffilename << filestem << "_surfaces_" << setfill('0') << setw(4)
	     << timestep << ends;
    outfile.open(surffilename.str().c_str());
    //    cout << " ... dumping file " << filename.str() << endl;
    outfile.setf(ios::left);
    int nr = N_r/4;
    outfile << "# Data at time " << time << ", radius " << psi->r(nr) << endl;
    outfile << "# " << setw(14) << "theta" << setw(16) << "phi"  << setw(16) 
	    << "psi_ana"  << endl;
    outfile << "#==================================================" << endl;
    outfile.setf(ios::right);
    //    cout << " Dimensions of gridfuction in dumper : (" << N_r << "," 
    //	 << N_theta << "," << N_phi << ") " << endl;
    for (int j = 1; j < N_theta-1; j++) 
      for (int k = 1; k < N_phi-1; k++) {
	outfile << setprecision(8) << setw(16) << psi->theta(j) 
		<< setw(16) << psi->phi(k) 
		<< setw(16) << (*psi)[nr][j][k] 
		<< setw(16) << (*psiana)[nr][j][k] << endl; 
      }
    outfile.close();
    
  }
};
 
#endif /* DUMPER_H */

