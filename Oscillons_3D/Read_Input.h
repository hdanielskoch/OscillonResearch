// Tell emacs that this is -*-c++-*- mode
#include <omp.h>
#include "nr3.h"
//================================================================
//
// Source File for input
//
//================================================================
int Read_Input(int argc, char * argv[], 
         int & N_r, int & N_theta, int & N_phi, 
         double & Courant,
         double & aMax, double & r_max,
         double & t_dump, int & dump_step,
         double & a, double & lambda, 
         double & m, double alphaF, char * dump_output, 
         int & type)
{
  int error = 0;
  if (argc != 2) {
    cerr << "Syntax: Scalar Input_file" << endl;
    error = 1;
  } else {
    ifstream infile;

    //First argument is the input file
    infile.open(argv[1]);

    if(!infile){
      cerr << "Can't open " << argv[1] << " for input." << endl;
      error = 2;
      return error;
    } else {
      cout << "Reading input from file " << argv[1] << endl;
    }
    char buf[100],c;
    char run_case[8];

    //Gets all values from input file and writes them to variables  (Very well done TWB)
    infile.get(buf,100,'='); infile.get(c); infile >> N_r;
    infile.get(buf,100,'='); infile.get(c); infile >> N_theta;
    infile.get(buf,100,'='); infile.get(c); infile >> N_phi;
    infile.get(buf,100,'='); infile.get(c); infile >> Courant;
    infile.get(buf,100,'='); infile.get(c); infile >> aMax;
    infile.get(buf,100,'='); infile.get(c); infile >> r_max;
    infile.get(buf,100,'='); infile.get(c); infile >> t_dump;
    infile.get(buf,100,'='); infile.get(c); infile >> dump_step;
    infile.get(buf,100,'='); infile.get(c); infile >> a;
    infile.get(buf,100,'='); infile.get(c); infile >> lambda;
    infile.get(buf,100,'='); infile.get(c); infile >> m;
    infile.get(buf,100,'='); infile.get(c); infile >> alphaF;
    infile.get(buf,100,'='); infile.get(c); infile >> type;
    infile.get(buf,100,'='); infile.get(c); infile.get(c); 
    infile.get(dump_output,64);
    
    if(infile.eof()) {
      cerr << "Error reading input file " << argv[1] << endl;
      error = 3;
    }

    //Prints parameters of grid
    cout << "   Size of grid : (" << N_r << "," << N_theta << "," 
   << N_phi << ") " << endl;
    cout << "   Using Courant factor : " << Courant << endl;
    cout << "   Maximum scale factor : " << aMax << endl;
    cout << "   Grid extends to r_max : " << r_max << endl;
    cout << "   Dumping after times : " << t_dump << endl;
    cout << "   Dumping output every " << dump_step << " timesteps " << endl;
   cout << " and a = " << a << endl;
   cout << "lambda = " << lambda << endl;
   cout << "Effective Mass = " << m << endl;
   cout << "alpha factor = " << alphaF << endl;
    cout << "   Writing dumps into files '" << dump_output << "'" << endl;
    
    //Type specifies PIRK or Lax-Friedrichs HRSC
    if (type == 1) 
      cout << "   Running with PIRK " << endl;
    else if (type == 2)
      cout << "   Running with Lax-Friedrichs HRSC " << endl;
    else
      error = 4;
  }
  return error;
}
