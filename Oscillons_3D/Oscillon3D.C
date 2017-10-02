#include <omp.h>
#include "nr3.h"
//#include <sys/time.h>
 
//================================================================
//
// Code for Scalar Wave Code in Spherical Coordinates
//
//================================================================
#include "Read_Input.h"
#include "Oscillon.h"
#include "dumper.h"
//================================================================
//
// Main code
//
//================================================================
 
int main(int argc, char * argv[])
{
  //================================================================
  // get rid of old output
  //================================================================
  //cout << " Erase old output files? [ 1 = yes ] ";
  int answer = 1;
  //cin >> answer;

  cout << "\n Using this many threads: " << thread_num << "\n\n";

  //Time it
  cout << "Beginning to solve!" << endl;
  double start = omp_get_wtime();


  if (answer == 1) system("rm -f output_*\n");
  int error = 0;
  //================================================================
  // Read Input
  //================================================================
  int N_r = 0;          // Size of grid/number of grid piints
  int N_theta = 0;
  int N_phi = 0;
  double Courant;
  double aMax;     //Max scale factor
  double r_max;     //Max radius of grid
  double t_dump;    //Time at which grid writes data to a file
  int dump_step;    //Time interval at which grid writes data to a file
  int type;       //Determines if function should use PIRK or not
  double psi_init;  //Amplitude of initial psi function
  double A;  // Initial data: psi = psi_init * exp(-(r/A)^2) ????
  double a; //Scale factor
  double lambda, g, m, alphaF; //Oscillon parameters
  char dump_output[64];
  error = Read_Input(argc,argv,N_r,N_theta,N_phi,Courant,
         aMax,r_max,t_dump,dump_step,a,lambda,m,alphaF,dump_output,type);
  cout << alphaF << endl;

  if (error) return error;
  //================================================================
  // Allocate dumper Class
  //================================================================
  cout << " allocating dumper ..." << endl;
  dumper dump(dump_step,dump_output);
  //================================================================
  // Allocate Oscillon Class
  //================================================================
  cout << " allocating Oscillon ..." << endl;
  Oscillon oscillon(N_r,N_theta,N_phi,r_max,Courant,&dump,type);
  //================================================================
  // Initialize
  //================================================================
  error = oscillon.Initialize(lambda, m, alphaF);
  //================================================================
  // Integrate to time t_max
  //================================================================
  //oscillon.Integrate(aMax);

  //Finish Timing
  double end = omp_get_wtime();
  double solveTime  = end-start;
  cout << "Solution found in " << solveTime << " seconds" << endl;
  
  return error;
}
