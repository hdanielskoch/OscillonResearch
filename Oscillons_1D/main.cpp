//Wave model main function
#include <omp.h>
#include "Wave.h"
//#include <sys/time.h>

using namespace std;


int main(int argc, char* argv[])
{
		
	//Read in file
	string fileName = string(argv[1]);

	Wave wave(fileName);

	//Start solving the Problem and time it
    cout << "Beginning to solve!" << endl;
    double start = omp_get_wtime();
    

	wave.init();

	wave.propagate();

	double end = omp_get_wtime();


    double solveTime  = end-start;
    cout << "Solution found in " << solveTime << " seconds" << endl;


	wave.deconstructor();

	return 0;
}














