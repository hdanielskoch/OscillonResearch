#include <omp.h>
#include "nr3.h"
#include "gridfunction.h"
#include "dumper.h"
#include "tridag.h"
#include <string>
#include <iostream>
#include <sstream> 
#include "ludcmp.h"

int main(){

	int N=256;
	VecDoub a(N),b(N),c(N),s(N),r(N),phi(N),dPhi(N);

	Doub pi = acos(-1.0);

	Doub rOut = 16.0;
	Doub rSurf = 1.0;
	Doub rho = 1.0;
	Doub source = 4.0*pi*rho;

	Doub dr = rOut / N;
	Doub dr2 = dr*dr;

	for(int i=0; i<N; i++){

		r[i] = (0.5 + i) * dr;
		phi[i] = 0.0;

	}

	Doub residual = 0.0;
	Doub residNorm = 1.0;
	Doub tol = 0.00000001;
	int numIts = 0;

	//Iterate until norm is below residual
	while((residNorm > tol) && (numIts < 20)){

		numIts ++;
		//Reset norm
	    residNorm = 0;

		for(int i=0; i<N; i++){

			if(i==0){
				residual = (phi[i+1] - phi[i]) / (dr2) + 2.0 / r[i] *  (phi[i+1] - phi[i]) / (2.0 * dr) - source;
			} else if (i == (N-1)){
				//r[i+1]*phi[i+1] = r[i]*phi[i] = const at boundary, plug in to laplace operator
				residual = ((r[i]/(r[i]+dr) - 2.0)/dr2 + 1.0 / (dr*(r[i] + dr))) * phi[i] + (1.0/dr2 -1.0/(r[i]*dr)) * phi[i-1];

			} else{
				//Inside rSurf
				if(r[i] < rSurf){
					residual = (phi[i+1] - 2.0*phi[i] + phi[i-1]) / (dr2) + 2.0 / r[i] *  (phi[i+1] - phi[i-1]) / (2.0 * dr) - source;
				} else {
					residual = (phi[i+1] - 2.0*phi[i] + phi[i-1]) / (dr2) + 2.0 / r[i] *  (phi[i+1] - phi[i-1]) / (2.0 * dr);
	            }
			}
			s[i] = -1.0*residual;
			//cout << residual << endl;
			residNorm += abs(residual);

		}
		//cout << endl;

		for(int i=0; i<N; i++){

			if (i==0){
				b[i] = -1.0 / dr2 - 2.0 /(r[i] * 2.0 * dr);
				c[i] = 1.0 / dr2 + 2.0 / (r[i] * 2.0 * dr);
			} else if(i==(N-1)){
				b[i] = (r[i]/(r[i]+dr) - 2.0)/dr2 + 1.0 / (dr*(r[i] + dr));
				a[i] = 1.0/dr2 -1.0/(r[i]*dr);
			}
			else{
				a[i] = 1.0/dr2 - 2.0 / (r[i]* 2.0 * dr);
				b[i] = -2.0 / dr2;
				c[i] = 1.0/dr2 + 2.0 / (r[i]* 2.0 * dr);
			}
			if(r[i] < rSurf){
				//s[i] = source;
			} else{
				//s[i] = 0.0;
			}
		}
		tridag(a,b,c,s,dPhi);

		for(int i=0; i<N; i++){
			phi[i] += dPhi[i];
			cout << r[i] << "    " << phi[i] << endl;
		}

	    cout << residNorm << endl;
} //End of while loop

	//Write values to a file
	cout << "Writing to file!" << endl;
    ofstream outfile;
    ostringstream rayfilename;
    rayfilename << "test" << setfill('0') << setw(4) << ends;
    outfile.open(rayfilename.str().c_str());
    //    cout << " ... dumping file " << filename.str() << endl;
    outfile.setf(ios::left);
    outfile << "# Data at time " << time << endl;
    outfile << "# " << setw(14) << "r" << setw(16) << "psi" << endl;
    outfile << "#==================================================" << endl;
    outfile.setf(ios::right);

    for (int i = 0; i < N; i++) {
      outfile << setprecision(8) << setw(16) << r[i] << setw(16) << phi[i]  << endl; 
    }
    outfile.close();

}