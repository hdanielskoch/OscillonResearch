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

	int N= 256;
	VecDoub a(N),b(N),c(N),s(N),r(N),phi(N),dPhi(N);

	Doub rOut = 10.0;

	Doub alphaF = 0.9803;
    Doub alphaC = sqrt(27.0/160.0);
    Doub alpha = alphaF * alphaC;

	Doub dr = rOut / N;
	Doub dr2 = dr*dr;

	//Initialize Guess
	for(int i=0; i<N; i++){
		r[i] = (0.5 + i) * dr;

		if(r[i] < (0.5 * rOut)){
			phi[i] = exp(-1.0*r[i]*r[i]/10.0); //Gaussian
			//phi[i] = exp(-1.0*r[i]) / alpha;
		} else {
			phi[i] = exp(-1.0*alpha*r[i]) / r[i];
		}
		//Gaussian
		//phi[i] = exp(-1.0*r[i]*r[i]/10.0);
	}

	Doub residual = 0.0;
	Doub residNorm = 1.0;
	Doub tol = 0.000000001;
	int numIts = 0;

	//Iterate until norm is below residual
	while((residNorm > tol) && (numIts < 20)){

		numIts ++;
		//Reset norm
	    residNorm = 0;

		for(int i=0; i<N; i++){

			if(i==0){
				residual = (phi[i+1] - phi[i]) / (dr2) + 2.0 / r[i] *  (phi[i+1] - phi[i]) / (2.0 * dr) 
				         - alpha*alpha * phi[i] + 3.0/4.0 * phi[i]*phi[i]*phi[i] - 5.0/8.0 * phi[i]*phi[i]*phi[i]*phi[i]*phi[i];
			
			//Since phi goes with 1/r, r[i+1]*phi[i+1] = r[i]*phi[i] = const at boundary, plug in to laplace operator with nonlinear terms
			} else if (i == (N-1)){
				residual = ((r[i]*exp(-1.0*alpha*dr)/(r[i]+dr) - 2.0)/dr2 + exp(-1.0*alpha*dr) / (dr*(r[i] + dr))) * phi[i] + (1.0/dr2 -1.0/((r[i])*dr)) * phi[i-1]
				           - alpha*alpha * phi[i] + 3.0/4.0 * phi[i]*phi[i]*phi[i] - 5.0/8.0 * phi[i]*phi[i]*phi[i]*phi[i]*phi[i];
					
				//residual = 0.0;//exp(-1.0*alpha * r[i]) * ( (alpha*alpha - alpha)/r[i] + (2.0*alpha - 1.0)/(r[i]*r[i]) + 2.0/(r[i]*r[i]*r[i]) )
							//- alpha*alpha * phi[i] + 3.0/4.0 * phi[i]*phi[i]*phi[i] - 5.0/8.0 * phi[i]*phi[i]*phi[i]*phi[i]*phi[i];

			} else{
				residual = (phi[i+1] - 2.0*phi[i] + phi[i-1]) / (dr2) + 2.0 / r[i] *  (phi[i+1] - phi[i-1]) / (2.0 * dr)
								- alpha*alpha * phi[i] + 3.0/4.0 * phi[i]*phi[i]*phi[i] - 5.0/8.0 * phi[i]*phi[i]*phi[i]*phi[i]*phi[i];
			}
			s[i] = -1.0*residual;
			residNorm += abs(residual);
		}

		for(int i=0; i<N; i++){
			if (i==0){
				b[i] = -1.0 / dr2 - 2.0 /(r[i] * 2.0 * dr) - alpha*alpha + 9.0/4.0 * phi[i]*phi[i] - 25.0/8.0 * phi[i]*phi[i]*phi[i]*phi[i];
				c[i] = 1.0 / dr2 + 2.0 / (r[i] * 2.0 * dr);

			//Same boundary conidtion as above, but with linearized non-linear terms
			} else if(i==(N-1)){
				b[i] = (r[i]*exp(-1.0*alpha*dr)/(r[i] + dr) - 2.0)/dr2 + exp(-1.0*alpha*dr)/(dr*(r[i] + dr)) - alpha*alpha + 9.0/4.0 * phi[i]*phi[i] - 25.0/8.0 * phi[i]*phi[i]*phi[i]*phi[i];
				a[i] = 1.0/dr2 - 2.0/(r[i] * 2.0 * dr);
				//b[i] = 1.0;//exp(-1.0*alpha * r[i]) * ( (alpha*alpha - alpha)/r[i] + (2.0*alpha - 1.0)/(r[i]*r[i]) + 2.0/(r[i]*r[i]*r[i]) )
							//- alpha*alpha + 9.0/4.0 * phi[i]*phi[i] - 25.0/8.0 * phi[i]*phi[i]*phi[i]*phi[i];
				//a[i] = 0.0;//(1.0/dr2 -1.0/(r[i]*dr));
			}
			else{
				a[i] = 1.0/dr2 - 2.0 / (r[i]* 2.0 * dr);
				b[i] = -2.0 / dr2 - alpha*alpha + 9.0/4.0 * phi[i]*phi[i] - 25.0/8.0 * phi[i]*phi[i]*phi[i]*phi[i];
				c[i] = 1.0/dr2 + 2.0 / (r[i]* 2.0 * dr);
			}
		}
		tridag(a,b,c,s,dPhi);

		for(int i=0; i<N; i++){
			phi[i] += dPhi[i];
			//cout << r[i] << "    " << phi[i] << endl;
		}

	    cout << "Norm of residual" << residNorm << endl;
	} //End of while loop

	//Write values to a file
	cout << "Writing to file!" << endl;
    ofstream outfile;
    ostringstream rayfilename;
    rayfilename << "test.out" << setfill('0') << setw(4) << ends;
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