//Relaxation.c

#include "nr3.h"
#include "tridag.h"


int main(){

  //Write values to files
	ofstream outfile_pos;
	outfile_pos.open("TriDag");
	outfile_pos.setf(ios::left);
	outfile_pos << setw(16) << "# p " << setw(16) << " Phi " << endl;
	outfile_pos << "#====================================================" << endl;

   //Number of grid cells for phi
	int M = 1000;

  //Define x boundaries using cell centered grid (never actually use r_B)
	Doub rhoInit = 0.0;
	Doub rhoEnd = 100.0;


  //Initiate phi which will be incremented by dPhi
  //Initialize radius
	VecDoub phi(M), r(M);

  //From 1D analytical solution
	Doub lambda = 0.0000025;
	Doub g = lambda / sqrt(0.2);
	//Amplitude of Cos term
    Doub alphaC = sqrt(27.0/160.0);
    Doub alphaF = 1.001;
    Doub alpha = alphaF * alphaC;
	Doub u = sqrt(1.0 - (alpha / alphaC) * (alpha / alphaC));
	Doub Psi_0 = sqrt(9 * lambda * (1.0 - u) / (10 * g*g));

  //Make guess at initial phi
	for(int i=0; i < M; i++){

		r[i] = (Doub)i / rhoEnd;

		phi[i] = sqrt(lambda) * Psi_0 * sqrt((1.0 + u) / (1.0 + u * cosh(2.0 * alpha * lambda * (r[i] / g))));
	}


  //phi[0] = phiInit;

  Doub dr = 100.0 / M; //cell size (0.1 for M=1000)
  Doub dr2 = dr*dr;

  //Initialize vectors for tridiagonal system
  VecDoub main_Diag(M),up_Diag(M-1),down_Diag(M-1);

  //Some tridiagonal matrix *y = s
  //s is the residual
  //dPsi is the increment of psi
  VecDoub s(M),dphi(M);

  //Set boundary conditions
  //topleft (at x=1)
  //Doub c_A = (4.0 + lambda*delta_x) / (4.0 - lambda*delta_x); //signs reversed

  //bottom right (at x=0)
  //Doub c_B;
  Doub residNorm = 0;
  Doub tolerance = 0.001;

  //Iterate until norm of residual is small enough
  while(residNorm > tolerance){

  	for(int i=0; i < M; i++){
	//Boundary conidtions for main diagonal
  		if (i==0) {

  			main_Diag[i] = (-1.0/dr2) + (1.0 / (r * dr))- alpha*alpha + 9.0/4.0 * phi[j]*phi[j] - 25.0/8.0 * phi[j]*phi[j]phi[j]*phi[j];;

  			up_Diag[i] = (1.0/dr2) + (1.0 / (r * dr));


  		} else if(i == (M-1)) {

  			main_Diag[i] = (-1.0/dr2) + (1.0 / (r * dr))- alpha*alpha + 9.0/4.0 * phi[j]*phi[j] - 25.0/8.0 * phi[j]*phi[j]phi[j]*phi[j];;

  			down_Diag[i] = (1.0/dr2) - (1.0 / (r * dr));


  		} else {
			//Set diagonal matrices
  			Doub up_Diag_Value = (1.0/dr2) + (1.0 / (r * dr));
  			Doub down_Diag_Value = (1.0/dr2) - (1.0 / (r * dr));
  			Doub main_Diag_Value = (-2.0 / dr2) - alpha*alpha + 9.0/4.0 * phi[j]*phi[j] - 25.0/8.0 * phi[j]*phi[j]phi[j]*phi[j];

  			up_Diag[i] = up_Diag_Value;
  			down_Diag[i-1] = down_Diag_Value;
  			main_Diag[i] = main_Diag_Value;

  		}
  	}

	//Calculate residual
  	for(int j=0; j < M; j++){
	//On left boundary
  		if(j==0){

  			residual[j] = (phi[j+1] - phi[j]) / (dr2) + 2.0 / r *  (phi[j+1] - phi[j]) / dr 
  			- alpha*alpha * phi[j] + 3.0/4.0 * phi[j]*phi[j]*phi[j] - 5.0/8.0 * phi[j]*phi[j]phi[j]*phi[j]*phi[j];

  			else if(j==(M-1)){

  				residual[j] = (- phi[j] + phi[j-1]) / (dr2) + 2.0 / r *  (phi[j] - phi[j-1]) / (2*dr) 
  				- alpha*alpha * phi[j] + 3.0/4.0 * phi[j]*phi[j]*phi[j] - 5.0/8.0 * phi[j]*phi[j]phi[j]*phi[j]*phi[j];

  			}

  		} else{

  			residual[j] = (phi[j+1] - 2*phi[j] + phi[j-1]) / (dr2) + 2.0 / r *  (phi[j+1] - phi[j-1]) / (2*dr) 
  			- alpha*alpha * phi[j] + 3.0/4.0 * phi[j]*phi[j]*phi[j] - 5.0/8.0 * phi[j]*phi[j]phi[j]*phi[j]*phi[j];

  			s[j] = residual[j];
  		}

  	}

	//Use the tridag method to solve for y
  	tridag(down_Diag,main_Diag,up_Diag,s,dphi);

	//Set new phi values
  	for(int i=0; i < M; i++){

  		phi[i] += dphi[i];

  	}

  	//Calculate norm of residual
  	for(int i=0; i < M; i++){

  		residNorm += residual[j];

  	}

  } //End of while loop, residual is now below norm

  //Write values to a file
 



