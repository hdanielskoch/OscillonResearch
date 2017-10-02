//Wave.cpp
//NOT Dover

//All functions are performed in the class Wave
#include <omp.h>
#include "nr3.h"
#include "Wave.h"
#include "ran.h"
#include "gamma.h"
#include "deviates.h"
//#include "dumper.h"
#include <cstdio>


using namespace std;

Wave::Wave(string fileName)
{
	
	/*#pragma omp parallel
  	{
    int thread_Num = omp_get_num_threads();
    cout << "Num threads " << thread_Num << endl;
	}*/
  

    /*#pragma omp for
    for(int i = 0; i < 10; ++i)
    {
      cout << i << endl;
    }
  	}*/


	//Read in parameters from a file
	cout << "Reading in " << fileName << endl << endl << endl;

	ifstream file(fileName.c_str());

	//File could not be opened
    if (!file) {
		std::cout << "Failed to open file. Exiting." << std::endl;
		return;
    }

    string line;

    //Run through file and print inputs
    while(getline(file, line))
    {
    	if(line.compare("Minimum position:") == 0){
    		getline(file, line);

    		stringstream ss(line);
    		ss >> xMin;
    		cout << "xMin: " << xMin << endl;
    		getline(file, line);
    	}
    	if(line.compare("Maximum position:") == 0){
    		getline(file, line);

    		stringstream ss(line);
    		ss >> xMax;
    		cout << "xMax: " << xMax << endl;
    		getline(file, line);
    	}
    	if(line.compare("Number of position grid points:") == 0){
    		getline(file, line);

    		stringstream ss(line);
    		ss >> nx;
    		cout << "Numer of grid points: " << nx << endl;
    		getline(file, line);
    	}
    	if(line.compare("Maximum scale factor a:") == 0){
    		getline(file, line);

    		stringstream ss(line);
    		ss >> aMax;
    		cout << "Maximum a: " << aMax << endl;
    		getline(file, line);
    	}
    	if(line.compare("Courant condition:") == 0){
    		getline(file, line);

    		stringstream ss(line);
    		ss >> courantCond;
    		cout << "Courant Factor: " << courantCond << endl;
    		getline(file, line);
    	}
    	if(line.compare("Alpha Factor:") == 0){
    		getline(file, line);

    		stringstream ss(line);
    		ss >> alphaF;
    		cout << "Alpha = " << alphaF << " * alphaC " << endl;
    		getline(file, line);
    	}
    	if(line.compare("Effective Mass:") == 0){
    		getline(file, line);

    		stringstream ss(line);
    		ss >> m;
    		cout << "m = " << m << endl;
    		getline(file, line);
    	}
    	if(line.compare("lambda (dimensionless parameter):") == 0){
    		getline(file, line);

    		stringstream ss(line);
    		ss >> lambda;
    		cout << "lambda = " << lambda << endl;
    		getline(file, line);
    	}
    	if(line.compare("Maximum Grid Factor:") == 0){
    		getline(file, line);

    		stringstream ss(line);
    		ss >> gridFactMax;
    		cout << "Maximum Grid Factor = " << gridFactMax << endl;
    		getline(file, line);
    	}
    	if(line.compare("File Increment:") == 0){
    		getline(file, line);

    		stringstream ss(line);
    		ss >> fileInc;
    		cout << "File Increment = " << fileInc << endl;
    		getline(file, line);
    	}
    	if(line.compare("H implementation") == 0){
    		getline(file, line);

    		stringstream ss(line);
    		ss >> HImplement;
    		if(HImplement == 1){
    			cout << "Implementing H = a^-0.5 " << endl;
    		} else {
    			cout << "Implementing H = 1 / (3m_pl^2) * <p> " << endl;
    		}
    		getline(file, line);
    	}
    	
    }
    cout << endl;

	//Allocate data structures to hold memory up at a = gridFactMax
	//Doesn't initially fill all the grid points, but it has the capacity to

    x_Pos = new VecDoub(nx);

	v_O = new VecDoub(nx);
	v_N = new VecDoub(nx);
	v_T = new VecDoub(nx);

	u_O = new VecDoub(nx);
	u_N = new VecDoub(nx);
	u_T = new VecDoub(nx);

	u_I = new VecDoub(nx);
	u_A = new VecDoub(nx);
	uTemp = new VecDoub(nx);
	vTemp = new VecDoub(nx);

}



//Initialize u and v with initial displacement and velocity
//Choose an initial, stationary gaussian displacement
void Wave::init()
{

	//position increment
	xInc = xMax / Doub(nx);

	//time increment
	tInc = courantCond * xInc;

	cout << "deltaX = " << xInc << endl << "deltaT = " << tInc << endl;

	g = lambda / sqrt(0.2);

	//Changeable parameter
	Doub m_Mpl = 0.000005;

	H_i = sqrt(lambda / (10 * g*g)) * m * m_Mpl;

	xi = lambda / g; 

	//Mid point of analytical initial function
	Doub xMid1 = (xMax + xMin) / 8.0;

	Doub xMid2 = xMid1 + xMax;

	Doub xMid3 = xMid1 - xMax;

	//Initial data that uses analytical solution
	Doub alphaC = sqrt(27.0/160.0);
	alpha = alphaF * alphaC;
	Doub u = sqrt(1.0 - (alpha / alphaC) * (alpha / alphaC));
	Doub Psi_0 = sqrt(9 * lambda * (1.0 - u) / (10 * g*g) );



	//Appendix A Implementation of random fluctuations
	Doub L_0 = nx; //Size of initial box
	int N = 50; //Number of terms in sum

	Doub N_Init = -N/2.0 + 1.0;
	Doub N_Final = N/2.0;


	Doub psiSum; //Sum in initial data
	//Doub psiPSum;
	Doub k_N, w_N;
	Doub pi = 2 / cos(0);
	Doub psiH_I = 0.36; //Mean phi at 0??
	Doub aReal,bImag; //Real and Imaginary parts of alphaN

	//Generate alpha_N randomly

	//Phase
	Doub ranSeedPhase = clock();
	Ran ranPhase(ranSeedPhase);

	//Amplitudes are drawn from Gaussian distribution with variance lambda / 2
	Doub mu = 0.0;  //Generally taken to be 0
	Doub sigma = sqrt(lambda / 2.0);
	Doub ranSeedAmp = 3 * clock();
	Normaldev ranAmp(mu, sigma, ranSeedAmp);

	Doub phaseUnit = 0;
	Doub phase = 0;
	Doub ampUnit = 0;
	Doub amp = 0;


	int i;
	//int threadNum;
	//This loop can be done in any order
	//However all the x_Pos must be filled first becuase u_0 uses x_Pos
	/*#pragma omp parallel num_threads(1)
	{
	#pragma omp for collapse(1) private(i)*/
	


	//Iterate over all positions of initial (smallest) grid size
	for(i = 0; i < nx; i++){

		//Fill position vecDoub
		(*x_Pos)[i] = ((Doub)i + 0.5) * xInc;
		
		//Calculate sums for initial displacement and derivatives
		for(int n = N_Init; n < N_Final; n++){

			//Generate Random Phase
			phaseUnit = ranPhase.doub(); //Phase between 0 and 1
			phase = phaseUnit * 2 * pi;

			//Generate Random Amplitude
			ampUnit = ranAmp.dev(); //Amplitude of gaussian
			amp = ampUnit; //Multiply it by Amplitude of gaussian (about 1000)???

			//Determine real and imaginary parts of alphaN
			aReal = amp * cos(phase);
			bImag = amp * sin(phase);

			k_N = 2 * pi * n / L_0;

			w_N = sqrt(1.0 + (2 * sin(k_N * xInc/2.0) / xInc) * (2 * sin(k_N * xInc/2.0) / xInc));

			psiSum = sqrt(1.0 / (2.0 * w_N)) * (2.0 * aReal * cos(k_N * i * xInc) - 2.0 * bImag * sin(k_N * i * xInc));

			//psiPSum = (1.0 / imag) * sqrt(w_N / 2.0) * alpha_N * exp(imag * k_N * (*x_Pos)[i]) - cc;
		}


		//Calculate Psi at t=0
		(*u_O)[i] = psiH_I + (1.0 / sqrt(L_0)) * psiSum;
		cout << "Initial u: " << (*u_O)[i] << endl;


		//(*u_O)[i] = 0.35;

		//Gaussian displacement
		//(*u_O)[x] = exp(-1.0 * ((*x_Pos)[x]-xMid) * ((*x_Pos)[x]-xMid) / (sigma * sigma));
		
		/*//Start with Oscillon Solution
		(*u_O)[i] = sqrt(lambda) * Psi_0 * sqrt((1.0 + u) / (1.0 + u * cosh(2.0 * alpha * lambda * ((*x_Pos)[i] - xMid1) / g))) 
		 + sqrt(lambda) * Psi_0 * sqrt((1.0 + u) / (1.0 + u * cosh(2.0 * alpha * lambda * ((*x_Pos)[i] - xMid2) / g)))
		 + sqrt(lambda) * Psi_0 * sqrt((1.0 + u) / (1.0 + u * cosh(2.0 * alpha * lambda * ((*x_Pos)[i] - xMid3) / g)));
		*/
		
		//Write initial data to be used in analytical solution
		(*u_I)[i] = sqrt(lambda) * Psi_0 * sqrt((1.0 + u) / (1.0 + u * cosh(2.0 * alpha * lambda * ((*x_Pos)[i] - xMid1) / g))) 
		 + sqrt(lambda) * Psi_0 * sqrt((1.0 + u) / (1.0 + u * cosh(2.0 * alpha * lambda * ((*x_Pos)[i] - xMid2) / g)))
		 + sqrt(lambda) * Psi_0 * sqrt((1.0 + u) / (1.0 + u * cosh(2.0 * alpha * lambda * ((*x_Pos)[i] - xMid3) / g)));


		//Set vectors to 0 for all positions
		/*(*v_O)[i] = 0.0;
		(*u_N)[i] = 0.0;
		(*u_T)[i] = 0.0;
		(*v_N)[i] = 0.0;
		(*v_T)[i] = 0.0; 
		(*u_A)[i] = 0.0; //Analytical
		(*u_I)[i] = 0.0;*/

	}
	cout << "Successfully set up initial data! " << endl;

}

//Find u and v values at all positions over time
void Wave::propagate()
{

	//Count iterations
	int timeSteps = 0;


	//Initial values
	a_Fried = 1.0;
	a_Approx = 1.0;
	H = H_i;

	//Calculate conformal time
	dn = tInc / a;
	Doub t = 0;

	//Use aDouble to double grid resolution every time universe expands by a factor of 2
	gridFactor = 1;

	Doub realTime = 0;
	int numFiles = 0;

	//Doub rhoAvg = 0;

	//Iterate over all CONFORMAL time starting at second time step
	while(a_Approx < aMax){

		//rhoAvg = calcRhoAvg();


		//Evolve H(a) and a(n) with Append A
		//H = 1.0 / (3 * mPL * mPL) * rhoAvg;
		//a = a * (1.0 + a * H * dn);
		//H = H_i / (0.5 * H_i * t + 1.0)

		//Evolve H and a with 5.1 Implementation
		//Old evolution of a with H = a^-0.5
		a_Approx = (0.5 * H_i * t + 1.0) * (0.5 * H_i * t + 1.0);
		a_Fried = exp(H * t);
		//H = H_i / (0.5 * H_i * t + 1.0);

		//Calculuate H
		//setHubble(); 
		//Old method with phi_p
		if(HImplement == 1){

			H = H_i / (0.5 * H_i * t + 1.0);

		//Calculate H with Friedmann equation
		} else{
			Doub rhoAvg = calcRhoAvg();
			H = rhoAvg * m_Mpl * m_Mpl / (3 * sqrt(lambda));
		}


		//Compare H values
		cout << "H calculuated with rhoAvg: " << H << " Plebian H: " << H_i / (0.5 * H_i * t + 1.0) << endl;

		//Universe has expanded by a factor of 2
		if(a >= 2 * gridFactor){

			xInc = xInc / 2.0;
			gridFactor = 2 * gridFactor;
			cout << "Grid Factor is now: " << gridFactor << endl;

			//Reset x_Pos to accommodate double the grid points
			resetPos();
		}

		//Refine timestep
		tInc = courantCond * xInc;

		//Calculate conformal time
		dn = tInc / a_Approx;

		//Calculate real time
		realTime = realTime + tInc;

		//cout << "Time: " << t << " a: " << a << " tInc : " << tInc << " dn: " << dn << endl;

		//calculate s2 constant
		//s2 = (tInc / xInc) * (tInc / xInc);

		//Predict u and v
		predictor();

		//Correct u and v
		corrector();

		//analytical(realTime);

		//Write to file
		//Only write to file every (fileInc) files
		if(timeSteps % fileInc == 0){
		//if(timeSteps == 4828000){
		  //fileWrite(timeSteps,t);
		  numFiles ++;
		}

		timeSteps ++;
		t+= dn;
	       

		//Prevent from spitting out too many files
		if(numFiles > 1000){
			break;
		}

	}
	cout << "Time: " << t << " a_Approx: " << a_Approx << " a_Fried" << a_Fried << " tInc : " << tInc << " dn: " << dn << " H_i: " << H_i << endl;
	cout << "Time steps " << timeSteps << endl;
}

//Predicts u and v values with v_T and u_T values at n+1
void Wave::predictor()
{

	//a = (0.5 * H_i * t + 1.0) * (0.5 * H_i * t + 1.0);
	//a2 = a * a;
	//H = H_i / (0.5 * H_i * t + 1.0);

	Doub uRHS;
	Doub vRHS;
	Doub V_p;
	int i;

	//Run in parallel
	#pragma omp parallel num_threads(thread_Num)
	{
	#pragma omp for collapse(1) private(i,uRHS,vRHS,V_p)
	
	
	//Iterate through all space
	for(i = 0; i < (gridFactor * nx); i++){

		//Calculate potential;

		V_p = m*m * (*u_O)[i] - m*m * (*u_O)[i] * (*u_O)[i] * (*u_O)[i] + m*m * g*g/(lambda*lambda) * (*u_O)[i] * (*u_O)[i] * (*u_O)[i] * (*u_O)[i] * (*u_O)[i];
;
		
		//v = -du/dt
		uRHS = -1.0 * (*v_O)[i];

		//Predict u and next time t+dn
		(*u_T)[i] = (*u_O)[i] + (uRHS * dn);


		//Left boundary,
		if(i == 0) {

			//Taylor Series Inpterpolation
			//(*v_T)[x] = (*v_O)[x+1] * (1.0 - (xInc - dn) / xInc) + ((*v_O)[x] * (xInc - dn) / xInc);
			
			//vRHS = -1 HaV - d2u/dx2 + aV'(u)
			vRHS = -1.0 * H * a * (*v_O)[i] - ((*u_O)[i+1] - 2.0 * (*u_O)[i] + (*u_O)[nx-1]) / (xInc*xInc) + a*a * V_p;

			(*v_T)[i] = (*v_O)[i] + (vRHS * dn);

		//Right boundary, interpolate for v_T
		} else if(i == (gridFactor * nx-1)){
			
			//Taylor Series BC
			//(*v_T)[x] = (*v_O)[x-1] * (1.0 - (xInc - dn) / xInc) + ((*v_O)[x] * (xInc - dn) / xInc);

			//Periodic BC
			//vRHS = -1 HaV - d2u/dx2 + aV'(u)
			vRHS = -1.0 * H * a * (*v_O)[i] - ((*u_O)[0] - 2.0 * (*u_O)[i] + (*u_O)[i-1]) / (xInc*xInc) + a*a * V_p;

			(*v_T)[i] = (*v_O)[i] + (vRHS * dn);

		//Not on boundary
		} else {

				//vRHS = -1 HaV - d2u/dx2 + aV'(u)
				vRHS = -1.0 * H * a * (*v_O)[i] - ((*u_O)[i+1] - 2.0 * (*u_O)[i] + (*u_O)[i-1]) / (xInc*xInc) + a*a * V_p;

				(*v_T)[i] = (*v_O)[i] + (vRHS * dn);

		}

	}
	}

	
	
}

//Corrects v_T and u_T values to determine final u and v values at n+1
void Wave::corrector()
{
	/*a = (0.5 * H_i * t + 1.0) * (0.5 * H_i * t + 1.0);
	a2 = a * a;
	H = H_i / (0.5 * H_i * t + 1.0);*/

	//Doub uRHS;
	//Doub vRHS;
	Doub V_pUT;

	int i;

	//Run in parallel (uRHS is used section?)
	#pragma omp parallel num_threads(thread_Num)
	{
	#pragma omp for collapse(1) private(i,V_pUT)
	
	

	//Iterate through all space
	for(i= 0; i < (gridFactor * nx); i++){

		//Correct u
		/*uRHS = -0.5 * ((*v_O)[x] + (*v_T)[x]);

		(*u_N)[x] = (*u_O)[x] + uRHS * dn;*/

		//Condensed
		(*u_N)[i] = 0.5 * ((*u_O)[i] + (*u_T)[i]) - 0.5 * (*v_T)[i] * dn;

		//Calculate Potentials for u_O and u_T
		//V_pU = m*m * (*u_O)[x] - m*m * (*u_O)[x] * (*u_O)[x] * (*u_O)[x] + m * g*g/(lambda*lambda) * (*u_O)[x] * (*u_O)[x] * (*u_O)[x] * (*u_O)[x] * (*u_O)[x];
		V_pUT = m*m * (*u_T)[i] - m*m * (*u_T)[i] * (*u_T)[i] * (*u_T)[i] + m * g*g/(lambda*lambda) * (*u_T)[i] * (*u_T)[i] * (*u_T)[i] * (*u_T)[i] * (*u_T)[i];

		//Correct v 
		//Left boundary with 
		if(i == 0){
			//Taylor Expansion
			//(*v_N)[x] = (*v_O)[x+1] * (1.0 - (xInc - dn) / xInc) + ((*v_O)[x] * (xInc - dn) / xInc);

			//Periodic BC
			/*vRHS = - 0.5 * (H * a * ((*v_O)[x] + (*v_T)[x]) + ((*u_O)[x+1] - 2.0 * (*u_O)[x] + (*u_O)[nx-1] + (*u_T)[x+1] -2.0 * (*u_T)[x] + (*u_T)[nx-1]) / (xInc*xInc)
			 -1.0 * a*a * (V_pUT + V_pU));

			(*v_N)[x] = (*v_O)[x] + vRHS * dn;*/

			//Condensed (got rid of dn in second term)
			(*v_N)[i] = 0.5 * (*v_O)[i] + 0.5 * (*v_T)[i] - 0.5 * (H * a * (*v_T)[i] + ((*u_T)[i+1] -2.0 * (*u_T)[i] + (*u_T)[nx-1]) / (xInc*xInc) -1.0 * a*a * V_pUT) * dn;


		//Right boundary
		} else if(i == (gridFactor * nx-1)){
			//Taylor Expansion
			//(*v_N)[x] = (*v_O)[x-1] * (1.0 - (xInc - dn) / xInc) + ((*v_O)[x] * (xInc - dn) / xInc);
			
			//Periodic BC
			/*vRHS = - 0.5 * (H * a * ((*v_O)[x] + (*v_T)[x]) + ((*u_O)[0] - 2.0 * (*u_O)[x] + (*u_O)[x-1] + (*u_T)[0] -2.0 * (*u_T)[x] + (*u_T)[x-1]) / (xInc*xInc)
			 -1.0 * a*a * (V_pUT + V_pU));

			 (*v_N)[x] = (*v_O)[x] + vRHS * dn;*/

			//Condensed 
			(*v_N)[i] = 0.5 * (*v_O)[i] + 0.5 * (*v_T)[i] - 0.5 * (H * a * (*v_T)[i] + ((*u_T)[0] -2.0 * (*u_T)[i] + (*u_T)[i-1]) / (xInc*xInc) -1.0*a*a * V_pUT) * dn;

		//Not on a boundary 
		} else {

			//vRHS = -0.5 * (Hv_0 + Hv_T + [d^2u_T/dx^2 + d^2u_O/dx^2] -a^2 * (V'(u_T) -V'(u_O))
			/*vRHS = - 0.5 * (H * a * ((*v_O)[x] + (*v_T)[x]) + ((*u_O)[x+1] - 2.0 * (*u_O)[x] + (*u_O)[x-1] + (*u_T)[x+1] -2.0 * (*u_T)[x] + (*u_T)[x-1]) / (xInc*xInc)
			 -1.0 * a*a * (V_pUT + V_pU));
		

			(*v_N)[x] = (*v_O)[x] + vRHS * dn;*/

			//Condensed 
			(*v_N)[i] = 0.5 * (*v_O)[i] + 0.5 * (*v_T)[i] - 0.5 * (H * a * (*v_T)[i] + ((*u_T)[i+1] -2.0 * (*u_T)[i] + (*u_T)[i-1]) / (xInc*xInc) -1.0*a*a * V_pUT) * dn;

		}
	}
	}

	
	
	#pragma omp parallel num_threads(thread_Num)
	{
	#pragma omp for collapse(1) private(i)
	//#pragma omp parallel for collapse(1) private(i) num_threads(4)
	
	//Set old u and v equal to new u and v
	for(i = 0; i < (gridFactor * nx); i++){
		(*u_O)[i] = (*u_N)[i];
		(*v_O)[i] = (*v_N)[i];

	}
	}
	

}

void Wave::analytical(Doub realTime)
{
	int i;
	Doub w = m * sqrt(1.0 - (lambda * lambda * alpha * alpha/ (g*g)));

	string sol = "uAna";

//Run in parallel (uRHS is used section?)
	#pragma omp parallel num_threads(thread_Num)
	{
	#pragma omp for collapse(1) private(i)	

	//Iterate through all space
	for(i= 0; i < (gridFactor * nx); i++){

		(*u_A)[i] = (*u_I)[i] * cos(w * realTime);
	}
	}
}

Doub Wave::calcRhoAvg()
{	

	//Calculate Total energy from each position
	Doub rhoTotal = 0;
	Doub potential;
	int i;
	//#pragma omp parallel num_threads(4)
	//{
	//#pragma omp for collapse(1) private(i)

	for(i = 0; i < (gridFactor * nx); i++){
		
		potential = m*m * ( 0.5 * (*u_O)[i] * (*u_O)[i] - 1.0 / 4.0 * (*u_O)[i] * (*u_O)[i] * (*u_O)[i] * (*u_O)[i] 
					+ g*g / (6.0 * lambda * lambda) * (*u_O)[i] * (*u_O)[i] * (*u_O)[i] * (*u_O)[i] * (*u_O)[i] * (*u_O)[i]);

		cout << "Potential: " << potential << endl;
		cout << "1st order" << (*u_O)[i] << endl;
		
		//rho = 1/2 * m*m (d/dt(phi))^2 + V(phi) + 0.5 * (grad(phi))^2
		rhoTotal += 0.5 * m*m * (*v_O)[i] * (*v_O)[i] + potential 
					+ ((*u_O)[i+1] * (*u_O)[i+1] - 2.0 * (*u_O)[i+1] * (*u_O)[i-1] + (*u_O)[i-1]) / (16 * xInc * xInc);
		
	//}
	}

	Doub rhoAvg = rhoTotal / xInc;
	cout << rhoAvg << endl << endl << endl << endl;

	return rhoAvg;
}

void Wave::setHubble()
{
	Doub rhoAvg = 0.0;
	Doub rhoTotal = 0.0;
	Doub vPotential,phiActual,phiNext,vActual;
	Doub plankMass = m / (0.000005);

	int i;

	#pragma omp parallel num_threads(thread_Num)
	{
	#pragma omp for collapse(1) private(i,vPotential, phiActual, phiNext, vActual, rhoTotal)	
	//Calculate spatially averaged energy density
	for( i = 0; i < (gridFactor * nx); i++){

		phiActual = (*u_O)[i] / sqrt(lambda);

		vPotential = 0.5 * phiActual * phiActual - 0.25 * phiActual * phiActual * phiActual * phiActual
		+ g * g / (6.0 * lambda * lambda) * phiActual * phiActual * phiActual * phiActual * phiActual * phiActual;

		if(i==0){
			cout << "Potential: " << vPotential << endl;
		}

		phiNext = (*u_O)[i+1] / sqrt(lambda);

		vActual = (*v_O)[i] / sqrt(lambda);

		//On boundary (B.C. for gradient of u)
		if(i== (gridFactor * nx -1 )){

			// v = -du/dt
			rhoTotal = rhoTotal - 0.5 * m * vActual + m * vPotential + m / (2.0 * a * a) * ((*u_O)[0] / sqrt(lambda) - phiActual) / 2.0 * (phiNext - phiActual) / 2.0;
		
		//Not on boundary
		} else {
			rhoTotal = rhoTotal - 0.5 * m * vActual + m * vPotential + m / (2.0 * a * a) * (phiNext - phiActual) / 2.0 * (phiNext - phiActual) / 2.0;
		}
	}
	}

	//Calculate H
	rhoAvg = rhoTotal / xMax;

	H = rhoAvg / (3.0 * plankMass * plankMass);



}

void Wave::resetPos()
{
	int i;
	uTemp->resize(gridFactor * nx);
	vTemp->resize(gridFactor * nx);

	//Can't parallelize this because order matters!
	//Also, this only occurs 4 times at most (usually only once in my tests so far)

	for(i = 0; i < (gridFactor * nx); i++){

		//For each new grid point create another grid point with the same value right next door
		//Write vector starting at right end moving towards left end
		//(*u_O)[(prevGridFact * nx) -1 -i/2] = (*u_O)[nx * gridFactor -1];
		//(*u_O)[(prevGridFact * nx) -1-i] = (*u_O)[nx * gridFactor -2];

		//Interpolate to grid points in between known values

		//Even i
		if(i % 2 == 0){

			//Left Boundary
			if(i == 0){

				(*uTemp)[i] = 0.75 * (*u_O)[i/2] + 0.25 * (*u_O)[nx-1];
				(*vTemp)[i] = 0.75 * (*u_O)[i/2] + 0.25 * (*u_O)[nx-1];
			//Even and not on boundary
			} else {

				(*uTemp)[i] = 0.75 * (*u_O)[i/2] + 0.25 * (*u_O)[i/2 -1];
				(*vTemp)[i] = 0.75 * (*u_O)[i/2] + 0.25 * (*u_O)[i/2 -1];
			}
		//Odd i
		} else {

			//Right Boundary
			if(i == (gridFactor * nx -1)){
				cout << "Right Boundary: " << 0.75 * (*u_O)[i/2] + 0.25 * (*u_O)[0] << endl;
				(*uTemp)[i] = 0.75 * (*u_O)[i/2] + 0.25 * (*u_O)[0];
				(*vTemp)[i] = 0.75 * (*u_O)[i/2] + 0.25 * (*u_O)[0];
			//Odd and not on boundary
			} else {

				(*uTemp)[i] = 0.75 * (*u_O)[i/2] + 0.25 * (*u_O)[i/2 +1];
				(*vTemp)[i] = 0.75 * (*u_O)[i/2] + 0.25 * (*u_O)[i/2 +1];
			}
		}		
	}

	//Write temporary u and v to actual u and v

	//Resize vectors
	u_O->resize(gridFactor * nx);
	v_O->resize(gridFactor * nx);
	u_T->resize(gridFactor * nx);
	v_T->resize(gridFactor * nx);
	u_N->resize(gridFactor * nx);
	v_N->resize(gridFactor * nx);
	u_A->resize(gridFactor * nx);
	v_O->resize(gridFactor * nx);
	x_Pos->resize(gridFactor * nx);

	for(int i = 0; i < (gridFactor * nx); i++){

		(*x_Pos)[i] = ((Doub)i + 0.5) * xInc;
		(*u_O)[i] = (*uTemp)[i];
		(*v_O)[i] = (*vTemp)[i];
	}

}

void Wave::fileWrite(int timeStep, Doub t)
{

	//Write u to a file at various time steps
	int width = 16;

    ofstream outfile_u;
    ostringstream fileName_u;
    fileName_u << "u_dis_" << setfill('0') << setw(6) << timeStep << ".out" << ends;

    outfile_u.open(fileName_u.str().c_str());
    outfile_u.setf(ios::left);
    outfile_u << setw(width) << "#      x     " << setw(width) << " u(x,t) " << setw(width) << " u_Analytical at scale factor (a): " << a << endl;
    outfile_u << "#====================================================" << endl;

    //Iterate over all positions
    for(int i = 0; i < (gridFactor * nx); i++){

    	//Print sqrt(lambda) * Phi
		outfile_u << setprecision(8) << setw(width) << (*x_Pos)[i] << setw(width) << (*u_O)[i] << endl; //setw(width) << (*u_A)[i] << endl;
	}
	
}

void Wave::deconstructor()
{
	//Delete allocated memory for vectors
	delete u_N;
	delete x_Pos;
	delete u_O;
	delete u_T;
	delete v_N;
	delete v_O;
	delete v_T;
	delete u_A;
	delete u_I;
	delete uTemp;
	delete vTemp;
}

Wave::~Wave(){}


























