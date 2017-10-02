// Tell emacs that this is -*-c++-*- mode
//================================================
// Class containing all the integration stuff...
//================================================
#ifndef SW_H
#define SW_H

#include <omp.h>
#include "nr3.h"
#include "gridfunction.h"
#include "dumper.h"
#include "tridag.h"
#include <string>
#include <iostream>
#include <sstream> 
#include "ludcmp.h"

//const static int thread_num = 1;

//
class Oscillon {
private:

  //Initialize grid functions for u and v old, new, 1, and u analytical solutions
  gf3d u_old;
  gf3d v_old;
  gf3d u_1, v_1;
  gf3d u_new, v_new;
  gf3d u_ana;

  //Initialize L terms
  gf3d L1_old, L1_1;
  gf3d L2_old, L2_1, L2_new;

  //Initialize vector of values
  VecDoub r, theta, x, sintheta, phi; //x is cos(theta)

  //Number of grid points in each coord
  int N_r, N_theta, N_phi;

  int gridFactor; //Resolution of grid

  //Increments and time
  Doub dr, dx, dphi, t, dt, a, dn;

  //Parameters
  Doub lambda, g, m, alphaF, m_mpl, H_i;

  //Create dumper class
  dumper *dump;

  // for initial data and analytical solution
  Doub psi_init, A;  
  Doub PI;
  //
  //================================================
  // Methods
  //================================================
  //
public:
  //================================================
  // Constructor
  //
  // Important note concerning grid: input are number of gridpoints
  // covering PHYSICAL space, not including "ghost zones".  All gridfunctions
  // have one ghost zone on each side, so we increase all gridnumbers by two.
  //
  //================================================
  Oscillon(int N_r_i, int N_theta_i, int N_phi_i, 
    Doub r_max, Doub Courant, dumper * dump_i, int type) : 

    //Creates dumper class with two extra points for ghost zones
    dump(dump_i), t(0.0),r(N_r_i+2),theta(N_theta_i+2),sintheta(N_theta_i+2),
    x(N_theta_i+2), phi(N_phi_i+2)
  {
    cout << " Constructing ScalarWave ... " << endl;
    PI = acos(-1.0); 

    //================================================
    // Allocate and set up coordinate vectors
    //================================================

    //Extra two points for ghost zones
    N_r     = N_r_i     + 2;
    N_theta = N_theta_i + 2;
    N_phi   = N_phi_i   + 2;
    Setup_Grid(r_max,Courant);  
    //================================================
    // Set up 3D grid functions
    //================================================
    u_old.setup(N_r, N_theta, N_phi, &r, &theta, &x, &phi);
    u_new.setup(N_r, N_theta, N_phi, &r, &theta, &x, &phi);
    u_1.setup(N_r, N_theta, N_phi, &r, &theta, &x, &phi);

    v_old.setup(N_r, N_theta, N_phi, &r, &theta, &x, &phi);
    v_new.setup(N_r, N_theta, N_phi, &r, &theta, &x, &phi);
    v_1.setup(N_r, N_theta, N_phi, &r, &theta, &x, &phi);

    u_ana.setup(N_r, N_theta, N_phi, &r, &theta, &x, &phi);

    L1_old.setup(N_r, N_theta, N_phi, &r, &theta, &x, &phi);
    L1_1.setup(N_r, N_theta, N_phi, &r, &theta, &x, &phi);

    L2_old.setup(N_r, N_theta, N_phi, &r, &theta, &x, &phi);    
    L2_new.setup(N_r, N_theta, N_phi, &r, &theta, &x, &phi);
    L2_1.setup(N_r, N_theta, N_phi, &r, &theta, &x, &phi);

  };
  //================================================
  // Destructor
  //================================================
  ~Oscillon() {}
  //================================================
  // Set up cell_centered grid (for now equidistant, but can be changed!)
  //================================================
  int Setup_Grid(double r_max, double Courant) {
    cout << " setting up grid ..." << endl;
    //
    // radial 
    //
    dr = r_max/double(N_r-2);


    //Parallelize
    int i;
    #pragma omp parallel num_threads(thread_num)
    {
    #pragma omp for collapse(1) private(i)
    //Iterate over all r
    for (int i = 0; i < N_r; i++)

      //Create cell-centered grid with ghost points on ends
      r[i] = (i - 0.5)*dr;
    }
    //
    // theta: keep track for minimum dtheta
    //
    Doub dtheta_min = PI;
    // dx = 2.0/double(N_theta-2);
    // for (int j = 1; j < N_theta-1; j++) {
    //   x[j] = 1.0 - (j - 0.5)*dx;
    //   theta[j] = acos(x[j]);
    //   sintheta[j] = sin(theta[j]);
    //   if (j > 1)
    //  if ( (theta[j] - theta[j-1]) < dtheta_min ) dtheta_min = (theta[j] - theta[j-1]);
    // }
    // x[0] = x[1];
    // theta[0] = - theta[1];
    // sintheta[0] = - sintheta[1];
    // x[N_theta-1] = x[N_theta-2];
    // theta[N_theta-1] = PI + theta[1];
    // sintheta[N_theta-1] = sintheta[0];

    //Increment of theta
    Doub dtheta = PI/double(N_theta - 2);

    //Parallelize
    int j;
    #pragma omp parallel num_threads(thread_num)
    {
    #pragma omp for collapse(1) private(j)

    //Not sure if I can parallelize this yet???
    //iterate over all theta
    for (int j = 0; j < N_theta; j++) {

      //Initialize to actual values
      theta[j] = (j - 0.5)*dtheta;

      sintheta[j] = sin(theta[j]);
      x[j] = cos(theta[j]);
    }
    }
    


    dtheta_min = dtheta;
    // for (int j = 0; j < N_theta; j++) {
    //   cout << " x[" << j << "] = " << x[j] << "; theta = " 
    //     << theta[j] << "; sintheta = " << sintheta[j] << endl;
    // }
    //
    // phi
    //

    //Phi increment
    dphi = 2.0*PI/double(N_phi-2);

    //Parallelize
    int k;
    #pragma omp parallel num_threads(thread_num)
    {
    #pragma omp for collapse(1) private(k)
    //Iterate over all phi
    for (int k = 0; k < N_phi; k++)
      phi[k] = (k - 0.5)*dphi;

    }
    //
    // dt 
    //
    // find minimum physical grid spacing
    Doub dtheta_x = r[1] * dtheta_min;  //r_0.5 * theta increment = dr * dtheta ??

    Doub dphi_x = r[1] * sintheta[1] * dphi; //dr * dtheta * dphi
    Doub dx_min = (dtheta_x < dphi_x) ? dtheta_x : dphi_x; //???
    dt = Courant * dx_min;
    cout << " for dr = " << dr << " and dt = " << dt << endl;
    return 0;
  };
  //================================================
  // Initialize
  //================================================
  int Initialize(double lambda_In, double m_In, double alphaF_In) 
  { 
    cout << " Setting up initial data ..." << endl;

    m = m_In;
    lambda = lambda_In;
    alphaF = alphaF_In;
    cout << alphaF << endl;

    //g = lambda / sqrt(0.2);
    g = 5.0;
    cout << "g: " << g << endl;

    //Changeable parameter
    Doub m_Mpl = 0.000005;

    H_i = sqrt(lambda / (10 * g*g)) * m * m_Mpl;
    cout << "H_i: " << H_i << endl;
    
    //Get analytical solution
    AnalyticalSolution();

    //Used just to determine if functions run or there is a problem
    int error = 0;

    //Parallelize
    int i,j,k;
    #pragma omp parallel num_threads(thread_num)
    {
    #pragma omp for collapse(3) private(k,j,i)

    //Iterate through all phi
    for (k = 0; k < N_phi; k++) 

      //Iterate through all theta
      for (j = 0; j < N_theta; j++) 

         //Iterate through all r
         for (i = 0; i < N_r; i++) {

              //Use analytical function values as initial u
             u_old[i][j][k] = u_ana[i][j][k];
             v_old[i][j][k] = 0.0;
         }
    }


    t = 0.0;
    cout << " ... done!" << endl;
    return error;
  }
  //================================================
  // Integrate
  //================================================

  //Beef of the code
  int Integrate(double aMax)
  {

    int timestep = 0;

    //Initial values
      Doub aTild = 1.0; //Prediction of a
      Doub aDot, aDotTild; //first derivative and firstDeriv Prediction 
      Doub a_Old = 1.0; //Stores the previous a
      Doub H = H_i;

      //Calculate conformal time (don't actually need?)
      dn = dt / a;

      //Use aDouble to double grid resolution every time universe expands by a factor of 2
      gridFactor = 1;

      Doub realTime = 0.0;
      int numFiles = 0;

    //Write time, u numerical(initial) and u analytical to a file
    dump->dump(a,timestep,&u_old,&u_ana);

    //Iterate until maximum time is reached
    cout << "aMax: " << endl;

    while (a < aMax) {

      //Predictor
      aDot = H_i * pow(a,1.5);

      aTild = a + aDot * dn;

      //Corrector
      aDotTild = H_i * pow(aTild,1.5);

      a = a + 0.5 * (aDot + aDotTild) * dn;

      //H = H_i * a ^ -0.5
      H = H_i / sqrt(a);

      dn = dt / a;

      t += dn;
      timestep++; 

      //Get numerical solution for each time
      PIRK_TimeStep();

      //Get analytical solution for each time
      AnalyticalSolution();

      //Write numerical and analytical solution to a file if time is a multiple of desired timestep
      dump->cond_dump(a,timestep,&u_old,&u_ana);
    }
  }



private:
  //================================================
  // Carries out one PIRK time step
  //================================================
  int PIRK_TimeStep()
  {
    //
    // Predictor step
    //
    // update u

    //Gets L1 = v = du/dt values (I think we can make this faster)
    Compute_L1(u_old,v_old,L1_old);

    //Parallelize
    int i,j,k;
    #pragma omp parallel num_threads(thread_num)
    {
    #pragma omp for collapse(3) private(i,j,k)
    //Iterate through r, t, p (without ghost points)
    for (i = 1; i < N_r-1; i++)    
      for (j = 1; j < N_theta-1; j++)
         for (k = 1; k < N_phi-1; k++) { 

            //Predictor step (u_T = u_O + u_RHS * dT)   
           u_1[i][j][k] = u_old[i][j][k] + dn*L1_old[i][j][k];
    };
    }

    //Fill ghost points for predicted u
    u_1.fill_ghosts();
  
    //Sets u on outer boundaries (at N_r - 1)
    OuterBoundary(u_1,u_old);

    // update v, using updated values of u

    //Find d2u/dt2 for u_old and u_T
    Compute_L2(u_old,L2_old);
    Compute_L2(u_1,L2_1);


    #pragma omp parallel num_threads(thread_num)
    {
    #pragma omp for collapse(3) private(i,j,k)
    //Predictor for v
    //Iterate over r, t, p (not including ghost points)
    for (i = 1; i < N_r-1; i++)    
      for (j = 1; j < N_theta-1; j++)
         for (k = 1; k < N_phi-1; k++) { 

          //v_T = v_O + v_RHS * dt 
           v_1[i][j][k] = v_old[i][j][k] + 0.5*dn*(L2_old[i][j][k] + L2_1[i][j][k]);
    };
    }

    //Fill ghost points for v prediction
    v_1.fill_ghosts();
    OuterBoundary(v_1,v_old);

    //
    // Corrector step
    //
    // update u

    //Set v = L1_1
    Compute_L1(u_1,v_1,L1_1);

    #pragma omp parallel num_threads(thread_num)
    {
    #pragma omp for collapse(3) private(i,j,k)
    //Iterate over all r, t, p
    for (int i = 1; i < N_r-1; i++)    
      for (int j = 1; j < N_theta-1; j++)
         for (int k = 1; k < N_phi-1; k++) {    

            //u_N = u_O + u_RHS * dt
            u_new[i][j][k] = u_old[i][j][k] + 0.5*dn*(L1_old[i][j][k] + L1_1[i][j][k]);
     };
    }

    //Fill ghost points and calculate outer boundary
    u_new.fill_ghosts();
    OuterBoundary(u_new,u_old);

    // update v

    //Find d2u/dt2 for u_N
    Compute_L2(u_new,L2_new);

    #pragma omp parallel num_threads(thread_num)
    {
    #pragma omp for collapse(3) private(i,j,k)

    //Iterate over r,t,p
    for (i = 1; i < N_r-1; i++)    
      for (j = 1; j < N_theta-1; j++)
         for (k = 1; k < N_phi-1; k++) {   

          //v_N = v_O + v_RHS * dt 
           v_new[i][j][k] = v_old[i][j][k] + 0.5*dn*(L2_old[i][j][k] + L2_new[i][j][k]);
    };
    }

    //Fill ghost points and calculate outer boundary
    v_new.fill_ghosts();
    OuterBoundary(v_new,v_old);

    //
    // Copy results
    //

    #pragma omp parallel num_threads(thread_num)
    {
    #pragma omp for collapse(3) private(i,j,k)

    //Set old u and v to new u and v
    for (i = 0; i < N_r; i++) 
      for (j = 0; j < N_theta; j++)
         for (k = 0; k < N_phi; k++) {
           u_old[i][j][k] = u_new[i][j][k];
           v_old[i][j][k] = v_new[i][j][k];
  };
  }
}
  //================================================
  // Computes L1 terms for PIRK time steps
  //================================================

  //Seems like all we need is v here?
  int Compute_L1(gf3d & u, gf3d & v, gf3d  & L1)
  {
    int i = 0;
    int j = 0;
    int k = 0;

    #pragma omp parallel num_threads(thread_num)
    {
    #pragma omp for private(i,j,k)
    //Iterate through all r, t, p (not including ghost points)
    for (i = 1; i < N_r-1; i++)    
      for (j = 1; j < N_theta-1; j++)
        for (k = 1; k < N_phi-1; k++) { 
           //k = pos % (N_phi-1);
           //j = pos / (N_phi-1);
            //L1 = v   
           L1[i][j][k] = v[i][j][k];
       };
    }
    return 0;
  }
  //================================================
  // Computes L2 terms for PIRK time steps
  //================================================
  int Compute_L2(gf3d & u, gf3d & L2)
  {
    //
    // interior of grid
    //

    int i,j,k;
    Doub r_l, sintheta_l, x_l, u_ddr, u_dr, u_ddtheta, u_dtheta, u_ddphi, V_Prime;
    #pragma omp parallel
    {
    #pragma omp for collapse(3) private(i,j,k,r_l,sintheta_l,x_l,u_ddr,u_dr,u_ddtheta,u_dtheta,u_ddphi,V_Prime)
    //Iterate through all r, t, p (not including ghost points)
    for (i = 1; i < N_r-1; i++)    
      for (j = 1; j < N_theta-1; j++)
         for (k = 1; k < N_phi-1; k++) {  
            r_l = r[i];
            x_l = x[j];
            sintheta_l = sintheta[j];
            V_Prime = u.potentialPrime(i,j,k,g,lambda);


            Doub u_ddr, u_dr, u_ddtheta, u_dtheta, u_ddphi;
            u.derivs(i, j, k, u_ddr, u_dr, u_ddtheta, u_dtheta, u_ddphi);

            L2[i][j][k]  = u_ddr + 2.0*u_dr / r_l
            + (u_ddtheta + x_l / sintheta_l * u_dtheta
              + u_ddphi / (sintheta_l * sintheta_l)) / (r_l*r_l)
                + a*a* V_Prime;
       }; 
    }
    //
    // outer boundary will be handled in Timestep...
    //
  }

  //================================================
  // Updates outer boundary for psi, using psiold
  //================================================
  int OuterBoundary(gf3d & psi, gf3d & psiold)
  {
    Doub Courant = dn/dr;
    // use linear interpolation to find psi at intersection
    // of radially outgoing characteristic with previous timeslice

    int i,j,k;
    Doub psi_last, r_last;
    #pragma omp parallel num_threads(thread_num)
    {
    #pragma omp for private(i,j,k, psi_last, r_last)
    //Iterate over theta and phi
    for (j = 0; j < N_theta; j++) 
      for (k = 0; k < N_phi; k++) {   

        //k = pos % N_phi;
        //j = pos / N_phi;

        //Do I need to account for conformal time in Courant condition??????
        //????? right now it changes with time
         //Linear combination of psi at N_r-2 and N_r-1 weighted with courant cond
         psi_last = Courant*psiold[N_r-2][j][k] + (1.0 - Courant)*psiold[N_r-1][j][k];

         //Lin combo of r at N_r-2 and N_r-1 weighted with courant cond
         r_last   = Courant*r[N_r-2]            + (1.0 - Courant)*r[N_r-1];

         //Set psi on boundary ?????
         psi[N_r-1][j][k] = (r_last/r[N_r-1]) * psi_last;
      }
    }
  }
  //================================================
  // Find analytical solution for wave at current time
  //================================================
  int AnalyticalSolution()
  {
    //3D Oscillon Initial Solution

    //Amplitude of Cos term
    g = 7.0;
    alphaF = 0.9;
    Doub alphaC = sqrt(27.0/160.0);
    cout << alphaF << endl;
    Doub alpha = alphaF * alphaC;
    //Doub w = sqrt(m*m - (lambda*lambda / (g*g) * alpha*alpha));
    //Doub amp = m / sqrt(lambda) * lambda / g * PhiCrit;

    //From 1D analytical solution
    Doub u = sqrt(1.0 - (alpha / alphaC) * (alpha / alphaC));
    //Doub Psi_0 = sqrt(9.0 * lambda * (1.0 - u) / (10 * g*g));
    Doub phiCrit = sqrt(9.0/10.0);
    Doub phi_0 = 1.1 * phiCrit * sqrt(1.0 - u);
    //Doub w = m*m * (1.0 - (lambda/g) * (lambda/g) * alpha*alpha);

    VecDoub main_Diag(N_r), up_Diag(N_r), down_Diag(N_r);
    VecDoub phiInit(N_r), dphiInit(N_r),s(N_r);

    VecDoub p(N_r); // change of variables from r to p
    Doub dp = dr / sqrt(g);
    Doub rhoSurf = 8.0;


  
  //Don't use ghost points for rho
  for(int i=0; i < N_r; i++){

  	p[i] = (i + 0.5) * dp;


    //Make guess at initial phi using phi_0
    if (p[i] < rhoSurf){
      phiInit[i] = 1.0-p[i]*p[i]/500.0;
    } else{
      phiInit[i] = 6.976/p[i];
    }
    cout << p[i] << "   " << phiInit[i] << endl;
    

    //Initialize dphiInit
    dphiInit[i] = 0.0;

  }

  Doub dp2 = dp*dp;
  Doub rho = 1.0;
  Doub pi = acos(-1.0);
  Doub source = 4*pi*rho;

  Doub residNorm = 1.0;
  Doub residual = 0.0;
  Doub tolerance = 0.0001;
  int numIts = 0;

  //Set phi0
  //phiInit[0] = 1.01*sqrt(9.0/10.0);
  //dphiInit[0] = 0;

  //Iterate until norm of residual is small enough
  while((residNorm > tolerance) && (numIts < 100)){

  	 numIts ++;

     //Reset norm
    residNorm = 0;

  	  //Calculate residual
    for(int i=0; i < N_r; i++){
  	   //On left boundary
      if(i==0){
        //Derivative = 0 at p=0
        residual = (phiInit[i+1] - phiInit[i]) / (dp2) + 2.0 / p[i] *  (phiInit[i+1] - phiInit[i]) / (2.0 * dp)
                   -alpha*alpha * phiInit[i] + 3.0/4.0 * phiInit[i]*phiInit[i]*phiInit[i] - 5.0/8.0 * phiInit[i]*phiInit[i]*phiInit[i]*phiInit[i]*phiInit[i];

        //phi goes as 1/p * exp(-alpha*p) and we calculate derivatives from this
      } else if(i==(N_r-1)){

          residual = 0.0;//(-1.0 * phiInit[i] + phiInit[i-1]) / (dp2) + 2.0 / p[i] *  (phiInit[i] - phiInit[i-1]) / (2.0 * dp);
          //- alpha*alpha * phiInit[i] + 3.0/4.0 * phiInit[i]*phiInit[i]*phiInit[i] - 5.0/8.0 * phiInit[i]*phiInit[i]*phiInit[i]*phiInit[i]*phiInit[i];


          //exp(-1.0*alpha * p[i]) / (alpha*alpha * p[i]*p[i]) - alpha*alpha / p[i] * exp(-1.0*alpha) 
            //          + 3.0/4.0 * 1.0/(p[i]*p[i]*p[i]) * exp(-3.0*alpha*p[i]) -5.0/8.0 * 1.0/(p[i]*p[i]*p[i]*p[i]*p[i])*exp(-5.0*alpha*p[i]);


          //Not on a boundary
      } else {

          residual = (phiInit[i+1] - 2.0*phiInit[i] + phiInit[i-1]) / (dp2) + 2.0 / p[i] *  (phiInit[i+1] - phiInit[i-1]) / (2.0 * dp)
              - alpha*alpha * phiInit[i] + 3.0/4.0 * phiInit[i]*phiInit[i]*phiInit[i] - 5.0/8.0 * phiInit[i]*phiInit[i]*phiInit[i]*phiInit[i]*phiInit[i];
           

      	  /* //poisson = 1
      	   if(p[i] < rSurf){

      		    residual = (phiInit[i+1] - 2.0*phiInit[i] + phiInit[i-1]) / (dp2) + 2.0 / p[i] *  (phiInit[i+1] - phiInit[i-1]) / (2.0 * dp) - source;
      	   //poisson = 0 
      	   } else {

        	   residual = (phiInit[i+1] - 2.0*phiInit[i] + phiInit[i-1]) / (dp2) + 2.0 / p[i] *  (phiInit[i+1] - phiInit[i-1]) / (2.0 * dp);
              //- alpha*alpha * phiInit[i] + 3.0/4.0 * phiInit[i]*phiInit[i]*phiInit[i] - 5.0/8.0 * phiInit[i]*phiInit[i]*phiInit[i]*phiInit[i]*phiInit[i];
    	     }*/

      }
      s[i] = -1.0*residual;
      residNorm += abs(residual);

    }

    cout << "residNorm" << residNorm << endl;


    for(int i=0; i < N_r; i++){
  //Boundary conidtions for main diagonal
      if (i==0) {
        main_Diag[i] = (-1.0/dp2) - (1.0 / (p[i] * dp)) - alpha*alpha + 9.0/4.0 * phiInit[i]*phiInit[i] - 25.0/8.0 * phiInit[i]*phiInit[i]*phiInit[i]*phiInit[i];
        up_Diag[i] = (1.0/dp2) + (1.0 / (p[i] * dp));

      } else if(i == (N_r - 1)) {
        main_Diag[i] = 1.0;//(-1.0/dp2) + (1.0 / (p[i] * dp)) - alpha*alpha + 9.0/4.0 * phiInit[i]*phiInit[i] - 25.0/8.0 * phiInit[i]*phiInit[i]*phiInit[i]*phiInit[i];
        //exp(-1.0*alpha * p[i]) / (alpha*alpha * p[i]*p[i]) - alpha*alpha / p[i] * exp(-1.0*alpha) 
                  //    + 9.0/4.0 * 1.0/(p[i]*p[i]*p[i]) * exp(-3.0*alpha*p[i]) -25.0/8.0 * 1.0/(p[i]*p[i]*p[i]*p[i]*p[i])*exp(-5.0*alpha*p[i]);
        //(-1.0/dp2) + (1.0 / (p[i] * dp))- alpha*alpha + 9.0/4.0 * phiInit[i]*phiInit[i] - 25.0/8.0 * phiInit[i]*phiInit[i]*phiInit[i]*phiInit[i];
        down_Diag[i] = 0.0;//(1.0/dp2) - (1.0 / (p[i] * dp));

      } else {
      	//Set diagonal vectors
        up_Diag[i] = (1.0/dp2) + (1.0 / (p[i] * dp));
        down_Diag[i] = (1.0/dp2) - (1.0 / (p[i] * dp));
        main_Diag[i] = (-2.0 / dp2) - alpha*alpha + 9.0/4.0 * phiInit[i]*phiInit[i] - 25.0/8.0 * phiInit[i]*phiInit[i]*phiInit[i]*phiInit[i];
      }
    }

  	//Use the tridag method to solve for y
    tridag(down_Diag,main_Diag,up_Diag,s,dphiInit);

  	//Set new phi values
    for(int i=0; i < N_r; i++){
      phiInit[i] += dphiInit[i];
    }


     /*//Write residual values to a file
    ofstream outfile;
    ostringstream rayfilename;

    //Convert numIts to a string to print out file
    stringstream ss;
    ss << numIts;
    string sNumIts = ss.str();
    
    rayfilename << "_residual_" << sNumIts << setfill('0') << setw(4) << ends;
    outfile.open(rayfilename.str().c_str());
    //    cout << " ... dumping file " << filename.str() << endl;
    outfile.setf(ios::left);
    outfile << "# Data at time " << time << endl;
    outfile << "# " << setw(14) << "p" << setw(16) << "residual" << endl;
    outfile << "#==================================================" << endl;
    outfile.setf(ios::right);

    for (int i = 0; i < N_r; i++) {
      outfile << setprecision(8) << setw(16) << p[i] << setw(16) << s[i]  << endl; 
    }

    outfile.close();*/

    cout << "Number of Iterations " << numIts << endl;

  } //End of while loop, residual is now below norm

  for(int i=0; i < N_r; i++){

  	//cout << phiInit[i] << endl;
  }
  cout << "Writing to file!" << endl;

  //Write values to a file
    ofstream outfile;
    ostringstream rayfilename;
    rayfilename << "_raysInit_" << setfill('0') << setw(4) << ends;
    outfile.open(rayfilename.str().c_str());
    //    cout << " ... dumping file " << filename.str() << endl;
    outfile.setf(ios::left);
    outfile << "# Data at time " << time << endl;
    outfile << "# " << setw(14) << "r" << setw(16) << "psi" << endl;
    outfile << "#==================================================" << endl;
    outfile.setf(ios::right);

    for (int i = 0; i < N_r; i++) {
      outfile << setprecision(8) << setw(16) << p[i] << setw(16) << phiInit[i]  << endl; 
    }
    outfile.close();





  int i,j,k;
    #pragma omp parallel num_threads(thread_num)
    {

      #pragma omp for collapse(3) private(i,j,k)

      //Iterate over all r
      for (i = 0; i < N_r; i++) {

        //Iterate over all theta and phi
         for (j = 0; j < N_theta; j++) {
             for (k = 0; k < N_phi; k++) {

                  u_ana[i][j][k] = phiInit[i] / sqrt(g);

             }
          }
      }
    } //end openmp
  };
  //================================================
  // Functions for scalar wave
  //================================================

  //Create initial analytical psi function
  //Why is there an x?
  inline Doub F(Doub x, Doub psi0, Doub A) {
    return psi0 * x * exp(- x*x/(A*A));
  };

  //Create initial psi function derivative

  inline Doub dFdx(Doub x, Doub psi0, Doub A) {

    //Why no a^2 in -2x^2 term??
    return psi0 * (1.0 - 2.0 * x * x) * exp(- x*x/(A*A));
  }; 
};


#endif  /* SW_H */
