// Tell emacs that this is -*-c++-*- mode
//================================================
// Class containing types for grid functions 
//
// Based on Mat3DDoub in numerical recipes, but includes
// functionality for spherical coordinate systems.
//
// Explain grid...
//================================================
#ifndef GF_H
#define GF_H

const static int thread_num = 1;

class gf3d {
 private:

  //Numer of grid points
  int nr;
  int nt;
  int np;
  Doub ***v;
  VecDoub *r_p, *theta_p, *x_p, *phi_p;
  Doub delta_phi;  // NOTE: assumes equidistant grid in phi!!
 public:
  //================================================
  // Constructor
  //================================================

  //Sets number of grid points to 0
  gf3d(): nr(0), nt(0), np(0), v(NULL)
    //  , r(NULL), theta(NULL), x(NULL), phi(NULL) 
  { // cout << " in default gridfunction constructor " << endl; 
  };
  //================================================
  // Set up grid functions
  //================================================

  //n,m,k are numbers of grid points in r, t, p directions
  int setup(int n, int m, int k, 
      VecDoub * r_i, VecDoub * theta_i, VecDoub * x_i, VecDoub * phi_i) 
  {
    //Assigns number of grid points from input to class variables
    nr = n;
    nt = m;
    np = k;

    v = new Doub**[n]; //????

    //Assigns r, theta and phi inputs to class vectors 
    r_p = r_i;
    theta_p = theta_i;
    x_p = x_i;
    phi_p = phi_i;
    int i,j;

    //????
    v[0] = new Doub*[n*m]; //Holds r and theta values in a 1D, specially ordered vector
    v[0][0] = new Doub[n*m*k]; //Holds r, theta, and phi values in a 2D, specially ordered matrix

    //Iterate through theta values for first row of matrix (r=0)
    //Neglects theta=0 (singular)


    //Can't parallelize normally
    for(j=1; j<m; j++){

      v[0][j] = v[0][j-1] + k;
    }

    //Iterate through r values for 1D vector and 2D vector for first column (theta =0)
    //Neglects r=0 (singular)

    //Can't parallelize normally
    for(i=1; i<n; i++) {
        v[i] = v[i-1] + m;
        v[i][0] = v[i-1][0] + m*k;

        //Iterate through theta values again for each r
        for(j=1; j<m; j++){
           v[i][j] = v[i][j-1] + k;
        }
    }
    delta_phi = (*phi_p)[1] - (*phi_p)[0];

    // sanity check if delta phi is constant
    if (delta_phi != (*phi_p)[2] - (*phi_p)[1]) 
      cout << " ERROR in gf3d: dphi not constant ! " << endl;
  }
  //End of Setup function (This indentention is funky)

  inline Doub** operator[](const int i) //subscripting: pointer to row i
  {
    return v[i];
  }
  inline const Doub* const * operator[](const int i) const //????
  {
    return v[i];
  }
  //================================================
  // copy operator
  //================================================
  gf3d & operator=(const gf3d & rhs) {
    cout << " in copy operator! " << endl;

    //Get dimensions of grid
    nr = rhs.dim1();
    nt = rhs.dim2();
    np = rhs.dim3();
    int n = nr;
    int m = nt;
    int k = np;

    //Delete all v vectors/matrices
    if (v!=NULL) { delete [] (v[0][0]); delete [] (v[0]); delete [] (v); }

    //Initialize v vectors/matrices
    v = new Doub**[n];
    v[0] = new Doub*[n*m];
    v[0][0] = new Doub[n*m*k];
    int i,j;

    //Fill v vector/matrix as done above in setup

    //Can't parallelize normally
    for(j=1; j<m; j++) v[0][j] = v[0][j-1] + k;
        for(i=1; i<n; i++) {

          v[i] = v[i-1] + m;
          v[i][0] = v[i-1][0] + m*k;

          for(j=1; j<m; j++) v[i][j] = v[i][j-1] + k;
    }

    //Set pointers for r,t,p,x
    r_p = rhs.r_pointer();
    theta_p = rhs.theta_pointer();
    x_p = rhs.x_pointer();
    phi_p = rhs.phi_pointer();

    delta_phi = rhs.dphi();

    //Parallelize
    /*#pragma omp parallel num_threads(thread_num)
    {
    #pragma omp for collapse(3) private(i,j,k)*/

    //For each value in the grid set rhs = v
    for (i = 0; i < nr; i++) 
      for (j = 0; j < nt; j++)
         for (k = 0; k < np; k++)
             v[i][j][k] = rhs[i][j][k];
    //}
    return *this; // Prof Baum....We're dissapointed
  }

  //Inline functions for dimensions of grid
  inline int dim1() const { return nr; }
  inline int dim2() const { return nt; }
  inline int dim3() const { return np; }

  //Inline functions for pointers to r,t,p,x
  inline VecDoub * r_pointer() const { return r_p;}
  inline VecDoub * theta_pointer() const { return theta_p;}
  inline VecDoub * x_pointer() const { return x_p;}
  inline VecDoub * phi_pointer() const { return phi_p;}

  //Inline functions for getting a value at an r,t,x,p grid index
  inline Doub r(int i) { return (*r_p)[i]; }
  inline Doub theta(int j) { return (*theta_p)[j]; }
  inline Doub x(int j) { return (*x_p)[j]; }
  inline Doub phi(int k) { return (*phi_p)[k]; }
  inline Doub dphi() const { return delta_phi; }


  //================================================
  // Fill inner ghost zones, assuming a spherical grid.
  //================================================
  int fill_ghosts() {
    //
    // start with lower r:
    // 
    //Parallelize
    int i,j,k,J,K;
    /*#pragma omp parallel num_threads(thread_num)
    {
    #pragma omp for collapse(2) private(i,j,J,K)*/
    //Iterate through all theta and phi points
    for (j = 1; j < nt - 1; j++)
      for (k = 1; k < np - 1; k++) {

        //Set super indices to be the value (theta - pi, phi - pi) across from the position
         J = nt - j - 1;
         K = 1 + (k + (np - 4)/2) % (np - 2);

         //Set first r value (r=0) to the r value above it at the super index
         v[0][j][k] = v[1][J][K];
      }
    //}


    //
    // now both lower and upper theta:
    //

    //Parallelize

    /*#pragma omp parallel num_threads(thread_num)
    {
    #pragma omp for collapse(2) private(i,k,K)*/
    //Iterate over all r, phi values (on z axis (theta=0))
    for (i = 0; i < nr; i++)
      for (k = 1; k < np - 1; k++) {

         //Set super index for phi value across from position
         K = 1 + (k + (np - 4)/2) % (np - 2);
         v[i][0][k]    = v[i][1][K];
         v[i][nt-1][k] = v[i][nt-2][K];
      }
    //}
    //
    // finally both lower and upper phi:
    //

    //Why do we need to do this??? We can set phi to be 0 in the Laplace Eq...
    //Ohhh, its just because we need to make sure phi values are continues at 0-2pi boundary
    //Hopefully this equation is azimuthally symmetric though??

    //Parallelize
    /*#pragma omp parallel num_threads(thread_num)
    {
    #pragma omp for collapse(2) private(i,j)*/
    //Iterate over all r and theta values
    for (int i = 0; i < nr; i++)
      for (int j = 0; j < nt; j++) {

         //Not using super indices
         v[i][j][0]    = v[i][j][np-2];
         v[i][j][np-1] = v[i][j][1];
      } 
    //}
  };

  //================================================
  // Compute derivative of potential for given coordinate
  //================================================
  int potentialPrime(int i,int j,int k, Doub g, Doub lambda)
  {
    //potential' = u - u^3 + g^2/lambda^2 * u^5
    return v[i][j][k] - v[i][j][k]*v[i][j][k]*v[i][j][k] + (g*g) / (lambda*lambda) * v[i][j][k]*v[i][j][k]*v[i][j][k]*v[i][j][k]*v[i][j][k];
  }


  //================================================
  // Derivative operators
  //================================================
  //Uses finite differences

  //
  // *MARK- function modifies all derivatives at once to eliminate calls to operator
  //
  inline int derivs(int i, int j, int k, Doub & ddr, Doub & dr, Doub & ddtheta, Doub & dtheta, Doub & ddphi) {
    const Doub v_ccc = v[i][j][k];
    const Doub v_pci = v[i + 1][j][k];
    const Doub v_mci = v[i - 1][j][k];

    const Doub r_lcc = (*r_p)[i];
    const Doub r_lpc = (*r_p)[i + 1];
    const Doub r_lmc = (*r_p)[i - 1];

    dr = (v_pci - v_mci) / (r_lpc - r_lmc);
    ddr = ((v_pci - v_ccc) / (r_lpc - r_lcc) - (v_ccc - v_mci)/(r_lcc - r_lmc)) * 2.0 / (r_lpc - r_lmc);

    const Doub v_pcj = v[i][j + 1][k];
    const Doub v_mcj = v[i][j - 1][k];

    const Doub theta_lcc = (*theta_p)[j];
    const Doub theta_lpc = (*theta_p)[j + 1];
    const Doub theta_lmc = (*theta_p)[j - 1];

    dtheta = (v_pcj - v_mcj) / (theta_lpc - theta_lmc);
    ddtheta = ((v_pcj - v_ccc) / (theta_lpc - theta_lcc) - (v_ccc - v_mcj) / (theta_lcc - theta_lmc)) * 2.0 / (theta_lpc - theta_lmc);
    
    const Doub v_pck = v[i][j][k + 1];
    const Doub v_mck = v[i][j][k - 1];

    ddphi = (v_pck - 2.0*v_ccc + v_mck) / (delta_phi * delta_phi);
    
    return 0;
  }


  //1st Derivative with respect to r
  inline Doub dr(int i, int j, int k) {
    return (v[i+1][j][k] - v[i-1][j][k])/((*r_p)[i+1] - (*r_p)[i-1]);
  }

   //2nd Derivative with respect to r
  inline Doub ddr(int i, int j, int k) {

    //Need to find the difference between each r value because of log(r) scale
    //In the case of constant delta_r, this agrees with cartesian finite difference model
    return ( (v[i+1][j][k] - v[i][j][k])/((*r_p)[i+1] - (*r_p)[i]) 
        - (v[i][j][k] - v[i-1][j][k])/((*r_p)[i] - (*r_p)[i-1]) ) *
      2.0 / ((*r_p)[i+1] - (*r_p)[i-1]);
  }

  //1st Derivative with respect to theta
  inline Doub dtheta(int i, int j, int k) {
    return (v[i][j+1][k] - v[i][j-1][k])/((*theta_p)[j+1] - (*theta_p)[j-1]);
  }

  //2nd Derivative with respect to theta
  inline Doub ddtheta(int i, int j, int k) {

    //Similar to ddr
    return ( (v[i][j+1][k] - v[i][j][k])/((*theta_p)[j+1] - (*theta_p)[j]) 
        - (v[i][j][k] - v[i][j-1][k])/((*theta_p)[j] - (*theta_p)[j-1]) ) *
      2.0 / ((*theta_p)[j+1] - (*theta_p)[j-1]);
  }
  // inline Doub dx(int i, int j, int k) {
  //   return (v[i][j+1][k] - v[i][j-1][k])/(2.0*delta_x);
  // }
  // inline Doub ddx(int i, int j, int k) {
  //   return (v[i][j+1][k] - 2.0*v[i][j][k] + v[i][j-1][k])/(delta_x*delta_x);
  // }

  //1st Derivative with respect to phi
  //Assume phi scale is linearly modelled ==> we can use regular finite difference model
  inline Doub dphi(int i, int j, int k) {
    return (v[i][j][k+1] - v[i][j][k-1])/(2.0*delta_phi);
  }

  //2nd Derivative with respect to phi
  inline Doub ddphi(int i, int j, int k) {
    return (v[i][j][k+1] - 2.0*v[i][j][k] + v[i][j][k-1])/(delta_phi*delta_phi);
  }
  //================================================
  // Destructor
  //================================================

  //Deletes grid
  ~gf3d()
    {
      if (v != NULL) {
  // cout << " deleting gf3d... " << endl; 
        delete[] (v[0][0]);
        delete[] (v[0]);
        delete[] (v);
      }
    }
};

#endif  /* GF_H */
