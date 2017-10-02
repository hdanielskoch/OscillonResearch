//Spins.h

#ifndef Spins_H
#define	Spins_H

#include "nr3.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

const static int thread_Num = 16;


using namespace std;

class Wave {
public:
	//Constructs grid functions
	Wave(string);

	void init();

	//Initializes 3 vectors for both u and v and loops through predictor functions
	void propagate();

	void predictor();

	void corrector();

	void analytical(Doub);

	void fileWrite(int, Doub);

	void setHubble();

	Doub calcRhoAvg();

	void resetPos();

	void deconstructor();

	~Wave();

private:
	//Initialize grid functions
	VecDoub *v_O;
	VecDoub *v_N;
	VecDoub *v_T;

	VecDoub *u_O;
	VecDoub *u_N;
	VecDoub *u_T;
	VecDoub *u_A;
	VecDoub *u_I;
	VecDoub *uTemp;
	VecDoub *vTemp;

	VecDoub *x_Pos;

	//Parameters that are read in
	Doub xMin;
	Doub xMax;
	int nx;
	Doub aMax;
	Doub courantCond;
	Doub mPL;
	Doub alphaF;
	Doub m;
	Doub m_Mpl;
	Doub lambda;
	int gridFactMax;
	int fileInc;
	int HImplement;
	Doub psiInit;

	//Calculated parameters
	Doub xInc;
	Doub dt;
	Doub dn;
	Doub xi;
	Doub s2;
	Doub H_i;
	Doub m2; //m^2
	Doub a2; //a^-2
	Doub g;
	Doub g2; //a^2
	Doub a;
	Doub H;
	Doub alpha;
	int gridFactor;
};

#endif	/* Spins_H */