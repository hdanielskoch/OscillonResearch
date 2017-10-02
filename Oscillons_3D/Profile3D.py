#Integrate Lane Emden equation

import math
#
# import some math utilities
#
from pylab import *
#
# import plotting utilities
#
import matplotlib.pyplot as plt
import scipy.optimize as opt

#define global variables
phi = 0.0
kappa = 0.0
alpha = 0.0
#
# compute changes in phi and theta
#
def dPhi(kappa):
    dphidr = kappa
    
    return dphidr

def dkappa(r,phi,kappa,alpha):
    dkappadr = -2.0*kappa/r + alpha*alpha*phi -3.0/4.0*phi**3 + 5.0/8.0*phi**5
    
    return dkappadr

def rootFunc(r):
    return kappa + phi*alpha + phi/r
#
# find theta and phi as function of xi
# integrate from xi=0, theta=1, phi=0 to theta = 1
#
def oscillon(r_0, phiF, kappa_0, alphaF, rMax, steps):

    global alpha
    global phi
    global kappa

    phiCrit = math.sqrt(9.0/10.0)
    alphaCrit = math.sqrt(27.0/160.0)
    alpha = alphaF*alphaCrit
    phi_0 = phiF*phiCrit

    #alpha = math.sqrt(3.0/8.0*phi_0*phi_0 - 5.0/24.0*phi_0**4)
    
    dr = rMax / steps
    phi = phi_0
    kappa = kappa_0
    r = r_0
    r_list = [ r ]
    phi_list = [ phi ]
    kappa_list = [ kappa ]


#  integrate until r=rMax
    while phiF < 0.999:

        phi = phiF*phiCrit
        kappa = kappa_0
        r = r_0

        while (r < rMax):
            phi = phi + dPhi(kappa)*dr

            kappa = kappa + dkappa(r,phi,kappa,alpha)*dr*dr

            r = r + dr
            r_list.append(r)
            phi_list.append(phi)
            kappa_list.append(kappa)

    #run rootfinding method
        print ("phiF ", phiF, "funct : ", kappa + phi*alpha + phi/rMax)
        phiF += 0.01
    #rootfind to find what r would be at desired kappa value
    #rNew = opt.bisect(rootFunc,rLeft,rRight) 
    #print("r calculated: ",rNew)

    tol = 0.01


#
# ... and plot results, including labels and legend
# 
    print(r)
    print(phi)
    print(kappa)
    plt.plot(r_list,phi_list, label="Phi vs r")
    #plt.plot(xi_list,theta_Ana_list, label="analytical")
    plt.xlabel('r')
    plt.ylabel('Phi')
    plt.legend()
    savefig("phiShooting.pdf")
    show()


