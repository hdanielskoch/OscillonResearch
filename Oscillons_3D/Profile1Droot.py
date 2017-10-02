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

r = 0.0
phi = 0.0
alphaF = 0.0
phi_list = [r]
r_list = [phi]

def dPhi(kappa):
    dphidr = kappa
    
    return dphidr

def dkappa(r,phi,kappa,alpha):
    dkappadr = alpha*alpha*phi -3.0/4.0*phi**3 + 5.0/8.0*phi**5
    
    return dkappadr

#Function of phiF which is the factor that multiplies phi_0
#returns the value that optimized to be 0
def oscillon(phiF):

    global alphaF

    phiCrit = math.sqrt(9.0/10.0)
    alphaCrit = math.sqrt(27.0/160.0)
    alpha = alphaF*alphaCrit
    phi_0 = phiF*phiCrit
    rMax = 100.0
    steps = 1000.0
    kappa_0 = 0.0
    r_0 = 0.0001


    #alpha = math.sqrt(3.0/8.0*phi_0*phi_0 - 5.0/24.0*phi_0**4)
    
    dr = rMax / steps
    phi = phi_0
    kappa = kappa_0
    r = r_0
    global r_list
    global phi_list
    r_list = [ r ]
    phi_list = [ phi ]
    kappa_list = [ kappa ]


    #integrate until r=rMax
    while (r < rMax):
        phi = phi + dPhi(kappa)*dr

        kappa = kappa + dkappa(r,phi,kappa,alpha)*dr*dr

        r = r + dr
        r_list.append(r)
        phi_list.append(phi)
        kappa_list.append(kappa)

    #run rootfinding method
    #print ("phiF ", phiF, "funct : ", kappa + phi*alpha + phi/rMax)


    return kappa + phi*alpha + phi/r

def main(alphaF_in):

    global alphaF
    alphaF = alphaF_in
  
    phiCrit = sqrt(9.0/10.0);
    tol = 0.000000001

    #//Initialize bounds for bisect with phiF
    low = 0.5;
    high = 1.0;
    #for i in frange(0.99,0.99,0.01):
        #alphaF = i
    #rootfind to find what r would be at desired kappa value
    phi_0 = opt.bisect(oscillon,low,high,xtol=tol) 
    print phi_0
#
# ... and plot results, including labels and legend
# 
    #global r_list
    #global phi_list
    #print(phi)
    #print(kappa)
    plt.plot(r_list,phi_list, label="Phi vs r")
    #plt.plot(xi_list,theta_Ana_list, label="analytical")
    plt.xlabel('r')
    plt.ylabel('Phi')
    plt.legend()
    #savefig("phiShooting.pdf")
    show()


