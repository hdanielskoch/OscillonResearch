#Integrate Lane Emden equation


#
# import some math utilities
#
from pylab import *
#
# import plotting utilities
#
import matplotlib.pyplot as plt
#
# compute changes in phi and theta
#
def dPhi(kappa):
    dphidr = kappa
    
    return dphidr

def dkappa(r,phi,kappa,alpha):
    dkappadr = -2.0*kappa/r + alpha*alpha*phi -3.0/4.0*phi**3+5.0/8.0**5
    
    return dkappadr
#
# find theta and phi as function of xi
# integrate from xi=0, theta=1, phi=0 to theta = 1
#
def oscillon(r_0, phi_0, kappa_0, alpha, rMax, steps):
#
# get started...
#
    
    dr = rMax / steps
    phi = phi_0
    kappa = kappa_0
    r = r_0
    r_list = [ r ]
    phi_list = [ phi ]
    kappa_list = [ kappa ]
    print(dr)

#
# ... integrate ...
#
    while (r < rMax):
        phi = phi + dPhi(kappa)*dr

        kappa = kappa + dkappa(r,phi,kappa,alpha)*dr*dr

        print("r: ", r, "phi: ", phi, "kappa: ", kappa)
        r = r + dr
        r_list.append(r)
        phi_list.append(phi)
        kappa_list.append(kappa)

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


