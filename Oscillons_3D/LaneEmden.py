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
def derivPhi(theta, xi, n):
    dphidxi = - theta**n * xi*xi
    
    return dphidxi

def derivTheta(phi, xi):
    dthetadxi = phi / (xi*xi)
    
    return dthetadxi
#
# find theta and phi as function of xi
# integrate from xi=0, theta=1, phi=0 to theta = 1
#
def laneEmden(theta_0, phi_0, xi_0, theta_final, n, steps):
#
# get started...
#
    delta_xi = theta_0 / steps
    xi = xi_0
    theta = theta_0
    phi = phi_0
    xi_list = [ xi ]
    phi_list = [ phi ]
    theta_list = [ theta ]
    theta_Ana_list = [theta]

#
# ... integrate ...
#
    while (theta > theta_final):
        dphidxi = derivPhi(theta,xi,n)
        phi = phi + dphidxi * delta_xi #might be  minus

        dthetadxi = derivTheta(phi, xi)
        theta = theta + dthetadxi * delta_xi

        print("xi: ", xi, "phi: ", phi, "theta: ", theta)
        r = r + dr
        absPhi = abs(phi)
        xi_list.append(xi)
        phi_list.append(absPhi)
        theta_list.append(theta**n)
        theta_Ana_list.append(math.sin(xi) / xi)

#
# ... and plot results, including labels and legend
# 
    print(xi)
    print(absPhi)
    print(theta)
    print(xi**3 /(3*absPhi))
    plt.plot(xi_list,theta_list, label="numerical")
    #plt.plot(xi_list,theta_Ana_list, label="analytical")
    plt.xlabel('xi')
    plt.ylabel('p/p_c')
    plt.legend()
    savefig("p|p_c vs Xi.pdf")
    #show()


