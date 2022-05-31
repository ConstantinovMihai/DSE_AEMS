import Performance_Analysis as eq
import numpy as np
import matplotlib.pyplot as plt

"""
Propeller Modeling Parameters
:param: Gp - Propeller Weight [g]
:param: Hp - Propeller Pitch [m]
:param: Bp - Blade Number [m]
:param: Dp - Propeller Diameter [m] 
:param: Cp - Blade Average Chord Length [m]
:param: N - Rotations per minute [m]
"""

"C O N S T A N T S"
eta = 0.85 #Correction factor for downwash
labda = 0.75 #correction coefficient lift
dzeta = 0.55 #correction coefficient weight
e = 0.83 #oswald efficiency factor
Cfd = 0.015 #zero-lift drag coefficient
alpha0 = 0 #Zero-lift angle of attack
K0 = 6.11 #Cl-alpha [rad^-1]

totalweight=2.8*9.81

def propellerThrust(labda,dzeta,K0,eta,Hp,alpha0,N,Dp,**kwargs):
    A = 5
    thrust = eq.thrustCoefficient(labda,dzeta,2,K0,eta,Hp,Dp,alpha0,A) * 1.225 * (N/60)**2 * Dp**4
    return thrust

def propellerTorque(Cfd, K0, e,eta,Hp,Dp,alpha0, labda, dzeta, N):
    A = 5
    C_d = eq.dragCoefficient(Cfd, A, K0, e, eta, Hp, Dp, alpha0)
    moment = eq.momentCoefficient(C_d, A, labda, dzeta, 2) * 1.225 * (N/60)**2 * Dp**5
    return moment

def check_propellers(propeller_matrix, total_mass, labda,dzeta,K0,eta,alpha0):
    """
    Calculates whether propeller options can provide enough thrust for T/W = 2
    propeller_matrix should have a row for each option: ["prop_name", diameter (m), pitch (m), max rpm]
    """
    required_thrust = total_mass * 9.81 / 2
    result_mat = []
    for row in propeller_matrix:
        name = row[0]
        T = propellerThrust(labda,dzeta,2,K0,eta,row[2],alpha0,row[3],row[1])
        if T > required_thrust:
            result_mat.append([name, "yes"])
        else:
            result_mat.append([name, "no"])

    return result_mat # return matrix with [[prop name, yes/no],[prop name, yes/no],[prop name, yes/no]]

# test_mat = [["APC 6×4.1SF", 0.1524, 0.10414, 20000],["T-Motor SW 13x5", 0.3302, 0.127, 9600],\
#             ["APC 11x12E", 0.2794, 0.3048, 13636.36364], ["Mezjlik 14x4.5", 0.3356,0.1143, 12900],["T-Motor SW 11x4.2",0.2794,0.10668,11000],\
#             ["T-Motor SW 15x5.6",0.381,0.14224,],["T-Motor CF 14x4.8",0.3556,0.12192,],["APC 10x4.6SF", 0.254, 0.11684, 6500],["APC 11x12WE", 0.254, 0.3048, 15000]]


testmat=np.array([["APC 11x4.6SF", 0.2794, 0.11684, 15000],["APC 11x12E", 0.2794, 0.3048, 13636.36364],["T-Motor SW 11x4.2",0.2794,0.10668,11000]])

def propellerEfficiency(labda,dzeta,K0,eta,alpha0,e,Cfd,propellerMatrix,**kwargs):
    Dp=propellerMatrix[:,1]
    Dp=Dp.astype(np.float)
    Hp=propellerMatrix[:,2]
    Hp=Hp.astype(np.float)
    names=propellerMatrix[:,0]
    for i in range(len(Dp)):
        N=np.arange(10,12000,10)
        thrust=propellerThrust(labda,dzeta,K0,eta,Hp[i],alpha0,N,Dp[i])
        efficiency=propellerThrust(labda,dzeta,K0,eta,Hp[i],alpha0,N,Dp[i])/( propellerTorque(Cfd, K0, e,eta,Hp[i],Dp[i],alpha0, labda, dzeta, N) * (N * 0.1047198 ) )
        plt.plot(efficiency,thrust, label=names[i])
    plt.xlim(5,16)
    plt.ylim(0,0.03)
    plt.show()

    return
propellerEfficiency(labda,dzeta,K0,eta,alpha0,e,Cfd,testmat)