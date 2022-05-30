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
eta = 0.75 #Correction factor for downwash
labda = 0.75 #correction coefficient lift
dzeta = 0.65 #correction coefficient weight
e = 0.83 #oswald efficiency factor
Cfd = 0.015 #zero-lift drag coefficient
alpha0 = 0 #Zero-lift angle of attack
K0 = 6.11 #Cl-alpha [rad^-1]

totalweight=2.7*9.81

def thrustPropeller(labda,dzeta,Bp,K0,eta,Hp,alpha0,height,Temp,N,Dp,**kwargs):
    A = 6.5
    thrust = eq.thrustCoefficient(labda,dzeta,Bp,K0,eta,Hp,Dp,alpha0,A) * eq.rho(height, Temp) * (N/60)**2 * Dp**4
    return thrust

def momentPropeller(Cfd, K0, e,eta,Hp,Dp,alpha0, labda, dzeta, Bp):
    moments={}
    speed = {}
    for i in range(len(Hp)):
        A = 5
        C_d = eq.dragCoefficient(Cfd, A, K0, e, eta, Hp[i], Dp[i], alpha0)
        N=np.arange(0,eq.speed(totalweight,1.225,Dp[i], eq.thrustCoefficient(labda,dzeta,Bp[i],K0,eta,Hp[i],Dp[i],alpha0,A), Bp[i]),10)
        moment = eq.momentCoefficient(C_d, A, labda, dzeta, Bp[i]) * 1.225 * (N/60)**2 * Dp[i]**5
        print(Hp[i]/0.0254,Dp[i]/0.0254,Bp[i]/0.0254)
        moments[f"Propeller {Dp[i]/0.0254}x{Hp[i]/0.0254}"]=moment
        speed[ f"Propeller {Dp[ i ] / 0.0254}x{Hp[ i ] / 0.0254}" ] = N

    plt.plot(list(speed.values()), list(moments.values()))
    plt.legend()  # To draw legend
    plt.show()
    return

Dp = np.array([14.0, 15.0, 16.0, 16.0, 17.0, 18.0, 20.0])*0.0254
Hp = np.array([4.8, 5, 5.4, 6.1, 5.8, 6, 6.2])*0.0254
Bp = np.array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])

momentPropeller(Cfd, K0, e,eta,Hp,Dp,alpha0, labda, dzeta, Bp)
