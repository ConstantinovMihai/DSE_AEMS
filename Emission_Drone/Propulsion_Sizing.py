import Performance_Analysis.py as eq
import numpy as np
import matplotlib.pyplot as plt

"""
Propeller Modeling Parameters
:param: Gp - Propeller Weight
:param: Hp - Propeller Pitch
:param: Bp - Blade Number
:param: Dp - Propeller Diameter 
:param: Cp - Blade Average Chord Length
:param: N - Rotations per minute
"""

"C O N S T A N T S"
 #Aspect ratio Dp/Cp
eta = 0.9 #Correction factor for downwash
labda = 0.8 #correction coefficient lift
dzeta = 0.65 #correction coefficient weight
e = 0.8 #oswald efficiency factor
Cfd = 0.015 #zero-lift drag coefficient
alpha0 = 0 #Zero-lift angle of attack
K0 = 6.11



def thrustPropeller(labda,dzeta,Bp,K0,eta,Hp,alpha0,height,Temp,N,Dp, Cp,**kwargs):
    A = Dp / Cp
    thrust = eq.thrustCoefficient(labda,dzeta,Bp,K0,eta,Hp,Dp,alpha0,A) * eq.rho(height, Temp) * (N/60)**2 * Dp**4
    return thrust

def momentPropeller(height, Temp, Cfd,N, Cp, K0, e,eta,Hp,Dp,alpha0, labda, dzeta, Bp):
    A = Dp / Cp
    C_d = eq.dragCoefficient(Cfd, A, K0, e,eta,Hp,Dp,alpha0)
    moment = eq.momentCoefficient(C_d, A, labda, dzeta, Bp) * eq.rho(height, Temp) *  (N/60)**2 * Dp**5
    return moment

def motor_U_I(M, N, KV0, Um0, Im0, Rm):
    U = eq.motorVoltage(M, KV0, Um0, Im0, N, Rm)
    I = eq.motorCurrent( M, KV0, Um0, Im0, Rm)

    print(U, I)
    return U, I

def check_propellers(propeller_matrix, total_mass, labda,dzeta,K0,eta,alpha0,height,Temp, Cp,**kwargs):
    """
    Calculates whether propeller options can provide enough thrust for T/W = 2
    propeller_matrix should have a row for each option: [prop_name, diameter (m), pitch (m), max rpm, number of blades]
    """
    required_thrust = total_mass * 9.81 / 2
    result_mat = []
    for row in propeller_matrix:
        name = row[0]
        T = thrustPropeller(labda,dzeta,row[4],K0,eta,row[2],alpha0,height,Temp,row[3],row[1],Cp)
        if T > required_thrust:
            result_mat.append([name, "yes"])
        else:
            result_mat.append([name, "no"])

    return result_mat # return matrix with [[prop name, yes/no],[prop name, yes/no],[prop name, yes/no]]