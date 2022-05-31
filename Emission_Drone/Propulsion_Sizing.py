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
    """
    :param: labda - [-]
    :param: dzeta - [-]
    :param: N - [rpm]
    :param: K0 - [-]
    :param: eta - [-]
    :param: Hp - [m]
    :param: Dp - [m]
    :param: alpha0 - [rad]
    """
    A = 5
    thrust = eq.thrustCoefficient(labda,dzeta,2,K0,eta,Hp,Dp,alpha0,A) * 1.225 * (N/60)**2 * Dp**4
    return thrust

def propellerTorque(Cfd, K0, e,eta,Hp,Dp,alpha0, labda, dzeta, N):
    """
    :param: labda - [-]
    :param: dzeta - [-]
    :param: N - [rpm]
    :param: K0 - [-]
    :param: eta - [-]
    :param: Hp - [m]
    :param: Dp - [m]
    :param: alpha0 - [rad]
    :param: Cfd - [-]
    """
    A = 5
    C_d = eq.dragCoefficient(Cfd, A, K0, e, eta, Hp, Dp, alpha0)
    moment = eq.momentCoefficient(C_d, A, labda, dzeta, 2) * 1.225 * (N/60)**2 * Dp**5
    return moment

def motor_U_I(M, N, KV0, Um0, Im0, Rm):
    U = eq.motorVoltage(M, KV0, Um0, Im0, N, Rm)
    I = eq.motorCurrent( M, KV0, Um0, Im0, Rm)

    print(U, I)
    return U, I


def check_propellers(propeller_matrix, total_mass, labda,dzeta,K0,eta,alpha0, **kwargs):
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


def flight_time(battery_w,  battery_cap, frame_w, no_propellers, prop_eff):
    return ((battery_cap*battery_w)/(frame_w+battery_w))*prop_eff*((frame_w+battery_w)/no_propellers)

test_mat = [["APC 6×4.1SF", 0.1524, 0.10414, 20000],["T-Motor SW 13x5", 0.3302, 0.127, 9600],\
            ["APC 11x12E", 0.2794, 0.3048, 13636.36364]]

testmat=np.array([["APC 11x4.6SF", 0.2794, 0.11684, 15000],["APC 11x12E", 0.2794, 0.3048, 13636.36364],["T-Motor SW 11x4.2",0.2794,0.10668,11000],["DJI Mavic 3", 0.239,0.135,13000]])

def propellerEfficiency(labda,dzeta,K0,eta,alpha0,e,Cfd,propellerMatrix,**kwargs):
    Dp=propellerMatrix[:,1]
    Dp=Dp.astype(np.float)
    Hp=propellerMatrix[:,2]
    Hp=Hp.astype(np.float)
    names=propellerMatrix[:,0]
    for i in range(len(Dp)):
        N=np.arange(10,12010,10)
        thrust = propellerThrust(labda,dzeta,K0,eta,Hp[i],alpha0,N,Dp[i])

        torque = propellerTorque(Cfd, K0, e,eta,Hp[i],Dp[i],alpha0, labda, dzeta, N)

        # efficiency=(thrust)/( torque * N * (2*np.pi/60)  )
        efficiency = [[thrust[i] for i in range(len(N))], [thrust[i] / (torque[i] * N[i] * (2*np.pi/60)) for i in range(len(N))]]

        plt.plot(efficiency[0],efficiency[1], label=names[i])

    plt.xlim(6, 16)
    plt.ylim(0, 0.1)
    plt.xlabel("Thrust in N")
    plt.ylabel("Propeller Efficiency in N/W")
    plt.legend()
    plt.show()

    return
propellerEfficiency(labda,dzeta,K0,eta,alpha0,e,Cfd,testmat)
print(check_propellers(test_mat, 2.8, labda,dzeta,K0,eta,alpha0, Cp,**kwargs))

def motor_efficiency(M, rpm, KV0, Um0, Im0, Rm):
    U, I = motor_U_I(M, rpm, KV0, Um0, Im0, Rm)
    efficiency = M * rpm / (U*I)
    return efficiency

def compare_motor_efficiencies(motor_matrix, Hp, Dp, labda, dzeta, K0, eta, alpha0, Cfd, e):
    """
    motor_matrix should have a row for each motor with: ['name', KV0, Um0, Im0, Rm]
    """
    rpm_range = np.arange(0,6000,400)
    for motor in motor_matrix:
        name = motor[0]
        KV0 = motor[1]
        Um0 = motor[2]
        Im0 = motor[3]
        Rm = motor[4]
        thrusts = []
        efficiencies = []
        for rpm in rpm_range:
            thrust = propellerThrust(labda,dzeta,K0,eta,Hp,alpha0,rpm,Dp)
            thrusts.append(thrust)
            torque = propellerTorque(Cfd, K0, e,eta,Hp,Dp,alpha0, labda, dzeta, rpm)
            eff = motor_efficiency(torque, rpm, KV0, Um0, Im0, Rm)
            efficiencies.append(eff)
        plt.scatter(thrusts, efficiencies, label=name)
    plt.legend()
    plt.show()

test_motor_matrix  = [['multistar 2306 2150kv', 2150, ]]
