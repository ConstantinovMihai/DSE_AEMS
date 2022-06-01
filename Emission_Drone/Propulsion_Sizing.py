import Performance_Analysis as eq
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
"""
Propeller Modeling Parameters
:param: Gp - Propeller Weight [g]
:param: Hp - Propeller Pitch [m]
:param: Bp - Blade Number [m]
:param: Dp - Propeller Diameter [m] 
:param: Cp - Blade Average Chord Length [m]
:param: N - Rotations per minute [m]
"""


def matrix_gen(SheetNumber):
    df = pd.read_excel(r'.\Propellers.xlsx', sheet_name = SheetNumber)
    a = df.to_numpy()
    return a
propellerMatrix=matrix_gen(0)


"C O N S T A N T S"
eta = 0.85 #Correction factor for downwash
labda = 0.75 #correction coefficient lift
dzeta = 0.5 #correction coefficient weight
e = 0.83 #oswald efficiency factor
Cfd = 0.015 #zero-lift drag coefficient
alpha0 = 0 #Zero-lift angle of attack
K0 = 6.11 #Cl-alpha [rad^-1]

totalweight=3.028*9.81


def propellerThrust(propellerMatrix, labda, dzeta, K0, eta, alpha0, total_weight, **kwargs):
    """
    :param: labda - [-]
    :param: dzeta - [-]
    :param: N - [rpm]
    :param: K0 - [-]
    :param: eta - [-]
    :param: Hp - [m]
    :param: Dp - [m]
    :param: alpha0 - [rad]
    :param: total_weight - [kg]
    """

    A = 5  # set the aspect ratio
    for row in propellerMatrix:  # compute the performance of each prop
        names = row[0]
        Dp = row[1]
        Hp = row[2]
        maxrpm = row[3]
        N = np.arange(10, maxrpm, 10)
        thrust = eq.thrustCoefficient(labda, dzeta, 2, K0, eta, Hp, Dp, alpha0, A) * 1.225 * (N/60)**2 * Dp**4

        # prune the thrust values based on maximum (climb) and then plot
        if thrust[-1] > (total_weight*9.81) / 2:
            plt.plot(N, thrust, label=names)

    plt.xlabel("RPM")
    plt.ylabel("Thrust in N")
    plt.ylim(0, 17)
    plt.axhline((total_weight*9.81)/4, linestyle="dotted", label="Thover")
    plt.axhline((total_weight*9.81)/2, linestyle="dashed", label="Tmax")
    plt.legend()
    plt.show()

    return


#propellerThrust(propellerMatrix[1:28, 1:5], labda, dzeta, K0, eta, alpha0, 3.028)


def propellerTorque(propellerMatrix,Cfd, K0, e,eta,Hp,Dp,alpha0, labda, dzeta, N):
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
    for row in propellerMatrix:
        names=row[0]
        Dp=row[1]
        Hp=row[2]
        maxrpm=row[3]
        N=np.arange(10,maxrpm,10)
        C_d = eq.dragCoefficient(Cfd, A, K0, e, eta, Hp, Dp, alpha0)
        moment = eq.momentCoefficient(C_d, A, labda, dzeta, 2) * 1.225 * (N/60)**2 * Dp**5
        plt.plot(N, moment, label=names)
    plt.xlabel("RPM")
    plt.ylabel("Moment in Nm")
    plt.legend()
    plt.show()
    return

# def motor_U_I(M, KV0, Um0, Im0, N, Rm):
#     U = eq.motorVoltage(M, KV0, Um0, Im0, N, Rm)
#     I = eq.motorCurrent( M, KV0, Um0, Im0, Rm)
#     print(U, I)
#     return U, I

def motor_U_I(M, N, KT, Im0, Rm):
    """
    Calculates motor voltage U [V] and current I [A].
    Parameters: motor torque M = propeller torque [Nm], motor speed N = propeller speed [rpm]
    Torque constant KT [Nm/A] = 9.55((Um0 − Im0*Rm) / (KV0*Um0)), Nominal no-load current Im0 [A],
    resistance Rm [Ohm]
    """
    U = (M / KT + Im0) * Rm + KT * N / 9.55
    I = M / KT + Im0
    return U, I

# def check_propellers(propeller_matrix, total_mass, labda,dzeta,K0,eta,alpha0, **kwargs):
#     """
#     Calculates whether propeller options can provide enough thrust for T/W = 2
#     propeller_matrix should have a row for each option: ["prop_name", diameter (m), pitch (m), max rpm]
#     """
#     required_thrust = total_mass * 9.81 / 2
#     result_mat = []
#     for row in propeller_matrix:
#         name = row[0]
#         T = propellerThrust(labda,dzeta,K0,eta,row[2],alpha0,0.9*row[3],row[1])
#         if T > required_thrust:
#             result_mat.append([name, "yes"])
#         else:
#             result_mat.append([name, "no"])
#
#     return print(result_mat) # return matrix with [[prop name, yes/no],[prop name, yes/no],[prop name, yes/no]]


# check_propellers(usefullmatrix[1:28,1:5],3.028,labda,dzeta,K0,eta,alpha0)



# test_mat = [["APC 6×4.1SF", 0.1524, 0.10414, 20000],["T-Motor SW 13x5", 0.3302, 0.127, 9600],\
#             ["APC 11x12E", 0.2794, 0.3048, 13636.36364], ["Mezjlik 14x4.5", 0.3356,0.1143, 12900],["T-Motor SW 11x4.2",0.2794,0.10668,11000],\
#             ["T-Motor SW 15x5.6",0.381,0.14224,],["T-Motor CF 14x4.8",0.3556,0.12192,],["APC 10x4.6SF", 0.254, 0.11684, 6500],["APC 11x12WE", 0.254, 0.3048, 15000]]


# def flight_time(battery_w,  battery_cap, frame_w, no_propellers, prop_eff):
#     return ((battery_cap*battery_w)/(frame_w+battery_w))*prop_eff*((frame_w+battery_w)/no_propellers)

test_mat = [["APC 6×4.1SF", 0.1524, 0.10414, 20000],["T-Motor SW 13x5", 0.3302, 0.127, 9600],\
            ["APC 11x12E", 0.2794, 0.3048, 13636.36364]]

testmat=np.array([["APC 6×4.1SF", 0.1778, 0.10414, 20000],["APC 11x4.6SF", 0.2794, 0.11684, 15000],["APC 11x12E", 0.2794, 0.3048, 13636.36364],["T-Motor SW 11x4.2",0.2794,0.10668,11000],["DJI Mavic 3", 0.239,0.135,13000]])

def propellerEfficiency(labda,dzeta,K0,eta,alpha0,e,Cfd,propellerMatrix,**kwargs):
    Dp=propellerMatrix[:,1]
    Dp=Dp.astype(np.float)
    Hp=propellerMatrix[:,2]
    Hp=Hp.astype(np.float)
    names=propellerMatrix[:,0]
    for i in range(len(Dp)):
        N=np.arange(10,20010,10)
        thrust = propellerThrust(labda,dzeta,K0,eta,Hp[i],alpha0,N,Dp[i])

        torque = propellerTorque(Cfd, K0, e,eta,Hp[i],Dp[i],alpha0, labda, dzeta, N)

        # efficiency=(thrust)/( torque * N * (2*np.pi/60)  )
        efficiency = [[thrust[i] for i in range(len(N))], [thrust[i] / (torque[i] * N[i] * (2*np.pi/60)) for i in range(len(N))]]

        #plt.plot(efficiency[0],efficiency[1], label=names[i])
        plt.plot(N,thrust,label=names[i]+"Lel")
    plt.xlabel("Thrust in N")
    plt.ylabel("Propeller Efficiency in N/W")
    plt.legend()
    plt.show()

    return
#N = np.arange(10,14000,10)
#plt.plot(N,propellerThrust(labda,dzeta,K0,eta,0.3048,alpha0,N,0.2794))
#plt.show()
#plt.plot(N,propellerTorque(Cfd, K0, e,eta,0.3048,0.2794,alpha0, labda, dzeta, N))
#plt.show()


#propellerEfficiency(labda,dzeta,K0,eta,alpha0,e,Cfd,testmat)
def calc_propellerThrust(labda,dzeta,K0,eta,Hp,alpha0,N,Dp,**kwargs):
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

def calc_propellerTorque(Cfd, K0, e,eta,Hp,Dp,alpha0, labda, dzeta, N):
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

def motor_efficiency(M, rpm, KT, Im0, Rm):
    U, I = motor_U_I(M, rpm, KT, Im0, Rm)
    N = rpm * 0.10472
    print('U=',U,'I=', I)
    efficiency = M * N / (U*I)
    return efficiency

def compare_motor_efficiencies(motor_matrix, Hp, Dp, labda, dzeta, K0, eta, alpha0, Cfd, e):
    """
    motor_matrix should have a row for each motor with: ['name', KT, Im0, Rm]
    """
    rpm_range = np.arange(0,10000,400) # Input max rpm here.
    for motor in motor_matrix:
        name = motor[0]
        KT = motor[1]
        Im0 = motor[2]
        Rm = motor[3]
        thrusts = []
        efficiencies = []
        for rpm in rpm_range:
            thrust = calc_propellerThrust(labda,dzeta,K0,eta,Hp,alpha0,rpm,Dp)
            thrusts.append(thrust)
            torque = calc_propellerTorque(Cfd, K0, e,eta,Hp,Dp,alpha0, labda, dzeta, rpm)
            eff = motor_efficiency(torque, rpm, KT, Im0, Rm)
            efficiencies.append(eff)
            print('Thrust=', thrust,'torque=', torque)
        plt.scatter(thrusts, efficiencies, label=name)
    plt.xlabel('Thrust [N]')
    plt.ylabel('Motor efficicency')
    plt.legend()
    plt.show()

# test_motor_mat = [['EC frameless DT 50 S', 0.0666, 0.162, 0.583],['ECX FLAT 32 L', 0.0215, 0.180, 0.445]]
# compare_motor_efficiencies(test_motor_mat, 0.11684, 0.2794, labda, dzeta, K0, eta, alpha0, Cfd, e)


def Battery_endurance(Cb, Cmin, Ub, Im, Um, Re, Bmass):
    """
    Calculates endurance for battery in minutes
    Battery parameters: Battery capacity Cb [Ah], battery minimum capacity Cmin (set at 0.2Cb),
    battery voltage Ub [V]
    Set parameters (set for hovering thrust): motor current Im [A], motor voltage Um [V], ESC resistance Re [ohm]
    """
    sigma = Um + Im * Re / Ub
    Ie = sigma * Im
    Iother = 1
    Ib = 4 * Ie + Iother
    hovering_endurance = (Cb - Cmin) * (60/1000) / Ib
    print('hovering_endurance: ',hovering_endurance,' minutes - mass:', Bmass, 'g')
    return hovering_endurance
