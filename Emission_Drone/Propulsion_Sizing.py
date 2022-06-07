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

propellerMatrix = matrix_gen(0)
sixinchprops = matrix_gen(2)
seveninchprops = matrix_gen(3)
eightinchprops = matrix_gen(4)
nineinchprops = matrix_gen(5)
teninchprops = matrix_gen(6)
eleveninchprops = matrix_gen(7)

"C O N S T A N T S"
eta = 0.85 #Correction factor for downwash
labda = 0.75 #correction coefficient lift
dzeta = 0.5 #correction coefficient weight
e = 0.83 #oswald efficiency factor
Cfd = 0.015 #zero-lift drag coefficient
alpha0 = 0 #Zero-lift angle of attack
K0 = 6.11 #Cl-alpha [rad^-1]

totalweight=3.028*9.81


def propellerThrust(labda, dzeta, K0, eta, alpha0,Dp,Hp,N,**kwargs):
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
    A = 5  # set the aspect ratio
    thrust = eq.thrustCoefficient(labda, dzeta, 2, K0, eta, Hp, Dp, alpha0, A) * 1.225 * (N/60)**2 * Dp**4
    return thrust


def propellerTorque(Cfd, K0, e, eta, alpha0, labda, dzeta, Hp, Dp, N):
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
    torque = eq.momentCoefficient(C_d, A, labda, dzeta, 2) * 1.225 * (N/60)**2 * Dp**5
    return torque


def propellerComparison(propellerMatrix, labda, dzeta, K0, eta, alpha0,total_weight):

    viableoptions= {}
    result1mat=[]
    result2mat = [ ]
    torquevalues = []
    fig, axs = plt.subplots(2)
    for row in propellerMatrix:
        names = row[ 0 ]
        Dp = row[ 1 ]
        Hp = row[ 2 ]
        maxrpm = row[ 3 ]

        N = np.arange(100, maxrpm, 100)

        thrust = propellerThrust(labda, dzeta, K0, eta, alpha0,Dp,Hp,N)
        if thrust[ -1 ] > (total_weight * 9.81) / 2:
            axs[0].plot(N, thrust, label=names)

            found = [False, False]
            for idx, val in enumerate(thrust):
                if val >=  (total_weight * 9.81) / 4 and val <= 1.05 * (total_weight * 9.81) / 4 and not found[0]:
                    found[0] = True
                    rpmNominalThrust = N[idx]
                    result1mat.append([ names, Dp, Hp, maxrpm, rpmNominalThrust ])
                    axs[ 0 ].axvline(rpmNominalThrust, linestyle="dotted", c=np.random.rand(3, ),
                                     label=names + "Nom Thrust RPM")
                if val >=  (total_weight * 9.81) / 2 and val <= 1.05 * (total_weight * 9.81) / 2 and not found[1]:
                    found[1] = True
                    rpmMaxThrust = N[ idx ]
                    result2mat.append([ names, Dp, Hp, maxrpm, rpmMaxThrust ])
                    axs[ 0 ].axvline(rpmMaxThrust, linestyle="dotted",c=np.random.rand(3,),  label=names+"Max Thrust RPM")
        axs[ 0 ].set_xlabel("RPM")
        axs[ 0 ].set_ylabel("Thrust in N")
        axs[ 0 ].set_ylim(0, 17)
    for i in range(len(result1mat)):
        result1mat[i]=result1mat[i]+[result2mat[i][4]]

    print(result1mat)
    for entry in result1mat:
        names1 = entry[ 0 ]
        Dp1 = entry[ 1 ]
        Hp1 = entry[ 2 ]
        maxRpm = entry[ 3 ]
        nomRpm = entry[4]
        maxThrustRpm = entry[5]
        N1 = np.arange(100, maxRpm, 100)
        torque = propellerTorque(Cfd, K0, e, eta, alpha0, labda, dzeta, Hp1, Dp1, N1)
        axs[ 1 ].plot(N1, torque, label=names1)

        axs[ 1 ].set_xlabel("RPM")
        axs[ 1 ].set_ylabel("Torque in Nm")
        axs[ 1 ].axvline(entry[4], linestyle="dotted",c=np.random.rand(3,),  label=names1+"Nominal Rpm")
        axs[ 1 ].legend(fontsize=7, loc=2)

        nomtorque = propellerTorque(Cfd, K0, e, eta, alpha0, labda, dzeta, Hp1, Dp1, nomRpm)
        maxtorque = propellerTorque(Cfd, K0, e, eta, alpha0, labda, dzeta, Hp1, Dp1, maxThrustRpm)
        torquevalues.append([names1, Dp1, Hp1, maxRpm, nomRpm,maxThrustRpm, nomtorque,maxtorque ])

    axs[ 0 ].axhline((total_weight * 9.81) / 4, linestyle="dotted", label="Thover")
    axs[ 0 ].axhline((total_weight * 9.81) / 2, linestyle="dashed", label="Tmax")
    axs[ 0 ].legend(fontsize=7, loc=2)
    plt.show()



    return print(torquevalues)


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
    efficiency = M * N / (U*I)
    return efficiency


def compare_motor_efficiencies(motor_matrix, Hp, Dp, labda, dzeta, K0, eta, alpha0, Cfd, e):
    """
    Makes 2 graphs to decide on a motor: Thrust vs efficiency and mass vs efficiency at hovering
    motor_matrix should have a row for each motor with: ['name', KT, Im0, Rm, mass]
    torque constant KT [Nm/A], Nominal no-load current Im0 [A], motor resistance Rm [ohm], motor mass [g]
    """
    rpm_range = np.arange(2200,10000,100) # Input max rpm here.

    # Lists for efficiencies at 7.4N vs mass graph
    masses = []
    efficiencies74N = []
    names = []

    # thrust vs efficiency graphs
    for motor in motor_matrix:
        name = motor[0]
        names.append(name)
        KT = motor[1]
        Im0 = motor[2]
        Rm = motor[3]
        masses.append(motor[4])

        thrusts = []
        efficiencies = []
        for rpm in rpm_range:
            thrust = calc_propellerThrust(labda,dzeta,K0,eta,Hp,alpha0,rpm,Dp)
            thrusts.append(thrust)
            torque = calc_propellerTorque(Cfd, K0, e,eta,Hp,Dp,alpha0, labda, dzeta, rpm)
            eff = motor_efficiency(torque, rpm, KT, Im0, Rm)
            efficiencies.append(eff)
            if rpm == 5800:
                print(name, thrust)
                efficiencies74N.append(eff)
        plt.scatter(thrusts, efficiencies, label=name)
    plt.axvline(x=totalweight/4, color='b', linestyle='dotted', label='Hovering thrust')
    plt.axvline(x=totalweight/2, color='r', linestyle='dotted', label='T/W = 2')
    plt.xlabel('Thrust [N]')
    plt.ylabel('Motor efficiency')
    plt.legend()
    plt.show()

    # graph efficiencies at 7.4N vs mass
    print(names, masses, efficiencies74N)
    plt.scatter(masses, efficiencies74N)
    plt.xlabel('masses [g]')
    plt.ylabel('Motor efficiency at 7.4N')
    for i, label in enumerate(names):
        plt.annotate(label, (masses[i], efficiencies74N[i]))
    plt.show()


def Battery_endurance(Tb, Im, Um):
    """
    Calculates endurance for battery in minutes
    Battery parameters: Battery capacity Cb [mAh], battery minimum capacity Cmin (set at 0.2Cb),
    battery voltage Ub [V]
    Set parameters (set for hovering thrust): motor current Im [A], motor voltage Um [V], ESC resistance Re [ohm]
    hovering time Tb [min]
    """
    Ub_list = np.arange(10,40, 0.5)
    Cb_list = []
    for Ub in Ub_list:
        sigma = Um / Ub
        Ie = sigma * Im
        Iother = 1
        Ib = 4 * Ie + Iother
        Cb_effective =  Tb * Ib * (1000/60)
        Cb_actual = Cb_effective * (10/8)
        Cb_list.append(Cb_actual)
    plt.scatter(Ub_list, Cb_list)
    plt.xlabel('Battery Voltage [V]')
    plt.ylabel('Battery capacity [mAh]')
    plt.show()
    return


def matrix_gen_options(SheetNumber):
    df = pd.read_excel(r'.\Propulsion_options.xlsx', sheet_name = SheetNumber)
    a = df.to_numpy()
    return a


"""
Using/testing the fuctions (Uncomment what you want to use):
"""
# N = np.arange(10,14000,10)
# plt.plot(N,propellerThrust(labda,dzeta,K0,eta,0.3048,alpha0,N,0.2794))
# plt.show()
# plt.plot(N,propellerTorque(Cfd, K0, e,eta,0.3048,0.2794,alpha0, labda, dzeta, N))
# plt.show()
# propellerEfficiency(labda,dzeta,K0,eta,alpha0,e,Cfd,testmat)
# propellerComparison(propellerMatrix[1:28, 1:5], labda, dzeta, K0, eta, alpha0,3.028)
# propellerComparison(eleveninchprops, labda, dzeta, K0, eta, alpha0,3.028)
motor_matrix = matrix_gen_options(2)[2:17,1:6]
compare_motor_efficiencies(motor_matrix, 0.1397, 0.2794, labda, dzeta, K0, eta, alpha0, Cfd, e)
Im, Um = motor_U_I(0.15044989095978484, 5800, 0.0295, 0.277, 0.447)
Battery_endurance(25, Im, Um)

"""
Current choices:
"""
# Choice iteration 1 propeller: 11 inch x 5.5 inch (22.9631112g * 4 = 91.8524448g)
# Choice iteration 1 motor: Maxon 608131 (113g * 4 = 452g)
# Choice iteration 1 ESC: Maxon 438725 (12g * 4 = 48g)
# Choice iteration 1 Battery: 1350g
# Total mass it 1: 3389.85g
# torque at hovering: 0.15044989095978484
# 5800 rpm