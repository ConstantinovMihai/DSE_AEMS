import numpy as np
import matplotlib.pyplot as plt
import math as m

#Altitude Parameters

h=10 #meters
T=25 #Celsius

#Basic Parameters

G= 14.7 #Newton
nt= 4 #

#Component Parameters

""" P R O P E L L E R """
Dp = 29 *0.0254  #Propeller Diameter [Inch to m]
Hp = 9.5 * 0.0254 #Propeller Pitch [Inch to m]
Bp = 2 #Blade Number [ nr of blades per propeller]

""" M O T O R"""
KV0 = 890 #
ImMax = 19 #Maximum Continuous Current [A]
Im0 = 0.5 #Nominal No-Load Current [A]
Um0 = 10 #Nominal No-Load Voltage [V]
Rm = 0.101 #Motor Resistance [Ohm]

""" E S C """
IeMax = 30 #Maximum Current  [A]
Re = 0.008 #Resistance [Ohm]

""" B A T T E R Y """
Cb = 16000 #Battery Capacity [mAh]
Rb = 0.01 #Battery Resistance [Ohm]
Ub = 12 #Battery Voltage [V]
Kb = 45 #Maximum Discharge Rate [C]
Iother = 0.9 #Current of other things like flight controller [A]

#Other Parameters

A = 5
eta = 0.85
labda = 0.75
dzeta = 0.5
e = 0.83
Cfd = 0.015
alpha0 = 0
K0 = 6.11
C1 = 3
C2 =1.5
Cmin = 0.2 * Cb
nr = 4

def rho(h,T):
    rho= ( ( 273*101325*( 1-0.0065* (h/(273+T)) )**5.2561 )/(101325*( 273+T)) ) * 1.293
    return rho

def thrustCoefficient(labda,dzeta,Bp,K0,eta,Hp,Dp,alpha0,A):

    C_T=0.25 * m.pi**3 *labda * dzeta**2 * Bp * K0 * (eta * m.atan((Hp/ (m.pi*Dp))) -alpha0 )/ (m.pi * A + K0)

    return C_T
def dragCoefficient(Cfd, A, K0, e,eta,Hp,Dp,alpha0):

    C_d= Cfd+(m.pi * A * K0**2)/(e) * (eta*m.atan((Hp/(np.pi*Dp))) -alpha0)**2/(m.pi * A + K0)**2
    return C_d

def momentCoefficient(C_d, A, labda, dzeta, Bp):

    C_M= 1/(8 * A) * m.pi**2 *C_d * labda * dzeta**2 * Bp**2
    return C_M

def speed(G,rho,Dp,Ct, nr):

    N = 60 * m.sqrt(G/(rho * Dp**4 * Ct * nr))
    return N

def torque(rho, Dp, CM, N):
    M= rho * Dp**5 * CM* (N/60)**2
    return M

def motorVoltage(M, KV0, Um0, Im0, N, Rm):

    Um = ( (M * KV0 * Um0)/(9.55 * (Um0-Im0*Rm)) + Im0) * Rm + (Um0 - Im0*Rm)/(KV0*Um0) * N
    return Um

def motorCurrent( M, KV0, Um0, Im0, Rm):
    Im = (M * KV0 * Um0) / (9.55 * (Um0 - Im0 * Rm)) + Im0
    return Im


def backElectromotive(Um0, Im0, Rm, KV0):

    KE= (Um0-Im0*Rm)/(KV0*Um0)
    return KE

def nominalNoLoadSpeed(KV0,Um0):
    Nm0 = KV0 * Um0
    return Nm0

def noLoadSpeed(KV0,Um):
    Nhat0 = KV0 * Um
    return Nhat0

def nominalPowerConsumption(Um0, Im0, Rm):
    PFe0 = Um0 * Im0 - Im0**2 * Rm
    return PFe0

def actualNoLoadCurrent(PFe0,KE, Nhat0, Nm0):
    Ihat0 = (PFe0) / (KE * Nhat0) * (Nhat0 / Nm0)**1.3
    return Ihat0

def torqueConstant(KE):
    KT = 9.55 * KE
    return KT

def equivalentDCVoltage(Um,Im,Re):
    Ueo = Um + Im * Re
    return Ueo

def throttleCommand(Ueo,Ub):
    sigma = Ueo / Ub
    return sigma


def inputCurrent(sigma,Im):
    Ie=sigma*Im
    return Ie


def ESCVoltage(Ub, nr, Ie, Iother, Rb):
    Ue=Ub-(nr*Ie+Iother)*Rb
    return Ue


def batteryCurrent(nr,Ie, Iother):
    Ib = nr * Ie + Iother
    return Ib


def hoveringTime(Cb,Cmin,Ib):
    Thover = 0.06 * (Cb - Cmin) / Ib
    return Thover




#Calculating Aerodynamic Coefficients
rho=rho(h,T)
Cd=dragCoefficient(Cfd, A, K0, e,eta,Hp,Dp,alpha0)
CT=thrustCoefficient(labda,dzeta,Bp,K0,eta,Hp,Dp,alpha0,A)
CM=momentCoefficient(Cd, A, labda, dzeta, Bp)

#Calculating Motor Parameters
N= speed(G,rho,Dp,CT, nr)
M = torque(rho,Dp, CM,N)
Um = motorVoltage(M, KV0, Um0, Im0, N, Rm)
Im = motorCurrent( M, KV0, Um0, Im0, Rm)
KE = backElectromotive(Um0, Im0, Rm, KV0)
Nm0 = nominalNoLoadSpeed(KV0,Um0)
Nhat0 = noLoadSpeed(KV0,Um)
PFe0 = nominalPowerConsumption(Um0, Im0, Rm)
Ihat0 = actualNoLoadCurrent(PFe0,KE, Nhat0, Nm0)
KT = torqueConstant(KE)

#ESC Parameters
Ueo = equivalentDCVoltage(Um,Im,Re)
sigma = throttleCommand(Ueo,Ub)
Ie = inputCurrent(sigma,Im)
Ue = ESCVoltage(Ub, nr, Ie, Iother, Rb)

#Battery Parameters
Ib = batteryCurrent(nr,Ie, Iother)
Thover = hoveringTime(Cb,Cmin,Ib)

print(f"Your propeller has size {Dp}x{Hp} [m] or {Dp/0.0254}x{Hp/0.0254} [inch]")
print(f"With a battery capacity of {Cb} [mAh], you can hover for {Thover} [min].")
print(Ib)