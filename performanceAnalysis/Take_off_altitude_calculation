import math
import matplotlib.pyplot as plt

### INPUT ###

#h = 1000

#############


p0 = 101325 # Pa
rho0 = 1.225 # kg/m^3
Temp0 = 288.15 # K
g0 = 9.80665 # m/s^2
alpha = -0.0065 # K/m
R = 287.04

A = 5
epsilon = 0.85
lambdav = 0.75
xi = 0.5
e = 0.83
Cfd = 0.015
alpha0 = 0
K0 = 6.11
Bp = 4
Hp = 0.1524
Dp = 0.40
W = 5.85*g0*2
N = 6000 # [RPM]

T = 150
h = 0

h_lst = []
T_lst = []
W_lst = []
W2_lst = []
W2 = 5.85*g0

for h in range(0, 15000, 1):
    if h <= 11000:
        Temp = Temp0 + alpha * h
        p = p0 * (1 + alpha * ( h / Temp0 ) )**5.2561
        rho = p / (R * Temp)

    if h >= 11000 and h <= 20000:
        Temp11 = Temp0 + alpha * 11000
        p11 = p0 * (1 + alpha * ( 11000 / Temp0 ) ) ** 5.2561
        h11 = 11000
        Temp = 216.65
        p = p11 * math.exp(-g0 * (h - h11 )/(R*Temp11))
        rho = p / (R * Temp)

    CT = 0.25*math.pi**3*lambdav*xi**2*Bp*K0*(epsilon*math.atan(Hp/(math.pi*Dp))-alpha0)/(math.pi*A+K0)
    T = 4 * CT * rho * (N/60)**2 * Dp ** 4

    h_lst.append(h)
    T_lst.append(T)
    W_lst.append(W)
    W2_lst.append(W2)
    
plt.plot(h_lst, T_lst, label = "Lift")
plt.plot(h_lst, W_lst, label = "L/MTOW = 2.0")
plt.plot(h_lst, W2_lst, label = "L/MTOW = 1.0")
plt.xlabel("Altitude [m] of take-off ")
plt.ylabel("Force [N]")
plt.grid(True)
plt.legend()
plt.axis([0, 15000, 0, 220])
plt.show()

run = True
i = 0
while run == True:
    i += 1
    Tlocal = T_lst[i]
    Wlocal = W_lst[i]
    hlocal = h_lst[i]
    if Tlocal <= Wlocal:
        run = False
        print("Max altitude with L/MTOW = 2: ", hlocal, "[m]")

running = True
j = 0
while running == True:
    j += 1
    Tlocal2 = T_lst[j]
    Wlocal2 = W2_lst[j]
    h2local = h_lst[j]
    if Tlocal2 <= Wlocal2:
        running = False
        print("Max altitude with L/MTOW = 1: ", h2local, "[m]")
