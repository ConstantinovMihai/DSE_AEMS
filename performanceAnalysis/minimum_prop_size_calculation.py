import math
rp0 = 0.70 # [m]
rho = 1.225 # [kgm^-3]
mass = 5.85 # [kg]
arms = 4
nt = arms

# Typical values
A = 6.5
epsilon = 0.9
lambdav = 0.8
ksi = 0.55
alpha0 = -(math.pi/36)
K_0 = 6.11
B_p = 4 # [-] Blade number
H_p = 0.1143 # propeller pitch

Ni = 4000 # [RPM]

rp_lst = [rp0]
deltarp = 50

while deltarp > 0.01:
    rp = rp_lst[-1]
    TincltwoLW = 2*mass*9.81/4 # N, includes 2.0 L/W factor required
    CT = 0.25*math.pi**3*lambdav*ksi**2*B_p*K_0*(epsilon*math.atan(H_p/(math.pi*2*rp))-alpha0)/(math.pi*A+K_0)
    rp_new = 0.5 * (TincltwoLW/(CT * rho * (Ni/60)**2))**0.25
    rp_lst.append(rp_new)
    deltarp = 100*abs(rp_lst[-1] - rp_lst[-2])/rp_lst[-2]

print("Minimum propeller radius = ", rp_lst[-1], " [m]")
