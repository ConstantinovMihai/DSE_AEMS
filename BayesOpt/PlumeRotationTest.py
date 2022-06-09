import numpy as np
import matplotlib.pyplot as plt
import math

import matplotlib.animation as ma


def gaussianPlumeRot(recieverPosition , sourcePosition, time, h, H, aircraftSpeed, windVector):  # x,y coordinates relative to source, z in real coordinates
    # variables = [aType, aPath, aEvent, thrust, temperature, humidity];
    # constants = [alpha, beta]
    # wind = [uX, uY, uZ]

    distance = recieverPosition - sourcePosition
    alpha = 0.02
    wind = aircraftSpeed*math.exp(-alpha*distance[0])
    wind += windVector
    wind[1] = -wind[1]
    sigma_x = 0.22*distance[0]*(1+0.0001*distance[0])**(-0.5)
    sigma_y = sigma_x
    sigma_z = 0.2*distance[0]
    Q = 10000
    Q1 = Q/((2*math.pi)**(1.5)*sigma_x*sigma_y*sigma_z)
    Qx = math.exp(-((distance[0])*wind[0] - (distance[1])*wind[1] - (wind[0]**2 + wind[1]**2)*(time))**2/(2*(sigma_x**2)*(wind[0]**2 + wind[1]**2)))
    Qy = math.exp(-(distance[0]*wind[1] + distance[1]*wind[0])**2/(2*(sigma_y**2)*(wind[0]**2 + wind[1]**2)))
    Qz = math.exp(-(distance[2] + h)**2/(2*sigma_z**2)) + math.exp(-(distance[2] - h + 2*H)**2/(2*sigma_z**2))
    C = Q1*Qx*Qy*Qz
    return C**0.1

minX = 0.1
maxX = 400
minY = -100
maxY = 100
dX = 2
dY = 2

aircraftSpeed = np.array([100, 0, 0])
windSpeed = np.array([0, 20, 0])


rx, ry = np.arange(minX, maxX, dX), np.arange(minY, maxY, dY)
gx, gy = np.meshgrid(rx, ry)

# Generate Initial Samples
X_2D = np.c_[gx.ravel(), gy.ravel()]
print(X_2D.reshape(gx.shape).shape)
time = [0,1,2,3]
Y_2D = []
Y_2DTime = []
for t in time:
    for points in X_2D:
        positionReciever = np.array([points[0], points[1], 5])
        positionSource = np.array([0, 0, 0])
        Y_2D.append(gaussianPlumeRot(positionReciever,positionSource, t, 2, 5,  aircraftSpeed, windSpeed))
    #Y_2D = np.array(Y_2D)
    Y_2DTime.append(Y_2D)
Y_2DTime = np.array(Y_2DTime)
print(Y_2DTime.shape)




"""
plt.imshow(Y_2D.reshape(gx.shape).T, interpolation="nearest",vmin= 0, vmax=1, origin="upper")
plt.xlabel("Y axis [m]")
plt.ylabel("X axis [m]")
plt.show()"""
