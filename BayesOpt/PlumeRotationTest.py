import numpy as np
import matplotlib.pyplot as plt
import math

"""
def gaussianPlume(x, y, aircraftSpeed, windSpeed):  # x,y coordinates relative to source, z in real coordinates
    # variables = [aType, aPath, aEvent, thrust, temperature, humidity];
    # constants = [alpha, beta]
    # wind = [uX, uY, uZ]
    y = y + abs(y)*windSpeed/aircraftSpeed
    speed = aircraftSpeed
    z = 5
    H = 5
    Q = 100
    const = np.array([0.08, 0.0001, 0.5, 0.06, 0.0015, 0.5])
    sigma_y = const[0] * x / np.power((1 + const[1] * x), const[2])
    sigma_z = const[3] * x / np.power((1 + const[4] * x), const[5])
    C = Q / speed * (1 / (2 * np.pi * sigma_y * sigma_z)) * math.exp(-0.5 * (y / sigma_y) ** 2) * (
            math.exp(-0.5 * ((z - H) / sigma_z) ** 2) + math.exp(-0.5 * ((z + H) / sigma_z) ** 2))

    return C**0.1

"""

def gaussianPlume(recieverPosition , sourcePosition, time, h, H, aircraftSpeed, windVector):  # x,y coordinates relative to source, z in real coordinates
    # variables = [aType, aPath, aEvent, thrust, temperature, humidity];
    # constants = [alpha, beta]
    # wind = [uX, uY, uZ]

    distance = recieverPosition - sourcePosition
    alpha = 1
    wind = aircraftSpeed*math.exp(-alpha*distance[0])
    wind += windVector
    sigma_x = 0.22*distance[0]*(1+0.0001*distance[0])**(-0.5)
    sigma_y = sigma_x
    sigma_z = 0.2*distance[0]
    Q = 100
    Q1 = Q/((2*math.pi)**(1.5)*sigma_x*sigma_y*sigma_z)
    Qx = math.exp(-((distance[0])*wind[0] - (distance[1])*wind[1] - (wind[0]**2 + wind[1]**2)*(time))**2/(2*(sigma_x**2)**(wind[0]**2 + wind[1]**2)))
    Qy = math.exp(-((distance[0])*wind[1] + (distance[1])*wind[0])**2/(2*(sigma_y**2)**(wind[0]**2 + wind[1]**2)))
    Qz = math.exp(-(distance[2] + h)**2/(2*sigma_z**2)) + math.exp(-(distance[2] + h + 2*H)**2/(2*sigma_z**2))
    C = Q1*Qx*Qy*Qz
    return C

minX = 0.1
maxX = 100
minY = -50
maxY = 50
dX = 0.5
dY = 0.5

aircraftSpeed = np.array([200, 0, 0])
windSpeed = np.array([3, 5, 1])


rx, ry = np.arange(minX, maxX, dX), np.arange(minY, maxY, dY)
gx, gy = np.meshgrid(rx, ry)

# Generate Initial Samples
X_2D = np.c_[gx.ravel(), gy.ravel()]

Y_2D = []
for points in X_2D:
    positionReciever = np.array([points[0], points[1], 5])
    positionSource = np.array([0, 0, 0])
    Y_2D.append(gaussianPlume(positionReciever,positionSource, 2, 5, 5,  aircraftSpeed, windSpeed))
Y_2D = np.array(Y_2D)
plt.imshow(Y_2D.reshape(gx.shape).T, interpolation="nearest",vmin= 0.3, vmax=1, origin="upper")
plt.show()

