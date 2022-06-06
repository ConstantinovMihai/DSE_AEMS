import numpy as np
import matplotlib.pyplot as plt

def gaussianPlume(x, y):  # x,y coordinates relative to source, z in real coordinates
    # variables = [aType, aPath, aEvent, thrust, temperature, humidity];
    # constants = [alpha, beta]
    # wind = [uX, uY, uZ]
    z = 5
    H = 5
    speed = 10
    Q = 100
    const = np.array([0.08, 0.0001, 0.5, 0.06, 0.0015, 0.5])
    sigma_y = const[0] * x / np.power((1 + const[1] * x), const[2])
    sigma_z = const[3] * x / np.power((1 + const[4] * x), const[5])
    C = Q / speed * (1 / (2 * np.pi * sigma_y * sigma_z)) * math.exp(-0.5 * (y / sigma_y) ** 2) * (
            math.exp(-0.5 * ((z - H) / sigma_z) ** 2) + math.exp(-0.5 * ((z + H) / sigma_z) ** 2))

    return C**0.1

minX = 0.1
maxX = 100
minY = -20
maxY = 20
dX = 1
dY = 1


rx, ry = np.arange(minX, maxX, dX), np.arange(minY, maxY, dY)
gx, gy = np.meshgrid(rx, ry)

# Generate Initial Samples
X_2D = np.c_[gx.ravel(), gy.ravel()]

Y_2D = []
for points in X_2D:
    Y_2D.append(gaussianPlume(points[0], points[1]))
Y_2D = np.array(Y_2D)
plt.imshow(Y_2D.reshape(gx.shape).T)
plt.show()


