import numpy as np
import matplotlib.pyplot as plt
import math

np.random.seed(2443)


def gaussianPlume(x, y, z, H, variables, constants, aircraftSpeed, wind): #x,y coordinates relative to source, z in real coordinates
    # variables = [aType, aPath, aEvent, thrust, temperature, humidity];
    # constants = [alpha, beta]
    # wind = [uX, uY, uZ]

    speed = aircraftSpeed + wind
    sigma_y = constants[0][0] * x / (1 + constants[0][1] * x) ** constants[0][2]  # gaussian diffusion for y coordinate
    sigma_z = constants[0][3] * x / (1 + constants[0][4] * x) ** constants[0][5]  # gaussian diffusion for z coordinate
    Q = (constants[1][0] * variables[0] + constants[1][1] * variables[1] + constants[1][2] * variables[2] //
         + constants[1][3] * variables[3] + constants[1][4] * variables[4] //
         + constants[1][5] * variables[5]) ** constants[1][6]  # source
    C = Q / speed * (1 / (2 * np.pi * sigma_y * sigma_z)) * math.exp(-0.5 * (y / sigma_y) ** 2) * (
                math.exp(-0.5 * ((z - H) / sigma_z) ** 2) + math.exp(-0.5 * ((z + H) / sigma_z) ** 2))
    return C

def main():
    # random constant generation
    alpha = []
    for i in range(6):
        alpha.append(np.random.uniform(0,1))
    beta = []
    for i in range(7):
        beta.append(np.random.uniform(0,1))

    # initialise coordinate system
    X = np.arange(0.1, 10, 0.1)
    Y = np.arange(-5, 5, 0.1)
    Z = 10
    H = 5
    sourceLocation = [0, 0]  # only x and y as H gives Z

    # initialise variables
    windSpeed = 10
    variables = [1000, 1000, 1000, 1000, 1000, 1000]
    aircraftSpeed = 100

    fDomain = []
    for x in X:
        #aircraftSpeed -= 1
        for y in Y:
            fDomain.append(
                gaussianPlume(x - sourceLocation[0], y - sourceLocation[1], Z, H, variables, [alpha, beta], aircraftSpeed, windSpeed))

    fDomain = np.array(fDomain)
    fDomain = np.reshape(fDomain, (99, 100))

    print(fDomain.shape)
    print(fDomain)
    plt.imshow(fDomain, interpolation="nearest", origin="upper")
    plt.gca().invert_yaxis()
    plt.colorbar()
    plt.show()

main()
