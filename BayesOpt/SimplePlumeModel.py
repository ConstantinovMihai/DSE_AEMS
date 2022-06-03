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
   # sigma_y = 0.01
    sigma_z = 10000
    Q = (constants[1][0] * variables[0] + constants[1][1] * variables[1] + constants[1][2] * variables[2] //
         + constants[1][3] * variables[3] + constants[1][4] * variables[4] //
         + constants[1][5] * variables[5]) ** constants[1][6]  # source
    Q = 100
    C = Q / speed * (1 / (2 * np.pi * sigma_y * sigma_z)) * math.exp(-0.5 * (y / sigma_y) ** 2) * (
                math.exp(-0.5 * ((z - H) / sigma_z) ** 2) + math.exp(-0.5 * ((z + H) / sigma_z) ** 2))
    return C



def gaussianPlumeModel(x : np.array, y : np.array, z : np.array, H : np.array, vel : np.array, wind : np.array, Q : np.array):
    """
    Returns a 2d np.array which represent the Gaussian plume model

    Args:
        x (np.array) - contains the x coordinates of the receiver wrt to the source 
        y (np.array) - contains the y coordinates of the receiver wrt to the source 
        z (np.array) - contains the z - height of the receiver wrt to the source 
        H (np.array) - contains the Height coordinates of the source wrt the ground
        vel (np.array) - velocity of the aircraft
        wind (np.array) - velocity of the wind
        Q (np.array) - intensity of the source
    """

    def generateSigmas(x : np.array, const = np.array([0.08, 0.0001, 0.5, 0.06, 0.0015, 0.5])):
        """
        Returns the standard deviation for the atmospheric model
        Args:
            x (np.array) - contains the x coordinates of the receiver wrt to the source 
            const (np.array) - contains the semi empirical constants (see de Vischer book) 
        """
        # sanity check that the const array has the correct dimension
        assert(len(const) == 6)

        sigma_y = const[0] * x / np.power((1 + const[1] * x), const[2])
        sigma_z = const[3] * x / np.power((1 + const[4] * x), const[5])

        return sigma_y, sigma_z


    # generate the standard deviations for the model
    sigma_y, sigma_z = generateSigmas(x)
    
    # advection in the x direction
    C_x = Q / vel * (1 / (2 * np.pi * sigma_y * sigma_z))

    # diffusion on the y direction
    phi_y = np.exp(-0.5 * (y / sigma_y) ** 2)
    
    # diffusion in the z direction
    phi_z = (np.exp(-0.5 * ((z - H) / sigma_z) ** 2) + np.exp(-0.5 * ((z + H) / sigma_z) ** 2))

    # the concentration vector (final output of the gaussian plume model)
    C = C_x * phi_y * phi_z
    
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
    windSpeed = 100
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

if __name__ == "__main__":
    main()
