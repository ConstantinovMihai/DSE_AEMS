import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
from numpy.linalg import cholesky, det
from scipy.linalg import solve_triangular
from scipy.optimize import minimize
from matplotlib import animation, cm
import math
from sqlalchemy import Constraint
from domains import Point
from domains import Domain
from domains import Rect
import matplotlib.pyplot as plt
import time

#TODO Implement Proper Total Distance Travelled Calculation
#TODO Implement Averages for ScatterPlotGenerate
#TODO SumGaussianPlume giving wrong values


def gaussianPlumeInstant(receiverPosition: np.array, sourcePosition: np.array, time: float, h: float, H:float, aircraftSpeed: np.array, windVector: np.array, Q: float =10000):  # x,y coordinates relative to source, z in real coordinates
    """
    Function called by sumGaussianPlume
    Calculates the concentration as a result of the plume at a specific instant in time
    :param receiverPosition: specific position of measurement of value e.g. drone
    :param sourcePosition: specific position of the source of the emission e.g. aircraft
    :param time: time of event
    :param h: height of the source with respect to the ground
    :param H: float
    :param aircraftSpeed: np.array[xVelocity, yVelocity, zVelocity]
    :param windVector: np.array[xVelocity, yVelocity, zVelocity]
    :param Q: float
    :return:
    """
    #calculate distance between reciever and position
    distance = receiverPosition - sourcePosition
    #aircraft plume dispersion
    alpha = 0.02
    wind = aircraftSpeed*math.exp(-alpha*distance[0])
    #create overall dispersion
    wind += windVector
    #y axis wind flipped due to weird flipped behavior
    wind[1] = -wind[1]
    #dispersion coefficients
    sigma_x = 0.22*distance[0]*(1+0.0001*distance[0])**(-0.5)
    sigma_y = sigma_x
    sigma_z = 0.2*distance[0]
    #Generate Gaussian Plume Model based on Kazakh Paper
    Q = 10000
    Q1 = Q/((2*math.pi)**(1.5)*sigma_x*sigma_y*sigma_z)
    Qx = math.exp(-((distance[0])*wind[0] - (distance[1])*wind[1] - (wind[0]**2 + wind[1]**2)*(time))**2/(2*(sigma_x**2)*(wind[0]**2 + wind[1]**2)))
    Qy = math.exp(-(distance[0]*wind[1] + distance[1]*wind[0])**2/(2*(sigma_y**2)*(wind[0]**2 + wind[1]**2)))
    Qz = math.exp(-(distance[2] + h)**2/(2*sigma_z**2)) + math.exp(-(distance[2] - h + 2*H)**2/(2*sigma_z**2))
    C = Q1*Qx*Qy*Qz
    #Return concentration - scaled by power of 0.1 for easier analysisof plot
    return C**0.1

def sumGaussianPlume(receiverPosition: np.array, aircraftPositionEvent: np.array, aircraftEventTime: float, timeStep: float, hEvent: np.array, aircraftSpeedEvent: np.array, windVector:np.array, QEvent: np.array ,H:float = 0):
    """
    Calculates the total Gaussian Plume experienced by the reciever due to the source moving in time
    :param receiverPosition: constant position of reciever over time
    :param aircraftPositionEvent: array containing position of aircraft over time contained as a matrix
    :param timeStep: total time step of the process
    :param h: height of the source with respect to the ground (considered as changing over time due to takeoff/landing)
    :param H: mixing layer Height (set to 0)
    :param aircraftSpeedEvent: speed of aircraft throughtout event
    :param windVector: speed of wind set as a constant vector throughtout event
    :param QEvent: source of emissions over time
    :return:
    """
    C = 0
    time = aircraftEventTime
    i = 0
    #sum over time
    for aircraftPosition in aircraftPositionEvent:
        #H is set to zero
        C += timeStep * gaussianPlumeInstant(receiverPosition, aircraftPosition, time, hEvent[i], 0, aircraftSpeedEvent[i], windVector, QEvent[i])
        time -= timeStep #as you move through the aircraft event the time between the measurement and the emission decreases
        i += 1
    return C

def kernel(X1, X2, l=1.0, sigma_f=1.0):
    """
    Isotropic squared exponential kernel.

    Args:
        X1: Array of m points (m x d).
        X2: Array of n points (n x d).

    Returns:
        (m x n) matrix.
    """
    #implementation of RBF kernel
    sqdist = np.sum(X1 ** 2, 1).reshape(-1, 1) + np.sum(X2 ** 2, 1) - 2 * np.dot(X1, X2.T)
    kernelRBF = sigma_f ** 2 * np.exp(-0.5 / l ** 2 * sqdist) #at 200 iter RMSE 0.152
    dist = np.sqrt(sqdist)
    #implementation of Matern kernel
    kernelMatern = sigma_f ** 2 * (1+ np.sqrt(5)*dist/l + 5*dist**2/(3*l**2))*np.exp(-5*dist/l) #at 200 iter RMSE 0.07784
    return kernelRBF

def posterior(X_s, X_train, Y_train, l=1.0, sigma_f=1.0, sigma_y=1e-8):
    """
    Computes the suffifient statistics of the posterior distribution
    from m training data X_train and Y_train and n new inputs X_s.

    Args:
        X_s: New input locations (n x d).
        X_train: Training locations (m x d).
        Y_train: Training targets (m x 1).
        l: Kernel length parameter.
        sigma_f: Kernel vertical variation parameter.
        sigma_y: Noise parameter.

    Returns:
        Posterior mean vector (n x d) and covariance matrix (n x n).
    """
    K = kernel(X_train, X_train, l, sigma_f) + sigma_y ** 2 * np.eye(len(X_train))
    K_s = kernel(X_train, X_s, l, sigma_f)
    K_ss = kernel(X_s, X_s, l, sigma_f) + 1e-8 * np.eye(len(X_s))
    K_inv = inv(K)

    # Equation (7)
    mu_s = K_s.T.dot(K_inv).dot(Y_train)

    # Equation (8)
    cov_s = K_ss - K_s.T.dot(K_inv).dot(K_s)

    return mu_s, cov_s


def nll_fn(X_train, Y_train, noise, naive=False):
    """
    Returns a function that computes the negative log marginal
    likelihood for training data X_train and Y_train and given
    noise level.

    Args:
        X_train: training locations (m x d).
        Y_train: training targets (m x 1).
        noise: known noise level of Y_train.
        naive: if True use a naive implementation of Eq. (11), if
               False use a numerically more stable implementation.

    Returns:
        Minimization objective.
    """

    Y_train = Y_train.ravel()

    def nll_naive(theta):
        # Naive implementation of Eq. (11). Works well for the examples
        # in this article but is numerically less stable compared to
        # the implementation in nll_stable below.
        K = kernel(X_train, X_train, l=theta[0], sigma_f=theta[1]) + \
            noise ** 2 * np.eye(len(X_train))
        return 0.5 * np.log(det(K)) + \
               0.5 * Y_train.dot(inv(K).dot(Y_train)) + \
               0.5 * len(X_train) * np.log(2 * np.pi)

    def nll_stable(theta):
        # Numerically more stable implementation of Eq. (11) as described
        # in http://www.gaussianprocess.org/gpml/chapters/RW2.pdf, Section
        # 2.2, Algorithm 2.1.

        K = kernel(X_train, X_train, l=theta[0], sigma_f=theta[1]) + \
            noise ** 2 * np.eye(len(X_train))
        L = cholesky(K)

        S1 = solve_triangular(L, Y_train, lower=True)
        S2 = solve_triangular(L.T, S1, lower=False)

        return np.sum(np.log(np.diagonal(L))) + \
               0.5 * Y_train.dot(S2) + \
               0.5 * len(X_train) * np.log(2 * np.pi)

    if naive:
        return nll_naive
    else:
        return nll_stable


# Distance-based Upper Confidence Bound
def UCB(X_2D, X_2D_train, Y_2D_train, noise_2D, kappa):
    """
    DUCB acquisition function
    :param X_2D: Total coordinates
    :param X_2D_train: coordinates sampled
    :param Y_2D_train: value of coordinates sampled
    :param noise_2D: noise of data
    :param kappa: #Exploration/Exploitation trade-off
    :return: total value of ucb value
    """
    #Assign mean and variance to each point in the domain based on concentration
    res = minimize(nll_fn(X_2D_train, Y_2D_train, noise_2D), [1, 1],
                   bounds=((1e-5, None), (1e-5, None)),
                   method='L-BFGS-B')
    mu_s, cov_s = posterior(X_2D, X_2D_train, Y_2D_train, *res.x, sigma_y=noise_2D)
    mu_global  = np.mean(mu_s)
    # return the UCB metric formula
    # weird normalisation in mu_global takes place to help intuitively understand UCB
    return mu_s/mu_global + kappa * np.sqrt(np.diag(cov_s))

def DUCB(X_2D, X_2D_train, Y_2D_train_concentrations, noise_2D, kappa, gamma, domain : Domain):
    """
    Compute total DUCB acquisition function for all concentrations
    :param X_2D: Total coordinates
    :param X_2D_train: coordinates sampled
    :param Y_2D_train_concentrations: value of coordinates sampled with all concentrations
    :param noise_2D: noise of data
    :param kappa: #Exploration/Exploitation
    :param domain: Domain for distance function
    :return: total DUCB function
    """
    # compute total acquistionFunc
    acquisitionFunc = 0
    # go through the acquisitionFunc for each concentration
    for Y_2D_train in Y_2D_train_concentrations:
        tempAcquisitionFunc = UCB(X_2D, X_2D_train, Y_2D_train, noise_2D, kappa)
        #normalise acquisitionFunction for each concentration
        meanAcquisitionFunc = np.mean(tempAcquisitionFunc)
        acquisitionFunc += tempAcquisitionFunc/meanAcquisitionFunc

    #compute distance to travel to all locations within domain
    lastMeasurement = X_2D_train[-1]
    distance = []
    for point in X_2D:
        p = Point(point[0], point[1], point[2])
        l = Point(lastMeasurement[0], lastMeasurement[1], lastMeasurement[2])
        distance.append(domain.computeDistance(p, l))
    distance = np.array(distance)
    # return the DUCB metric formula
    return acquisitionFunc + gamma*distance

def RMSE(mu, mu_s):
    """
    compute error of total model from real values
    :param mu: true mean of model
    :param mu_s: estimated mean of model
    :return: root mean square error
    """
    rmse = 0
    n = len(mu)
    for i in range(n):
        rmse += (mu[i]-mu_s[i])**2
    return math.sqrt(rmse/n)

#To be fixed
def distanceTravelled(X_2D_train : np.array):
    """
    Computes total distance travelled
    :param X_2D_train: locations sampled
    :return: euclidian distance
    """
    dist = 0
    X0 = X_2D_train[0][0]
    X1 = X_2D_train[0][1]
    X2 = X_2D_train[0][2]
    for X in X_2D_train:
        dist += math.sqrt((X[0] - X0)**2 + (X[1]-X1)**2 + (X[2]-X2)**2)
        X0 = X[0]
        X1 = X[1]
        X2 = X[2]
    return dist



def generateMesh(xDomain: np.array, yDomain: np.array, zDomain: np.array, constr1: np.array, constr2: np.array):
    """
    Generates discretised mesh based on domains and constraints
    :param xDomain: np.array([xMin, xMax, dX])
    :param yDomain: np.array([yMin, yMax, dY])
    :param zDomain: np.array([zMin, zMax, dZ])
    :param constr1: np.array([xPoint1, yPoint1, zPoint1]) set top lower left most corner of domain
    :param constr2: np.array([xPoint2, yPoint2, zPoint2]) set back top right most corner of damin
    :return: discretised domain points
    """
    rx, ry, rz = np.arange(xDomain[0], xDomain[1], xDomain[2]), np.arange(yDomain[0], yDomain[1], yDomain[2]), np.arange(zDomain[0], zDomain[1], zDomain[2])
    gx, gy, gz = np.meshgrid(rx, ry, rz)
    tempX_3D = np.c_[gx.ravel(), gy.ravel(), gz.ravel()]
    status = False
    for location in tempX_3D:
        if not(constr1[0] <= location[0] <= constr2[0] and constr1[1] <= location[1] <= constr2[1] and constr1[2] <= location[2] <= constr2[2]):
            if status:
                X_3D = np.vstack((X_3D, location))
            else:
                X_3D = np.array(location)
                status = True

    return X_3D


def ScatterPlotGenerate3D(X_3D: np.array):
    """
    Generates scatter plot
    :param X_3D: domain
    :return: plot
    """
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X_3D[:, 0], X_3D[:, 1], X_3D[:, 2])
    ax.set_xlabel('x label')
    ax.set_ylabel('y label')
    ax.set_ylabel('z label')
    ax.set_xlim3d(min(X_3D[:, 0]), max(X_3D[:, 0]))
    ax.set_ylim3d(min(X_3D[:, 1]), max(X_3D[:, 1]))
    ax.set_zlim3d(min(X_3D[:, 2]), max(X_3D[:, 2]))
    plt.show()


def randomGenerate(X_3D: np.array, nRand: int):
    """
    Returns random number of samples
    :param X_3D: discretised mesh
    :param nRand: number of random points needed
    :return: return sample point
    """
    return X_3D[np.random.randint(len(X_3D), size=nRand)]

def generateAircraftTakeOffEvent(startLocation : np.array, totalTime: float, timeStep: float):
    """
    Simulate aircraft takeOff Event
    :param startLocation: select location from which the aircraft takesoff from
    :param totalTime: total time of interest during measurement
    :param timeStep: time step of measurements
    :return: aircraft position over time, aircraft speed over time, height of emission source over time
    """
    a = 6 #aircraft acceleration set as constant based on paper "Simplified Algorithm to Model Aircraft Acceleration During Takeoff"
    takeOffVelocity = 70 #~737 TakeOff Velocity
    takeOffGradient = 3 #Above FAA legislation for two engine aircraft
    velocity = 0
    nStep = int(totalTime/timeStep)
    h = 5
    aircraftPosition = startLocation
    aircraftPositionEvent = np.array([])
    hEvent = []
    aircraftSpeedEvent = []
    for i in range(nStep):
        #compute velocity over time
        velocity += a*timeStep
        aircraftSpeedEvent.append([velocity, 0, 0])
        #compute position over time
        aircraftPosition -= np.array([velocity, 0, 0])
        aircraftPositionEvent = np.concatenate((aircraftPositionEvent, aircraftPosition))
        #compute height over time
        if velocity > takeOffVelocity:
            h += velocity*timeStep*takeOffGradient/100
        hEvent.append(h)
    return aircraftPositionEvent.reshape(nStep, 3), np.array(aircraftSpeedEvent), np.array(hEvent)


def pollutionAircraftEvent(aircraftStartLocation: np.array, receiverLocation : np.array, totalTime: float, timeStep: float):
    """
    Calculates the total pollution experienced by a point over time and integrates the generation of aircraft events
    :param aircraftStartLocation:
    :param receiverLocation:
    :param totalTime:
    :param timeStep:
    :return: list with polluton levels over time
    """
    Q = np.array([10000]*int(totalTime/timeStep))
    windVector = np.array([0,5,0])
    aircraftPositionEvent, aircraftSpeedEvent, hEvent = generateAircraftTakeOffEvent(aircraftStartLocation, totalTime, timeStep)
    C = np.array([])
    timeArray = np.linspace(1, totalTime, int(totalTime/timeStep))
    for time in timeArray:
        print(time)
        tempC = np.array([sumGaussianPlume(receiverLocation, aircraftPositionEvent, time, timeStep, hEvent, aircraftSpeedEvent, windVector, Q)])
        C=np.concatenate((C, tempC))
    return C
aircraftStartLocation = np.array([200.0, 0, 0])
receiverLocation = np.array([300,0,0])
totalTime = 10
timeStep = 1
print(pollutionAircraftEvent(aircraftStartLocation, receiverLocation, totalTime, timeStep))
"""
startLocation = np.array([100.0, 0, 0])
totalTime = 15
timeStep = 0.5
print(generateAircraftTakeOffEvent(startLocation, totalTime, timeStep))
minX, maxX, dX = 0.1, 40, 2
minY, maxY, dY = -20, 20, 2
minZ, maxZ, dZ = 0, 30, 2
constr1 = [0,-10,10]
constr2 = [20,10,20]
X_3D = generateMesh([minX, maxX, dX],[minY, maxY, dY], [minZ, maxZ, dZ], constr1, constr2)
#ScatterPlotGenerate3D(X_3D)
print(randomGenerate(X_3D, 5))

aircraftPositionEvent = np.array([[0,0,0], [2,0,0]])
aircraftSpeedEvent = np.array([[100,0,0], [100,0,0]])
timeEvent = 1
aircraftSpeedEvent = np.array([[3,5,0], [3,5,0]])
QEvent = np.array([[10000],[10000]])
hEvent = np.array([[0.5], [0.5]])
C = sumGaussianPlume([10,0,0], aircraftPositionEvent, timeEvent, 0.5, hEvent,0, aircraftSpeedEvent, np.array([0,10,0]), QEvent)
print(C)
"""

def main3D():
    minX, maxX, dX = 0.1, 40, 2
    minY, maxY, dY = -20, 20, 2
    minZ, maxZ, dZ = 0, 30, 2
    constr1 = [0,-5,10]
    constr2 = [5,5,20]
    constr = Rect(Point(constr1[0],constr1[1], constr1[2]), Point(constr2[0], constr2[1], constr2[2]))
    domain = Domain(Point(minX, minY, minZ), Point(maxX, maxY, maxZ), constr)
    rx, ry, rz = np.arange(minX, maxX, dX), np.arange(minY, maxY, dY), np.arange(minZ, maxZ, dZ)
    gx, gy, gz = np.meshgrid(rx, ry, rz)
    X_3D = np.c_[gx.ravel(), gy.ravel(), gz.ravel()]
    # Explore/Exploit TradeOff
    # kappa = 15 #exploration/exploitation constant
    kappa = 3000
    #gamma = -0.1  # cost-to-evaluate
    gamma = 0
    initialSamples = 20 # random initial samples
    nIter = 10  # number of points selected by BO algorithm
    noise_3D = 0.01  # Needs a small noise otherwise kernel can become positive semi-definite which leads to minimise() not working

    X_3D_train = np.array([[np.random.uniform(minX, maxX), np.random.uniform(minY, maxY), np.random.uniform(minZ, maxZ)]])
    for i in range(initialSamples):
        X_3D_train = np.vstack((X_3D_train, [np.random.uniform(minX, maxX), np.random.uniform(minY, maxY), np.random.uniform(minZ, maxZ)]))

    Y_3D_train = []
    for points in X_3D_train:
        Y_3D_train.append(gaussianPlume(points[0], points[1], points[2]) + noise_3D * np.random.randn())
    Y_3D_train = np.array(Y_3D_train)

    Y_3D = []
    for points in X_3D:
        Y_3D.append(gaussianPlume(points[0], points[1], points[2]) + noise_3D * np.random.randn())
    Y_3D = np.array(Y_3D)

    for i in range(nIter):
        print("sampling number: ", i)
        sampleLocation = proposeLocation(X_3D, X_3D_train, Y_3D_train, noise_3D, kappa, gamma, domain)
        print("sample location: ", sampleLocation)
        X_3D_train = np.vstack((X_3D_train, [sampleLocation[0], sampleLocation[1], sampleLocation[2]]))
        Y_3D_train = np.hstack(
            (Y_3D_train, gaussianPlume(sampleLocation[0], sampleLocation[1], sampleLocation[2]) + noise_3D * np.random.randn()))


    print("Total Distance Travelled: ", distanceTravelled(X_3D_train[initialSamples:]))
    res = minimize(nll_fn(X_3D_train, Y_3D_train, noise_3D), [1, 1],
                   bounds=((1e-5, None), (1e-5, None)),
                   method='L-BFGS-B')

    mu_s, cov_s = posterior(X_3D, X_3D_train, Y_3D_train, *res.x, sigma_y=noise_3D)
    """
    print(RMSE(Y_3D, mu_s))
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(gx, gy, gz, c=mu_s)
    plt.show()"""
    cutPlot(X_3D, X_3D_train, Y_3D_train, noise_3D, gx, gy, gaussianPlume)
    plt.show()

#main3D()
