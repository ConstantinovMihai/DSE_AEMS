import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
from numpy.linalg import cholesky, det
from scipy.linalg import solve_triangular
from scipy.optimize import minimize
from matplotlib import animation, cm
import math
from copy import deepcopy

from sqlalchemy import Constraint

from domains import Point
from domains import Domain
from domains import Rect

#new imports for 3D plot
import matplotlib.pyplot as plt

np.random.seed(2314)


def gaussianPlume(x, y, z):  # x,y coordinates relative to source, z in real coordinates
    # variables = [aType, aPath, aEvent, thrust, temperature, humidity];
    # constants = [alpha, beta]
    # wind = [uX, uY, uZ]
    H = 2
    speed = 10
    Q = 100
    const = np.array([0.08, 0.0001, 0.5, 0.06, 0.0015, 0.5])
    sigma_y = const[0] * x / np.power((1 + const[1] * x), const[2])
    sigma_z = const[3] * x / np.power((1 + const[4] * x), const[5])
    C = Q / speed * (1 / (2 * np.pi * sigma_y * sigma_z)) * math.exp(-0.5 * (y / sigma_y) ** 2) * (
            math.exp(-0.5 * ((z - H) / sigma_z) ** 2) + math.exp(-0.5 * ((z + H) / sigma_z) ** 2))

    return C**0.15


def kernel(X1, X2, l=1.0, sigma_f=1.0):
    """
    Isotropic squared exponential kernel.

    Args:
        X1: Array of m points (m x d).
        X2: Array of n points (n x d).

    Returns:
        (m x n) matrix.
    """
    sqdist = np.sum(X1 ** 2, 1).reshape(-1, 1) + np.sum(X2 ** 2, 1) - 2 * np.dot(X1, X2.T)
    kernelRBF = sigma_f ** 2 * np.exp(-0.5 / l ** 2 * sqdist) #at 200 iter RMSE 0.152
    dist = np.sqrt(sqdist)
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
def DUCB(X_2D, X_2D_train, Y_2D_train, noise_2D, kappa, gamma, position, domain : Domain):
    """
    DUCB acquisition function

    Args:
        :param:
        :param:

    Returns
    """
    res = minimize(nll_fn(X_2D_train, Y_2D_train, noise_2D), [1, 1],
                   bounds=((1e-5, None), (1e-5, None)),
                   method='L-BFGS-B')
    mu_s, cov_s = posterior(X_2D, X_2D_train, Y_2D_train, *res.x, sigma_y=noise_2D)
    mu_global  = np.mean(mu_s)

    lastMeasurement = position
    distance = []
    for point in X_2D:
        p = Point(point[0], point[1], point[2])
        l = Point(lastMeasurement[0], lastMeasurement[1], lastMeasurement[2])
        distance.append(domain.computeDistance(p, l))
    distance = np.array(distance)
    # return the DUCB metric formula
    # weird normalisation in mu_global takes place
    return mu_s/mu_global + kappa * np.sqrt(np.diag(cov_s)) + gamma*distance

def proposeLocation(X_2D, X_2D_train, Y_2D_train, noise_2D, kappa, gamma, workers_loc, domain: Domain):
    """
    Proposes location for the next batch of experiments
    Args:
        :param: workers_loc (np.array) - current locations of the drones
        :param: domain (Domain) - the domain class

    Returns a numpy array containing the proposed locations
    """
    # total number of workers
    nb_workers = len(workers_loc)

    # initially all workers are idle
    idle = np.ones(nb_workers)

    # copy the training data, which will be appended with expected values in
    # order to generate predictions for all points in the workers' batch
    x_train = deepcopy(X_2D_train)
    y_train = deepcopy(Y_2D_train)

    # proposed sampling locations
    sample_locs = np.zeros(nb_workers)
    workerLocation = np.array([])
    # while there are still idle workers, continue iterating
    while np.count_nonzero(idle):
        acquisitionFuncs = np.zeros((nb_workers, len(X_2D)))
        max_acq_func = np.zeros(nb_workers)

        for worker_idx in range(nb_workers):
            if idle[worker_idx]:
                acquisitionFuncs[worker_idx] = DUCB(X_2D, x_train, y_train, noise_2D, kappa, gamma, workers_loc[worker_idx], domain)
                max_acq_func[worker_idx] = np.argmax(acquisitionFuncs[worker_idx])
                #acquisitionFuncs[worker_idx] = DUCB(X_2D, x_train, y_train, noise_2D, kappa, gamma,
                                                    #workers_loc[worker_idx], domain)
                #max_acq_func[worker_idx] = np.argmax(DUCB(X_2D, x_train, y_train, noise_2D, kappa, gamma, workers_loc[worker_idx], domain))
                # if the new acquisition function is the highest one so far encountered

        # select the next sampling location
        assigned_worker = np.argmax(max_acq_func)
        idle[assigned_worker] = 0
        sample_locs[assigned_worker] = int(np.argmax(acquisitionFuncs[assigned_worker]))
        x_proposed = X_2D[int(sample_locs[assigned_worker])]
        workerLocation =
        # update the X_train, Y_train arrays
        # Assign mean and variance to each point in the domain based on concentrations
        res = minimize(nll_fn(X_2D_train, Y_2D_train, noise_2D), [1, 1],
                       bounds=((1e-5, None), (1e-5, None)),
                       method='L-BFGS-B')

        mu_s, cov_s = posterior(X_2D, X_2D_train, Y_2D_train, *res.x, sigma_y=noise_2D)
        # the expected value for the measurement of the proposed point
        y_expected = mu_s[int(sample_locs[assigned_worker])]

        x_train = np.vstack((x_train, [x_proposed[0], x_proposed[1], x_proposed[2]]))
        y_train = np.hstack((y_train, y_expected))

    # update the locations of X_2D_train
    X_2D_train = deepcopy(x_train)
    return x_proposed

"""
def synTS(X_2D, X_2D_train, Y_2D_train, domain, noise_2D, kappa, nb_iter: int, gamma, workers_loc, start_loc,
          total_time, workers_history: list, dt) -> None:
    
    Propose the next locations based on the synchronous Thompson Sampling Algorithm

    Args:
        :param: prior - the prior matrix
        :param: domain (Domain) - the domain class
        :param: workers_history - tracks the entire history of the locations of the workers
        :param: nb_iter (int) - total number of iterations
        :param: nb_workers - total number of workers
        :param: measurements - contains the measured quantities
        :param: workers_loc - contains the locations of the workers
    
    for i in range(nb_iter):
        windVector = np.array([0, 5, 0])
        aircraftParameters = 1
        workers_loc = proposeLocation(X_2D, X_2D_train, Y_2D_train, windVector, aircraftParameters, total_time,
                                      noise_2D, kappa, gamma, workers_loc, domain)
        print(f"The workers are at {workers_loc}")
        # perform measurements at proposed locations
        workers_history.append(workers_loc)
        for sampleLocation in workers_loc:
            Y_2D_train = np.hstack((Y_2D_train, pollutionAircraftEvent(start_loc, sampleLocation, total_time, dt)))

def proposeLocation(X_2D, X_2D_trainDummy, Y_2D_trainDummy, noise_2D, kappa, gamma, workerPosition, domain : Domain):
    #Proposes a sampling location based on the acquisition function


    idleWorkers = np.zeros(len(workerPosition))
    for i in range(len(workerPosition)):
        sampleLocationList = []
        acquisitionMaxList = []
        res = minimize(nll_fn(X_2D_trainDummy, Y_2D_trainDummy, noise_2D), [1, 1],
                       bounds=((1e-5, None), (1e-5, None)),
                       method='L-BFGS-B')

        mu_s, cov_s = posterior(X_2D, X_2D_trainDummy, Y_2D_trainDummy, *res.x, sigma_y=noise_2D)
        if idleWorkers[i] == 0:
            acquisitionFunc = DUCB(X_2D, X_2D_trainDummy, Y_2D_trainDummy, noise_2D, kappa, gamma, workerPosition[i],domain)
            sampleLocation = np.argmax(acquisitionFunc)
            sampleLocationList.append(sampleLocation)
            acquisitionMax = np.max(acquisitionFunc)
            acquisitionMaxList.append(acquisitionMax)
        else:
            sampleLocationList.append(0)
            acquisitionMaxList.append(0)
        activeWorkerIndex = np.argmax(acquisitionMax)
        print(activeWorkerIndex)
        workerPosition[activeWorkerIndex] = X_2D[sampleLocationList[activeWorkerIndex]]
        X_2D_trainDummy = np.vstack((X_2D_trainDummy, X_2D[sampleLocationList[activeWorkerIndex]]))
        Y_2D_trainDummy = np.append(Y_2D_trainDummy, mu_s[sampleLocationList[activeWorkerIndex]])
        idleWorkers[activeWorkerIndex] = 1
    print(workerPosition)
    return workerPosition
"""
def RMSE(mu, mu_s):
    rmse = 0
    n = len(mu)
    for i in range(n):
        rmse += (mu[i]-mu_s[i])**2
    return math.sqrt(rmse/n)

def randomGenerate(X_3D: np.array, nRand: int):
    """
    Returns random number of samples
    :param X_3D: discretised mesh
    :param nRand: number of random points needed
    :return: return sample point
    """
    return X_3D[np.random.randint(len(X_3D), size=nRand)]

def generateMesh3D(xDomain: np.array, yDomain: np.array, zDomain: np.array, constr1: np.array, constr2: np.array):
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

def main3D():
    minX, maxX, dX = 0.1, 40, 5
    minY, maxY, dY = -20, 20, 5
    minZ, maxZ, dZ = 0, 10, 5
    constr1 = [0, -10, 0]
    constr2 = [20, 10, 10]
    X_3D = generateMesh3D([minX, maxX, dX], [minY, maxY, dY], [minZ, maxZ, dZ], constr1, constr2)
    constr = Rect(Point(constr1[0], constr1[1], constr1[2]), Point(constr2[0], constr2[1], constr2[2]))
    x_3dDomain = Domain(Point(minX, minY, minZ), Point(maxX, maxY, maxZ), constr)
    # Explore/Exploit TradeOff
    # kappa = 15 #exploration/exploitation constant
    kappa = 3000
    #gamma = -0.1  # cost-to-evaluate
    gamma = 0
    initialSamples = 10 # random initial samples
    nIter = 20  # number of points selected by BO algorithm
    noise_3D = 0.01  # Needs a small noise otherwise kernel can become positive semi-definite which leads to minimise() not working

    nWorkers = 2
    workerPosition = randomGenerate(X_3D, nWorkers)
    X_3D_train = workerPosition
    Y_3D = []
    for points in X_3D:
        Y_3D.append(gaussianPlume(points[0], points[1], points[2]))
    Y_3D = np.array(Y_3D)

    #print(workerPosition)

    Y_3D_train = []
    for points in X_3D_train:
        Y_3D_train.append(gaussianPlume(points[0], points[1], points[2]) + noise_3D * np.random.randn())
    Y_3D_train = np.array(Y_3D_train)


    RMSEArray = np.array([])
    for i in range(nIter):
        print("sampling number: ", i)
        workerPosition = proposeLocation(X_3D, X_3D_train, Y_3D_train, noise_3D, kappa, gamma, workerPosition, x_3dDomain)
        print(workerPosition)
        sampleLocation = workerPosition
        print("sample location: ", sampleLocation)

        X_3D_train = np.vstack((X_3D_train, [sampleLocation[0], sampleLocation[1], sampleLocation[2]]))
        Y_3D_train = np.hstack(
            (Y_3D_train, gaussianPlume(sampleLocation[0], sampleLocation[1], sampleLocation[2]) + noise_3D * np.random.randn()))
        res = minimize(nll_fn(X_3D_train, Y_3D_train, noise_3D), [1, 1],
                      bounds=((1e-5, None), (1e-5, None)),
                      method='L-BFGS-B')

        mu_s, cov_s = posterior(X_3D, X_3D_train, Y_3D_train, *res.x, sigma_y=noise_3D)

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
main3D()