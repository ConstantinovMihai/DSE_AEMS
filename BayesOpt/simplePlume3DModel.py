import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
from numpy.linalg import cholesky, det
from scipy.linalg import solve_triangular
from scipy.optimize import minimize
from matplotlib import animation, cm
import math
#new imports for 3D plot
import matplotlib.pyplot as plt


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


def plot_gp_2D(gx, gy, mu, X_train, Y_train, title, i):
    ax = plt.gcf().add_subplot(1, 2, i, projection='3d')
    ax.plot_surface(gx, gy, mu.reshape(gx.shape), cmap=cm.coolwarm, linewidth=0, alpha=0.2, antialiased=False)
    ax.scatter(X_train[:, 0], X_train[:, 1], Y_train, c=Y_train, cmap=cm.coolwarm)
    ax.set_title(title)


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


def simpleFunc(x, y):
    return x + y + x * y ** 2


# Upper Confidence Bound
def UCB(X_2D, X_2D_train, Y_2D_train, noise_2D, kappa, gamma):
    res = minimize(nll_fn(X_2D_train, Y_2D_train, noise_2D), [1, 1],
                   bounds=((1e-5, None), (1e-5, None)),
                   method='L-BFGS-B')
    mu_s, cov_s = posterior(X_2D, X_2D_train, Y_2D_train, *res.x, sigma_y=noise_2D)
    mu_global  = np.mean(mu_s)
    return mu_s/mu_global + kappa * np.sqrt(np.diag(cov_s))

# Distance-based Upper Confidence Bound
def DUCB(X_2D, X_2D_train, Y_2D_train, noise_2D, kappa, gamma):
    res = minimize(nll_fn(X_2D_train, Y_2D_train, noise_2D), [1, 1],
                   bounds=((1e-5, None), (1e-5, None)),
                   method='L-BFGS-B')
    mu_s, cov_s = posterior(X_2D, X_2D_train, Y_2D_train, *res.x, sigma_y=noise_2D)
    mu_global  = np.mean(mu_s)
    lastMeasurement = X_2D_train[-1]
    distance = []
    for point in X_2D:
        distance.append((lastMeasurement[0] - point[0])**2 + (lastMeasurement[1] - point[1])**2 + (lastMeasurement[2] - point[2])**2)
    distance = np.sqrt(np.array(distance))
    print("average variance: ",np.mean(np.sqrt(np.diag(cov_s))))
    return mu_s/mu_global + kappa * np.sqrt(np.diag(cov_s)) + gamma*distance


def proposeLocation(X_2D, X_2D_train, Y_2D_train, noise_2D, kappa, gamma):
    acquisitionFunc = DUCB(X_2D, X_2D_train, Y_2D_train, noise_2D, kappa, gamma)
    sampleLocation = np.argmax(acquisitionFunc)
    return X_2D[sampleLocation]

def parametrisationPlot(X_2D, X_2D_train, Y_2D_train, noise_2D, gx, gy):
    # Plotting Results
    plt.figure(figsize=(14, 7))

    #mu_s, cov_s = posterior(X_2D, X_2D_train, Y_2D_train, sigma_y=noise_2D)
    #plot_gp_2D(gx, gy, mu_s, X_2D_train, Y_2D_train,
    #           f'Before parameter optimization: l={1.00} sigma_f={1.00}', 1)

    res = minimize(nll_fn(X_2D_train, Y_2D_train, noise_2D), [1, 1],
                   bounds=((1e-5, None), (1e-5, None)),
                   method='L-BFGS-B')
    mu_s, cov_s = posterior(X_2D, X_2D_train, Y_2D_train, *res.x, sigma_y=noise_2D)
    plot_gp_2D(gx, gy, mu_s, X_2D_train, Y_2D_train,
               f'After parameter optimization: l={res.x[0]:.2f} sigma_f={res.x[1]:.2f}', 1)
    plt.show()

def RMSE(mu, mu_s):
    rmse = 0
    n = len(mu)
    for i in range(n):
        rmse += (mu[i]-mu_s[i])**2
    return math.sqrt(rmse/n)

def cutPlot(X_2D, X_2D_train, Y_2D_train, noise_2D, gx, gy, func):
    mu = []
    for X in X_2D:
        mu.append(func(X[0], X[1]))
    mu = np.array(mu)
    res = minimize(nll_fn(X_2D_train, Y_2D_train, noise_2D), [1, 1],
                   bounds=((1e-5, None), (1e-5, None)),
                   method='L-BFGS-B')

    mu_s, cov_s = posterior(X_2D, X_2D_train, Y_2D_train, *res.x, sigma_y=noise_2D)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    print("Root mean square error is: ",RMSE(mu, mu_s))
    ax1.imshow(mu.reshape(gx.shape).T, interpolation="nearest",vmin= 0, vmax=1, origin="upper") #Only C between 0 and 1 coloured for interpretation
    # plt.gca().invert_yaxis()
    ax2.imshow(mu_s.reshape(gx.shape).T, interpolation="nearest",vmin=0, vmax=1, origin="upper")
        # plt.gca().invert_yaxis()
    plt.show()

def distanceTravelled(X_2D_train):
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

def randomGenerate(minX, maxX, minY, maxY, minZ, maxZ, restrictions):
    np.random.uniform(minX, maxX), np.random.uniform(minY, maxY), np.random.uniform(minZ, maxZ)


def main3D():
    minX = 0.1
    maxX = 40
    minY = -20
    maxY = 20
    minZ = 0
    maxZ = 30
    dX = 2
    dY = 2
    dZ = 2
    restrictions = [[0, -5, 10], [5, 5, 20]]
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
        sampleLocation = proposeLocation(X_3D, X_3D_train, Y_3D_train, noise_3D, kappa, gamma)
        print("sample location: ", sampleLocation)
        X_3D_train = np.vstack((X_3D_train, [sampleLocation[0], sampleLocation[1], sampleLocation[2]]))
        Y_3D_train = np.hstack(
            (Y_3D_train, gaussianPlume(sampleLocation[0], sampleLocation[1], sampleLocation[2]) + noise_3D * np.random.randn()))


    print("Total Distance Travelled: ", distanceTravelled(X_3D_train[initialSamples:]))
    res = minimize(nll_fn(X_3D_train, Y_3D_train, noise_3D), [1, 1],
                   bounds=((1e-5, None), (1e-5, None)),
                   method='L-BFGS-B')

    mu_s, cov_s = posterior(X_3D, X_3D_train, Y_3D_train, *res.x, sigma_y=noise_3D)

    print(RMSE(Y_3D, mu_s))
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(gx, gy, gz, c=mu_s)
    plt.show()

main3D()