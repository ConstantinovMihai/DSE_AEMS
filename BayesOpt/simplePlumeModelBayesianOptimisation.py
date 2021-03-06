import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
from numpy.linalg import cholesky, det
from scipy.linalg import solve_triangular
from scipy.optimize import minimize
from matplotlib import animation, cm
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable

#np.random.seed(2314)
np.random.seed(2319)
def gaussianPlumeInstant(receiverPosition: np.array, time: float,windVector: np.array):  # x,y coordinates relative to source, z in real coordinates
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
    distance = receiverPosition
    #aircraft plume dispersion
    alpha = 0.02
    aircraftSpeed = np.array([200,0,0])
    wind = aircraftSpeed*math.exp(-alpha*distance[0])
    #create overall dispersion
    wind += windVector
    h = 5
    H = 0
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
    if C<=0: #floating point errors sometimes give negative values
        return 0
    else:
        return C

def sumGaussianPlume(receiverPosition: np.array,   windVector:np.array):
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
    :return: concentration of plume as a function over time
    """
    aircraftEventTime = 1
    timeStep = 0.1
    C = 0
    time = aircraftEventTime
    i = 0
    #sum over time
    while True:
        #H is set to zero
        C += timeStep * gaussianPlumeInstant(receiverPosition,time, windVector)
        time -= timeStep #as you move through the aircraft event the time between the measurement and the emission decreases
        i += 1
        if time < 0:
            return C


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
    print(title)
    ax1 = plt.gcf().add_subplot(1, 2, 1, projection='3d')
    ax1.plot_surface(gx, gy, mu.reshape(gx.shape), cmap=cm.coolwarm, linewidth=0, alpha=0.2, antialiased=False)
    ax1.scatter(X_train[:, 0], X_train[:, 1], Y_train, c=Y_train, cmap=cm.coolwarm)
    ax1.set_xlabel("x coordinate [m]")
    ax1.set_ylabel("y coordinate [m]")
    ax1.set_zlabel("\n Proportional concentration \n magnitude[-]")
    ax1.view_init(0,270)
    ax1 = plt.gcf().add_subplot(1, 2, 2, projection='3d')
    ax1.plot_surface(gx, gy, mu.reshape(gx.shape), cmap=cm.coolwarm, linewidth=0, alpha=0.2, antialiased=False)
    ax1.scatter(X_train[:, 0], X_train[:, 1], Y_train, c=Y_train, cmap=cm.coolwarm)
    ax1.set_xlabel("x coordinate [m]")
    ax1.set_ylabel("y coordinate [m]")
    #ax1.set_zlabel("Proportional concentration magnitude[-]")
    ax1.view_init(elev=90, azim=270)



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
        distance.append((lastMeasurement[0] - point[0])**2 + (lastMeasurement[1] - point[1])**2)
    distance = np.sqrt(np.array(distance))
    print("average variance: ",np.mean(np.sqrt(np.diag(cov_s))))
    return mu_s/mu_global + kappa * np.sqrt(np.diag(cov_s)) + gamma*distance


def proposeLocation(X_2D, X_2D_train, Y_2D_train, noise_2D, kappa, gamma):
    acquisitionFunc = DUCB(X_2D, X_2D_train, Y_2D_train, noise_2D, kappa, gamma)
    sampleLocation = np.argmax(acquisitionFunc)
    return X_2D[sampleLocation]

def parametrisationPlot(X_2D, X_2D_train, Y_2D_train, noise_2D, gx, gy):
    # Plotting Results
    plt.figure(figsize=(14, 14))

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
    mu = mu.reshape(gx.shape).T
    mu_s = mu_s.reshape(gx.shape).T
    ax1.set_title('True Values')
    ax1.set_xlabel('distance [m]')
    ax1.set_ylabel('distance [m]')
    ax1.imshow(mu, interpolation="none",vmin= 0, vmax=1, origin="upper") #Only C between 0 and 1 coloured for interpretation)
    ax2.set_title('Intelligent Sampling')
    im = ax2.imshow(mu_s, interpolation="none",vmin=0, vmax=1, origin="upper")
    #i, j = np.unravel_index(mu_s.argmin(), mu_s.shape)

    ax2.scatter(X_2D_train[:,1]+20, X_2D_train[:,0]-20, color='r')
    ax2.set_xlabel('distance [m]')
    ax2.set_ylabel('distance [m]')
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
        # plt.gca().invert_yaxis()
    plt.show()

def partialCutPlot(X_2D, X_2D_train, Y_2D_train, noise_2D, gx, gy, func):
    res = minimize(nll_fn(X_2D_train, Y_2D_train, noise_2D), [1, 1],
                   bounds=((1e-5, None), (1e-5, None)),
                   method='L-BFGS-B')

    mu_s, cov_s = posterior(X_2D, X_2D_train, Y_2D_train, *res.x, sigma_y=noise_2D)

    mu_s = mu_s.reshape(gx.shape).T
    fig, (ax) = plt.subplots(1, 1)
    ax.set_title('Intelligent Sampling')
    im = ax.imshow(mu_s, interpolation="none",vmin=0, vmax=1, origin="upper")
    #i, j = np.unravel_index(mu_s.argmin(), mu_s.shape)

    ax.scatter((X_2D_train[:,1])+20, (X_2D_train[:,0])-20, color='r')
    ax.set_xlabel('distance [m]')
    ax.set_ylabel('distance [m]')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
        # plt.gca().invert_yaxis()
    plt.show()

def distanceTravelled(X_2D_train):
    dist = 0
    X0 = X_2D_train[0][0]
    X1 = X_2D_train[0][1]
    for X in X_2D_train:
        dist += math.sqrt((X[0] - X0)**2 + (X[1]-X1)**2)
        X0 = X[0]
        X1 = X[1]
    return dist


# Determine Domain Size and Width
minX = 20
maxX = 80
minY = -20
maxY = 20
dX = 1
dY = 1

# Explore/Exploit TradeOff
#kappa = 15 #exploration/exploitation constant
kappa = 15
#gamma = -0.01 #cost-to-evaluate
gamma = 0
initialSamples = 1 #random initial samples
nIter = 7
#number of points selected by BO algorithm
noise_2D = 0.001  # Needs a small noise otherwise kernel can become positive semi-definite which leads to minimise() not working

rx, ry = np.arange(minX, maxX, dX), np.arange(minY, maxY, dY)
gx, gy = np.meshgrid(rx, ry)

# Generate Initial Samples
X_2D = np.c_[gx.ravel(), gy.ravel()]

X_2D_train = np.array([[np.random.uniform(minX, maxX), np.random.uniform(minY, maxY)]])
for i in range(initialSamples):
    X_2D_train = np.vstack((X_2D_train, [np.random.uniform(minX, maxX), np.random.uniform(minY, maxY)]))
# pls cut the initial 10% of the plot, the singularity is hard to reconstruct qnd physically unrealistic
#I think there are more important things to do rn such as the acuqisition then i can do that you want it now?
# not necessarily, but it would be smth nice to have and not very long to implement i think
# i think all that needs to be done is to not let GP sample the top 10% of the map, and plot only the bottom 90% of both graphs
Y_2D_train = []
for points in X_2D_train:
    Y_2D_train.append(gaussianPlume(points[0], points[1]))
Y_2D_train = np.array(Y_2D_train)

# Selection of new sampling locations
for i in range(nIter):
    #print("sampling number: ", i)
    for i in range(3):
        sampleLocation = proposeLocation(X_2D, X_2D_train, Y_2D_train, noise_2D, kappa, gamma)
        #print("sample location: ", sampleLocation)
        X_2D_train = np.vstack((X_2D_train, [sampleLocation[0], sampleLocation[1]]))
        Y_2D_train = np.hstack((Y_2D_train, gaussianPlume(sampleLocation[0], sampleLocation[1])))
cutPlot(X_2D, X_2D_train, Y_2D_train, noise_2D, gx, gy, gaussianPlume)
plt.show()
print("Total Distance Travelled: ", distanceTravelled(X_2D_train))
#parametrisationPlot(X_2D, X_2D_train, Y_2D_train, noise_2D, gx, gy)
