"""
Implements the noise model and simulates the pressure level for a receiver
Measuring noise from a point source passing by
"""

from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import random
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF
import gp_util
import NoiseBO
from numpy.linalg import inv
from numpy.linalg import cholesky, det
from scipy.linalg import solve_triangular
from scipy.optimize import minimize


from domains import Point
from domains import Domain
from domains import Rect

from sklearn.base import clone
from skopt import gp_minimize
from skopt.learning import GaussianProcessRegressor
from skopt.learning.gaussian_process.kernels import ConstantKernel, Matern

class Vector():
    def __init__(self, x, y, z) -> None:
        self.x = x
        self.y = y
        self.z = z

    def norm(self) -> float:
        """
        Computes the Euclidian Norm of a vector
        """
        return np.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)

    def subtract(self, other):
        """
        Implements vectorial substraction operation

        Args:
            other (Vector) - The vector to be subtracted
        Returns a Vector containing the result of the substraction
        """
        res = deepcopy(self)
        res.x -= other.x
        res.y -= other.y
        res.z -= other.z
        return res
    
    def print(self) -> None:
        """
        prints the components of the vector
        """
        print(f"x,y,z: {self.x} , {self.y}, {self.z}")

    def setX(self, x : float) -> None:
        """Sets the x coordinate of the vector

        Args:
            x (float): the new value of the x coordinate
        """
        self.x = x
    
    def setY(self, y : float) -> None:
        """Sets the y coordinate of the vector

        Args:
            y (float): the new value of the y coordinate
        """
        self.y = y

    def setZ(self, z : float) -> None:
        """Sets the z coordinate of the vector

        Args:
            z (float): the new value of the z coordinate
        """
        self.z = z


def noiseModel(pwl : float, r : float, alpha : float) -> float:
    """
    Computes the SPL (eq 3.63) 

    Args:
        pwl (float): power of the source
        r (float): distance from the source to the receiver
        alpha (float): atmospheric attenuation factor
    """
    return pwl - 10.8 - 20 * np.log10(r) - r * alpha


def simulatePassOver(d : float, vel : float, pwl : float, start :float, alpha : float, tot_time :float, dt : float = 1, plot = False):
    """
    Simulates the pass over of an aircraft on the runway

    Args:
        d (float): smallest distance from the a/c to the receiver (m)
        vel (float) : velocity of the a/c (m/s)
        start (float) : the starting distance on x coord from the source (start of sim) (m)
        alpha (float): atmospheric attenuation factor (dB/m)
        pwl (float) : power of the source (W)
        tot_time (float) : total simulation time (s)
        dt (float) : time step of the simulation (s)
        plot (bool) : if set true then noise vs time graph will be displayed

    Returns a np.array containing the noise corresponding to each step
    """
    # the origin of the SoC is on the path of the a/c
    
    # coordinates of a/c
    ac = Vector(start, 0, 0)
    # coordinates of receiver
    rec = Vector(0, d, 0.10)

    n_steps = int(tot_time / dt)
    noise = np.zeros(n_steps)
    times = np.arange(0, tot_time, dt)
    x_coords = np.arange(start, tot_time, dt * vel)

    for idx, x in enumerate(x_coords):
        ac.setX(x)
        diff = ac.subtract(rec)
        r = diff.norm()
        noise[idx] = noiseModel(pwl, r, alpha) 

    if plot:
        plt.plot(times, noise)
        plt.ylabel("Noise (dB)")
        plt.xlabel("Distance (m)")
        plt.show()

    return noise


def rationalQuadraticKernel(X1, X2, l=1.0, sigma_f=1.0, alpha = 2):
    """
    Isotropic squared exponential kernel.
    
    Args:
        X1: Array of m points (m x d).
        X2: Array of n points (n x d).

    Returns:
        (m x n) matrix.
    """
    sqdist = np.sum(X1**2, 1).reshape(-1, 1) + np.sum(X2**2, 1) - 2 * np.dot(X1, X2.T)
    return sigma_f**2 * np.power(1 + 0.5 / (l**2 * alpha) * sqdist, -alpha)


def posteriorQuadraticKernel(X_s, X_train, Y_train, l=1.0, sigma_f=1.0, alpha = 2, sigma_y=1e-8):
    """
    Computes the suffifient statistics of the posterior distribution 
    from m training data X_train and Y_train and n new inputs X_s.
    
    Args:
        X_s: New input locations (n x d).
        X_train: Training locations (m x d).
        Y_train: Training targets (m x 1).
        l: Kernel length parameter.
        alpha: kernel parameter
        sigma_f: Kernel vertical variation parameter.
        sigma_y: Noise parameter.
    
    Returns:
        Posterior mean vector (n x d) and covariance matrix (n x n).
    """
    K = rationalQuadraticKernel(X_train, X_train, l, sigma_f, alpha) + sigma_y**2 * np.eye(len(X_train))
    K_s = rationalQuadraticKernel(X_train, X_s, l, sigma_f, alpha)
    K_ss = rationalQuadraticKernel(X_s, X_s, l, sigma_f, alpha) + 1e-8 * np.eye(len(X_s))
    K_inv = inv(K)
    
    # Equation (7)
    mu_s = K_s.T.dot(K_inv).dot(Y_train)

    # Equation (8)
    cov_s = K_ss - K_s.T.dot(K_inv).dot(K_s)
    
    return mu_s, cov_s


def nll_fnQuadratic(X_train, Y_train, noise, naive=True):
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
        K = rationalQuadraticKernel(X_train, X_train, l=theta[0], sigma_f=theta[1]) + \
            noise**2 * np.eye(len(X_train))
        return 0.5 * np.log(det(K)) + \
               0.5 * Y_train.dot(inv(K).dot(Y_train)) + \
               0.5 * len(X_train) * np.log(2*np.pi)
        
    def nll_stable(theta):
        # Numerically more stable implementation of Eq. (11) as described
        # in http://www.gaussianprocess.org/gpml/chapters/RW2.pdf, Section
        # 2.2, Algorithm 2.1.
        
        K = rationalQuadraticKernel(X_train, X_train, l=theta[0], sigma_f=theta[1]) + \
            noise**2 * np.eye(len(X_train))
        L = cholesky(K)
        
        S1 = solve_triangular(L, Y_train, lower=True)
        S2 = solve_triangular(L.T, S1, lower=False)
        
        return np.sum(np.log(np.diagonal(L))) + \
               0.5 * Y_train.dot(S2) + \
               0.5 * len(X_train) * np.log(2*np.pi)

    if naive:
        return nll_naive
    else:
        return nll_stable

def noiseSim():
    minX = -600
    maxX = -600+1501
    minY = 0
    maxY = 20
    minZ = 0
    maxZ = 20
    dX = 1
    dY = 20
    dZ = 20
    constr = Rect(Point(0,-5,10), Point(0,-5,10))
    domain = Domain(Point(minX, minY, minZ), Point(maxX, maxY, maxZ), constr)
    rx, ry, rz = np.arange(minX, maxX, dX), np.arange(minY, maxY, dY), np.arange(minZ, maxZ, dZ)
    gx, gy, gz = np.meshgrid(rx, ry, rz)
    X_3D = np.c_[gx.ravel(), gy.ravel(), gz.ravel()]

    noise_data = simulatePassOver(100, 50, 150, -600, 0.0085, 1500)
    
    # Explore/Exploit TradeOff
    # kappa = 15 #exploration/exploitation constant
    kappa = 3000
    #gamma = -0.1  # cost-to-evaluate
    gamma = 0
    initialSamples = 20 # random initial samples
    nIter = 20  # number of points selected by BO algorithm
    noise_3D = 0.01  # Needs a small noise otherwise kernel can become positive semi-definite which leads to minimise() not working

    X_3D_train = np.array([[np.random.uniform(minX, maxX), np.random.uniform(minY, maxY), np.random.uniform(minZ, maxZ)]])
    for i in range(initialSamples):
        X_3D_train = np.vstack((X_3D_train, [np.random.uniform(minX, maxX), np.random.uniform(minY, maxY), np.random.uniform(minZ, maxZ)]))

    Y_3D_train = []
    for points in X_3D_train:
        Y_3D_train.append(noise_data[int(points[0])])
    Y_3D_train = np.array(Y_3D_train)

    Y_3D = []
    for points in X_3D:
        Y_3D.append(noise_data[int(points[0])] + noise_3D * np.random.randn())
    Y_3D = np.array(Y_3D)

    for i in range(nIter):
        print("sampling number: ", i)
        sampleLocation = NoiseBO.proposeLocation(X_3D, X_3D_train, Y_3D_train, noise_3D, kappa, gamma, domain)
        print("sample location: ", sampleLocation)
        X_3D_train = np.vstack((X_3D_train, [sampleLocation[0], sampleLocation[1], sampleLocation[2]]))
        Y_3D_train = np.hstack(
            (Y_3D_train, noise_data[int(sampleLocation[0])] + noise_3D * np.random.randn()))


if __name__ == "__main__":
    noiseSim()