"""
This file contains the implementation of the Monte Carlo class
Which performs the simulations 
"""

import numpy as np
import scipy.stats as ss
from copy import deepcopy
from WeightClass import Weights

class MonteCarlo(Weights):
    """
    Class containing the methods to perform the Monte Carlo sensitivity analysis
    It adds some noise in the matrix of weights and checks which design is predicted as optimum
    """
    def __init__(self):
        """
        Constructor for the MonteCarlo Class
        """
        Weights.__init__(self)


    def generateRandomNumbers(self, nb_instances, deviation=4):
        """ 
        Generate random numbers by drawing numbers from multinomial distributions
        whose probabilities are drawned from normal distribution
        :param: deviation - sets the standard deviation - the lower the deviation value, the higher the variance
        :param: nb_instances - the number of random numbers to be generated
        Returns a numpy array containg nb_instances numbers drawn from this distribution
        """          
        # defines a Truncated Normal Distribution
        # scale = 3 as with this standard distribution the sum is roughly equal to 1
        # / deviation adjust the standard deviation
        x = np.arange(-10, 11)
        xU, xL = x + 0.5, x - 0.5 
        prob = ss.norm.cdf(xU, scale = 3) - ss.norm.cdf(xL, scale = 3)
        prob = prob / prob.sum() # normalize the probabilities so their sum is 1
        nums = np.random.choice(x, size = nb_instances, p = prob) / deviation 
        return nums
    

    def rectifyValues(self, values):
        """
        Utility function which rectifies the negative values
        In values matrix to 0
        """
        values[values < 0] = 0
      
        return values


    def generateNoisyWeights(self, deviation=4, weight=True):
        """
        Generates a matrix of weights containg noisy data, 
        used to verify the robustness of the trade-off
        :param: deviation - sets the standard deviation - the lower the deviation value, the higher the variance
        :param: weight - if true it will assign noise also to the overall weights
        Returns a matrix containing noises
        """
        # deep copy the initial matrix, then add noise
        d = deepcopy(self)

        if weight:
            noise = self.generateRandomNumbers(len(d.overall_weights), deviation)
            d.overall_weights = d.overall_weights + noise
        

        # unfortunately we could not find a way to reshape 
        # either reshape the list into np.array without errors/deprecation warnings
        # list comprehension for noise would be possible though, it is something 
        # TODO: get rid of the matrix element by element iteration
        for score in d.scores:
            for row in score:
                noise = self.generateRandomNumbers(len(row), deviation)
                row += noise
                row = self.rectifyValues(row)

        self.rectifyValues(d.overall_weights)
        return d


    def iterateWeights(self, deviation=4, nb_experiments=100, weight=True):
        """
        Main method of the class, generates nb_experiments noisy matrices and checks 
        whether or not they predict the same final design option as the "clean" matrix
        :param: nb_experiments - the number of Monte Carlo simulations to be performed
        :param: deviation - sets the standard deviation - the lower the deviation value, the higher the variance
        :param: weight - if true it will assign noise also to the overall weights
        Returns the vector of frequencies
        """
        options = np.zeros(3)
        for i in range(nb_experiments):
            # generate noisy matrix check the final score it predicts
            d = self.generateNoisyWeights(deviation, weight)
            option = d.finalScore()
            options[option] += 1

        # self.generateVisualisation(options, nb_experiments)
        print(options / sum(options))
        return options / sum(options)
        

    def iterateDeviations(self, nb_experiments=100, plot=False, weight=True):
        """
        Iterates among different values for the deviation in order to explore the sensibility of the model
        :param: plot - if True, it will generate the histogram
        :param: weight - if true it will assign noise also to the overall weights
        To the chosen probability distribution
        """  
        deviations = np.arange(1, 7, 0.2)
        freqs = np.zeros((len(deviations), 3))

        for idx, deviation in enumerate(deviations):
            freq = self.iterateWeights(deviation, nb_experiments, weight)
            freqs[idx] = freq
        # 3 / deviation is the the standard deviation for
        # the Truncated Normal Distribution 
        self.generateVisualisation(3 / deviations, freqs)

    