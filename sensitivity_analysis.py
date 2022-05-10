""""
Sensitivity Analysis For The Trade-Off
"""

import numpy as np
from random import random
import re
import scipy.stats as ss
import matplotlib.pyplot as plt
from copy import deepcopy


class Weights:
    def __init__(self) -> None:
        self.overall_weights = []
        # stores the matrices for the trade off criteria
        self.scores = []
        self.winner = -1
        self.partial_scores = []


    # each sub-system is a matrix of weights, and their weighted sum represents the final score
    def readData(self, filename):
        """
        Reads the data from filename and populates the matrices attributes of Weights
        """
        file = open(filename, 'r')
        lines  = file.readlines()
        for line in lines:
            if line.startswith("\"Overall\""):
                self.overall_weights = np.array([float(s) for s in re.findall(r'-?\d+\.?\d*', line)])
            else:
                line = line.split(":")[1].split("],")
                #self.weights_crit = [float(s) for s in re.findall(r'-?\d+\.?\d*', line[0])]
                crits = []       
                for el in line:
                    crit = [float(s) for s in re.findall(r'-?\d+\.?\d*', el)]
                    crits.append(np.array(crit))
                self.scores.append(crits)
        # convert the self.scores from a list to an array
        self.scores = self.scores
           

    def partialScores(self):
        """
        Computes the partial score for each sub-trade off
        """
        # iterate among each criteria and 
        # perform the weighted average sum
        self.partial_scores = []
        for crit in self.scores:
            print(f"crit is {crit} ")
            print(f"np.average {np.average(crit[1:], axis=1, weights=crit[0])}")
            self.partial_scores.append(np.average(crit[1:], axis=1, weights=crit[0]))


    def finalScore(self):
        """
        Computes the final score and decides a winner
        """
        self.partialScores()
        self.winner = np.argmax(np.average(self.partial_scores, axis=0, weights=self.overall_weights))
        print(self.winner)



class OneAtATime(Weights):
    """
    Contains all the methods needed to perform the one at a time method, which assigns different values
    and scores to one element at a time in order to ensure the trade off is robust
    """
    def __init__(self):
        Weights.__init__(self)


    def iterateChanges(self, increments, vector):
        """
        auxiliary method used to iterate among the different combinations for the partial weights in a vector
        """
         # for the overall weights:
        for idx, weight in enumerate(vector):
            # change the overall weight to analyse its effect
            for increment in increments:
                vector[idx] = max(weight + increment, 0)
                self.finalScore()
            # after all the iterations, retrieve the all value for the weight
            vector[idx] = weight

    
    def sensitivityOverallWeights(self, increments):
        """
        explores the effect of changing the weights assigned to the overall criteria 
        """
        self.iterateChanges(increments, self.overall_weights)


    def sensitivityPartialWeightsScores(self, increments):
        """
        Explores the effect of changing weights for partial scores
        """
        for crit in self.scores:
            for weight in crit:
                self.iterateChanges(increments, weight)
    

    def createIncrements(self):
        """
        create the increments array
        """
        increments = np.arange(-2, 2, 0.5)
        increments  = np.delete(increments, np.where(increments == 0))
        return increments


    def perform(self):
        """
        Implements the one at a time technique
        """
        # check if varying the weights does not affect the outcome
        # by varying the outcome
        increments = self.createIncrements()
        self.sensitivityOverallWeights(increments)
        self.sensitivityPartialWeightsScores(increments)
       


class MonteCarlo(Weights):
    """
    Class containing the methods to perform the Monte Carlo sensitivity analysis
    It adds some noise in the matrix of weights and checks which design is predicted as optimum
    """
    def __init__(self):
        Weights.__init__(self)
    

    def generateRandomNumbers(self, nb_instances):
        """ 
        Generate random numbers by drawing numbers from multinomial distributions
        whose probabilities are drawned from normal distribution
        Returns nb_instances numbers drawn from this distribution
        """          
        x = np.arange(-10, 11)
        xU, xL = x + 0.5, x - 0.5 
        prob = ss.norm.cdf(xU, scale = 3) - ss.norm.cdf(xL, scale = 3)
        prob = prob / prob.sum() # normalize the probabilities so their sum is 1
        nums = np.random.choice(x, size = nb_instances, p = prob) / 4
        return nums
    

    def generateNoisyWeights(self):
        """
        Generates a matrix of weights containg noisy data, 
        used to verify the robustness of the trade-off
        """
        # deep copy the initial matrix, then add noise
        d = deepcopy(self)

        noise = self.generateRandomNumbers(len(d.overall_weights))
        d.overall_weights = d.overall_weights + noise
       #  d.overall_weights[d.overall_weights < 0] = 0

       
        noise = self.generateRandomNumbers(len(d.scores[0]) * len(d.scores))
        # print(f" (len(d.scores), len(d.scores[0])), { len(d.scores)}   {len(d.scores[0])}")
        # print(f"d.scores is \n{d.scores} ")
        # d.scores = d.scores + np.reshape(noise, (len(d.scores), len(d.scores[0])))
        print(f"modified d.scores are \n{d.scores}")
        #d.scores[d.scores < 0] = 0.01
        return d


    def iterateWeights(self, nb_experiments):
        """
        Main method of the class, generates nb_experiments noisy matrices and checks 
        whether or not they predict the same final design option as the "clean" matrix
        """
        for i in range(nb_experiments):
            # generate noisy matrix check the final score it predicts
            d = self.generateNoisyWeights()
            d.finalScore()

p = MonteCarlo()
p.readData(r"C:\Users\mihai\OneDrive\Desktop\DSE\trade_off.txt")
p.iterateWeights(100)


# p = OneAtATime()
# p.readData(r"C:\Users\mihai\OneDrive\Desktop\DSE\trade_off.txt")
# p.oneAtATime()


# TODO: implement Monte Carlo Technique
# TODO: implement Visualisation for the Monte Carlo Techniques