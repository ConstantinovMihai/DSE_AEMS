""""
Sensitivity Analysis For The Trade-Off
"""

from email.headerregistry import HeaderRegistry
import numpy as np
import re
import scipy.stats as ss
import matplotlib.pyplot as plt
from copy import deepcopy
import os
import pandas as pd

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
                crits = []       
                for el in line:
                    crit = [float(s) for s in re.findall(r'-?\d+\.?\d*', el)]
                    crits.append(np.array(crit))
                self.scores.append(crits)
        # convert the self.scores from a list to an array
        # self.scores = np.array(self.scores)
           

    def partialScores(self):
        """
        Computes the partial score for each sub-trade off
        """
        # iterate among each criteria and 
        # perform the weighted average sum
        self.partial_scores = []
        for crit in self.scores:
            weights = crit[0]
            # dirty trick to avoid the situation of the sum of weights being 0
            # overall it should not affect the final result
            # it just assumes all the weights are equal
            # indeed a bit biased towards a particular output, but in the long run it 
            # should not affect much the final result
            if not sum(crit[0]):
                weights = np.ones(len(crit[0])) 
            self.partial_scores.append(np.average(crit[1:], axis=1, weights=weights))


    def finalScore(self, print_it=False):
        """
        Computes the final score and decides a winner
        :param:  print (bool) - if set true it prints the index of the winner
        """
        self.partialScores()
        self.winner = np.argmax(np.average(self.partial_scores, axis=0, weights=self.overall_weights))

        if print_it:
            print_it(self.winner)

        return self.winner

        
    def visualiseResults(self, freq, nb_experiments=1, plot=False, x=[ "Balloon","VTOL", "N-copter"], print_it=False, name=""):
        """
        Functionality to be able to visualise the results of the sensitivity analysis 
        :param: nb_experiments - the number of Monte Carlo simulations to be performed
        :param: Freq - a vector with the frequencies of various options as simulated by the Monte Carlo method
        :param: plot - if True, it will generate the histogram
        :param: x - the different options whose frequencies has to be visualised
        :param: name - name under which to save the histograms
        Generates a histogram with the options' frequencies values
        """
        freq = freq / nb_experiments * 100
        if print_it:
            print(f"Options1 - {round(freq[0], 3)}")
            print(f"Options2 - {round(freq[1], 3)}")
            print(f"Options3 - {round(freq[2], 3)}") 
            print(f"Total iterations - {nb_experiments}")
        # plot the histogram
        if plot:
            df = pd.DataFrame({"Options":x,'Freq. (%)':freq})
            df.set_index("Options")['Freq. (%)'].plot.bar(rot=0)
            plt.savefig('Histogram' + name + '.png')
            plt.show()



class OneAtATime(Weights):
    """
    Contains all the methods needed to perform the one at a time method, which assigns different values
    and scores to one element at a time in order to ensure the trade off is robust
    """
    def __init__(self):
        """
        Constructor for the MonteCarlo Class
        """
        Weights.__init__(self)


    def iterateChanges(self, increments : np.array, vector : np.array):
        """
        auxiliary method used to iterate among the different combinations for the partial weights in a vector
        :param: increments (np.array) - increments added to the initial value to perform OAT
        Returns options array, containing the number of wins for each design option
        """
         # for the overall weights:
        options = np.zeros(3)
        for idx, weight in enumerate(vector):
            # change the overall weight to analyse its effect
            for increment in increments:
                vector[idx] = max(weight + increment, 0)
                decision = self.finalScore()
                # increment the corresponding value for the winner of the trade off
                options[decision] += 1
            # after all the iterations, retrieve the all value for the weight
            vector[idx] = weight

        return options

    
    def sensitivityOverallWeights(self, increments):
        """
        explores the effect of changing the weights assigned to the overall criteria 
        :param: increments (np.array) - increments added to the initial value to perform OAT
        """
        options = np.zeros(3)
        options += self.iterateChanges(increments, self.overall_weights)
        return options


    def sensitivityPartialWeightsScores(self, increments):
        """
        Explores the effect of changing weights for partial scores
        :param: increments (np.array) - increments added to the initial value to perform OAT
        """
        options = np.zeros(3)
        for crit in self.scores:
            for weight in crit:
                options += self.iterateChanges(increments, weight)
        return options


    def createIncrements(self, low_lim=-2, high_lim=2, step=0.25):
        """
        create the increments array
        :param: low_lim - the lower limit for the increments 
        :param: high_lim - the higher limit for the increments 
        """
        increments = np.arange(low_lim, high_lim, step)
        # remove the increment = 0 entry in the array
        increments  = np.delete(increments, np.where(increments == 0))
        return increments


    def perform(self, low_lim=-2, high_lim=2, step=0.25):
        """
        Implements the one at a time technique
        :param: low_lim - the lower limit for the increments 
        :param: high_lim - the higher limit for the increments 
        :param: step - the step explored by the the OAT
        """
        # check if varying the weights does not affect the outcome
        # by varying the outcome
        options = np.zeros(3)
        increments = self.createIncrements(low_lim, high_lim, step)
        
        options += self.sensitivityOverallWeights(increments)
        options += self.sensitivityPartialWeightsScores(increments)
        print(f"options {options}")

        self.visualiseResults(options)


    def iterateLimitsStep(self):
        steps = np.arange(0.1, 1, 0.05)
        lims = np.arange(0.5, 3, 0.1)
        for step in steps:
            print(f"step is {step}")
            for lim in lims:
                print(f"lim is {round(lim,1)}")
                self.perform(-lim, lim, step)
       


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
        values[values < 0] = 0.01
        return values


    def generateNoisyWeights(self, deviation=4):
        """
        Generates a matrix of weights containg noisy data, 
        used to verify the robustness of the trade-off
        :param: deviation - sets the standard deviation - the lower the deviation value, the higher the variance
        Returns a matrix containing noises
        """
        # deep copy the initial matrix, then add noise
        d = deepcopy(self)

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


    def iterateWeights(self, deviation=4, nb_experiments=100):
        """
        Main method of the class, generates nb_experiments noisy matrices and checks 
        whether or not they predict the same final design option as the "clean" matrix
        :param: nb_experiments - the number of Monte Carlo simulations to be performed
        :param: deviation - sets the standard deviation - the lower the deviation value, the higher the variance
        Returns the vector of frequencies
        """
        options = np.zeros(3)
        for i in range(nb_experiments):
            # generate noisy matrix check the final score it predicts
            d = self.generateNoisyWeights(deviation)
            option = d.finalScore()
            options[option] += 1

        self.visualiseResults(options, nb_experiments)
        print(options / sum(options))
        return options / sum(options)

    
    def generateVisualisation(self, x, y):
        """
        Generates visualisation of the various percentages of the trade off wins
        Based on the standard deviations used when generating the simulations
        :param: x (list)  - elements on the x axis
        :param: y (list of tuples) - elements on the y axis 
        """
        y1 = y[:,0]
        y2 = y[:,1]
        y3 = y[:,2]
        plt.plot(x, y1, color = 'r', label = "Balloon")
        plt.plot(x, y2, color = 'g', label = "VTOL")
        plt.plot(x, y3, color = 'b', label = "N-copter")
        plt.xlabel("Deviation")
        plt.ylabel("Frequency")
        plt.legend()
        plt.grid(True)
        plt.show()
        

    def iterateDeviations(self, nb_experiments=100, plot=False):
        """
        Iterates among different values for the deviation in order to explore the sensibility of the model
        :param: plot - if True, it will generate the histogram
        To the chosen probability distribution
        """  
        deviations = np.arange(0.5,5,0.25)
        freqs = np.zeros((len(deviations), 3))

        for idx, deviation in enumerate(deviations):
            freq = self.iterateWeights(deviation, nb_experiments)
            freqs[idx] = freq
        
        self.generateVisualisation(3 / deviations, freqs)

    
if __name__ == "__main__":
    p = MonteCarlo()
    currentDirectory = str(os.getcwd()) + "\\trade_off.txt"
    p.readData(currentDirectory)
    p.iterateDeviations(nb_experiments = 10000, plot=True)
   
