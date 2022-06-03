"""
File containing the implementation of the Weights class
Which stores the reading data, score calculation and winner determination methods
"""

import numpy as np
import re
import matplotlib.pyplot as plt


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
        
        # print(f"scores are {np.average(self.partial_scores, axis=0, weights=self.overall_weights)}")
        self.winner = np.argmax(np.average(self.partial_scores, axis=0, weights=self.overall_weights))
      
        if print_it:
            print(f'winner is {self.winner} ')

        return self.winner


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
        plt.ylabel("Frequency (%)")
        plt.legend()
        plt.grid(True)
        plt.show()
