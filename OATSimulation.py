"""
This file contains the implementation of the one at a time method
Which is one of the sensitivity analysis techniques 
"""

import numpy as np
from sqlalchemy import over
from WeightClass import Weights

class OneAtATime(Weights):
    """
    Contains all the methods needed to perform the one at a time method, which assigns different values
    and scores to one element at a time in order to ensure the trade off is robust
    """
    def __init__(self):
        """
        Constructor for the OAT Class
        """
        Weights.__init__(self)
      


    def iterateChanges(self, increments : np.array, vector : np.array, overall_weight : bool):
        """
        auxiliary method used to iterate among the different combinations for the partial weights in a vector
        :param: increments (np.array) - increments added to the initial value to perform OAT
        :param: vector (np.array) - contains the vector to iterate the changes on
        :param: over_weight (bool) - true if the changes are iterated among the overall weights
        Returns options array, containing the number of wins for each design option
        """
        # for the overall weights:
        options = np.zeros(3)
        for idx, weight in enumerate(vector):
            # change the overall weight to analyse its effect
            for increment in increments:
                # rectify the values
                new_score = weight + increment
                # check if the iteration is valid, if not, skip it
                skip = (new_score < 0) or (new_score > 4 and overall_weight) or (new_score > 5 and not overall_weight) 
                if skip:
                    continue

                vector[idx] = new_score
                decision = self.finalScore()
                # increment the corresponding value for the winner of the trade off
                options[decision] += 1
            # after all the iterations, retrieve the all value for the weight
            vector[idx] = weight

        return options

    
    def sensitivityOverallWeights(self, increments : np.array):
        """
        explores the effect of changing the weights assigned to the overall criteria 
        :param: increments (np.array) - increments added to the initial value to perform OAT
        Returns the nb of option for the partial weightsscores
        """
        options = np.zeros(3)
        options += self.iterateChanges(increments, self.overall_weights, overall_weight=True)
        return options


    def sensitivityPartialWeightsScores(self, increments : np.array):
        """
        Explores the effect of changing weights for partial scores
        :param: increments (np.array) - increments added to the initial value to perform OAT
        Returns the nb of option for the partial weightsscores
        """
        options = np.zeros(3)
        for crit in self.scores:
            for weight in crit:
                options += self.iterateChanges(increments, weight, overall_weight=False)
        return options


    def createIncrements(self, lim : float, step : float):
        """
        create the increments array
        :param: lim (float) - the higher limit for the increments 
        :param: step (float) - the step explored by the the OAT
        Returns an np.array containing the increments  (without value 0)
        """
        increments = np.arange(-abs(lim), abs(lim), step)
        # removes the increment = 0 entry in the array
        increments  = np.delete(increments, np.where(increments == 0))
        return increments


    def perform(self, lim : float, step : float, weight : bool):
        """
        Implements the one at a time technique
        :param: lim (float) - the higher limit for the increments 
        :param: step (float) - the step explored by the the OAT
        :param: weight - if true it will assign noise also to the overall weights
        Returns the number of "victories" for each design option, as a numpy array
        """
        # check if varying the weights does not affect the outcome
        # by varying the outcome
        # vector containing the nb of wins for each design option for this particular time step
        options = np.zeros(3)
        
        increments = self.createIncrements(lim, step)
        
        if weight:
            options += self.sensitivityOverallWeights(increments)
        
        options += self.sensitivityPartialWeightsScores(increments)

        # normalise the frequency values
        options = options / np.sum(options) * 100
        
        return options
        

    def iterateLimitsStep(self, lim : float, step : float, weight : bool):
        """
        The main loop for the One-At-A-Time technique, 
        iterates among different steps values
        :param: lim (float) - the higher limit for the increments 
        :param: step (float) - the step explored by the the OAT
        :param: weight - if true it will assign noise also to the overall weights
        And limits values and perform the OAT
        """
        # changes the step for the OAT, in order to check if the scoring system is not too course
        steps = np.arange(0.1, 0.15, 0.05)

        # different limits around which OAT will perform: a smaller value for lims corresponds
        # to a smaller scoring/weighting space around the initially assigned value being explored
        lims = np.arange(0.1, lim, step)
        
        # freqs (np.arrays of tuples) - contains the percentages of "victories" for each design option 
        freqs = np.zeros((len(lims), 3))

        for step in steps:
            print(f"step is {step}")
            for idx, lim in enumerate(lims):
                print(f"lim is {round(lim,2)}")
                freqs[idx] = self.perform(lim, step, weight)
              
                
             
            self.generateVisualisation(lims, freqs)
