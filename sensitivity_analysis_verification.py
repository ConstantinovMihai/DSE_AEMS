"""
Performs sensitivity analysis  verification 
"""

import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from random import random

# first, try some dummy value to perform sensitivity analysis for a very simple problem
# such as a rock thrown at an angle

# inputs: angle thetha, initial velocity V0, gravity
# quantities of interest: max height H, range R
# H = V^2 Sin^2(θ) / 2g

red_patch = mpatches.Patch(color='red', label='The red data')
blue_patch = mpatches.Patch(color='blue', label='The blue data')
orange_patch = mpatches.Patch(color ='orange', label = 'The orange data')
colors = ['red', 'blue', 'orange']
plt.legend(handles=[red_patch, blue_patch])

class Params:
    def __init__(self, vel, angle, g0 = 10) -> None:
        self.v = vel
        self.theta = angle
        self.g0 = g0


def rangeFormula(inputs : list):
    """
    simple dummy function to compute the range of a 
    ball thrown at an angle v (m/s), angle theta (deg) under 
    a gravity field g (m/s2)
    """
    v = inputs[0]
    theta = inputs[1]
    g = inputs[2]
    return v * v * np.sin(theta) * np.sin(theta) / (2 * g)


def oatAnalysis(inputs : list):
    """
    Performs the One-at-a-time sensitivity analysis technique
    :param: inputs : contains the list of in
    returns: a list of lists of ranges
    """
    for el, param in enumerate(inputs):
        x_val = []
        values = []
        aux = copy.deepcopy(param)
        inputs_aux = copy.deepcopy(inputs)
        low_val = 100
        high_val = 100 
        for i in range(-low_val,high_val):
            # vary the input from -10% initial values to 10% initial value
            inputs_aux[el] = aux + aux * (i/(5*high_val))
            # store the computed data of the new range
            x_val.append(1 + i/(5*high_val) )
            values.append(rangeFormula(inputs_aux))
        
        # visualise the spreadness of the values
        plt.legend()
        plt.xlabel("ratio of variation of one parameter")
        plt.ylabel("range")
        plt.scatter(x_val, values)
    plt.show()


def monteCarloSensAnalysis(inputs : list, nb_exp : int):
    """
    Performs the One-at-a-time sensitivity analysis technique
    :param: inputs : contains the list of in
    :param: nb_exp : number of experiments
    returns: a list of lists of ranges
    """
    x_val = []
    ranges = []
    for i in range(nb_exp):
        val = random()
        x_val.append(-0.1 + val / 5) 

if __name__ == "__main__":
    inputs = [50, 20/180 * np.pi, 20]
    oatAnalysis(inputs)

