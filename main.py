""""
Sensitivity Analysis For The Trade-Off
Uses the Weights class as a data structure to
store the data
And the OAT and MonteCarlo classes to perform the simulations
"""

import string
import os

from OATSimulation import OneAtATime
from MonteCarloSimulation import MonteCarlo


def runOATSimulation(filename : string, lim : float = 2.1, step : float = 0.1, weight = True):
    """
    Performs the OAT method and produces the frequency graph
    :param: filename (string) - the location of the file with data (scores and weights)
    :param: weight - if true it will assign noise also to the overall weights
    """
    oat = OneAtATime(filename)
    oat.finalScore()
    oat.iterateLimitsStep(lim, step, weight)


def runMonteCarloSimulation(filename : string, weight : bool, nb_experiments : int = 10000, plot : bool = True):
    """
    Performs the Monte Carlo method and produces the frequency graph
    :param: filename (string) - the location of the file with data (scores and weights)
    :param: weight (bool) - if true it will assign noise only to the overall weights
    :param: nb_experiments (float) - the number of experiments to be performed at each monte carlo simulation
    :param: plot (bool) - if set true the graph will be shown
    """
    p = MonteCarlo()
    p.readData(filename)
    p.iterateDeviations(nb_experiments, plot, weight)


if __name__ == "__main__":
    currentDirectory = str(os.getcwd()) + "\\tradeOffVerification.txt"
    step = 0.1
    runOATSimulation(currentDirectory, 2 + step, step, weight=True)
    #runMonteCarloSimulation(currentDirectory, weight=False, nb_experiments=2000)
    print("done")
    
    
