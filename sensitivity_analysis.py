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


def runOATSimulation(filename : string, lim : float = 2.1, step : float = 0.1):
    """
    Performs the OAT method and produces the frequency graph
    :param: filename (string) - the location of the file with data (scores and weights)
    """
    oat = OneAtATime()
    oat.readData(filename)
    oat.iterateLimitsStep(lim, step)


def runMonteCarloSimulation(filename : string, nb_experiments : int = 10000, plot : bool = True):
    """
    Performs the Monte Carlo method and produces the frequency graph
    :param: filename (string) - the location of the file with data (scores and weights)
    :param: nb_experiments (float) - the number of experiments to be performed at each monte carlo simulation
    :param: plot (bool) - if set true the graph will be shown
    """
    p = MonteCarlo()
    p.readData(filename)
    p.iterateDeviations(nb_experiments, plot)


if __name__ == "__main__":
    currentDirectory = str(os.getcwd()) + "\\trade_off.txt"

    runOATSimulation(currentDirectory, 2.1, 0.1)

    # runMonteCarloSimulation(currentDirectory, nb_experiments=5000)

    
    print("done")
    
    
