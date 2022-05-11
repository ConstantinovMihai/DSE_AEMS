import numpy as np

crit = [np.array([5.1, 0, 4.3, 2.4]), np.array([2.4, 4, 4, 2]), np.array([3, 2.2, 2, 2]), np.array([3, 2.2, 2, 4])]

crit = np.array(crit)

print(f"np.average {np.average(crit[1:], axis=1, weights=crit[0])}")