import numpy as np

def exponential_connect(weight, n1, n2, autoconnect=True):
    weights = np.random.exponential(1, n1*n2)*weight
    weights = weights.reshape((n1, n2))
    if not(autoconnect):
        weights = weights - np.diag(np.diag(weights))
    return weights