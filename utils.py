import numpy as np
from neuron import h


def exponential_connect(weight, n1, n2, autoconnect=True):
    """
    :param weight: desired average of exponential distribution
    :param n1: number of neurons of network 1
    :param n2: number of neurons of network 2
    :param autoconnect: False if diagonal weights should be zero, True otherwise
    :return: weight matrix
    """
    weights = np.random.exponential(1, n1*n2)*weight
    weights = weights.reshape((n1, n2))
    if not autoconnect:
        weights = weights - np.diag(np.diag(weights))
    return weights


def constant_connect(weight, n1, n2, autoconnect=True):
    """
    :param weight: desired connection weight
    :param n1: number of neurons of network 1
    :param n2: number of neurons of network 2
    :param autoconnect: False if diagonal weights should be zero, True otherwise
    :return: weight matrix
    """
    weights = weight*np.ones(n1*n2)
    weights = weights.reshape((n1, n2))
    if not autoconnect:
        weights = weights - np.diag(np.diag(weights))
    return weights


def e_net_connect(net1, net2, threshold, delay, weights):
    net1_net2_sin = list()
    for net1_neuron_i, net1_neuron in enumerate(net1):
        net1_neuron.soma.push()
        for net2_neuron_i, net2_neuron in enumerate(net2):
            net1_net2_sin.append(h.NetCon(net1_neuron.soma(0.5)._ref_v, net2_neuron.synE, threshold, delay,
                                          weights[net1_neuron_i, net2_neuron_i]))
        h.pop_section()
    return net1_net2_sin


def e_ct_net_connect(net1, net2, threshold, delay, weights):
    net1_net2_sin = list()
    for net1_neuron_i, net1_neuron in enumerate(net1):
        net1_neuron.soma.push()
        for net2_neuron_i, net2_neuron in enumerate(net2):
            net1_net2_sin.append(h.NetCon(net1_neuron.soma(0.5)._ref_v, net2_neuron.synE_CT, threshold, delay,
                                          weights[net1_neuron_i, net2_neuron_i]))
        h.pop_section()
    return net1_net2_sin


def i_net_connect(net1, net2, threshold, delay, weights):
    net1_net2_sin = list()
    for net1_neuron_i, net1_neuron in enumerate(net1):
        net1_neuron.soma.push()
        for net2_neuron_i, net2_neuron in enumerate(net2):
            net1_net2_sin.append(h.NetCon(net1_neuron.soma(0.5)._ref_v, net2_neuron.synI, threshold, delay,
                                          weights[net1_neuron_i, net2_neuron_i]))
        h.pop_section()
    return net1_net2_sin