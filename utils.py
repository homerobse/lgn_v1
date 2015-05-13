import numpy as np
from neuron import h


def exponential_connect(weight, n1, n2, selfconnect=True):
    """
    :param weight: desired average of exponential distribution
    :param n1: number of neurons of network 1
    :param n2: number of neurons of network 2
    :param selfconnect: False if diagonal weights should be zero, True otherwise
    :return: exponentially distributed weight matrix (n1 x n2)
    """
    weights = np.random.exponential(1, n1*n2)*weight
    weights = weights.reshape((n1, n2))
    if not selfconnect:
        weights = weights - np.diag(np.diag(weights))
    return weights


def constant_connect(weight, n1, n2, selfconnect=True):
    """
    :param weight: desired connection weight
    :param n1: number of neurons of network 1
    :param n2: number of neurons of network 2
    :param selfconnect: False if diagonal weights should be zero, True otherwise
    :return: constant value weight matrix (n1 x n2)
    """
    weights = weight*np.ones(n1*n2)
    weights = weights.reshape((n1, n2))
    if not selfconnect:
        weights = weights - np.diag(np.diag(weights))
    return weights


def e_net_connect(net1, net2, threshold, delay, weights):
    net1_net2_syn = list()
    for net1_neuron_i, net1_neuron in enumerate(net1):
        net1_neuron.soma.push()
        for net2_neuron_i, net2_neuron in enumerate(net2):
            net1_net2_syn.append(h.NetCon(net1_neuron.soma(0.5)._ref_v, net2_neuron.synE, threshold, delay,
                                          weights[net1_neuron_i, net2_neuron_i]))
        h.pop_section()
    return net1_net2_syn


def e_net_connect_delay_dist(net1, net2, threshold, delay_distbtn, weights):
    net1_net2_syn = list()
    for net1_neuron_i, net1_neuron in enumerate(net1):
        net1_neuron.soma.push()
        for net2_neuron_i, net2_neuron in enumerate(net2):
            net1_net2_syn.append(h.NetCon(net1_neuron.soma(0.5)._ref_v, net2_neuron.synE, threshold,
                                          delay_distbtn[net1_neuron_i], weights[net1_neuron_i, net2_neuron_i]))
        h.pop_section()
    return net1_net2_syn


def e_ct_net_connect_delay_dist(net1, net2, threshold, delay_distbtn, weights):
    net1_net2_syn = list()
    for net1_neuron_i, net1_neuron in enumerate(net1):
        net1_neuron.soma.push()
        for net2_neuron_i, net2_neuron in enumerate(net2):
            net1_net2_syn.append(h.NetCon(net1_neuron.soma(0.5)._ref_v, net2_neuron.synE_CT, threshold,
                                          delay_distbtn[net1_neuron_i], weights[net1_neuron_i, net2_neuron_i]))
        h.pop_section()
    return net1_net2_syn


def e_ct_net_connect(net1, net2, threshold, delay, weights):
    net1_net2_syn = list()
    for net1_neuron_i, net1_neuron in enumerate(net1):
        net1_neuron.soma.push()
        for net2_neuron_i, net2_neuron in enumerate(net2):
            net1_net2_syn.append(h.NetCon(net1_neuron.soma(0.5)._ref_v, net2_neuron.synE_CT, threshold, delay,
                                          weights[net1_neuron_i, net2_neuron_i]))
        h.pop_section()
    return net1_net2_syn


def i_net_connect(net1, net2, threshold, delay, weights):
    net1_net2_syn = list()
    for net1_neuron_i, net1_neuron in enumerate(net1):
        net1_neuron.soma.push()
        for net2_neuron_i, net2_neuron in enumerate(net2):
            net1_net2_syn.append(h.NetCon(net1_neuron.soma(0.5)._ref_v, net2_neuron.synI, threshold, delay,
                                          weights[net1_neuron_i, net2_neuron_i]))
        h.pop_section()
    return net1_net2_syn
