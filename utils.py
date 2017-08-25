import numpy as np
from neuron import h
from cells import Pyrcell
from random import random


def createNetwork(n_neurons=4):

    network = h.List()
    network_rec = h.List()

    for i in range(n_neurons):
        p = Pyrcell()
        network.append(p)
        network_rec.append(h.Vector())
        network_rec[i].record(network[i].soma(0.5)._ref_v)

    return network, network_rec


def createNetworkL6(n_neurons=4):

    # network = h.List()
    # network_rec = h.List()
    #
    # for i in range(n_neurons):
    #     p = L6cell()
    #     network.append(p)
    #     network_rec.append(h.Vector())
    #     network_rec[i].record(network[i].soma(0.5)._ref_v)
    return createNetwork(n_neurons)


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

def partial_e_net_connect(net1, net2, net1_prop_b, net1_prop_e, net2_prop_b, net2_prop_e, threshold, delay, weights):
    """
    Partially connect neurons in two networks with excitatory synapse
    :param net1: First network list (h.List()) of neurons
    :param net2: Second network list (h.List()) of neurons
    :param net1_prop_b: Proportion in network 1 by which we begin to connect neurons
    :param net1_prop_e: Proportion in network 1 by which we end to connect neurons
    :param net2_prop_b: Proportion in network 2 by which we begin to connect neurons
    :param net2_prop_e: Proportion in network 2 by which we end to connect neurons
    :param threshold: voltage threshold that generates spike in neuron in net1
    :param delay: time between spike in net1 and PSP in net2 (ms)
    :param weights: matrix of connection weights (strength of connection)
    :return:
    """
    net1_net2_syn = list()
    len_net1 = len(net1)
    len_net2 = len(net2)
    for neuron_i in range(int(len_net1*net1_prop_b), int(len_net1*net1_prop_e)):
        net1[neuron_i].soma.push()
        for neuron_j in range(int(len_net2*net2_prop_b), int(len_net2*net2_prop_e)):
            net1_net2_syn.append(h.NetCon(net1[neuron_i].soma(0.5)._ref_v, net2[neuron_j].synE,
                                           threshold, delay, weights[neuron_i, neuron_j]))
        h.pop_section()


def e_net_connect(net1, net2, threshold, delay, weights, prob):
    """
    Connects two networks with an excitatory synapse
    :param net1: First network list (h.List()) of neurons
    :param net2: Second network list (h.List()) of neurons
    :param threshold: voltage threshold that generates spike in neuron in net1
    :param delay: time between spike in net1 and PSP in net2 (ms)
    :param weights: matrix of connection weights (strength of connection)
    :param prob: connection probability
    :return: list of synapses
    """
    net1_net2_syn = list()
    for net1_neuron_i, net1_neuron in enumerate(net1):
        net1_neuron.soma.push()
        for net2_neuron_i, net2_neuron in enumerate(net2):
            if random() < prob:
                net1_net2_syn.append(h.NetCon(net1_neuron.soma(0.5)._ref_v, net2_neuron.synE, threshold, delay,
                                              weights[net1_neuron_i, net2_neuron_i]))
        h.pop_section()
    return net1_net2_syn


def e_net_connect_delay_dist(net1, net2, threshold, delay_distbtn, weights, prob):
    net1_net2_syn = list()
    for net1_neuron_i, net1_neuron in enumerate(net1):
        net1_neuron.soma.push()
        for net2_neuron_i, net2_neuron in enumerate(net2):
            if random() < prob:
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


def i_net_connect(net1, net2, threshold, delay, weights, prob):
    net1_net2_syn = list()
    for net1_neuron_i, net1_neuron in enumerate(net1):
        net1_neuron.soma.push()
        for net2_neuron_i, net2_neuron in enumerate(net2):
            if random() < prob:
                net1_net2_syn.append(h.NetCon(net1_neuron.soma(0.5)._ref_v, net2_neuron.synI, threshold, delay,
                                              weights[net1_neuron_i, net2_neuron_i]))
        h.pop_section()
    return net1_net2_syn
