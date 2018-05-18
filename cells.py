from neuron import h

class Pyrcell:
    "Pyramidal cell"

    def __init__(self):

        self.soma = h.Section(name='soma', cell=self)
        self.soma.L = 30
        self.soma.nseg = 1
        self.soma.diam = 1

        self.soma.Ra = 100
        self.soma.cm = 1

        #self.soma.insert("hhvitor2")
        self.soma.insert("hh_original")

        # self.soma.gna_hhvitor2 = 0.039
        # self.soma.gkd_hhvitor2 = 0.006
        # self.soma.gl_hhvitor2 = 0
        # self.soma.gm_hhvitor2 = 0
        # self.soma.gt_hhvitor2 = 0
        # self.soma.gleak_hhvitor2 = 0.0000273

        self.synE = h.Exp2Syn(0.5, sec=self.soma)
        self.synE.e = 0
        self.synE.tau1 = 0.1  #  Christoph Borgers, Nancy Kopell (2013) "Synchronization in Networks of Excitatory and Inhibitory Neurons with Sparse, Random Connectivity" pg 514
        self.synE.tau2 = 2    #  Christoph Borgers, Nancy Kopell (2013) "Synchronization in Networks of Excitatory and Inhibitory Neurons with Sparse, Random Connectivity" pg 517

        self.synE_CT = h.Exp2Syn(0.5, sec=self.soma)
        self.synE_CT.e = 0
        self.synE_CT.tau1 = 5               # using time constants for NMDA synapses reported at Wulfram Gerstner's book http://neuronaldynamics.epfl.ch/online/Ch3.S1.html
        self.synE_CT.tau2 = 40              # TODO check this comparing with Li, Guido and Bickford 2003

        self.synI = h.Exp2Syn(0.5, sec=self.soma)
        self.synI.e = -100

        self.synI.tau1 = 0.1  #  Christoph Borgers, Nancy Kopell (2013) "Synchronization in Networks of Excitatory and Inhibitory Neurons with Sparse, Random Connectivity" pg 514
        self.synI.tau2 = 10   #  Christoph Borgers, Nancy Kopell (2013) "Synchronization in Networks of Excitatory and Inhibitory Neurons with Sparse, Random Connectivity" pg 517

        self.stm = h.IClamp(0.5, sec=self.soma)
        self.stm.amp = 0
        self.stm.dur = 1000
        self.stm.delay = 100


# class L6cell:
#     "Layer 6 cell"
#
#     def __init__(self):
#
#         self.soma = h.Section(name='soma', cell=self)
#         self.soma.L = 30
#         self.soma.nseg = 1
#         self.soma.diam = 1
#
#         self.soma.Ra = 100
#         self.soma.cm = 1
#
#         #self.soma.insert("hhvitor2")
#         self.soma.insert("hh_original")
#
#         self.soma.gna_hhvitor2 = 0.039
#         self.soma.gkd_hhvitor2 = 0.006
#         self.soma.gl_hhvitor2 = 0
#         self.soma.gm_hhvitor2 = 0
#         self.soma.gt_hhvitor2 = 0
#         self.soma.gleak_hhvitor2 = 0.0000273
#
#         self.synE = h.Exp2Syn(0.5, sec=self.soma)
#         self.synE.tau1 = 10
#         self.synE.tau2 = 20
#
#         self.synI = h.Exp2Syn(0.5, sec=self.soma)
#         self.synI.e = -100
#
#         self.synI.tau1 = 5
#         self.synI.tau2 = 40
#
#         self.stm = h.IClamp(0.5, sec=self.soma)
#         self.stm.amp = 0
#         self.stm.dur = 1000
#         self.stm.delay = 100