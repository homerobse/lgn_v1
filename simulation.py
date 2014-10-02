try:
    import sys
    import DLFCN
    sys.setdlopenflags(DLFCN.RTLD_NOW | DLFCN.RTLD_GLOBAL)
except:
    pass

from neuron import h
import numpy as np
from plotting import plot_all


h.load_file("nrngui.hoc")  # load standard run system

h.dt = 1
global dt
dt = h.dt
h.tstop = 500
global nsamples
nsamples = h.tstop*h.dt
h.v_init = -67


class Pyrcell:
    "Pyramidal cell"

    def __init__(self):
       
        self.soma = h.Section(name='soma', cell=self)
        self.soma.L = 30
        self.soma.nseg = 1
        self.soma.diam = 1

        self.soma.Ra = 100
        self.soma.cm = 1

        #self.soma.insert("hh")
        self.soma.insert("hhvitor2")
        
        self.soma.gna_hhvitor2 = 0.039
        self.soma.gkd_hhvitor2 = 0.006
        self.soma.gl_hhvitor2 = 0
        self.soma.gm_hhvitor2 = 0
        self.soma.gt_hhvitor2 = 0
        self.soma.gleak_hhvitor2 = 0.0000273
   
        self.synE = h.Exp2Syn(0.5, sec=self.soma)
        self.synE.tau1 = 1
        self.synE.tau2 = 3

        self.synI = h.Exp2Syn(0.5, sec=self.soma)
        self.synI.e = -100

        self.synI.tau1 = 4
        self.synI.tau2 = 2

        self.stm = h.IClamp(0.5, sec=self.soma)
        self.stm.amp = 0
        self.stm.dur = 1000
        self.stm.delay = 100
        
        
class L6cell:
    "Layer 6 cell"

    def __init__(self):
       
        self.soma = h.Section(name='soma', cell=self)
        self.soma.L=30
        self.soma.nseg=1
        self.soma.diam=1

        self.soma.Ra = 100
        self.soma.cm = 1

        #self.soma.insert("hh")
        self.soma.insert("hhvitor2")
        
        self.soma.gna_hhvitor2 = 0.039
        self.soma.gkd_hhvitor2 = 0.006
        self.soma.gl_hhvitor2 = 0
        self.soma.gm_hhvitor2 = 0
        self.soma.gt_hhvitor2 = 0
        self.soma.gleak_hhvitor2 = 0.0000273
   
        self.synE = h.Exp2Syn(0.5, sec=self.soma)
        self.synE.tau1 = 10
        self.synE.tau2 = 20

        self.synI = h.Exp2Syn(0.5, sec=self.soma)
        self.synI.e = -100
        self.synI.tau1 = 5
        self.synI.tau2 = 40
        
        #self.synI.tau1=1
        #self.synI.tau2=2        
        
        self.stm = h.IClamp(0.5, sec=self.soma)
        self.stm.amp = 0
        self.stm.dur = 1000
        self.stm.delay = 100

def createnetwork(Nneurons=4):

    network = h.List()
    network_rec = h.List()
 
    for Nindex in range(Nneurons):
        p = Pyrcell()
        network.append(p)
        network_rec.append(h.Vector())
        network_rec[Nindex].record(network[Nindex].soma(0.5)._ref_v)
        
    return network, network_rec
    
def createnetworkL6(Nneurons=4):

    network = h.List()
    network_rec = h.List()
 
    for Nindex in range(Nneurons):
        p = L6cell()
        network.append(p)
        network_rec.append(h.Vector())
        network_rec[Nindex].record(network[Nindex].soma(0.5)._ref_v)
        
    return network, network_rec


def simulate(nruns, with_V1_L4, with_V1_L6, with_TRN,
             input, input_glut_threshold, input_glut_delay, input_glut_weight, input_gaba_threshold, input_gaba_delay, input_gaba_weight,
             Nneurons, NGABAn, NneuronsL6, NGABA_L6, NneuronsL4, NGABAL4, NGABA_trn,
             delay_distbtn_E_L6_LGN, delay_E_L4_E_LGN, delay_E_LGN_I_L4,
             W_E_LGN_TRN, W_E_L6_TRN, W_E_L4_E_L6, W_E_LGN_E_L4, W_E_L4_E_LGN, W_E_L6_E_LGN, W_E_LGN_I_L4, W_E_L4_TRN,
             connect_E_LGN_E_L4, connect_E_LGN_I_L4, connect_E_L4_E_LGN, connect_E_LGN_E_L6, connect_E_L6_E_LGN, connect_E_L4_TRN, connect_E_L6_TRN):

    print "Simulating %d runs" % (nruns)

    for nsim in range(nruns):

        # creating LGN network
        GABAneurons, GABAneurons_rec = createnetwork(NGABAn)
        GABAneurons_W = np.random.exponential(1, NGABAn*NGABAn)*4/100000.
        GABAneurons_W = GABAneurons_W.reshape((NGABAn, NGABAn))
        GABAneurons_W = GABAneurons_W - np.diag(np.diag(GABAneurons_W))

        Glutneurons, Glutneurons_rec = createnetwork(Nneurons)
        Glutneurons_W = np.random.exponential(1, Nneurons*Nneurons)*4/100000.
        Glutneurons_W = Glutneurons_W.reshape((Nneurons, Nneurons))
        Glutneurons_W = Glutneurons_W - np.diag(np.diag(Glutneurons_W))

        #create connections in network 1 (LGN)
        GlutGlut_sin = list()
        for neuron_i in range(len(Glutneurons)):
            Glutneurons[neuron_i].soma.push()
            for neuron_j in range(len(Glutneurons)):
                GlutGlut_sin.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v, Glutneurons[neuron_j].synE,
                                             0, 1, Glutneurons_W[neuron_i, neuron_j]))
            h.pop_section()


        GABAGABA_sin = list()
        for neuron_i in range(len(GABAneurons)):
            GABAneurons[neuron_i].soma.push()
            for neuron_j in range(len(GABAneurons)):
                GABAGABA_sin.append(h.NetCon(GABAneurons[neuron_i].soma(0.5)._ref_v,GABAneurons[neuron_j].synI,
                                             0, 1, GABAneurons_W[neuron_i, neuron_j]))

        GABAGlut_sin = list()
        GlutGABA_sin = list()
        for Gneuron_i in range(len(GABAneurons)):
            for neuron_i in range(len(Glutneurons)):
                GABAneurons[Gneuron_i].soma.push()
                GABAGlut_sin.append(h.NetCon(GABAneurons[Gneuron_i].soma(0.5)._ref_v, Glutneurons[neuron_i].synI,
                                             0, 1, 1./100000))
                h.pop_section()

                Glutneurons[neuron_i].soma.push()
                delayGlutGABA = 1
                #GlutGABA_sin.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v,GABAneurons[Gneuron_i].synE,0,delayGlutGABA,1./100000))
                GlutGABA_sin.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v,GABAneurons[Gneuron_i].synE,
                                             0, delayGlutGABA, 1./100000000)) # this is to turn off the connection (E to I)
                h.pop_section()


        if with_V1_L4:
            #create connections in network 2  (V1 superficial)

            GABAneurons2, GABAneurons_rec2 = createnetwork(NGABAL4)
            GABAneurons_W2 = np.random.exponential(1, NGABAL4*NGABAL4)*4/100000.
            GABAneurons_W2 = GABAneurons_W2.reshape((NGABAL4, NGABAL4))
            GABAneurons_W2 = GABAneurons_W2 - np.diag(np.diag(GABAneurons_W2))

            Glutneurons2, Glutneurons_rec2 = createnetwork(NneuronsL4)
            Glutneurons_W2 = np.random.exponential(1, NneuronsL4*NneuronsL4)*4/100000.
            Glutneurons_W2 = Glutneurons_W2.reshape((NneuronsL4, NneuronsL4))
            Glutneurons_W2 = Glutneurons_W2 - np.diag(np.diag(Glutneurons_W2))

            GlutGlut_sin2 = list()
            for neuron_i in range(len(Glutneurons2)):
                Glutneurons2[neuron_i].soma.push()
                for neuron_j in range(len(Glutneurons2)):
                    GlutGlut_sin2.append(h.NetCon(Glutneurons2[neuron_i].soma(0.5)._ref_v,  Glutneurons2[neuron_j].synE,
                                                  0, 1, Glutneurons_W2[neuron_i, neuron_j]))
                h.pop_section()


            GABAGABA_sin2 = list()
            for neuron_i in range(len(GABAneurons2)):
                GABAneurons2[neuron_i].soma.push()
                for neuron_j in range(len(GABAneurons2)):
                    GABAGABA_sin2.append(h.NetCon(GABAneurons2[neuron_i].soma(0.5)._ref_v, GABAneurons2[neuron_j].synI,
                                                  0, 1, GABAneurons_W2[neuron_i, neuron_j]))

            GABAGlut_sin2 = list()
            GlutGABA_sin2 = list()
            for Gneuron_i2 in range(len(GABAneurons2)):
                for neuron_i2 in range(len(Glutneurons2)):

                    GABAneurons2[Gneuron_i2].soma.push()
                    GABAGlut_sin2.append(h.NetCon(GABAneurons2[Gneuron_i2].soma(0.5)._ref_v, Glutneurons2[neuron_i2].synI,
                                                  0, 1, 6./10000))
                    h.pop_section()

                    Glutneurons2[neuron_i2].soma.push()
                    GlutGABA_sin2.append(h.NetCon(Glutneurons2[neuron_i2].soma(0.5)._ref_v, GABAneurons2[Gneuron_i2].synE,
                                                  0, 1, 1./100000))
                    h.pop_section()

            if connect_E_LGN_E_L4:
                #extrinsic connections
                #connections from Glutamatergic neurons of network 1 (LGN) to network 2 (V1)

                Glutnt1nt2_sin = list()
                delayGlutnt1nt2 = 3
                for neuron_i in range(len(Glutneurons)):
                    Glutneurons[neuron_i].soma.push()
                    for neuron_j in range(len(Glutneurons2)):
                        Glutnt1nt2_sin.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v, Glutneurons2[neuron_j].synE,
                                                       0, delayGlutnt1nt2, W_E_LGN_E_L4[neuron_i, neuron_j]))
                    h.pop_section()

            if connect_E_L4_E_LGN:
                #connections from Glutamatergic neurons of network 2 (V1) to network 1 (LGN)
                Glutnt2nt1_sin = list()
                for neuron_i in range(len(Glutneurons2)*1/4, len(Glutneurons2)):
                    Glutneurons2[neuron_i].soma.push()
                    for neuron_j in range(len(Glutneurons)*1/4, len(Glutneurons)):
                        Glutnt2nt1_sin.append(h.NetCon(Glutneurons2[neuron_i].soma(0.5)._ref_v, Glutneurons[neuron_j].synE,
                                                       0, delay_E_L4_E_LGN, W_E_L4_E_LGN[neuron_i, neuron_j]))
                    h.pop_section()

            # Population 1) 15 LGN E cells connect to 15 V1 L4 E cells
            # Population 2) 5 LGN E cells connect to 5 V1 L4 I cells
            #
            # Population 1 and population 2 are different
            #
            # Hirsch et al., 1998
            if connect_E_LGN_I_L4:
                # connections from Glutamatergic neurons of network (LGN) to network V1 L4
                sin_E_LGN_I_L4 = list()
                len_LGN = len(Glutneurons)
                for neuron_i in range(0, len_LGN*1/4):
                    Glutneurons[neuron_i].soma.push()
                    for neuron_j in range(0, len(GABAneurons2)*1/4):
                        sin_E_LGN_I_L4.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v, GABAneurons2[neuron_j].synE,
                                                       0, delay_E_LGN_I_L4, W_E_LGN_I_L4[neuron_i, neuron_j]))
                    h.pop_section()

        if with_V1_L6:
            #create connections in network 2  (V1 L6)

            GABAneuronsL6, GABAneurons_recL6 = createnetworkL6(NGABA_L6)
            GABAneurons_WL6 = np.random.exponential(1, NGABA_L6*NGABA_L6)*4/100000.
            GABAneurons_WL6 = GABAneurons_WL6.reshape((NGABA_L6, NGABA_L6))
            GABAneurons_WL6 = GABAneurons_WL6 - np.diag(np.diag(GABAneurons_WL6))

            GlutneuronsL6, Glutneurons_recL6 = createnetworkL6(NneuronsL6)
            Glutneurons_WL6 = np.random.exponential(1, NneuronsL6*NneuronsL6)*4/100000.
            Glutneurons_WL6 = Glutneurons_W2.reshape((NneuronsL6, NneuronsL6))
            Glutneurons_WL6 = Glutneurons_W2 - np.diag(np.diag(Glutneurons_WL6))

            GlutGlutL6_sin = list()
            for neuron_i in range(len(GlutneuronsL6)):
                GlutneuronsL6[neuron_i].soma.push()
                for neuron_j in range(len(GlutneuronsL6)):
                    GlutGlutL6_sin.append(h.NetCon(GlutneuronsL6[neuron_i].soma(0.5)._ref_v, GlutneuronsL6[neuron_j].synE,
                                                   0, 1, Glutneurons_WL6[neuron_i, neuron_j]))
                h.pop_section()

            GABAGABAL6_sin = list()
            for neuron_i in range(len(GABAneuronsL6)):
                GABAneuronsL6[neuron_i].soma.push()
                for neuron_j in range(len(GABAneurons2)):
                    GABAGABAL6_sin.append(h.NetCon(GABAneuronsL6[neuron_i].soma(0.5)._ref_v, GABAneuronsL6[neuron_j].synI,
                                                   0, 1, GABAneurons_WL6[neuron_i, neuron_j]))
                    # TODO check if we need a pop_section here

            GABAGlutL6_sin = list()
            GlutGABAL6_sin = list()
            for Gneuron_i2 in range(len(GABAneuronsL6)):
                for neuron_i2 in range(len(GlutneuronsL6)):

                    GABAneuronsL6[Gneuron_i2].soma.push()
                    GABAGlutL6_sin.append(h.NetCon(GABAneuronsL6[Gneuron_i2].soma(0.5)._ref_v, GlutneuronsL6[neuron_i2].synI,
                                                   0, 1, 6./10000))
                    h.pop_section()

                    GlutneuronsL6[neuron_i2].soma.push()
                    GlutGABAL6_sin.append(h.NetCon(GlutneuronsL6[neuron_i2].soma(0.5)._ref_v, GABAneuronsL6[Gneuron_i2].synE,
                                                   0, 1, 1./100000))
                    h.pop_section()

            #connections from V1 input layer to L6

            GlutGlut_L4L6_sin = list()
            for neuron_i in range(len(Glutneurons2)):
                Glutneurons2[neuron_i].soma.push()
                for neuron_j in range(len(GlutneuronsL6)):
                    GlutGlut_L4L6_sin.append(h.NetCon(Glutneurons2[neuron_i].soma(0.5)._ref_v, GlutneuronsL6[neuron_j].synE,
                                                      0, 1, W_E_L4_E_L6[neuron_i, neuron_j]))
                h.pop_section()

            if connect_E_L6_E_LGN:
                #connections from Glutamatergic neurons of network 2 (V1) to network 1 (LGN)

                k = 0
                GlutL6nt1_sin = list()
                for neuron_i in range(len(GlutneuronsL6)):
                    GlutneuronsL6[neuron_i].soma.push()
                    for neuron_j in range(len(Glutneurons)):
                        GlutL6nt1_sin.append(h.NetCon(GlutneuronsL6[neuron_i].soma(0.5)._ref_v, Glutneurons[neuron_j].synE,
                                                      0, delay_distbtn_E_L6_LGN[k], W_E_L6_E_LGN[neuron_i, neuron_j]))
                        k += 1
                    h.pop_section()

        if with_TRN:
            #create TRN neurons (inhibitory only)
            GABAneurons_trn, GABAneurons_trn_rec = createnetwork(NGABA_trn)
            GABAneurons_trnW = np.random.exponential(1, NGABA_trn*NGABA_trn)*4/100000.
            GABAneurons_trnW = GABAneurons_trnW.reshape((NGABA_trn, NGABA_trn))
            GABAneurons_trnW = GABAneurons_trnW - np.diag(np.diag(GABAneurons_trnW))

            GABAGABA_trn_sin = list()
            for neuron_i in range(len(GABAneurons_trn)):
                GABAneurons_trn[neuron_i].soma.push()
                for neuron_j in range(len(GABAneurons_trn)):
                    GABAGABA_trn_sin.append(h.NetCon(GABAneurons_trn[neuron_i].soma(0.5)._ref_v,GABAneurons_trn[neuron_j].synI,
                                                     0, 1, GABAneurons_trnW[neuron_i, neuron_j]))
                    # h.pop_section() #  TODO tem que ter esse pop_section??

            #connections from Glutamatergic neurons of network 2 (V1) to TRN
            if with_V1_L4 and connect_E_L4_TRN:
                GlutGABAtneurons_sin2 = list()
                for neuron_i in range(len(Glutneurons2)):
                    Glutneurons2[neuron_i].soma.push()
                    for neuron_j in range(len(GABAneurons_trn)):
                        GlutGABAtneurons_sin2.append(h.NetCon(Glutneurons2[neuron_i].soma(0.5)._ref_v, GABAneurons_trn[neuron_j].synE, 0, 5,
                                                              W_E_L4_TRN[neuron_i, neuron_j]))
                    h.pop_section()

            if with_V1_L6 and connect_E_L6_TRN:

                GlutGABAtneurons_sinL6trn = list()
                for neuron_i in range(len(GlutneuronsL6)):
                    GlutneuronsL6[neuron_i].soma.push()
                    for neuron_j in range(len(GABAneurons_trn)):
                        GlutGABAtneurons_sinL6trn.append(h.NetCon(GlutneuronsL6[neuron_i].soma(0.5)._ref_v, GABAneurons_trn[neuron_j].synE, 0, 5,
                                                                  W_E_L6_TRN[neuron_i, neuron_j]))
                    h.pop_section()

            #connections from Glutamatergic neurons of network 1 (LGN) to TRN
            GlutGABAtneurons_sin1 = list()
            delayGlutGABAtneurons = 1
            for neuron_i in range(len(Glutneurons)):
                Glutneurons[neuron_i].soma.push()
                for neuron_j in range(len(GABAneurons_trn)):
                    GlutGABAtneurons_sin1.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v, GABAneurons_trn[neuron_j].synE, 0, delayGlutGABAtneurons,
                                                          W_E_LGN_TRN[neuron_i, neuron_j]))
                h.pop_section()

            #connectinos from GABAergic neurons of TRN to network 1 (LGN)
            GABAGlutneurons_W_trn_l6 = 1/1000000.*np.random.exponential(1, NGABA_trn*Nneurons)
            GABAGlutneurons_W_trn_l6 = GABAGlutneurons_W_trn_l6.reshape(NGABA_trn, Nneurons)

            GABAGlutneurons_sin = list()
            for neuron_i in range(len(GABAneurons_trn)):
                GABAneurons_trn[neuron_i].soma.push()
                for neuron_j in range(len(Glutneurons)):
                    GABAGlutneurons_sin.append(h.NetCon(GABAneurons_trn[neuron_i].soma(0.5)._ref_v,Glutneurons[neuron_j].synE,0,1,
                                                        GABAGlutneurons_W_trn_l6[neuron_i, neuron_j]))
                h.pop_section()


        #generate inputs to network 1
        netStim = list()
        for stim_i in range(0, input['nstims']):
            netStim.append(h.NetStimPois(input['input']))
            netStim[stim_i].start = 0
            netStim[stim_i].mean = input['stimrate']  # 100 = 10 Hz, 10 = 100 Hz, 1 = 1000Hz, 5 = 200 Hz, 6 = 150 Hz
            netStim[stim_i].number = 0

        #stim = h.NetCon(netStim[0],Glutneurons[0].synE,0.5,0,4./(100000))
        #stim_rec = h.Vector()
        #stim.record(stim_rec)
        #stim8 = h.NetCon(netStim[0],Glutneurons[1].synE,0.5,0,4./(100000))
        #stim9 = h.NetCon(netStim[0],Glutneurons[2].synE,0.5,0,4./(100000))
        #stim10 = h.NetCon(netStim[0],Glutneurons[3].synE,0.5,0,4./(100000))
        #stim11 = h.NetCon(netStim[0],Glutneurons[4].synE,0.5,0,4./(100000))
        #
        #
        #
        #stim3 = h.NetCon(netStim[0],GABAneurons[0].synE,0.5,0,4./(100000))
        #stim4 = h.NetCon(netStim[0],GABAneurons[1].synE,0.5,0,4./(100000))
        #stim5 = h.NetCon(netStim[0],GABAneurons[2].synE,0.5,0,4./(100000))
        #stim6 = h.NetCon(netStim[0],GABAneurons[3].synE,0.5,0,4./(100000))
        #stim7 = h.NetCon(netStim[0],GABAneurons[4].synE,0.5,0,4./(100000))

        stim = h.NetCon(netStim[0], Glutneurons[0].synE, input_glut_threshold, input_glut_delay, input_glut_weight)
        stim_rec = h.Vector()
        stim.record(stim_rec)
        stim2 = h.NetCon(netStim[1], Glutneurons[1].synE, input_glut_threshold, input_glut_delay, input_glut_weight)
        #stim3 = h.NetCon(netStim[0], Glutneurons[1].synE, input_glut_threshold, input_glut_delay, input_glut_weight)
        stim8 = h.NetCon(netStim[2], Glutneurons[2].synE, input_glut_threshold, input_glut_delay, input_glut_weight)
        stim9 = h.NetCon(netStim[3], Glutneurons[3].synE, input_glut_threshold, input_glut_delay, input_glut_weight)
        #stim10 = h.NetCon(netStim[0], Glutneurons[4].synE, input_glut_threshold, input_glut_delay, input_glut_weight)


        stim3 = h.NetCon(netStim[0], GABAneurons[0].synE, input_gaba_threshold, input_gaba_delay, input_gaba_weight)
        stim4 = h.NetCon(netStim[1], GABAneurons[1].synE, input_gaba_threshold, input_gaba_delay, input_gaba_weight)
        stim5 = h.NetCon(netStim[2], GABAneurons[2].synE, input_gaba_threshold, input_gaba_delay, input_gaba_weight)
        #stim6 = h.NetCon(netStim[0], GABAneurons[3].synE, input_gaba_threshold, input_gaba_delay, input_gaba_weight)
        #stim7 = h.NetCon(netStim[0], GABAneurons[4].synE, input_gaba_threshold, input_gaba_delay, input_gaba_weight)


        timeaxis = h.Vector()
        timeaxis.record(h._ref_t)

        h.run()

        # all
        #meanLGN, meanTRN, meanV1input, meanV1output = plot_all(timeaxis, stim_rec, with_V1_L4, with_V1_L6, with_TRN,
                                                               #Glutneurons_rec2, GABAneurons_trn_rec, Glutneurons_recL6, Glutneurons_rec,
                                                               #GABAneurons_recL6, GABAneurons, GABAneurons2, GABAneurons_rec, GABAneurons_rec2,
                                                               #Nneurons, NneuronsL6, NneuronsL4, NGABA_trn, NGABA_L6)
        #only LGN TRN
        # meanLGN, meanTRN = plot_all(timeaxis, stim_rec,  with_V1_L4, with_V1_L6, with_TRN,
        #                                GABAneurons_trn_rec, Glutneurons_rec,
        #                                GABAneurons, GABAneurons_rec,
        #                                Nneurons, NGABA_trn)

        # LGN, TRN, L4, no L6
        meanLGN, meanTRN, meanV1input = plot_all(timeaxis, stim_rec, with_V1_L4, with_V1_L6, with_TRN,
                                                               Glutneurons_rec2, GABAneurons_trn_rec, Glutneurons_rec,
                                                               GABAneurons, GABAneurons2, GABAneurons_rec, GABAneurons_rec2,
                                                               Nneurons, NneuronsL4, NGABA_trn)


        ofname = "../data_files/" "sim" + str(nsim+96) + ".txt"

        n = len(timeaxis)
        indx = np.arange(1, n, 40)

        w = np.array(meanLGN)
        u = np.array(meanTRN)
        #x = np.array(meanV1input)
        #y = np.array(meanV1output)
        z = np.array(timeaxis)
        #np.savetxt(ofname, (w[indx], u[indx], x[indx], y[indx], z[indx]))
        np.savetxt(ofname, (w[indx], u[indx], z[indx]))

        print "Progress: %d runs simulated %d runs missing" % (nsim + 1, nruns - nsim - 1)


def detect_spikes(voltagesignals):
    nneurons = len(voltagesignals)
    for neuron_i in range(nneurons):
        auxsignal = np.asarray(voltagesignals[neuron_i])
        thrscross_ind = [i for (i, val) in enumerate(auxsignal) if val > 30]

        spikeidx = list()
        aux = 0
        for thrscross_i in range(len(thrscross_ind)):
            if thrscross_ind[thrscross_i] > aux:
                auxspan = range(thrscross_ind[thrscross_i], thrscross_ind[thrscross_i] + 15)
                spikeidx.append(auxspan.index(max(auxspan))+thrscross_ind[thrscross_i])
                aux = thrscross_ind[thrscross_i]+15
                
    return spikeidx