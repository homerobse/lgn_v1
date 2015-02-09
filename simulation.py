try:
    import sys
    import DLFCN
    sys.setdlopenflags(DLFCN.RTLD_NOW | DLFCN.RTLD_GLOBAL)
except:
    pass

from neuron import h
import numpy as np
import pdb

from plotting import plot_all
from utils import e_net_connect, i_net_connect, e_ct_net_connect

h.load_file("nrngui.hoc")  # load standard run system

h.dt = 1
global dt
dt = h.dt
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

        self.synE_CT = h.Exp2Syn(0.5, sec=self.soma)
        self.synE_CT.tau1 = 5               #  this was made to be slower than the synE - but artificial values were used because the synE time constants couldn't get smaller than 1
        self.synE_CT.tau2 = 15              #  TODO check this with Li, Guido and Bickford 2003

        self.synI = h.Exp2Syn(0.5, sec=self.soma)
        self.synI.e = -100

        self.synI.tau1 = 4
        self.synI.tau2 = 2

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
#         #self.soma.insert("hh")
#         self.soma.insert("hhvitor2")
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


def createnetwork(n_neurons=4):

    network = h.List()
    network_rec = h.List()
 
    for i in range(n_neurons):
        p = Pyrcell()
        network.append(p)
        network_rec.append(h.Vector())
        network_rec[i].record(network[i].soma(0.5)._ref_v)
        
    return network, network_rec


def createnetworkL6(n_neurons=4):

    # network = h.List()
    # network_rec = h.List()
    #
    # for i in range(n_neurons):
    #     p = L6cell()
    #     network.append(p)
    #     network_rec.append(h.Vector())
    #     network_rec[i].record(network[i].soma(0.5)._ref_v)
    return createnetwork(n_neurons)


def simulate(n_runs, total_time, with_V1_L4, with_V1_L6, with_TRN,
             input, input_glut_threshold, input_glut_delay, input_glut_weight, input_gaba_threshold, input_gaba_delay, input_gaba_weight,
             n_E_LGN, n_I_LGN, n_E_L6, n_I_L6, n_E_L4, n_I_L4, n_TRN,
             delay_distbtn_E_L6_LGN, delay_E_L4_E_LGN, delay_E_LGN_I_L4, delay_E_LGN_E_L4, delay_E_LGN_E_L6,
             delay_E_LGN_TRN, delay_E_L4_TRN, delay_E_L6_TRN, delay_E_LGN_I_LGN, delay_I_LGN_E_LGN,
             W_E_LGN_E_LGN, W_I_LGN_I_LGN, W_E_L4_E_L4, W_I_L4_I_L4, W_E_L6_E_L6, W_I_L6_I_L6, W_TRN_TRN,
             W_I_LGN_E_LGN, W_I_L4_E_L4, W_E_L4_I_L4, W_I_L6_E_L6, W_E_L6_I_L6,
             W_E_LGN_TRN, W_TRN_E_LGN, W_E_L6_TRN, W_E_L4_E_L6, W_E_LGN_E_L4, W_E_L4_E_LGN, W_E_L6_E_LGN, W_E_LGN_E_L6, W_E_LGN_I_L4, W_E_L4_TRN,
             connect_E_LGN_E_L4, connect_E_LGN_I_L4, connect_E_L4_E_LGN, connect_E_LGN_E_L6, connect_E_L6_E_LGN, connect_E_L4_TRN, connect_E_L6_TRN,
             connect_E_LGN_TRN, connect_TRN_E_LGN, connect_E_L4_E_L6):

    print "* * * Simulating %d runs * * *" % n_runs
    h.tstop = total_time
    for nsim in range(n_runs):

        # creating LGN network
        I_LGN, I_LGN_rec = createnetwork(n_I_LGN)
        E_LGN, E_LGN_rec = createnetwork(n_E_LGN)

        #create connections in network 1 (LGN)
        GlutGlut_sin = list()
        for neuron_i in range(len(E_LGN)):
            E_LGN[neuron_i].soma.push()
            for neuron_j in range(len(E_LGN)):
                GlutGlut_sin.append(h.NetCon(E_LGN[neuron_i].soma(0.5)._ref_v, E_LGN[neuron_j].synE,
                                             0, 1, W_E_LGN_E_LGN[neuron_i, neuron_j]))
            h.pop_section()


        GABAGABA_sin = list()
        for neuron_i in range(len(I_LGN)):
            I_LGN[neuron_i].soma.push()
            for neuron_j in range(len(I_LGN)):
                GABAGABA_sin.append(h.NetCon(I_LGN[neuron_i].soma(0.5)._ref_v,I_LGN[neuron_j].synI,
                                             0, 1, W_I_LGN_I_LGN[neuron_i, neuron_j]))

        GABAGlut_sin = list()
        GlutGABA_sin = list()
        for Gneuron_i in range(len(I_LGN)):
            for neuron_i in range(len(E_LGN)):
                I_LGN[Gneuron_i].soma.push()
                GABAGlut_sin.append(h.NetCon(I_LGN[Gneuron_i].soma(0.5)._ref_v, E_LGN[neuron_i].synI,
                                             0, delay_I_LGN_E_LGN, W_I_LGN_E_LGN[Gneuron_i, neuron_i]))
                h.pop_section()

                E_LGN[neuron_i].soma.push()
                #GlutGABA_sin.append(h.NetCon(E_LGN[neuron_i].soma(0.5)._ref_v,I_LGN[Gneuron_i].synE,0,delay_E_LGN_I_LGN,1./100000))
                GlutGABA_sin.append(h.NetCon(E_LGN[neuron_i].soma(0.5)._ref_v, I_LGN[Gneuron_i].synE,
                                             0, delay_E_LGN_I_LGN, 1./100000000)) # this is to turn off the connection (E to I)
                h.pop_section()

        E_L4, E_L4_rec = createnetwork(n_E_L4)
        I_L4, I_L4_rec = createnetwork(n_I_L4)
        if with_V1_L4:
            #create connections in network 2  (V1 superficial)

            e_net_connect(E_L4, E_L4, 0, 1, W_E_L4_E_L4)
            i_net_connect(I_L4, I_L4, 0, 1, W_I_L4_I_L4)

            e_net_connect(E_L4, I_L4, 0, 1, W_E_L4_I_L4)
            i_net_connect(I_L4, E_L4, 0, 1, W_I_L4_E_L4)

#ALL-to-ALL connectivity
            if connect_E_LGN_E_L4:
                #extrinsic connections
                #connections from Glutamatergic neurons of network 1 (LGN) to network 2 (V1)

                e_net_connect(E_LGN, E_L4, 0, delay_E_LGN_E_L4, W_E_LGN_E_L4)

            #topographic connectivity        
#            if connect_E_LGN_E_L4:
#                #extrinsic connections
#                #connections from Glutamatergic neurons of network 1 (LGN) to network 2 (V1)
#
#                Glutnt1nt2_sin = list()
#                for neuron_i in range(len(E_LGN)):
#                    E_LGN[neuron_i].soma.push()
#                    Glutnt1nt2_sin.append(h.NetCon(E_LGN[neuron_i].soma(0.5)._ref_v, E_L4[neuron_i].synE,
#                                         0, delay_E_LGN_E_L4, W_E_LGN_E_L4[neuron_i, neuron_i]))
#                    h.pop_section()     
#                    

            if connect_E_L4_E_LGN:
                #connections from Glutamatergic neurons of network 2 (V1) to network 1 (LGN)
                Glutnt2nt1_sin = list()
                for neuron_i in range(len(E_L4)*1/4, len(E_L4)):
                    E_L4[neuron_i].soma.push()
                    for neuron_j in range(len(E_LGN)*1/4, len(E_LGN)):
                        Glutnt2nt1_sin.append(h.NetCon(E_L4[neuron_i].soma(0.5)._ref_v, E_LGN[neuron_j].synE_CT,
                                                       0, delay_E_L4_E_LGN, W_E_L4_E_LGN[neuron_i, neuron_j]))
                    h.pop_section()

            # Population 1) 15 LGN E cells connect to 15 V1 L4 E cells
            # Population 2) 5 LGN E cells connect to 5 V1 L4 I cells
            #
            # Population 1 and population 2 are different
            #
            # Hirsch et al., 1998
            #All-to-ALL connectivity
            if connect_E_LGN_I_L4:
                # connections from Glutamatergic neurons of network (LGN) to network V1 L4
                sin_E_LGN_I_L4 = list()
                len_LGN = len(E_LGN)
                for neuron_i in range(0, len_LGN*1/4):
                    E_LGN[neuron_i].soma.push()
                    for neuron_j in range(0, len(I_L4)*1/4):
                        sin_E_LGN_I_L4.append(h.NetCon(E_LGN[neuron_i].soma(0.5)._ref_v, I_L4[neuron_j].synE,
                                                       0, delay_E_LGN_I_L4, W_E_LGN_I_L4[neuron_i, neuron_j]))
                    h.pop_section()
            
#            #topographic connectivity        
#            if connect_E_LGN_I_L4:
#                # connections from Glutamatergic neurons of network (LGN) to network V1 L4
#                sin_E_LGN_I_L4 = list()
#                len_LGN = len(E_LGN)
#                for neuron_i in range(0, len_LGN*1/4):
#                    E_LGN[neuron_i].soma.push()
#                    sin_E_LGN_I_L4.append(h.NetCon(E_LGN[neuron_i].soma(0.5)._ref_v, I_L4[neuron_i].synE,
#                                                       0, delay_E_LGN_I_L4, W_E_LGN_I_L4[neuron_i, neuron_i]))
#                    h.pop_section()                 

        I_L6, I_L6_rec = createnetworkL6(n_I_L6)
        E_L6, E_L6_rec = createnetworkL6(n_E_L6)
        if with_V1_L6:
            #create connections in network 2  (V1 L6)
            GlutGlutL6_sin = list()
            for neuron_i in range(len(E_L6)):
                E_L6[neuron_i].soma.push()
                for neuron_j in range(len(E_L6)):
                    GlutGlutL6_sin.append(h.NetCon(E_L6[neuron_i].soma(0.5)._ref_v, E_L6[neuron_j].synE,
                                                   0, 1, W_E_L6_E_L6[neuron_i, neuron_j]))
                h.pop_section()

            GABAGABAL6_sin = list()
            for neuron_i in range(len(I_L6)):
                I_L6[neuron_i].soma.push()
                for neuron_j in range(len(I_L4)):
                    GABAGABAL6_sin.append(h.NetCon(I_L6[neuron_i].soma(0.5)._ref_v, I_L6[neuron_j].synI,
                                                   0, 1, W_I_L6_I_L6[neuron_i, neuron_j]))
                h.pop_section()  # TODO we put this pop_secion here. There was no pop_section here before

            GABAGlutL6_sin = list()
            GlutGABAL6_sin = list()
            for Gneuron_i2 in range(len(I_L6)):
                for neuron_i2 in range(len(E_L6)):

                    I_L6[Gneuron_i2].soma.push()
                    GABAGlutL6_sin.append(h.NetCon(I_L6[Gneuron_i2].soma(0.5)._ref_v, E_L6[neuron_i2].synI,
                                                   0, 1, W_I_L6_E_L6[Gneuron_i2, neuron_i2]))
                    h.pop_section()

                    E_L6[neuron_i2].soma.push()
                    GlutGABAL6_sin.append(h.NetCon(E_L6[neuron_i2].soma(0.5)._ref_v, I_L6[Gneuron_i2].synE,
                                                   0, 1, W_E_L6_I_L6[neuron_i2, Gneuron_i2]))
                    h.pop_section()

            # connections from V1 input (L4) layer to L6
            if connect_E_L4_E_L6:
                e_net_connect(E_L4, E_L6, 0, 1, W_E_L4_E_L6)

#ALL-to-ALL connections of feedback
            if connect_E_L6_E_LGN:
                #connections from Glutamatergic neurons of network 2 (V1) to network 1 (LGN)
                k = 0
                GlutL6nt1_sin = list()
                for neuron_i in range(len(E_L6)):
                    E_L6[neuron_i].soma.push()
                    for neuron_j in range(len(E_LGN)):
                        GlutL6nt1_sin.append(h.NetCon(E_L6[neuron_i].soma(0.5)._ref_v, E_LGN[neuron_j].synE_CT,
                                                      0, delay_distbtn_E_L6_LGN[k], W_E_L6_E_LGN[neuron_i, neuron_j]))
                        k += 1
                    h.pop_section()
                    
#            if connect_E_L6_E_LGN:
#                #connections from Glutamatergic neurons of network 2 (V1) to network 1 (LGN)
#                k = 0
#                GlutL6nt1_sin = list()
#                for neuron_i in range(len(E_L6)):
#                    E_L6[neuron_i].soma.push()
#                    GlutL6nt1_sin.append(h.NetCon(E_L6[neuron_i].soma(0.5)._ref_v, E_LGN[neuron_i].synE_CT,
#                                        0, delay_distbtn_E_L6_LGN[k], W_E_L6_E_LGN[neuron_i, neuron_i]))
#                    k += 1
#                    h.pop_section()                    

            if connect_E_LGN_E_L6:
                #connections from Glutamatergic neurons of LGN network to V1 L6 network
                e_net_connect(E_LGN, E_L6, 0, delay_E_LGN_E_L6, W_E_LGN_E_L6)

        #create TRN neurons (inhibitory only)
        TRN, TRN_rec = createnetwork(n_TRN)
        if with_TRN:
            i_net_connect(TRN, TRN, 0, 1, W_TRN_TRN)

            #connections from Glutamatergic neurons of network V1 L4 to TRN
            if with_V1_L4 and connect_E_L4_TRN:
                e_ct_net_connect(E_L4, TRN, 0, delay_E_L4_TRN, W_E_L4_TRN)

            if with_V1_L6 and connect_E_L6_TRN:
                e_ct_net_connect(E_L6, TRN, 0, delay_E_L6_TRN, W_E_L6_TRN)

            #ALL-to-ALL
            if connect_E_LGN_TRN:
                #connections from Glutamatergic neurons of network 1 (LGN) to TRN
                e_net_connect(E_LGN, TRN, 0, delay_E_LGN_TRN, W_E_LGN_TRN)

            #topographic
#            if connect_E_LGN_TRN:
#                #connections from Glutamatergic neurons of LGN to TRN
#                GlutGABAtneurons_sin1 = list()
#                delayGlutGABAtneurons = 1
#                for neuron_i in range(len(E_LGN)):
#                    E_LGN[neuron_i].soma.push()
#                    GlutGABAtneurons_sin1.append(h.NetCon(E_LGN[neuron_i].soma(0.5)._ref_v, TRN[neuron_i].synE, 0, delayGlutGABAtneurons,
#                                                              W_E_LGN_TRN[neuron_i, neuron_i]))
#                     h.pop_section()
            #ALL-to-ALL
            if connect_TRN_E_LGN:
                #connectinos from GABAergic neurons of TRN to LGN
                i_net_connect(TRN, E_LGN, 0, 1, W_TRN_E_LGN)

#            if connect_TRN_E_LGN:
#                #connectinos from GABAergic neurons of TRN to network 1 (LGN)
#                GABAGlutneurons_sin = list()
#                for neuron_i in range(len(TRN)):
#                    TRN[neuron_i].soma.push()
#                    GABAGlutneurons_sin.append(h.NetCon(TRN[neuron_i].soma(0.5)._ref_v, E_LGN[neuron_i].synI,
#                                                            0, 1, W_TRN_E_LGN[neuron_i, neuron_i]))
#                    h.pop_section()

        #generate inputs to network 1
        netStim = list()
        for stim_i in range(0, input['nstims']):
            netStim.append(h.NetStimPois(input['input']))
            netStim[stim_i].start = 0
            netStim[stim_i].mean = input['stimrate']  # 100 = 10 Hz, 10 = 100 Hz, 1 = 1000Hz, 5 = 200 Hz, 6 = 150 Hz
            netStim[stim_i].number = 0

        #stim = h.NetCon(netStim[0],E_LGN[0].synE,0.5,0,4./(100000))
        #stim_rec = h.Vector()
        #stim.record(stim_rec)
        #stim8 = h.NetCon(netStim[0],E_LGN[1].synE,0.5,0,4./(100000))
        #stim9 = h.NetCon(netStim[0],E_LGN[2].synE,0.5,0,4./(100000))
        #stim10 = h.NetCon(netStim[0],E_LGN[3].synE,0.5,0,4./(100000))
        #stim11 = h.NetCon(netStim[0],E_LGN[4].synE,0.5,0,4./(100000))
        #
        #
        #
        #stim3 = h.NetCon(netStim[0],I_LGN[0].synE,0.5,0,4./(100000))
        #stim4 = h.NetCon(netStim[0],I_LGN[1].synE,0.5,0,4./(100000))
        #stim5 = h.NetCon(netStim[0],I_LGN[2].synE,0.5,0,4./(100000))
        #stim6 = h.NetCon(netStim[0],I_LGN[3].synE,0.5,0,4./(100000))
        #stim7 = h.NetCon(netStim[0],I_LGN[4].synE,0.5,0,4./(100000))

        stim = h.NetCon(netStim[0], E_LGN[0].synE, input_glut_threshold, input_glut_delay, input_glut_weight)
        stim_rec = h.Vector()
        stim.record(stim_rec)
        stim2 = h.NetCon(netStim[1], E_LGN[1].synE, input_glut_threshold, input_glut_delay, input_glut_weight)
        #stim3 = h.NetCon(netStim[0], E_LGN[1].synE_CT, input_glut_threshold, input_glut_delay, input_glut_weight)
        stim8 = h.NetCon(netStim[2], E_LGN[2].synE, input_glut_threshold, input_glut_delay, input_glut_weight)
        stim9 = h.NetCon(netStim[3], E_LGN[3].synE, input_glut_threshold, input_glut_delay, input_glut_weight)
        #stim10 = h.NetCon(netStim[0], E_LGN[4].synE_CT, input_glut_threshold, input_glut_delay, input_glut_weight)


        stim3 = h.NetCon(netStim[0], I_LGN[0].synE, input_gaba_threshold, input_gaba_delay, input_gaba_weight)
        stim4 = h.NetCon(netStim[1], I_LGN[1].synE, input_gaba_threshold, input_gaba_delay, input_gaba_weight)
        stim5 = h.NetCon(netStim[2], I_LGN[2].synE, input_gaba_threshold, input_gaba_delay, input_gaba_weight)
        stim6 = h.NetCon(netStim[3], I_LGN[3].synE, input_gaba_threshold, input_gaba_delay, input_gaba_weight)
        #stim7 = h.NetCon(netStim[0], I_LGN[4].synE, input_gaba_threshold, input_gaba_delay, input_gaba_weight)


        timeaxis = h.Vector()
        timeaxis.record(h._ref_t)

        h.run()

        # #all
        meanLGN, meanTRN, meanV1input, meanV1output = plot_all(timeaxis, stim_rec, with_V1_L4, with_V1_L6, with_TRN,
                                                                E_L4_rec, TRN_rec, E_L6_rec, E_LGN_rec,
                                                                I_L6_rec, I_LGN, I_L4, I_LGN_rec, I_L4_rec,
                                                                n_E_LGN, n_E_L6, n_E_L4, n_TRN, n_I_L6)

        ofname = "../data_files/" "sim" + str(nsim+0) + ".txt"

        n = len(timeaxis)
        indx = np.arange(1, n, 40)

        w = np.array(meanLGN)
        u = np.array(meanTRN)
        x = np.array(meanV1input)
        y = np.array(meanV1output)
        z = np.array(timeaxis)
        np.savetxt(ofname, (w[indx], u[indx], x[indx], y[indx], z[indx]))

        print "Progress: %d runs simulated %d runs missing" % (nsim + 1, n_runs - nsim - 1)


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