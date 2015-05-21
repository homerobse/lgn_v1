from neuron import h
import numpy as np

from plotting import plot_all
from utils import e_net_connect, i_net_connect, e_ct_net_connect, e_net_connect_delay_dist, e_ct_net_connect_delay_dist

h.load_file("nrngui.hoc")  # load standard run system
# h.dt = 1
# global dt
# dt = h.dt
# global nsamples
# nsamples = h.tstop*h.dt
# h.v_init = -67


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
             n_e_lgn, n_i_lgn, n_e_l6, n_i_l6, n_e_l4, n_I_L4, n_trn,
             delay_distbtn_E_L6_LGN, delay_E_L4_E_LGN, delay_E_LGN_I_L4, delay_E_LGN_E_L4, delay_E_LGN_E_L6,
             delay_E_LGN_TRN, delay_E_L4_TRN, delay_distbtn_E_L6_TRN, delay_E_LGN_I_LGN, delay_I_LGN_E_LGN, delay_E_LGN_I_L6,
             lgn_params, W_E_L4_E_L4, W_I_L4_I_L4, w_e_l6_e_l6, w_i_l6_i_l6, W_TRN_TRN,
             W_I_L4_E_L4, W_E_L4_I_L4, w_i_l6_e_l6, w_e_l6_i_l6,
             W_E_LGN_TRN, W_TRN_E_LGN, w_e_l6_trn, W_E_L4_E_L6, W_E_LGN_E_L4, W_E_L4_E_LGN,
             w_e_l6_e_lgn, W_E_LGN_E_L6, W_E_LGN_I_L6, W_E_LGN_I_L4, W_E_L4_TRN,
             connect_E_LGN_E_L4, connect_E_LGN_I_L4, connect_E_L4_E_LGN, connect_E_LGN_I_L6, connect_E_LGN_E_L6, connect_E_L6_E_LGN, connect_E_L4_TRN, connect_E_L6_TRN,
             connect_E_LGN_TRN, connect_TRN_E_LGN, connect_E_L4_E_L6):

    print "* * * Simulating %d runs * * *" % n_runs
    h.tstop = total_time
    for nsim in range(n_runs):

        # creating LGN network
        i_lgn, I_LGN_rec = createnetwork(n_i_lgn)
        e_lgn, E_LGN_rec = createnetwork(n_e_lgn)

        #create connections in network 1 (LGN)
        e_lgn_e_lgn_syn = e_net_connect(e_lgn, e_lgn, 0, 1, lgn_params['w_e_lgn_e_lgn'])
        i_lgn_i_lgn_syn = i_net_connect(i_lgn, i_lgn, 0, 1, lgn_params['w_e_lgn_e_lgn'])
        i_lgn_e_lgn_syn = i_net_connect(i_lgn, e_lgn, 0, 1, lgn_params['w_e_lgn_e_lgn'])
        e_lgn_i_lgn_syn = e_net_connect(e_lgn, i_lgn, 0, 1, lgn_params['w_e_lgn_e_lgn'])  # weight should be set to zero

        e_l4, E_L4_rec = createnetwork(n_e_l4)
        I_L4, I_L4_rec = createnetwork(n_I_L4)
        if with_V1_L4:
            #create connections in network 2  (V1 superficial)

            e_l4_e_l4_sin = e_net_connect(e_l4, e_l4, 0, 1, W_E_L4_E_L4)
            i_l4_i_l4_sin = i_net_connect(I_L4, I_L4, 0, 1, W_I_L4_I_L4)
            e_l4_i_l4_sin = e_net_connect(e_l4, I_L4, 0, 1, W_E_L4_I_L4)
            i_l4_e_l4_sin = i_net_connect(I_L4, e_l4, 0, 1, W_I_L4_E_L4)

            # Population 1) 15 LGN E cells connect to 15 V1 L4 E cells
            # Population 2) 5 LGN E cells connect to 5 V1 L4 I cells
            #
            # Population 1 and population 2 are different
            #
            # Hirsch et al., 1998

            #extrinsic connections
            #ALL-to-ALL connectivity
            if connect_E_LGN_E_L4:
                #connections from Glutamatergic neurons of network LGN to network V1 L4
                # e_lgn_e_l4_sin = e_net_connect(e_lgn, e_l4, 0, delay_E_LGN_E_L4, W_E_LGN_E_L4)
                sin_E_LGN_E_L4 = list()
                len_LGN = len(e_lgn)
                for neuron_i in range(len(e_lgn)*1/4, len(e_lgn)):
                    e_lgn[neuron_i].soma.push()
                    for neuron_j in range(len(e_l4)*1/4, len(e_l4)):
                        sin_E_LGN_E_L4.append(h.NetCon(e_lgn[neuron_i].soma(0.5)._ref_v, e_l4[neuron_j].synE,
                                                       0, delay_E_LGN_E_L4, W_E_LGN_E_L4[neuron_i, neuron_j]))
                    h.pop_section()


                    #topographic connectivity
                    # if connect_E_LGN_E_L4:
                    #     #extrinsic connections
                    #     #connections from Glutamatergic neurons of network 1 (LGN) to network 2 (V1)
                    #
                    #     Glutnt1nt2_sin = list()
                    #     for neuron_i in range(len(e_lgn)):
                    #         e_lgn[neuron_i].soma.push()
                    #         Glutnt1nt2_sin.append(h.NetCon(e_lgn[neuron_i].soma(0.5)._ref_v, e_l4[neuron_i].synE,
                    #                              0, delay_E_LGN_E_L4, W_E_LGN_E_L4[neuron_i, neuron_i]))
                    #         h.pop_section()


            if connect_E_L4_E_LGN:
                #connections from Glutamatergic neurons of network 2 (V1) to network 1 (LGN)
                Glutnt2nt1_sin = list()
                for neuron_i in range(len(e_l4)*1/4, len(e_l4)):
                    e_l4[neuron_i].soma.push()
                    for neuron_j in range(len(e_lgn)*1/4, len(e_lgn)):
                        Glutnt2nt1_sin.append(h.NetCon(e_l4[neuron_i].soma(0.5)._ref_v, e_lgn[neuron_j].synE_CT,
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
                len_LGN = len(e_lgn)
                for neuron_i in range(0, len_LGN*1/4):
                    e_lgn[neuron_i].soma.push()
                    for neuron_j in range(0, len(I_L4)*1/4):
                        sin_E_LGN_I_L4.append(h.NetCon(e_lgn[neuron_i].soma(0.5)._ref_v, I_L4[neuron_j].synE,
                                                       0, delay_E_LGN_I_L4, W_E_LGN_I_L4[neuron_i, neuron_j]))
                    h.pop_section()

                #            #topographic connectivity
                #            if connect_E_LGN_I_L4:
                #                # connections from Glutamatergic neurons of network (LGN) to network V1 L4
                #                sin_E_LGN_I_L4 = list()
                #                len_LGN = len(e_lgn)
                #                for neuron_i in range(0, len_LGN*1/4):
                #                    e_lgn[neuron_i].soma.push()
                #                    sin_E_LGN_I_L4.append(h.NetCon(e_lgn[neuron_i].soma(0.5)._ref_v, I_L4[neuron_i].synE,
                #                                                       0, delay_E_LGN_I_L4, W_E_LGN_I_L4[neuron_i, neuron_i]))
                #                    h.pop_section()

        i_l6, I_L6_rec = createnetworkL6(n_i_l6)
        e_l6, E_L6_rec = createnetworkL6(n_e_l6)
        if with_V1_L6:
            # create connections in network 2  (V1 L6)
            e_l6_e_l6_sin = e_net_connect(e_l6, e_l6, 0, 1, w_e_l6_e_l6)
            i_l6_i_l6_sin = i_net_connect(i_l6, i_l6, 0, 1, w_i_l6_i_l6)
            i_l6_e_l6_syn = i_net_connect(i_l6, e_l6, 0, 1, w_i_l6_e_l6)
            e_l6_i_l6_syn = e_net_connect(e_l6, i_l6, 0, 1, w_e_l6_i_l6)

            # connections from V1 input (L4) layer to L6
            if connect_E_L4_E_L6:
                e_l4_e_l6_sin = e_net_connect(e_l4, e_l6, 0, 1, W_E_L4_E_L6)

#ALL-to-ALL connections of feedback
            if connect_E_L6_E_LGN:
                e_l6_e_lgn_sin = e_ct_net_connect_delay_dist(e_l6, e_lgn, 0, delay_distbtn_E_L6_LGN, w_e_l6_e_lgn)

#            if connect_E_L6_E_LGN:
#                #connections from Glutamatergic neurons of network 2 (V1) to network 1 (LGN)
#                k = 0
#                GlutL6nt1_sin = list()
#                for neuron_i in range(len(e_l6)):
#                    e_l6[neuron_i].soma.push()
#                    GlutL6nt1_sin.append(h.NetCon(e_l6[neuron_i].soma(0.5)._ref_v, e_lgn[neuron_i].synE_CT,
#                                        0, delay_distbtn_E_L6_LGN[k], w_e_l6_e_lgn[neuron_i, neuron_i]))
#                    k += 1
#                    h.pop_section()

            # TODO: Connectivity as Hirsch
            if connect_E_LGN_E_L6:
                e_lgn_e_l6_syn = e_net_connect(e_lgn, e_l6, 0, delay_E_LGN_E_L6, W_E_LGN_E_L6)

            # TODO: Connectivity as Hirsch
            if connect_E_LGN_I_L6:
                e_lgn_i_l6_syn = e_net_connect(e_lgn, i_l6, 0, delay_E_LGN_I_L6, W_E_LGN_I_L6)

        #create trn neurons (inhibitory only)
        trn, TRN_rec = createnetwork(n_trn)
        if with_TRN:
            trn_trn_syn = i_net_connect(trn, trn, 0, 1, W_TRN_TRN)

            #connections from Glutamatergic neurons of network V1 L4 to trn
            if with_V1_L4 and connect_E_L4_TRN:
                e_l4_trn_syn = e_net_connect(e_l4, trn, 0, delay_E_L4_TRN, W_E_L4_TRN)

            if with_V1_L6 and connect_E_L6_TRN:
                e_l6_trn_syn = e_net_connect_delay_dist(e_l6, trn, 0, delay_distbtn_E_L6_TRN, w_e_l6_trn)

            #ALL-to-ALL
            if connect_E_LGN_TRN:
                #connections from Glutamatergic neurons of network 1 (LGN) to trn
                e_lgn_trn_syn = e_net_connect(e_lgn, trn, 0, delay_E_LGN_TRN, W_E_LGN_TRN)

            #topographic
#            if connect_E_LGN_TRN:
#                #connections from Glutamatergic neurons of LGN to trn
#                GlutGABAtneurons_sin1 = list()
#                delayGlutGABAtneurons = 1
#                for neuron_i in range(len(e_lgn)):
#                    e_lgn[neuron_i].soma.push()
#                    GlutGABAtneurons_sin1.append(h.NetCon(e_lgn[neuron_i].soma(0.5)._ref_v, trn[neuron_i].synE, 0, delayGlutGABAtneurons,
#                                                              W_E_LGN_TRN[neuron_i, neuron_i]))
#                     h.pop_section()
            #ALL-to-ALL
            if connect_TRN_E_LGN:
                trn_e_lgn_sin = i_net_connect(trn, e_lgn, 0, 1, W_TRN_E_LGN)

        #generate inputs to network 1
        netStim = list()
        for stim_i in range(0, input['nstims']):
            netStim.append(h.NetStimPois(input['input']))
            netStim[stim_i].start = 0
            netStim[stim_i].mean = input['stimrate']  # 100 = 10 Hz, 10 = 100 Hz, 1 = 1000Hz, 5 = 200 Hz, 6 = 150 Hz
            netStim[stim_i].number = 0

        #stim = h.NetCon(netStim[0],e_lgn[0].synE,0.5,0,4./(100000))
        #stim_rec = h.Vector()
        #stim.record(stim_rec)
        #stim8 = h.NetCon(netStim[0],e_lgn[1].synE,0.5,0,4./(100000))
        #stim9 = h.NetCon(netStim[0],e_lgn[2].synE,0.5,0,4./(100000))
        #stim10 = h.NetCon(netStim[0],e_lgn[3].synE,0.5,0,4./(100000))
        #stim11 = h.NetCon(netStim[0],e_lgn[4].synE,0.5,0,4./(100000))
        #
        #
        #
        #stim3 = h.NetCon(netStim[0],i_lgn[0].synE,0.5,0,4./(100000))
        #stim4 = h.NetCon(netStim[0],i_lgn[1].synE,0.5,0,4./(100000))
        #stim5 = h.NetCon(netStim[0],i_lgn[2].synE,0.5,0,4./(100000))
        #stim6 = h.NetCon(netStim[0],i_lgn[3].synE,0.5,0,4./(100000))
        #stim7 = h.NetCon(netStim[0],i_lgn[4].synE,0.5,0,4./(100000))

        stim = h.NetCon(netStim[0], e_lgn[0].synE, input_glut_threshold, input_glut_delay, input_glut_weight)
        stim_rec = h.Vector()
        stim.record(stim_rec)
        stim2 = h.NetCon(netStim[1], e_lgn[1].synE, input_glut_threshold, input_glut_delay, input_glut_weight)
        #stim3 = h.NetCon(netStim[0], e_lgn[1].synE_CT, input_glut_threshold, input_glut_delay, input_glut_weight)
        stim8 = h.NetCon(netStim[2], e_lgn[2].synE, input_glut_threshold, input_glut_delay, input_glut_weight)
        stim9 = h.NetCon(netStim[3], e_lgn[3].synE, input_glut_threshold, input_glut_delay, input_glut_weight)
        #stim10 = h.NetCon(netStim[0], e_lgn[4].synE_CT, input_glut_threshold, input_glut_delay, input_glut_weight)


        stim3 = h.NetCon(netStim[0], i_lgn[0].synE, input_gaba_threshold, input_gaba_delay, input_gaba_weight)
        stim4 = h.NetCon(netStim[1], i_lgn[1].synE, input_gaba_threshold, input_gaba_delay, input_gaba_weight)
        stim5 = h.NetCon(netStim[2], i_lgn[2].synE, input_gaba_threshold, input_gaba_delay, input_gaba_weight)
        stim6 = h.NetCon(netStim[3], i_lgn[3].synE, input_gaba_threshold, input_gaba_delay, input_gaba_weight)
        #stim7 = h.NetCon(netStim[0], i_lgn[4].synE, input_gaba_threshold, input_gaba_delay, input_gaba_weight)


        timeaxis = h.Vector()
        timeaxis.record(h._ref_t)

        h.run()

        #all
        meanLGN, meanTRN, meanV1input, meanV1output = plot_all(timeaxis, stim_rec, with_V1_L4, with_V1_L6, with_TRN,
                                                               E_L4_rec, TRN_rec, E_L6_rec, E_LGN_rec,
                                                               I_L6_rec, i_lgn, I_L4, I_LGN_rec, I_L4_rec,
                                                               n_e_lgn, n_e_l6, n_e_l4, n_trn, n_i_l6)

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
