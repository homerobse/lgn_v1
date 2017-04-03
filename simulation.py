from neuron import h
import numpy as np

from plotting import plot_all
from utils import e_net_connect, i_net_connect, e_ct_net_connect, e_net_connect_delay_dist, e_ct_net_connect_delay_dist
from utils import createNetwork, createNetworkL6

h.load_file("nrngui.hoc")  # load standard run system
# h.dt = 1
# global dt
# dt = h.dt
# global nsamples
# nsamples = h.tstop*h.dt
# h.v_init = -67


def simulate(n_runs, total_time, with_V1_L4, with_V1_L6, with_TRN, input, con_input_lgn,
             n_e_lgn, n_i_lgn, n_e_l6, n_i_l6, n_e_l4, n_i_l4, n_trn,
             delay_distbtn_E_L6_LGN, delay_E_L4_E_LGN, delay_E_LGN_I_L4, delay_E_LGN_E_L4, delay_E_LGN_E_L6,
             delay_E_LGN_TRN, delay_E_L4_TRN, delay_distbtn_E_L6_TRN, delay_E_LGN_I_LGN, delay_I_LGN_E_LGN, delay_E_LGN_I_L6,
             lgn_params, l4_params, l6_params, W_TRN_TRN, W_E_LGN_TRN, W_TRN_E_LGN, w_e_l6_trn, W_E_L4_E_L6,
             W_E_LGN_E_L4, W_E_L4_E_LGN, w_e_l6_e_lgn, W_E_LGN_E_L6, W_E_LGN_I_L6, W_E_LGN_I_L4, W_E_L4_TRN,
             connect_E_LGN_E_L4, connect_E_LGN_I_L4, connect_E_L4_E_LGN, connect_E_LGN_I_L6, connect_E_LGN_E_L6, connect_E_L6_E_LGN, connect_E_L4_TRN, connect_E_L6_TRN,
             connect_E_LGN_TRN, connect_TRN_E_LGN, connect_E_L4_E_L6):

    print "* * * Simulating %d runs * * *" % n_runs
    h.tstop = total_time
    for n_sim in range(n_runs):
        print "#%d: Constructing circuits..." % (n_sim + 1)

        # creating LGN network
        i_lgn, I_LGN_rec = createNetwork(n_i_lgn)
        e_lgn, E_LGN_rec = createNetwork(n_e_lgn)

        #create connections in network 1 (LGN)
        e_lgn_e_lgn_syn = e_net_connect(e_lgn, e_lgn, 0, 1, lgn_params['w_e_lgn_e_lgn'], 1)
        i_lgn_i_lgn_syn = i_net_connect(i_lgn, i_lgn, 0, 1, lgn_params['w_i_lgn_i_lgn'], 1)
        i_lgn_e_lgn_syn = i_net_connect(i_lgn, e_lgn, 0, delay_I_LGN_E_LGN, lgn_params['w_i_lgn_e_lgn'], 1)
        e_lgn_i_lgn_syn = e_net_connect(e_lgn, i_lgn, 0, delay_E_LGN_I_LGN, lgn_params['w_e_lgn_i_lgn'], 1)  # weight should be set to zero

        e_l4, E_L4_rec = createNetwork(n_e_l4)
        i_l4, I_L4_rec = createNetwork(n_i_l4)
        if with_V1_L4:
            #create connections in network 2  (V1 superficial)
            # using values different from Hauesler and Maass yet, in order to be able to generate gamma
            e_l4_e_l4_sin = e_net_connect(e_l4, e_l4, 0, 1, l4_params['w_e_l4_e_l4'], 0.3)  # 17% connectivity Hauesler and Maass 2007
            i_l4_i_l4_sin = i_net_connect(i_l4, i_l4, 0, 1, l4_params['w_i_l4_i_l4'], 0.7)  # 50% connectivity Hauesler and Maass 2007
            e_l4_i_l4_sin = e_net_connect(e_l4, i_l4, 0, 1, l4_params['w_e_l4_i_l4'], 0.3)  # 19% connectivity Hauesler and Maass 2007
            i_l4_e_l4_sin = i_net_connect(i_l4, e_l4, 0, 1, l4_params['w_i_l4_e_l4'], 0.7)  # 10% connectivity Hauesler and Maass 2007

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
                    for neuron_j in range(0, len(i_l4)*1/4):
                        sin_E_LGN_I_L4.append(h.NetCon(e_lgn[neuron_i].soma(0.5)._ref_v, i_l4[neuron_j].synE,
                                                       0, delay_E_LGN_I_L4, W_E_LGN_I_L4[neuron_i, neuron_j]))
                    h.pop_section()

                #            #topographic connectivity
                #            if connect_E_LGN_I_L4:
                #                # connections from Glutamatergic neurons of network (LGN) to network V1 L4
                #                sin_E_LGN_I_L4 = list()
                #                len_LGN = len(e_lgn)
                #                for neuron_i in range(0, len_LGN*1/4):
                #                    e_lgn[neuron_i].soma.push()
                #                    sin_E_LGN_I_L4.append(h.NetCon(e_lgn[neuron_i].soma(0.5)._ref_v, i_l4[neuron_i].synE,
                #                                                       0, delay_E_LGN_I_L4, W_E_LGN_I_L4[neuron_i, neuron_i]))
                #                    h.pop_section()

        i_l6, I_L6_rec = createNetworkL6(n_i_l6)
        e_l6, E_L6_rec = createNetworkL6(n_e_l6)
        if with_V1_L6:
            # create connections in network 2  (V1 L6)
            # using values different from Hauesler and Maass yet, in order to be able to generate gamma
            e_l6_e_l6_sin = e_net_connect(e_l6, e_l6, 0, 1, l6_params['w_e_l6_e_l6'], 0.3)  # 9% connectivity Hauesler and Maass 2007 (heuristic from L5)
            i_l6_i_l6_sin = i_net_connect(i_l6, i_l6, 0, 1, l6_params['w_i_l6_i_l6'], 0.7)  # 60% connectivity Hauesler and Maass 2007 (heuristic from L5)
            e_l6_i_l6_syn = e_net_connect(e_l6, i_l6, 0, 1, l6_params['w_e_l6_i_l6'], 0.3)  # 10% connectivity Hauesler and Maass 2007 (heuristic from L5)
            i_l6_e_l6_syn = i_net_connect(i_l6, e_l6, 0, 1, l6_params['w_i_l6_e_l6'], 0.7)  # 12% connectivity Hauesler and Maass 2007 (heuristic from L5)

            # connections from V1 input (L4) layer to L6
            if connect_E_L4_E_L6:
                e_l4_e_l6_sin = e_net_connect(e_l4, e_l6, 0, 1, W_E_L4_E_L6,1)

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
                e_lgn_e_l6_syn = e_net_connect(e_lgn, e_l6, 0, delay_E_LGN_E_L6, W_E_LGN_E_L6,1)

            # TODO: Connectivity as Hirsch
            if connect_E_LGN_I_L6:
                e_lgn_i_l6_syn = e_net_connect(e_lgn, i_l6, 0, delay_E_LGN_I_L6, W_E_LGN_I_L6,1)

        #create trn neurons (inhibitory only)
        trn, TRN_rec = createNetwork(n_trn)
        if with_TRN:
            trn_trn_syn = i_net_connect(trn, trn, 0, 1, W_TRN_TRN,1)

            #connections from Glutamatergic neurons of network V1 L4 to trn
            if with_V1_L4 and connect_E_L4_TRN:
                e_l4_trn_syn = e_net_connect(e_l4, trn, 0, delay_E_L4_TRN, W_E_L4_TRN,1)

            if with_V1_L6 and connect_E_L6_TRN:
                e_l6_trn_syn = e_net_connect_delay_dist(e_l6, trn, 0, delay_distbtn_E_L6_TRN, w_e_l6_trn,1)

            #ALL-to-ALL
            if connect_E_LGN_TRN:
                #connections from Glutamatergic neurons of network 1 (LGN) to trn
                e_lgn_trn_syn = e_net_connect(e_lgn, trn, 0, delay_E_LGN_TRN, W_E_LGN_TRN,1)

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
                trn_e_lgn_sin = i_net_connect(trn, e_lgn, 0, 1, W_TRN_E_LGN,1)

        # generate inputs to LGN
        netStim = list()
        i_stims = list()
        e_stims = list()
        stim_rec = h.Vector()

        for stim_i in range(0, input['nstims']):
            netStim.append(h.NetStimPois(input['input']))
            netStim[stim_i].start = 0
            netStim[stim_i].mean = input['stimrate']  # 100 = 10 Hz, 10 = 100 Hz, 1 = 1000Hz, 5 = 200 Hz, 6 = 150 Hz
            netStim[stim_i].number = 0
            if stim_i < n_i_lgn:
                i_stims.append(h.NetCon(netStim[stim_i], i_lgn[stim_i].synE, con_input_lgn['gaba_threshold'],
                                        con_input_lgn['gaba_delay'], con_input_lgn['gaba_weight']))
            e_stims.append(h.NetCon(netStim[stim_i], e_lgn[stim_i].synE, con_input_lgn['glut_threshold'],
                                    con_input_lgn['glut_delay'], con_input_lgn['glut_weight']))

        e_stims[0].record(stim_rec)  # measure poisson input #0 to LGN Cell #0

        timeaxis = h.Vector()
        timeaxis.record(h._ref_t)
        print "#%d: Running simulation..." % (n_sim + 1)
        h.run()

        meanLGN, meanTRN, meanV1input, meanV1output = plot_all(timeaxis, stim_rec, with_V1_L4, with_V1_L6, with_TRN,
                                                               E_L4_rec, TRN_rec, E_L6_rec, E_LGN_rec,
                                                               I_L6_rec, i_lgn, i_l4, I_LGN_rec, I_L4_rec,
                                                               n_e_lgn, n_e_l6, n_e_l4, n_trn, n_i_l6)

        ofname = "../data_files/" "sim" + str(n_sim+0) + ".txt"

        n = len(timeaxis)
        indx = np.arange(1, n, 40)

        w = np.array(meanLGN)
        u = np.array(meanTRN)
        x = np.array(meanV1input)
        y = np.array(meanV1output)
        z = np.array(timeaxis)
        np.savetxt(ofname, (w[indx], u[indx], x[indx], y[indx], z[indx]))

        print "Progress: %d runs simulated %d runs missing" % (n_sim + 1, n_runs - n_sim - 1)


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
