from neuron import h
import numpy as np

from plotting import plot_all
from utils import e_net_connect, partial_e_net_connect, i_net_connect, e_ct_net_connect, e_net_connect_delay_dist, e_ct_net_connect_delay_dist
from utils import create_network, create_network_L6

# h.load_file("nrngui.hoc")  # load standard run system
# h.dt = 1
# global dt
# dt = h.dt
# global nsamples
# nsamples = h.tstop*h.dt
# h.v_init = -67


def simulate(OUTPUT_DIR, n_runs, total_time, temperature, with_v1_l4, with_v1_l6, with_trn, input, con_input_lgn,
             n_e_lgn, n_i_lgn, n_e_l6, n_i_l6, n_e_l4, n_i_l4, n_trn,
             delay_distbtn_e_l6_lgn, delay_e_l4_e_lgn, delay_e_lgn_i_l4, delay_e_lgn_e_l4, delay_e_lgn_e_l6,
             delay_e_lgn_trn, delay_e_l4_trn, delay_distbtn_e_l6_trn, delay_e_lgn_i_l6,
             lgn_params, l4_params, l6_params, trn_params, w_e_lgn_trn, w_trn_e_lgn, w_e_l6_trn, w_e_l4_e_l6,
             w_e_lgn_e_l4, w_e_l4_e_lgn, w_e_l6_e_lgn, w_e_lgn_e_l6, w_e_lgn_i_l6, w_e_lgn_i_l4, w_e_l4_trn,
             connect_e_lgn_e_l4, connect_e_lgn_i_l4, connect_e_l4_e_lgn, connect_e_lgn_i_l6, connect_e_lgn_e_l6, connect_e_l6_e_lgn, connect_e_l4_trn, connect_e_l6_trn,
             connect_e_lgn_trn, connect_trn_e_lgn, connect_e_l4_e_l6):
    h.celsius = temperature
    print "* * * Simulating %d runs * * *" % n_runs
    h.tstop = total_time
    for n_sim in range(n_runs):
        print "#%d: Constructing circuits..." % (n_sim + 1)

        # creating LGN network
        i_lgn, i_lgn_rec = create_network(n_i_lgn)
        e_lgn, e_lgn_rec = create_network(n_e_lgn)

        # create connections in network 1 (LGN)
        e_lgn_e_lgn_syn = e_net_connect(e_lgn, e_lgn, 0, 1, lgn_params['w_e_lgn_e_lgn'], 1)
        i_lgn_i_lgn_syn = i_net_connect(i_lgn, i_lgn, 0, 1, lgn_params['w_i_lgn_i_lgn'], 1)
        i_lgn_e_lgn_syn = i_net_connect(i_lgn, e_lgn, 0, lgn_params['delay_i_e'], lgn_params['w_i_lgn_e_lgn'], 1)
        e_lgn_i_lgn_syn = e_net_connect(e_lgn, i_lgn, 0, lgn_params['delay_e_i'], lgn_params['w_e_lgn_i_lgn'], 1)  # weight should be set to zero

        e_l4, e_l4_rec = create_network(n_e_l4)
        i_l4, i_l4_rec = create_network(n_i_l4)
        if with_v1_l4:
            # create connections in V1 L4
            e_l4_e_l4_sin = e_net_connect(e_l4, e_l4, 0, 1, l4_params['w_e_l4_e_l4'], l4_params['p_e_e'])
            i_l4_i_l4_sin = i_net_connect(i_l4, i_l4, 0, 1, l4_params['w_i_l4_i_l4'], l4_params['p_i_i'])
            e_l4_i_l4_sin = e_net_connect(e_l4, i_l4, 0, 1, l4_params['w_e_l4_i_l4'], l4_params['p_e_i'])
            i_l4_e_l4_sin = i_net_connect(i_l4, e_l4, 0, 1, l4_params['w_i_l4_e_l4'], l4_params['p_i_e'])


            # extrinsic connections

            # Population 1) 15 LGN E cells connect to 15 V1 L4 E cells
            # Population 2) 5 LGN E cells connect to 5 V1 L4 I cells
            #
            # Population 1 and population 2 are different
            #
            # Hirsch et al., 1998
            if connect_e_lgn_e_l4:
                # connections from Glutamatergic neurons of network LGN to network V1 L4
                partial_e_net_connect(e_lgn, e_l4, 1./4, 1, 1./4, 1, 0, delay_e_lgn_e_l4, w_e_lgn_e_l4)
                # topographically_connect(e_lgn, e_l4, 0, 1, 0, delay_e_lgn_e_l4, w_e_lgn_e_l4)

            if connect_e_l4_e_lgn:
                # TODO: feedback connections are only of 3/4 of neurons?
                # connections from Glutamatergic neurons of network 2 (V1) to network 1 (LGN)
                partial_e_net_connect(e_l4, e_lgn, 1./4, 1, 1./4, 1, 0, delay_e_l4_e_lgn, w_e_l4_e_lgn)

            if connect_e_lgn_i_l4:
                # connections from Glutamatergic neurons of network (LGN) to network V1 L4
                partial_e_net_connect(e_lgn, i_l4, 0, 1./4, 0, 1./4, 0, delay_e_lgn_i_l4, w_e_lgn_i_l4)
                # topographically_connect(e_lgn, i_l4, 0, 1./4, 0, delay_e_lgn_i_l4, w_e_lgn_i_l4)

        i_l6, i_l6_rec = create_network_L6(n_i_l6)
        e_l6, e_l6_rec = create_network_L6(n_e_l6)
        if with_v1_l6:
            # create connections in V1 L6
            e_l6_e_l6_sin = e_net_connect(e_l6, e_l6, 0, 1, l6_params['w_e_l6_e_l6'], l6_params['p_e_e'])
            i_l6_i_l6_sin = i_net_connect(i_l6, i_l6, 0, 1, l6_params['w_i_l6_i_l6'], l6_params['p_i_i'])
            e_l6_i_l6_syn = e_net_connect(e_l6, i_l6, 0, 1, l6_params['w_e_l6_i_l6'], l6_params['p_e_i'])
            i_l6_e_l6_syn = i_net_connect(i_l6, e_l6, 0, 1, l6_params['w_i_l6_e_l6'], l6_params['p_i_e'])

            # connections from V1 input (L4) layer to L6
            if connect_e_l4_e_l6:
                e_l4_e_l6_sin = e_net_connect(e_l4, e_l6, 0, 1, w_e_l4_e_l6, 1)

# ALL-to-ALL connections of feedback
            if connect_e_l6_e_lgn:
                e_l6_e_lgn_sin = e_ct_net_connect_delay_dist(e_l6, e_lgn, 0, delay_distbtn_e_l6_lgn, w_e_l6_e_lgn)

            # TODO: Connectivity as Hirsch
            if connect_e_lgn_e_l6:
                e_lgn_e_l6_syn = e_net_connect(e_lgn, e_l6, 0, delay_e_lgn_e_l6, w_e_lgn_e_l6, 1)

            # TODO: Connectivity as Hirsch
            if connect_e_lgn_i_l6:
                e_lgn_i_l6_syn = e_net_connect(e_lgn, i_l6, 0, delay_e_lgn_i_l6, w_e_lgn_i_l6, 1)

        # create trn neurons (inhibitory only)
        trn, trn_rec = create_network(n_trn)
        if with_trn:
            trn_trn_syn = i_net_connect(trn, trn, 0, trn_params['delay_i_i'], trn_params['w_trn_trn'],
                                        trn_params['p_i_i'])

            # connections from Glutamatergic neurons of network V1 L4 to trn
            if with_v1_l4 and connect_e_l4_trn:
                e_l4_trn_syn = e_net_connect(e_l4, trn, 0, delay_e_l4_trn, w_e_l4_trn, 1)

            if with_v1_l6 and connect_e_l6_trn:
                e_l6_trn_syn = e_net_connect_delay_dist(e_l6, trn, 0, delay_distbtn_e_l6_trn, w_e_l6_trn, 1)

            if connect_e_lgn_trn:
                # connections from Glutamatergic neurons of LGN to TRN
                # ALL-to-ALL
                e_lgn_trn_syn = e_net_connect(e_lgn, trn, 0, delay_e_lgn_trn, w_e_lgn_trn, 1)

                # # topographic
                # topographically_connect(e_lgn, trn, 0, 1, 0, delay_e_lgn_trn, w_e_lgn_trn)

            if connect_trn_e_lgn:
                # ALL-to-ALL
                trn_e_lgn_sin = i_net_connect(trn, e_lgn, 0, 1, w_trn_e_lgn, 1)

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

        mean_lgn, mean_trn, mean_v1_l4, mean_v1_l6 = plot_all(timeaxis, stim_rec, with_v1_l4, with_v1_l6, with_trn,
                                                               e_lgn_rec, i_lgn_rec, trn_rec, e_l4_rec, i_l4_rec, e_l6_rec, i_l6_rec,
                                                               n_e_lgn, n_i_lgn, n_trn, n_e_l4, n_i_l4, n_e_l6, n_i_l6)

        ofname = OUTPUT_DIR + "/sim" + str(n_sim+0) + ".txt"

        n = len(timeaxis)
        indx = np.arange(1, n, 40)  # store one in every 40 values

        w = np.array(mean_lgn)
        u = np.array(mean_trn)
        x = np.array(mean_v1_l4)
        y = np.array(mean_v1_l6)
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
