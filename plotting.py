import matplotlib.pyplot as plt
from matplotlib.mlab import psd
import numpy as np
from neuron import h


def plot_all(timeaxis, stim_rec, with_V1_L4, with_V1_L6, with_TRN,
             e_lgn_rec, i_lgn_rec, trn_rec, e_l4_rec, i_l4_rec, e_l6_rec,
             i_l6_rec,
             n_e_lgn, n_i_lgn, n_trn, n_e_l4, n_i_l4, n_e_l6, n_i_l6):

    Fs = 1/2.5e-05

    # LGN LFP
    mean_lgn = np.mean(e_lgn_rec, axis=0)
    mean_lgn = mean_lgn - np.mean(mean_lgn)

    (Pxx, freqpsd) = psd(mean_lgn, 20000/2, Fs)  # args are signal, nfft, Fs
    plt.figure(8)
    plt.plot(freqpsd, Pxx)
    plt.xlim([0, 150])
    plt.title('PSD of LGN LFP')

    mean_trn = np.zeros(np.shape(np.mean(trn_rec, 0)))
    if with_TRN:
        # TRN LFP
        mean_trn = np.mean(trn_rec, axis=0)
        mean_trn = mean_trn - np.mean(mean_trn)

        (Pxx, freqpsd) = psd(mean_trn, 20000/2, Fs)  # args are signal, nfft, Fs
        plt.figure(9)
        plt.plot(freqpsd, Pxx)
        plt.xlim([0, 150])
        plt.title('PSD of TRN LFP')

        # plot trn neurons

        plt.figure(3)
        b = plt.subplot(211)

        plt.plot(timeaxis, mean_trn)
        plt.title('Average membrane potential of TRN cells' )

        a = plt.subplot(212)
        for neuron_i in range(0, n_trn):
            plt.plot(timeaxis, trn_rec[neuron_i])
            plt.title('2-rest inputs net 2')

        plt.ylim([-100, 50])
        plt.title('GABAergic TRN neurons')

        plt.xlim([0, h.tstop])

    mean_v1_l4 = np.zeros(np.shape(np.mean(e_l4_rec, axis=0)))
    if with_V1_L4:
        # V1 L4 LFP
        mean_v1_l4 = np.mean(e_l4_rec, axis=0)
        mean_v1_l4 = mean_v1_l4 - np.mean(mean_v1_l4)
        (Pxx, freqpsd) = psd(mean_v1_l4, 20000 / 2, Fs)  # args are signal, nfft, Fs
        import pdb
        pdb.set_trace()
        plt.figure(6)
        plt.plot(freqpsd, Pxx)
        plt.xlim([0, 150])
        plt.title('PSD of V1 L4 LFP')

    mean_v1_l6 = np.zeros(np.shape(np.mean(e_l6_rec, 0)))
    if with_V1_L6:
        # V1 L6 LFP
        mean_v1_l6 = np.mean(e_l6_rec, axis=0)
        mean_v1_l6 = mean_v1_l6 - np.mean(mean_v1_l6)

        (Pxx, freqpsd) = psd(mean_v1_l6, 20000/2, Fs)  # args are signal, nfft, Fs
        plt.figure(10)
        plt.plot(freqpsd, Pxx)
        plt.xlim([0, 150])
        plt.title('PSD of V1 L6 LFP')

        # plot L6 neurons

        plt.figure(4)
        a = plt.subplot(311)

        plt.plot(timeaxis, mean_v1_l6)
        plt.title('Average membrane potential of V1 (L6) cells')

        a = plt.subplot(312)
        for neuron_i in range(2, n_e_l6):
            plt.plot(timeaxis, e_l6_rec[neuron_i])
        plt.title('Glut L6 neurons')

        a = plt.subplot(313)
        for neuron_i in range(0, n_i_l6):
            plt.plot(timeaxis, i_l6_rec[neuron_i])

        plt.ylim([-100, 50])
        plt.title('GABAergic L6 neurons')
        plt.xlim([0, h.tstop])

    plt.figure(1)
    a = plt.subplot(411)

    plt.plot(timeaxis, mean_lgn)
    plt.title('Average membrane potential of excitatory cells in LGN')

    a = plt.subplot(413)
    for neuron_i in range(4, n_e_lgn):
        plt.plot(timeaxis, e_lgn_rec[neuron_i])
        plt.title('Not receiving inputs')

    a = plt.subplot(412, sharex=a)
    for neuron_i in range(0, 4):
        plt.plot(timeaxis, e_lgn_rec[neuron_i])
        plt.title('Receiving inputs')

    plt.subplot(414, sharex=a)

    for neuron_i in range(n_i_lgn):
        plt.plot(timeaxis, i_lgn_rec[neuron_i])
    plt.ylim([-100, 50])
    plt.title('inhibitory cells in LGN')

    plt.xlim([0, h.tstop])

    if with_V1_L4:
        plt.figure(2)
        b = plt.subplot(411)

        plt.plot(timeaxis, mean_v1_l4)
        plt.title('average membrane potential of V1 (L4) cells')

        a = plt.subplot(413)
        for neuron_i in range(2, n_e_l4):
            plt.plot(timeaxis, e_l4_rec[neuron_i])
            plt.title('2-rest inputs net L4')

        a = plt.subplot(412, sharex=a)
        for neuron_i in range(0, 2):
            plt.plot(timeaxis, e_l4_rec[neuron_i])
            plt.title('0-2 inputs net L4')

        plt.subplot(414, sharex=a)

        for neuron_i in range(n_i_l4):
            plt.plot(timeaxis, i_l4_rec[neuron_i])
        plt.ylim([-100, 50])
        plt.title('GABAergic net in V1 (L4)')

        plt.xlim([0, h.tstop])

    # plot input
    plt.figure(5)
    plt.plot(stim_rec, np.ones([int(stim_rec.size()), 1]), 'ob')
    plt.xlim([timeaxis[0], timeaxis[-1]])
    plt.ylim([0., 2.])
    plt.title('Poisson input #0 to LGN Cell #0')

    # all
    return mean_lgn, mean_trn, mean_v1_l4, mean_v1_l6
