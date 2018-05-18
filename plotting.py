import operator
import matplotlib.pyplot as plt
from matplotlib.mlab import psd
import numpy as np
from neuron import h
from numpy.ma import mean


def plot_all(conn, output_dir, fname, n_sim, dt, timeaxis, stim_rec, with_V1_L4, with_V1_L6, with_TRN,
             e_lgn_rec, i_lgn_rec, trn_rec, e_l4_rec, i_l4_rec, e_l6_rec, i_l6_rec,
             n_e_lgn, n_i_lgn, n_trn, n_e_l4, n_i_l4, n_e_l6, n_i_l6):

    #  FFT resolution = Fs/nfft = 1.22 Hz, waveform frequency resolution = 1/500ms = 2 Hz
    Fs = 1e3/dt
    nfft = int(round(timeaxis[-1]/dt)+1)
    zero_padding = int(1<<(nfft-1).bit_length())  # set zero padding to be the first power of 2 bigger than nfft

    # LGN LFP
    mean_lgn = np.mean(e_lgn_rec, axis=0)
    mean_lgn = mean_lgn - np.mean(mean_lgn)
    plt.figure(1)
    a = plt.subplot(411)
    plt.plot(timeaxis, mean_lgn)
    plt.title('Average membrane potential of excitatory cells in LGN')

    a = plt.subplot(412, sharex=a)
    for neuron_i in range(0, 4):
        plt.plot(timeaxis, e_lgn_rec[neuron_i])
        plt.title('Receiving inputs')

    a = plt.subplot(413)
    for neuron_i in range(4, n_e_lgn):
        plt.plot(timeaxis, e_lgn_rec[neuron_i])
        plt.title('Not receiving inputs')


    plt.subplot(414, sharex=a)
    for neuron_i in range(n_i_lgn):
        plt.plot(timeaxis, i_lgn_rec[neuron_i])
    plt.ylim([-100, 50])
    plt.xlim([0, h.tstop])
    plt.title('inhibitory cells in LGN')

    (Pxx, freqpsd) = psd(mean_lgn, nfft, Fs, pad_to=zero_padding)  # args are signal, nfft, Fs
    plt.figure(8)
    plt.plot(freqpsd, Pxx)
    plt.xlim([0, 150])
    plt.title('PSD of LGN LFP')

    (lgn_index, lgn_mean_peak_power, lgn_peak_freq) = process_peak(Pxx, freqpsd)

    mean_trn = np.zeros(np.shape(np.mean(trn_rec, 0)))
    if with_TRN:
        # TRN LFP
        mean_trn = np.mean(trn_rec, axis=0)
        mean_trn = mean_trn - np.mean(mean_trn)

        (Pxx, freqpsd) = psd(mean_trn, nfft, Fs)  # args are signal, nfft, Fs
        (trn_index, trn_mean_peak_power, trn_peak_freq) = process_peak(Pxx, freqpsd)
        plt.figure(9)
        plt.plot(freqpsd, Pxx)
        plt.xlim([0, 150])
        plt.title('PSD of TRN LFP')

        # plot trn neurons

        plt.figure(3)
        b = plt.subplot(211)

        plt.plot(timeaxis, mean_trn)
        plt.title('Average membrane potential of TRN cells')

        a = plt.subplot(212)
        for neuron_i in range(0, n_trn):
            plt.plot(timeaxis, trn_rec[neuron_i])

        plt.ylim([-100, 50])
        plt.xlim([0, h.tstop])

    mean_v1_l4 = np.zeros(np.shape(np.mean(e_l4_rec, axis=0)))
    if with_V1_L4:
        # V1 L4 LFP
        mean_v1_l4 = np.mean(e_l4_rec, axis=0)
        mean_v1_l4 = mean_v1_l4 - np.mean(mean_v1_l4)
        (Pxx, freqpsd) = psd(mean_v1_l4, nfft, Fs)  # args are signal, nfft, Fs
        (l4_index, l4_mean_peak_power, l4_peak_freq) = process_peak(Pxx, freqpsd)
        plt.figure(6)
        plt.plot(freqpsd, Pxx)
        plt.xlim([0, 150])
        plt.title('PSD of V1 L4 LFP')

    mean_v1_l6 = np.zeros(np.shape(np.mean(e_l6_rec, 0)))
    if with_V1_L6:
        # V1 L6 LFP
        mean_v1_l6 = np.mean(e_l6_rec, axis=0)
        mean_v1_l6 = mean_v1_l6 - np.mean(mean_v1_l6)

        (Pxx, freqpsd) = psd(mean_v1_l6, nfft, Fs)  # args are signal, nfft, Fs
        (l6_index, l6_mean_peak_power, l6_peak_freq) = process_peak(Pxx, freqpsd)
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
        plt.xlim([0, h.tstop])
        plt.title('GABAergic L6 neurons')

    if with_V1_L4:
        plt.figure(2)
        b = plt.subplot(411)

        plt.plot(timeaxis, mean_v1_l4)
        plt.title('Average membrane potential of V1 (L4) cells')

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
        plt.title('GABAergic neurons in V1 (L4)')

        plt.xlim([0, h.tstop])

    plot_input(stim_rec, timeaxis)

    sweep_f = open(fname, 'a')
    sweep_f.write(str([[lgn_peak_freq, lgn_mean_peak_power], n_sim, conn]))
    sweep_f.write("\n")
    # # sweep_f.write("TRN\n")
    # # sweep_f.write([trn_peak_freq, trn_mean_peak_power])
    # # sweep_f.write("V1 L4\n")
    # # sweep_f.write([l4_peak_freq, l4_mean_peak_power])
    # # sweep_f.write("V1 L6\n")
    # # sweep_f.write([l6_peak_freq, l6_mean_peak_power])
    sweep_f.close()

    return mean_lgn, mean_trn, mean_v1_l4, mean_v1_l6


def plot_input(stim_rec, timeaxis):
    plt.figure(5)
    plt.plot(stim_rec, np.ones([int(stim_rec.size()), 1]), 'ob')
    plt.xlim([timeaxis[0], timeaxis[-1]])
    plt.ylim([0., 2.])
    plt.title('Poisson input #0 to LGN excitatory cell #0')


def process_peak(Pxx, freqpsd):
    """
    Calculates where is the peak (index and frequency) and how strong it is (power) for the Power Spectrum Density array
    :param Pxx: Power spectrum density of the signal
    :param freqpsd: Array with the corresponding frequencies of each point in the Pxx array
    :return:  (index, mean_peak_power, peak_freq)
    index: index of Pxx array where the peak is located
    mean_peak_power: mean between the peak power and two neighboring points # TODO: make this use Hz interval, instead of the neighboring points
    peak_freq: frequency of the peak of PSD
    """
    # index = detect_peaks(Pxx, mpd=4)
    index, value = max(enumerate(Pxx), key=operator.itemgetter(1))
    mean_peak_power = mean([Pxx[index - 1], Pxx[index], Pxx[index + 1]])  # using Fs = 1/2.5e-05 nfft = 2**15 and dt = 0.025 the difference in frequency between neighbour indexes is 1.25 Hz
    peak_freq = freqpsd[index]
    return (index, mean_peak_power, peak_freq)
