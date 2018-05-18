import os
import pdb
import numpy as np
import linecache
import matplotlib.pyplot as plt
from matplotlib.mlab import psd


def extract_lfp_and_time(sim_file):
    lgn_line = linecache.getline(sim_file, 1)  # first line
    lgn_lfp = np.array(map(float, lgn_line.split()))
    time_line = linecache.getline(sim_file, 5)  # fifth line
    time = np.array(map(float, time_line.split()))
    return time, lgn_lfp

def choose_plot_position(idx, c, connectivity):
    if idx<4: j = idx; i=0;
    elif idx<10: j = idx-4; i=1;
    elif idx<14: j = idx-10; i=2;
    else: j = 5; i=2;
    return i, j

sim_time = 500e-3
path = '/home/homerobse/lgn-v1/data_files'
# fname_dateprefix = '2018-04-02 00:39:27.361073_'
# fname_dateprefix = '2018-04-03 00:10:26.703816_'
# fname_dateprefix = '2018-04-03 01:05:34.703962_'
# fname_dateprefix = '2018-04-03 01:31:43.440352_'  # just 4 inputs
# fname_dateprefix = '2018-04-03 01:49:19.531271_'  # 40 inputs
fname_dateprefix = '2018-04-25 09:25:09.183849_'

weight = 4./100000
connectivity = [
                    [weight, 0, 0, 0],                  # E-E
                    [0, weight, 0, 0],                  #     I-I
                    [0, 0, weight, 0],                  #         I-E
                    [0, 0, 0, weight],                  #             E-I
                    [weight, weight, 0, 0],             # E-E I-I
                    [weight, 0, weight, 0],             # E-E     I-E
                    [weight, 0, 0, weight],             # E-E         E-I
                    [0, weight, weight, 0],             #     I-I I-E
                    [0, weight, 0, weight],             #     I-I     E-I
                    [0, 0, weight, weight],             #         I-E E-I
                    [weight, weight, weight, 0],        # E-E I-I I-E
                    [weight, weight, 0, weight],        # E-E I-I     E-I
                    [weight, 0, weight, weight],        # E-E     I-E E-I
                    [0, weight, weight, weight],        #     I-I I-E E-I
                    [weight, weight, weight, weight]    # E-E I-I I-E E-I
                ]

nruns = 10
dt = 0.025e-3
Fs = 1/dt/40  # only stored every 40th data point
nsamples = int(round(sim_time/dt/40))+1  # only stored every 40th data point
nfft = nsamples
zero_padding = int(1<<(nfft-1).bit_length())  # set zero padding to be the first power of 2 bigger than nfft

fig_lfp, axs_lfp = plt.subplots(nrows=3, ncols=6)
fig_psd, axs_psd = plt.subplots(nrows=3, ncols=6)
for idx, c in enumerate(connectivity):

    # get lfps, psds and time array
    lgn_pows = np.empty((nruns, zero_padding / 2 + 1))
    lgn_lfps = np.empty((nruns, nsamples))
    for nsim in np.arange(nruns):
        sim_file = os.path.join(path, fname_dateprefix + str(c) + '_sim-' + str(nsim) + '.txt')
        time, lgn_lfp = extract_lfp_and_time(sim_file)
        (lgn_pows[nsim, :], freqpsd) = psd(lgn_lfp, nfft, Fs, pad_to=zero_padding)  # args are signal, nfft, Fs
        # pdb.set_trace()
        lgn_lfps[nsim, :] = lgn_lfp

    i, j = choose_plot_position(idx, c, connectivity)

    # plot LFP
    ax = axs_lfp[i, j]
    mean_lfp = np.mean(lgn_lfps, axis=0)
    sem_lfp = np.std(lgn_lfps, axis=0) / nruns ** 0.5  # standard error of the mean
    plt.figure(1)
    ax.plot(time, mean_lfp)
    ax.fill_between(time, mean_lfp - sem_lfp, mean_lfp + sem_lfp, alpha=0.5, edgecolor='k', linestyle='dashdot')
    ax.set_title(str(c))
    fig_lfp.suptitle('LFP: E-E, I-I, I-E, E-I')

    axs_lfp[0, 4].set_visible(False)
    axs_lfp[0, 5].set_visible(False)
    axs_lfp[2, 4].set_visible(False)

    # plot PSD
    ax = axs_psd[i, j]
    mean_pxx = np.mean(lgn_pows, axis=0)
    sem_pxx = np.std(lgn_pows, axis=0) / nruns ** 0.5  # standard error of the mean
    ax.plot(freqpsd, mean_pxx, 'k-')
    ax.fill_between(freqpsd, mean_pxx - sem_pxx, mean_pxx + sem_pxx, alpha=0.5, edgecolor='k', linestyle='dashdot')
    ax.set_xlim([0, 150])
    ax.set_title(str(c))
    fig_psd.suptitle('PSD: E-E, I-I, I-E, E-I')

    axs_psd[0, 4].set_visible(False)
    axs_psd[0, 5].set_visible(False)
    axs_psd[2, 4].set_visible(False)

plt.show()
