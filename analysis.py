from matplotlib.pyplot import show, subplots
from numpy.ma import mean, std, arange

# fname = "../data_files/sweep_2018-03-13 22:42:30.233927.txt"
# fname = "../data_files/sweep_2018-03-14 00:18:53.949155.txt"
# fname = "../data_files/sweep_2018-03-14 00:21:41.664107.txt"
fname = "../data_files/sweep_2018-03-23 02:48:35.182757.txt"
n_sim = 20
n_sweeps = 15

with open(fname, 'r') as f:
    all = []
    for l in f.readlines():
        all.append(eval(l))

i_sweep = 0
mean_freqs = []
sem_freqs = []  # Standard error of the mean
mean_pows = []
sem_pows = []  # Standard error of the mean
while i_sweep < n_sweeps*n_sim:
    freqs = []
    pows = []
    for j_sim in xrange(n_sim):
        (peak_freq, mean_peak_power) = (all[i_sweep + j_sim][0][0], all[i_sweep + j_sim][0][1])
        freqs.append(peak_freq)
        pows.append(mean_peak_power)
    mean_freqs.append(mean(freqs))
    sem_freqs.append(std(freqs) / (n_sim**0.5))
    mean_pows.append(mean(pows))
    sem_pows.append(std(freqs) / (n_sim**0.5))
    i_sweep = i_sweep + n_sim

labels = ["|E-E            |", "|    I-I        |", "|        I-E    |", "|            E-I|", "|E-E I-I        |", "|E-E     I-E    |", "|E-E         E-I|", "|    I-I I-E    |", "|    I-I     E-I|", "|        I-E E-I|", "|E-E I-I I-E    |", "|E-E I-I     E-I|", "|E-E     I-E E-I|", "|    I-I I-E E-I|", "|E-E I-I I-E E-I|"]
fig, axs = subplots(nrows=2, ncols=1, sharex=True)
ax = axs[0]
xs = arange(n_sweeps)
ax.errorbar(xs, mean_pows, yerr=sem_pows, fmt='o', capsize=4)
ax.set_xticks(xs)
ax.set_xticklabels(labels)
ax.set_title('Power')
ax = axs[1]
ax.errorbar(xs, mean_freqs, yerr=sem_freqs, fmt='o', capsize=4)
ax.set_xticks(xs)
ax.set_xticklabels(labels)
ax.set_title('Frequencies')
show()
