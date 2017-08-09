from neuron import h
import matplotlib.pyplot as plt
import numpy as np

h.load_file("nrngui.hoc")  # load standard run system

soma = h.Section(name='soma')
print h.cas().name(), soma.name()
soma2 = h.Section(name='soma2')
print h.cas().name(), soma2.name()

soma.insert('hh_original')
soma2.insert('hh_original')

syn = h.Exp2Syn(0.5, sec=soma2)
syn.e = 0
syn.tau1 = 1
syn.tau2 = 3

net_con = h.NetCon(soma(0.5)._ref_v, syn, 0, 5, 0.07)

stim = h.IClamp(soma(0.5))
stim.delay = 1
stim.dur = 5
stim.amp = 20

rec_cell = h.Vector()
rec_cell.record(soma(0.5)._ref_v)
rec_cell2 = h.Vector()
rec_cell2.record(soma2(0.5)._ref_v)
rec_syn = h.Vector()
rec_syn.record(syn._ref_i)
rec_stim = h.Vector()
rec_stim.record(stim._ref_i)

timeaxis = h.Vector()
timeaxis.record(h._ref_t)

h.tstop = 80
h.run()

f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(timeaxis, rec_cell)
axarr[0].plot(timeaxis, rec_cell2)
axarr[0].plot(timeaxis, rec_stim)
axarr[0].set_title("IClamp stimulus and voltage traces")
axarr[0].set_ylabel('Voltage (mV)')
axarr[1].plot(timeaxis, rec_syn)
axarr[1].set_title("Synaptic current")
axarr[1].set_ylabel("Current (mA)")
plt.xlabel("Time (ms)")
plt.show()
