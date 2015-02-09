################################################################################
####### DIDN'T WORK - UNKNOWN PROBLEM TO CONNECT SOMA TO SOMA2

################################################################################

from neuron import h
import matplotlib.pyplot as plt
import numpy as np

h.load_file("nrngui.hoc")  # load standard run system

# soma = h.Section(cell='soma')
soma = h.Section(name='soma')
print h.cas().name(), soma.name()
# soma2 = h.Section(cell='asd')
soma2 = h.Section(name='soma2')
print h.cas().name(), soma2.name()

soma.insert('hh')
soma2.insert('hh')

syn = h.Exp2Syn(soma(0.5))
syn.tau1 = 1
syn.tau2 = 3

h.NetCon(soma(0.5)._ref_v, syn)

stim = h.IClamp(soma(0.5))
stim.delay = 1
stim.dur = 5
stim.amp = 20

stim2 = h.IClamp(soma2(0.5))
stim2.delay = 20
stim2.dur = 20
stim2.amp = 20

rec_cell = h.Vector()
rec_cell.record(soma(0.5)._ref_v)
rec_cell2 = h.Vector()
rec_cell2.record(soma2(0.5)._ref_v)
rec_stim = h.Vector()
rec_stim.record(stim._ref_i)
rec_stim2 = h.Vector()
rec_stim2.record(stim2._ref_i)


timeaxis = h.Vector()
timeaxis.record(h._ref_t)

print soma(0.5).v

h.tstop = 50
h.run()

print soma(0.5).v

plt.plot(timeaxis, rec_cell)
plt.plot(timeaxis, rec_cell2)
plt.plot(timeaxis, rec_stim)
plt.plot(timeaxis, rec_stim2)
plt.show()