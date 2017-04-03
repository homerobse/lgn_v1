from neuron import h
import matplotlib.pyplot as plt
import numpy as np

h.load_file("nrngui.hoc")  # load standard run system

soma = h.Section(name='soma')
soma2 = h.Section(name='soma2')
soma3 = h.Section(name='soma3')
soma4 = h.Section(name='soma4')

soma.insert('hh')
soma2.insert('hh')
soma3.insert('hh')
soma4.insert('hh')

syn2 = h.ExpSyn(0.5, sec=soma2)
syn4 = h.ExpSyn(0.5, sec=soma4)
# syn.e = 0
# syn.tau1 = 1
# syn.tau2 = 3


def connect(v1, syn2, v3, syn4):
    l = list()
    print h.cas().name()
    l.append(h.NetCon(v1, syn2, 0, 0, 0.5))
    soma3.push()
    print h.cas().name()
    l.append(h.NetCon(v3, syn4, 0, 0, 1))
    h.pop_section()
    print h.cas().name()
    return l

# net_con = h.NetCon(soma(0.5)._ref_v, syn, 30, 0, 1)
l = connect(soma(0.5)._ref_v, syn2, soma3(0.5)._ref_v, syn4)

stim = h.IClamp(soma(0.5))
stim.delay = 1
stim.dur = 5
stim.amp = 20

stim2 = h.IClamp(soma3(0.5))
stim2.delay = 20
stim2.dur = 20
stim2.amp = 20

rec_cell = h.Vector()
rec_cell.record(soma(0.5)._ref_v)
rec_cell2 = h.Vector()
rec_cell2.record(soma2(0.5)._ref_v)
rec_cell3 = h.Vector()
rec_cell3.record(soma3(0.5)._ref_v)
rec_cell4 = h.Vector()
rec_cell4.record(soma4(0.5)._ref_v)
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
plt.plot(timeaxis, rec_cell3)
plt.plot(timeaxis, rec_cell4)
plt.plot(timeaxis, rec_stim)
plt.plot(timeaxis, rec_stim2)
plt.show()