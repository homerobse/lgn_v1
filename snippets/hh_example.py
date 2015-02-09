from neuron import h
import matplotlib.pyplot as plt

h.load_file("nrngui.hoc")  # load standard run system

soma = h.Section()

soma.insert('hh')

stim = h.IClamp(soma(0.5))

stim.delay = 1
stim.dur = 7
stim.amp = 20

rec_cell = h.Vector()
rec_cell.record(soma(0.5)._ref_v)

rec_stim = h.Vector()
rec_stim.record(stim._ref_i)

timeaxis = h.Vector()
timeaxis.record(h._ref_t)

print soma(0.5).v

h.tstop = 10

h.run()

print soma(0.5).v

plt.plot(timeaxis, rec_cell)
plt.plot(timeaxis, rec_stim)
plt.show()