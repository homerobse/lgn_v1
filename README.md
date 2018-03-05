Models for oscillations in visual thalamocortical pathway

# Parameters in this model

Size of the network, input, connectivity and simulation parameters included in s_tests.py
##  Synapse parameters:
In cells.py there are these (and others)
```python
synE = h.Exp2Syn(0.5, sec=self.soma)
synE.tau1 = 1
synE.tau2 = 3
synE_CT = h.Exp2Syn(0.5, sec=self.soma)
synE_CT.tau1 = 5               #  this was made to be slower than the synE - but artificial values were used because the synE time constants couldn't get smaller than 1
synE_CT.tau2 = 15              #  TODO check this with Li, Guido and Bickford 2003
synI = h.Exp2Syn(0.5, sec=self.soma)
synI.e = -100
synI.tau1 = 4
synI.tau2 = 2
```

## Hodgkin & Huxley parameters
Look at hh_original.mod

# Building NEURON mechanisms and running the simulation
Run `$nrnivmodl hh_original.mod netstim_pois.mod` to build these NEURON mechanisms.

Afterwards, to run an example simulation, run `$nrngui -python s_tests.py`
