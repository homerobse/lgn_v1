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
Look at hhvitor2.mod

# Draft of written work

## Introduction
Explain gamma mechanisms and give orverview about Cortico-thalamic circuits

## Important things to be discussed in the paper:
* Show that if LGN has E-I and I-I connection it entrains in gamma
* Show that if we give oscillatory input it entrains LGN and cortex
* Show that for weak drive from the retina, fast time constants and strong feedback LGN entrains gamma
* Show that LGN and TRN oscillate if they have fast delays

## Conclusion:
1. Mechanistic explanation for the generation of gamma in the cortex and its absence in LGN
    1. Lack of E to I and I to I recurrence in the LGN
2. Cortical gamma does not entrain geniculate
    1. Corticogeniculate feedback synapses have a prolonged time constant
    2. Feedback is weak relative to feedforward drive from the retina