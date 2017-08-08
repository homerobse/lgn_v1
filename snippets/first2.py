# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 16:53:48 2014

@author: vitor
"""

try:
    import sys
    import DLFCN
    sys.setdlopenflags(DLFCN.RTLD_NOW | DLFCN.RTLD_GLOBAL)
except:
    pass

from neuron import h
import numpy as np
import matplotlib.pyplot as plt

h.load_file("nrngui.hoc")  # load standard run system

h.dt = 1
global dt
dt = h.dt
h.tstop = 5000
global nsamples
nsamples = h.tstop*h.dt
h.v_init = -67


class Pyrcell:
    "Pyramidal cell"

    def __init__(self):
       
        self.soma = h.Section(name='soma',cell=self)
        self.soma.L=30
        self.soma.nseg=1
        self.soma.diam=1

        self.soma.Ra=100
        self.soma.cm=1

        #self.soma.insert("hhvitor2")
        self.soma.insert("hhvitor2")
        
        self.soma.gna_hhvitor2=0.039
        self.soma.gkd_hhvitor2=0.006
        self.soma.gl_hhvitor2=0
        self.soma.gm_hhvitor2=0
        self.soma.gt_hhvitor2=0
        self.soma.gleak_hhvitor2 = 0.0000273
   
        self.synE = h.Exp2Syn(0.5,sec=self.soma)
        self.synE.tau1=1
        self.synE.tau2=3

        self.synI = h.Exp2Syn(0.5,sec=self.soma)
        self.synI.e=-100
        self.synI.tau1=4
        self.synI.tau2=2
        
        self.stm = h.IClamp(0.5,sec=self.soma)
        self.stm.amp=0
        self.stm.dur=1000
        self.stm.delay=100

def createnetwork(Nneurons=4):

    network=h.List()
    network_rec=h.List()
 
    for Nindex in range(Nneurons):
        network.append(Pyrcell())
        network_rec.append(h.Vector())
        network_rec[Nindex].record(network[Nindex].soma(0.5)._ref_v)
        
    return network,network_rec     
    
NGABAn = 1
GABAneurons,GABAneurons_rec = createnetwork(NGABAn)
    
Nneurons = 50
Glutneurons,Glutneurons_rec = createnetwork(Nneurons)
Glutneurons_W = np.random.exponential(1,Nneurons*Nneurons)/100000.
Glutneurons_W = Glutneurons_W.reshape((Nneurons,Nneurons))
Glutneurons_W = Glutneurons_W - np.diag(np.diag(Glutneurons_W))

GlutGlut_sin = list()
for neuron_i in range(len(Glutneurons)):
    
    Glutneurons[neuron_i].soma.push()
    
    for neuron_j in range(len(Glutneurons)):
        
        GlutGlut_sin.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v,Glutneurons[neuron_j].synE,0,5,Glutneurons_W[neuron_i,neuron_j]))

GABAGlut_sin = list()
GlutGABA_sin = list()
for Gneuron_i in range(len(GABAneurons)):
    for neuron_i in range(len(Glutneurons)):
        GABAneurons[Gneuron_i].soma.push()
        GABAGlut_sin.append(h.NetCon(GABAneurons[Gneuron_i].soma(0.5)._ref_v,Glutneurons[neuron_i].synI,0,7,3./10000))
        
        Glutneurons[neuron_i].soma.push()
        GlutGABA_sin.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v,GABAneurons[Gneuron_i].synE,0,5,1./50000))

nstims = 1
stimrate = 8
netStim = list()

for stim_i in range(0,nstims):
    netStim.append(h.NetStimPois(.5));
    netStim[stim_i].start=0;
    netStim[stim_i].mean = stimrate;
    netStim[stim_i].number=0;

stim = h.NetCon(netStim[0],Glutneurons[0].synE,0.5,0,1./(150000))
stim_rec = h.Vector()
stim.record(stim_rec)
stim2 = h.NetCon(netStim[0],Glutneurons[1].synE,0.5,0,1.1/(150000))
stim2 = h.NetCon(netStim[0],Glutneurons[2].synE,0.5,0,1.2/(150000))

timeaxis = h.Vector()
timeaxis.record(h._ref_t)

h.run() 

a = plt.subplot(312)
for neuron_i in range(3,Nneurons):
    plt.plot(timeaxis,Glutneurons_rec[neuron_i])

a = plt.subplot(311,sharex=a)
for neuron_i in range(3):
    plt.plot(timeaxis,Glutneurons_rec[neuron_i])    

plt.subplot(313,sharex=a)

for neuron_i in range(len(GABAneurons)):
    plt.plot(timeaxis,GABAneurons_rec[neuron_i])
plt.ylim([-100,50])

plt.xlim([0,h.tstop])

def detect_spikes(voltagesignals):
    nneurons = len(voltagesignals)
    spiketimes = list()
    for neuron_i in range(nneurons):
        auxsignal = np.asarray(voltagesignals[neuron_i])
        thrscross_ind = [i for (i, val) in enumerate(auxsignal) if val > 30]

        spikeidxAUX = list()
        aux = 0
        for thrscross_i in range(len(thrscross_ind)):
            if thrscross_ind[thrscross_i]>aux:
                auxspan = range(thrscross_ind[thrscross_i],thrscross_ind[thrscross_i]+15)
                auxspanvalues = list(auxsignal[auxspan])
                spikeidxAUX.append(auxspanvalues.index(max(auxspanvalues))+thrscross_ind[thrscross_i])
                aux = thrscross_ind[thrscross_i]+15
        
        spiketimes.append(spikeidxAUX)
                
    return spiketimes

def generateraster(spiketimes):
    spktimes = np.array([])
    neuronlabel = np.array([])
    for neuron_i in range(len(spiketimes)):
        spktimes = np.concatenate([spktimes,np.asarray(spiketimes[neuron_i])],0)
        neuronlabel = np.concatenate([neuronlabel,neuron_i*np.ones(len(spiketimes[neuron_i]))],0)

    return spktimes,neuronlabel
              
GABAspiketimes = detect_spikes(GABAneurons_rec)
GLUTspiketimes = detect_spikes(Glutneurons_rec)

spktimes,neuronlabel = generateraster(GLUTspiketimes)

timevect = np.array(timeaxis)
dttimevec = timevect[1]-timevect[0]

plt.figure(5)
plt.subplot(1,1,1,sharex=a)
plt.plot(spktimes*dttimevec,neuronlabel,'.k')


    
plt.show()

