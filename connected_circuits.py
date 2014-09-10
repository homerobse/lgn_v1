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
from matplotlib.mlab import psd

h.load_file("nrngui.hoc") # load standard run system

h.dt = 1
global dt
dt = h.dt
h.tstop = 500
global nsamples
nsamples = h.tstop*h.dt
h.v_init = -67

print "nsamples"
print nsamples

with_V1_L4 = True
with_TRN = False

class Pyrcell:
    "Pyramidal cell"

    def __init__(self):
       
        self.soma = h.Section(name='soma',cell=self)
        self.soma.L=30
        self.soma.nseg=1
        self.soma.diam=1

        self.soma.Ra = 100
        self.soma.cm = 1

        #self.soma.insert("hh")
        self.soma.insert("hhvitor2")
        
        self.soma.gna_hhvitor2=0.039
        self.soma.gkd_hhvitor2=0.006
        self.soma.gl_hhvitor2=0
        self.soma.gm_hhvitor2=0
        self.soma.gt_hhvitor2=0
        self.soma.gleak_hhvitor2 = 0.0000273
   
        self.synE = h.Exp2Syn(0.5, sec = self.soma)
        self.synE.tau1 = 1
        self.synE.tau2 = 3

        self.synI = h.Exp2Syn(0.5,sec=self.soma)
        self.synI.e=-100
        
        self.synI.tau1=1
        self.synI.tau2=2        
        
        self.stm = h.IClamp(0.5,sec=self.soma)
        self.stm.amp=0
        self.stm.dur=1000
        self.stm.delay=100
        

def createnetwork(Nneurons=4):

    network=h.List()
    network_rec=h.List()
 
    for Nindex in range(Nneurons):
        p = Pyrcell()
        network.append(p)
        network_rec.append(h.Vector())
        network_rec[Nindex].record(network[Nindex].soma(0.5)._ref_v)
        
    return network,network_rec     
    
NGABAn = 6
GABAneurons,GABAneurons_rec = createnetwork(NGABAn)
GABAneurons_W = np.random.exponential(1,NGABAn*NGABAn)*4/100000.
GABAneurons_W = GABAneurons_W.reshape((NGABAn,NGABAn))
GABAneurons_W = GABAneurons_W - np.diag(np.diag(GABAneurons_W))
    
Nneurons = 20
Glutneurons,Glutneurons_rec = createnetwork(Nneurons)
Glutneurons_W = np.random.exponential(1,Nneurons*Nneurons)*4/100000.
Glutneurons_W = Glutneurons_W.reshape((Nneurons,Nneurons))
Glutneurons_W = Glutneurons_W - np.diag(np.diag(Glutneurons_W))

#create connections in network 1 (LGN)
GlutGlut_sin = list()
for neuron_i in range(len(Glutneurons)):
    Glutneurons[neuron_i].soma.push()
    for neuron_j in range(len(Glutneurons)):
        GlutGlut_sin.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v,Glutneurons[neuron_j].synE,0,1,Glutneurons_W[neuron_i,neuron_j]))
    h.pop_section()
    

GABAGABA_sin = list()
for neuron_i in range(len(GABAneurons)):
    GABAneurons[neuron_i].soma.push()
    for neuron_j in range(len(GABAneurons)):
        GABAGABA_sin.append(h.NetCon(GABAneurons[neuron_i].soma(0.5)._ref_v,GABAneurons[neuron_j].synI,0,1,GABAneurons_W[neuron_i,neuron_j]))
    
GABAGlut_sin = list()
GlutGABA_sin = list()
for Gneuron_i in range(len(GABAneurons)):
    for neuron_i in range(len(Glutneurons)):
        GABAneurons[Gneuron_i].soma.push()
        GABAGlut_sin.append(h.NetCon(GABAneurons[Gneuron_i].soma(0.5)._ref_v,Glutneurons[neuron_i].synI,0,1,6./10000))
        h.pop_section()
        
        Glutneurons[neuron_i].soma.push()
        #GlutGABA_sin.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v,GABAneurons[Gneuron_i].synE,0,1,1./100000))
        GlutGABA_sin.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v,GABAneurons[Gneuron_i].synE,0,5,1./100000000)) #this is to turn off the connection (E to I)
        h.pop_section()

#create connections in network 2  (V1 superficial)      
        
NGABAn2 = 6
GABAneurons2,GABAneurons_rec2 = createnetwork(NGABAn2)
GABAneurons_W2 = np.random.exponential(1,NGABAn2*NGABAn2)*4/100000.
GABAneurons_W2 = GABAneurons_W2.reshape((NGABAn2,NGABAn2))
GABAneurons_W2 = GABAneurons_W2 - np.diag(np.diag(GABAneurons_W2))
    
Nneurons2 = 20
Glutneurons2,Glutneurons_rec2 = createnetwork(Nneurons2)
Glutneurons_W2 = np.random.exponential(1,Nneurons2*Nneurons2)*4/100000.
Glutneurons_W2 = Glutneurons_W2.reshape((Nneurons2,Nneurons2))
Glutneurons_W2 = Glutneurons_W2 - np.diag(np.diag(Glutneurons_W2))

GlutGlut_sin2 = list()
for neuron_i in range(len(Glutneurons2)):
    Glutneurons2[neuron_i].soma.push()
    for neuron_j in range(len(Glutneurons2)):
        GlutGlut_sin2.append(h.NetCon(Glutneurons2[neuron_i].soma(0.5)._ref_v,Glutneurons2[neuron_j].synE,0,1,Glutneurons_W2[neuron_i,neuron_j]))
    h.pop_section()
    

GABAGABA_sin2 = list()
for neuron_i in range(len(GABAneurons2)):
    GABAneurons2[neuron_i].soma.push()
    for neuron_j in range(len(GABAneurons2)):
        GABAGABA_sin2.append(h.NetCon(GABAneurons2[neuron_i].soma(0.5)._ref_v,GABAneurons2[neuron_j].synI,0,1,GABAneurons_W2[neuron_i,neuron_j]))
    
GABAGlut_sin2 = list()
GlutGABA_sin2 = list()
for Gneuron_i2 in range(len(GABAneurons2)):
    for neuron_i2 in range(len(Glutneurons2)):
        
        GABAneurons2[Gneuron_i2].soma.push()
        GABAGlut_sin2.append(h.NetCon(GABAneurons2[Gneuron_i2].soma(0.5)._ref_v,Glutneurons2[neuron_i2].synI,0,1,6./10000))
        h.pop_section()
        
        Glutneurons2[neuron_i2].soma.push()
        GlutGABA_sin2.append(h.NetCon(Glutneurons2[neuron_i2].soma(0.5)._ref_v,GABAneurons2[Gneuron_i2].synE,0,1,1./100000))
        h.pop_section()
#####
#extrinsic connections
#connections from Glutamatergic neurons of network 1 (LGN) to network 2 (V1)

GlutGlutneurons_W12 = np.random.exponential(1,Nneurons*Nneurons2)*4/100000.
GlutGlutneurons_W12 = GlutGlutneurons_W12.reshape((Nneurons,Nneurons2))
GlutGlutneurons_W12 = GlutGlutneurons_W12 - np.diag(np.diag(GlutGlutneurons_W12))

Glutnt1nt2_sin = list()
for neuron_i in range(len(Glutneurons)):
    Glutneurons[neuron_i].soma.push()
    for neuron_j in range(len(Glutneurons2)):
        Glutnt1nt2_sin.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v,Glutneurons2[neuron_j].synE,0,10,GlutGlutneurons_W12[neuron_i,neuron_j]))
    h.pop_section()

#connections from Glutamatergic neurons of network 2 (V1) to network 1 (LGN)

GlutGlutneurons_W21 = np.random.exponential(1,Nneurons*Nneurons2)*4/1000000.
#GlutGlutneurons_W21 = np.random.exponential(1,Nneurons*Nneurons2)*4/10000000000.             #turn off the connection
GlutGlutneurons_W21 = GlutGlutneurons_W21.reshape((Nneurons,Nneurons2))
GlutGlutneurons_W21 = GlutGlutneurons_W21 - np.diag(np.diag(GlutGlutneurons_W21))

Glutnt2nt1_sin = list()
for neuron_i in range(len(Glutneurons2)):
    Glutneurons2[neuron_i].soma.push()
    for neuron_j in range(len(Glutneurons)):
        Glutnt2nt1_sin.append(h.NetCon(Glutneurons2[neuron_i].soma(0.5)._ref_v,Glutneurons[neuron_j].synE,0,10,GlutGlutneurons_W21[neuron_i,neuron_j]))
    h.pop_section()


#generate inputs to network 1        
nstims = 5
stimrate = 1
netStim = list()
for stim_i in range(0,nstims):
    #netStim.append(h.NetStimPois(.5));
    netStim.append(h.NetStimPois(.1));
    netStim[stim_i].start=0;
    netStim[stim_i].mean = stimrate;
    netStim[stim_i].number=0;

#stim = h.NetCon(netStim[0],Glutneurons[0].synE,0.5,0,4./(100000))
stim = h.NetCon(netStim[0],Glutneurons[0].synE,0.5,0,4./(100000))
stim_rec = h.Vector()
stim.record(stim_rec)
#stim2 = h.NetCon(netStim[0],Glutneurons[1].synE,0.5,0,4./(100000))
stim2 = h.NetCon(netStim[1],Glutneurons[1].synE,0.5,0,4./(100000))
#stim8 = h.NetCon(netStim[2],Glutneurons[2].synE,0.5,0,4./(100000))
#stim9 = h.NetCon(netStim[3],Glutneurons[3].synE,0.5,0,4./(100000))
#stim10 = h.NetCon(netStim[4],Glutneurons[4].synE,0.5,0,4./(100000))



stim3 = h.NetCon(netStim[0],GABAneurons[0].synE,0.5,0,1./(100000))
#stim4 = h.NetCon(netStim[0],GABAneurons[1].synE,0.5,0,1./(100000))
#stim5 = h.NetCon(netStim[0],GABAneurons[2].synE,0.5,0,1./(100000))
#stim6 = h.NetCon(netStim[0],GABAneurons[3].synE,0.5,0,1./(100000))
#stim7 = h.NetCon(netStim[0],GABAneurons[4].synE,0.5,0,1./(100000))
    

timeaxis = h.Vector()
timeaxis.record(h._ref_t)

print np.shape(stim_rec.ref)

h.run() 



tmpmean = np.mean(Glutneurons_rec2,0)
tmpmean = tmpmean- np.mean(tmpmean)
tmpfft = np.fft.fft(tmpmean)
Fs =  1/2.5e-05


(Pxx, freqpsd) = psd(tmpmean, 20000/2, Fs) #args are signal, nfft, Fs
plt.figure()
plt.plot(freqpsd, Pxx)
plt.xlim([0,150])
plt.title('PSD of V1 LFP (PSD)')


####################################
#LGN LFP
tmpmean = np.mean(Glutneurons_rec,0)
tmpmean = tmpmean- np.mean(tmpmean)

(Pxx, freqpsd) = psd(tmpmean, 20000/2, Fs) #args are signal, nfft, Fs
plt.figure()
plt.plot(freqpsd, Pxx)
plt.xlim([0,150])
plt.title('PSD of LGN LFP (PSD)')




plt.figure()
a = plt.subplot(411)

#get fft of pyramidal neurons
#print "printings size"
#print np.shape(Glutneurons_rec)
tmpmean = np.mean(Glutneurons_rec,0)
#print np.shape(tmpmean)
#print "done printing size"
plt.plot(timeaxis,tmpmean)
plt.title('average membrane potential of excitatory cells in LGN')

a = plt.subplot(413)
for neuron_i in range(2,Nneurons):
    plt.plot(timeaxis,Glutneurons_rec[neuron_i])
    plt.title('Not receiving inputs')

a = plt.subplot(412,sharex=a)
for neuron_i in range(0,2):
    plt.plot(timeaxis,Glutneurons_rec[neuron_i])    
    plt.title('Receiving inputs')

plt.subplot(414,sharex=a)

for neuron_i in range(len(GABAneurons)):
    plt.plot(timeaxis,GABAneurons_rec[neuron_i])
plt.ylim([-100,50])
plt.title('inhibitory cells in LGN')

plt.xlim([0,h.tstop])


plt.figure()
b = plt.subplot(411)

#get fft of pyramidal neurons
#print "printings size"
#print np.shape(Glutneurons_rec)
tmpmean = np.mean(Glutneurons_rec2,0)
#print np.shape(tmpmean)
#print "done printing size"
plt.plot(timeaxis,tmpmean)
plt.title('average membrane potential of PY cells in V1' )

a = plt.subplot(413)
for neuron_i in range(2,Nneurons2):
    plt.plot(timeaxis,Glutneurons_rec2[neuron_i])
    plt.title('2-rest inputs net 2')

a = plt.subplot(412,sharex=a)
for neuron_i in range(0,2):
    plt.plot(timeaxis,Glutneurons_rec2[neuron_i])    
    plt.title('0-2 inputs net 2')

plt.subplot(414,sharex=a)

for neuron_i in range(len(GABAneurons2)):
    plt.plot(timeaxis,GABAneurons_rec2[neuron_i])
plt.ylim([-100,50])
plt.title('GABAergic net in V1')

plt.xlim([0,h.tstop])


def detect_spikes(voltagesignals):
    nneurons = len(voltagesignals)
    for neuron_i in range(nneurons):
        auxsignal = np.asarray(voltagesignals[neuron_i])
        thrscross_ind = [i for (i, val) in enumerate(auxsignal) if val > 30]

        spikeidx = list()
        aux = 0
        for thrscross_i in range(len(thrscross_ind)):
            if thrscross_ind[thrscross_i]>aux:
                auxspan = range(thrscross_ind[thrscross_i],thrscross_ind[thrscross_i]+15)
                spikeidx.append(auxspan.index(max(auxspan))+thrscross_ind[thrscross_i])
                aux = thrscross_ind[thrscross_i]+15
                
    return spikeidx
            
spikeidx = detect_spikes(GABAneurons_rec)
print spikeidx    
    
    
plt.show()

    
