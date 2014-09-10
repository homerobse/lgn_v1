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

with_V1_L4 = False
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

        self.synI.tau1=4
        self.synI.tau2=2       
        
        self.stm = h.IClamp(0.5,sec=self.soma)
        self.stm.amp=0
        self.stm.dur=1000
        self.stm.delay=100
        
        
class L6cell:
    "Layer 6 cell"

    def __init__(self):
       
        self.soma = h.Section(name='soma',cell=self)
        self.soma.L=30
        self.soma.nseg=1
        self.soma.diam=1

        self.soma.Ra=100
        self.soma.cm=1

        #self.soma.insert("hh")
        self.soma.insert("hhvitor2")
        
        self.soma.gna_hhvitor2=0.039
        self.soma.gkd_hhvitor2=0.006
        self.soma.gl_hhvitor2=0
        self.soma.gm_hhvitor2=0
        self.soma.gt_hhvitor2=0
        self.soma.gleak_hhvitor2 = 0.0000273
   
        self.synE = h.Exp2Syn(0.5,sec=self.soma)
        self.synE.tau1=10
        self.synE.tau2=20

        self.synI = h.Exp2Syn(0.5,sec=self.soma)
        self.synI.e=-100
        self.synI.tau1=5
        self.synI.tau2=40
        
        #self.synI.tau1=1
        #self.synI.tau2=2        
        
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
    
def createnetworkL6(Nneurons=4):

    network=h.List()
    network_rec=h.List()
 
    for Nindex in range(Nneurons):
        p = L6cell()
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
        GlutGABA_sin.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v,GABAneurons[Gneuron_i].synE,0,1,1./100000000)) #this is to turn off the connection (E to I)
        h.pop_section()
  
       


#generate inputs to network 1        
nstims = 5
stimrate = 1
               #100 = 10 Hz, 10 = 100 Hz, 1 = 1000Hz,5 = 200 Hz
netStim = list()
for stim_i in range(0,nstims):
    #netStim.append(h.NetStimPois(.5));
    netStim.append(h.NetStimPois(.5));
    netStim[stim_i].start=0;
    netStim[stim_i].mean = stimrate;
    netStim[stim_i].number=0;

#stim = h.NetCon(netStim[0],Glutneurons[0].synE,0.5,0,4./(100000))
#stim_rec = h.Vector()
#stim.record(stim_rec)
#stim8 = h.NetCon(netStim[0],Glutneurons[1].synE,0.5,0,4./(100000))
#stim9 = h.NetCon(netStim[0],Glutneurons[2].synE,0.5,0,4./(100000))
#stim10 = h.NetCon(netStim[0],Glutneurons[3].synE,0.5,0,4./(100000))
#stim11 = h.NetCon(netStim[0],Glutneurons[4].synE,0.5,0,4./(100000))
#
#
#
#stim3 = h.NetCon(netStim[0],GABAneurons[0].synE,0.5,0,4./(100000))
#stim4 = h.NetCon(netStim[0],GABAneurons[1].synE,0.5,0,4./(100000))
#stim5 = h.NetCon(netStim[0],GABAneurons[2].synE,0.5,0,4./(100000))
#stim6 = h.NetCon(netStim[0],GABAneurons[3].synE,0.5,0,4./(100000))
#stim7 = h.NetCon(netStim[0],GABAneurons[4].synE,0.5,0,4./(100000))



#stim = h.NetCon(netStim[0],Glutneurons[0].synE,0.5,0,4./(100000))
stim = h.NetCon(netStim[0],Glutneurons[0].synE,0.5,0,4./(100000))
stim_rec = h.Vector()
stim.record(stim_rec)
stim2 = h.NetCon(netStim[0],Glutneurons[1].synE,0.5,0,4./(100000))
#stim3 = h.NetCon(netStim[0],Glutneurons[2].synE,0.5,0,4./(100000))
#stim8 = h.NetCon(netStim[0],Glutneurons[3].synE,0.5,0,4./(100000))
#stim9 = h.NetCon(netStim[0],Glutneurons[4].synE,0.5,0,4./(100000))
#stim10 = h.NetCon(netStim[0],Glutneurons[5].synE,0.5,0,4./(100000))


stim3 = h.NetCon(netStim[0],GABAneurons[0].synE,0.5,0,1./(100000))
#stim4 = h.NetCon(netStim[0],GABAneurons[1].synE,0.5,0,1./(100000))
#stim5 = h.NetCon(netStim[0],GABAneurons[2].synE,0.5,0,1./(100000))
#stim6 = h.NetCon(netStim[0],GABAneurons[3].synE,0.5,0,1./(100000))
#stim7 = h.NetCon(netStim[0],GABAneurons[4].synE,0.5,0,1./(100000))
    

timeaxis = h.Vector()
timeaxis.record(h._ref_t)

print np.shape(stim_rec.ref)

h.run() 

#vitor
#x = rand(1000)
#y = sin(t*0.1*pi*2)
#z = x+y
#a = np.fft.fft(z)

#t= range(len(x))

#t = t/ fs

#fs = 1?
# f = fs*(asarray(range(len(a)))-1)/len(a)
#plot(f,abs(a))


####################################
#LGN LFP
tmpmean = np.mean(Glutneurons_rec,0)
tmpmean = tmpmean- np.mean(tmpmean)
tmpfft = np.fft.fft(tmpmean)

#plt.figure(7)
#plt.plot(freq[1:lenfft/2], abs(tmpfft[1:lenfft/2]))
#plt.xlim([0,150])
#plt.title('PSD of LGN LFP (FFT)')

Fs =  1/2.5e-05
(Pxx, freqpsd) = psd(tmpmean, 20000/2, Fs) #args are signal, nfft, Fs
plt.figure(8)
plt.plot(freqpsd, Pxx)
plt.xlim([0,150])
plt.title('PSD')
#plt.ylim([0. , 40.])


plt.figure(1)
a = plt.subplot(411)

#get fft of pyramidal neurons
#print "printings size"
#print np.shape(Glutneurons_rec)
tmpmean_LGN = np.mean(Glutneurons_rec,0)
#print np.shape(tmpmean)
#print "done printing size"
plt.plot(timeaxis,tmpmean_LGN)
plt.title('average membrane potential of excitatory cells')

a = plt.subplot(413)
for neuron_i in range(2,Nneurons):
    plt.plot(timeaxis,Glutneurons_rec[neuron_i])
    plt.title('Not receiving inputs')

a = plt.subplot(412,sharex=a)
for neuron_i in range(0,5):
    plt.plot(timeaxis,Glutneurons_rec[neuron_i])    
    plt.title('Receiving inputs')

plt.subplot(414,sharex=a)

for neuron_i in range(len(GABAneurons)):
    plt.plot(timeaxis,GABAneurons_rec[neuron_i])
plt.ylim([-100,50])
plt.title('inhibitory cells')

plt.xlim([0,h.tstop])


#plot input
plt.figure(5)
plt.plot(stim_rec, np.ones([stim_rec.size(),1]), 'ob')
plt.xlim([timeaxis[0], timeaxis[-1]])
plt.ylim([0. , 2.])


#nsim = 1
#ofname = "./data_files/" "sim" + str(nsim) + ".txt"

#np.savetxt('./data_files/output_file.txt', (tmpmean_TRN, tmpmean_V1, tmpmean_LGN, timeaxis))
#np.savetxt(ofname, (tmpmean_TRN, tmpmean_V1, tmpmean_LGN, timeaxis))


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

    
