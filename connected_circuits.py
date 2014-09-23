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
with_V1_L6 = False
with_TRN = False

connect_E_LGN_E_L4 = True
connect_E_L4_E_LGN = True
connect_E_LGN_E_L6 = False
connect_E_L6_E_LGN = False
connect_E_L4_TRN = False
connect_E_L6_TRN = False

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

nruns = 1

for nsim in range(nruns):

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
            delayGlutGABA = 5
            #GlutGABA_sin.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v,GABAneurons[Gneuron_i].synE,0,delayGlutGABA,1./100000))
            GlutGABA_sin.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v,GABAneurons[Gneuron_i].synE,0,delayGlutGABA,1./100000000)) #this is to turn off the connection (E to I)
            h.pop_section()


    if with_V1_L4:
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

        if connect_E_LGN_E_L4:
            #####
            #extrinsic connections
            #connections from Glutamatergic neurons of network 1 (LGN) to network 2 (V1)

            GlutGlutneurons_W12 = np.random.exponential(1,Nneurons*Nneurons2)*4/100000.
            GlutGlutneurons_W12 = GlutGlutneurons_W12.reshape((Nneurons,Nneurons2))

            Glutnt1nt2_sin = list()
            delayGlutnt1nt2 = 10
            for neuron_i in range(len(Glutneurons)):
                Glutneurons[neuron_i].soma.push()
                for neuron_j in range(len(Glutneurons2)):
                    Glutnt1nt2_sin.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v,Glutneurons2[neuron_j].synE,0,delayGlutnt1nt2,GlutGlutneurons_W12[neuron_i,neuron_j]))
                h.pop_section()

        if connect_E_L4_E_LGN:
            #connections from Glutamatergic neurons of network 2 (V1) to network 1 (LGN)

            GlutGlutneurons_W21 = np.random.exponential(1,Nneurons*Nneurons2)*4/1000000.
            #GlutGlutneurons_W21 = np.random.exponential(1,Nneurons*Nneurons2)*4/10000000000.  #turn off connection
            GlutGlutneurons_W21 = GlutGlutneurons_W21.reshape((Nneurons,Nneurons2))
            GlutGlutneurons_W21 = GlutGlutneurons_W21 - np.diag(np.diag(GlutGlutneurons_W21))  # BUG: the diagonal should NOT be subtracted because this is a connections between two different network, not a net with itself

            Glutnt2nt1_sin = list()
            delayGlutnt2nt1 = 10
            for neuron_i in range(len(Glutneurons2)):
                Glutneurons2[neuron_i].soma.push()
                for neuron_j in range(len(Glutneurons)):
                    Glutnt2nt1_sin.append(h.NetCon(Glutneurons2[neuron_i].soma(0.5)._ref_v,Glutneurons[neuron_j].synE,0,delayGlutnt2nt1,GlutGlutneurons_W21[neuron_i,neuron_j]))
                h.pop_section()


    if with_V1_L6:        
        #create connections in network 2  (V1 L6)      
                
        NGABAn2 = 6
        GABAneuronsL6,GABAneurons_recL6 = createnetworkL6(NGABAn2)
        GABAneurons_WL6 = np.random.exponential(1,NGABAn2*NGABAn2)*4/100000.
        GABAneurons_WL6 = GABAneurons_WL6.reshape((NGABAn2,NGABAn2))
        GABAneurons_WL6 = GABAneurons_WL6 - np.diag(np.diag(GABAneurons_WL6))
            
        Nneurons2 = 20
        GlutneuronsL6,Glutneurons_recL6 = createnetworkL6(Nneurons2)
        Glutneurons_WL6 = np.random.exponential(1,Nneurons2*Nneurons2)*4/100000.
        Glutneurons_WL6 = Glutneurons_W2.reshape((Nneurons2,Nneurons2))
        Glutneurons_WL6 = Glutneurons_W2 - np.diag(np.diag(Glutneurons_WL6))
        
        GlutGlutL6_sin = list()
        for neuron_i in range(len(GlutneuronsL6)):
            GlutneuronsL6[neuron_i].soma.push()
            for neuron_j in range(len(GlutneuronsL6)):
                GlutGlutL6_sin.append(h.NetCon(GlutneuronsL6[neuron_i].soma(0.5)._ref_v,GlutneuronsL6[neuron_j].synE,0,1,Glutneurons_WL6[neuron_i,neuron_j]))
            h.pop_section()
            

        GABAGABAL6_sin = list()
        for neuron_i in range(len(GABAneuronsL6)):
            GABAneuronsL6[neuron_i].soma.push()
            for neuron_j in range(len(GABAneurons2)):
                GABAGABAL6_sin.append(h.NetCon(GABAneuronsL6[neuron_i].soma(0.5)._ref_v,GABAneuronsL6[neuron_j].synI,0,1,GABAneurons_WL6[neuron_i,neuron_j]))
            
        GABAGlutL6_sin = list()
        GlutGABAL6_sin = list()
        for Gneuron_i2 in range(len(GABAneuronsL6)):
            for neuron_i2 in range(len(GlutneuronsL6)):
                
                GABAneuronsL6[Gneuron_i2].soma.push()
                GABAGlutL6_sin.append(h.NetCon(GABAneuronsL6[Gneuron_i2].soma(0.5)._ref_v,GlutneuronsL6[neuron_i2].synI,0,1,6./10000))
                h.pop_section()
                
                GlutneuronsL6[neuron_i2].soma.push()
                GlutGABAL6_sin.append(h.NetCon(GlutneuronsL6[neuron_i2].soma(0.5)._ref_v,GABAneuronsL6[Gneuron_i2].synE,0,1,1./100000))
                h.pop_section()   

        #connections from V1 input layer to L6

        GlutGlutneurons_WL4L6 = np.random.exponential(1,Nneurons2*Nneurons2)*4/100000.
        GlutGlutneurons_WL4L6 = GlutGlutneurons_WL4L6.reshape((Nneurons2,Nneurons2))

        GlutGlut_L4L6_sin = list()
        for neuron_i in range(len(Glutneurons2)):
            Glutneurons2[neuron_i].soma.push()
            for neuron_j in range(len(GlutneuronsL6)):
                GlutGlut_L4L6_sin.append(h.NetCon(Glutneurons2[neuron_i].soma(0.5)._ref_v,GlutneuronsL6[neuron_j].synE,0,1,GlutGlutneurons_WL4L6[neuron_i,neuron_j]))
            h.pop_section() 

        if connect_E_L6_E_LGN:
            #connections from Glutamatergic neurons of network 2 (V1) to network 1 (LGN)

            GlutGlutneurons_WL61 = np.random.exponential(1,Nneurons*Nneurons2)*4/1000000.
            #GlutGlutneurons_WL61 = np.random.exponential(1,Nneurons*Nneurons2)*4/10000000000.  #turn off connection
            GlutGlutneurons_WL61 = GlutGlutneurons_WL61.reshape((Nneurons,Nneurons2))

            GlutL6nt1_sin = list()
            delayGlutL6nt1 = 8
            for neuron_i in range(len(GlutneuronsL6)):
                GlutneuronsL6[neuron_i].soma.push()
                for neuron_j in range(len(Glutneurons)):
                    GlutL6nt1_sin.append(h.NetCon(GlutneuronsL6[neuron_i].soma(0.5)._ref_v,Glutneurons[neuron_j].synE,0, delayGlutL6nt1,GlutGlutneurons_WL61[neuron_i,neuron_j]))
                h.pop_section()

    if with_TRN:
        #create TRN neurons (inhibitory only)

        NGABA_trn = 10
        GABAneurons_trn, GABAneurons_trn_rec = createnetwork(NGABA_trn)
        GABAneurons_trnW = np.random.exponential(1,NGABA_trn*NGABA_trn)*4/100000.
        GABAneurons_trnW = GABAneurons_trnW.reshape((NGABA_trn,NGABA_trn))
        GABAneurons_trnW = GABAneurons_trnW - np.diag(np.diag(GABAneurons_trnW))

        GABAGABA_trn_sin = list()
        for neuron_i in range(len(GABAneurons_trn)):
            GABAneurons_trn[neuron_i].soma.push()
            for neuron_j in range(len(GABAneurons_trn)):
                GABAGABA_trn_sin.append(h.NetCon(GABAneurons_trn[neuron_i].soma(0.5)._ref_v,GABAneurons_trn[neuron_j].synI,0,1,GABAneurons_trnW[neuron_i,neuron_j]))
                #h.pop_section()

        #connections from Glutamatergic neurons of network 2 (V1) to TRN

        GlutGABAneurons_Wnet2trn = np.random.exponential(1,Nneurons2*NGABA_trn)*4/1000000.
        #GlutGABAneurons_Wnet2trn = np.random.exponential(1,Nneurons2*NGABA_trn)*4/100000000000.  #turn off connection
        GlutGABAneurons_Wnet2trn = GlutGABAneurons_Wnet2trn.reshape((Nneurons2,NGABA_trn))

        if with_V1_L4 and connect_E_L4_TRN:
            GlutGABAtneurons_sin2 = list()
            for neuron_i in range(len(Glutneurons2)):
                Glutneurons2[neuron_i].soma.push()
                for neuron_j in range(len(GABAneurons_trn)):
                    GlutGABAtneurons_sin2.append(h.NetCon(Glutneurons2[neuron_i].soma(0.5)._ref_v,GABAneurons_trn[neuron_j].synE,0,5,GlutGABAneurons_Wnet2trn[neuron_i,neuron_j]))
                h.pop_section()

        if with_V1_L6 and connect_E_L6_TRN:
            GlutGABAtneurons_sinL6trn = list()
            for neuron_i in range(len(GlutneuronsL6)):
                GlutneuronsL6[neuron_i].soma.push()
                for neuron_j in range(len(GABAneurons_trn)):
                    GlutGABAtneurons_sinL6trn.append(h.NetCon(GlutneuronsL6[neuron_i].soma(0.5)._ref_v,GABAneurons_trn[neuron_j].synE,0,5,GlutGABAneurons_Wnet2trn[neuron_i,neuron_j]))
                h.pop_section()

        #connections from Glutamatergic neurons of network 1 (LGN) to TRN

        GlutGABAtneurons_Wnet1trn = np.random.exponential(1,Nneurons*NGABA_trn)*4/1000000.
        GlutGABAtneurons_Wnet1trn = GlutGABAtneurons_Wnet1trn.reshape((Nneurons,NGABA_trn))

        GlutGABAtneurons_sin1 = list()
        delayGlutGABAtneurons = 1
        for neuron_i in range(len(Glutneurons)):
            Glutneurons[neuron_i].soma.push()
            for neuron_j in range(len(GABAneurons_trn)):
                GlutGABAtneurons_sin1.append(h.NetCon(Glutneurons[neuron_i].soma(0.5)._ref_v,GABAneurons_trn[neuron_j].synE,0,delayGlutGABAtneurons,GlutGABAtneurons_Wnet1trn[neuron_i,neuron_j]))
            h.pop_section()

        #connectinos from GABAergic nuerons of TRN to network 1 (LGN)
        GABAGlutneurons_W = np.random.exponential(1,NGABA_trn*Nneurons)*8/1000000.
        GABAGlutneurons_W = GABAGlutneurons_W.reshape((NGABA_trn,Nneurons))

        GABAGlutneurons_sin = list()
        for neuron_i in range(len(GABAneurons_trn)):
            GABAneurons_trn[neuron_i].soma.push()
            for neuron_j in range(len(Glutneurons)):
                GABAGlutneurons_sin.append(h.NetCon(GABAneurons_trn[neuron_i].soma(0.5)._ref_v,Glutneurons[neuron_j].synE,0,1,GABAGlutneurons_W[neuron_i,neuron_j]))
            h.pop_section()



    #generate inputs to network 1
    nstims = 5
    stimrate = 1
                   #100 = 10 Hz, 10 = 100 Hz, 1 = 1000Hz,5 = 200 Hz
    netStim = list()
    for stim_i in range(0,nstims):
        input = 0.1
        netStim.append(h.NetStimPois(input))
        netStim[stim_i].start=0
        netStim[stim_i].mean = stimrate
        netStim[stim_i].number=0

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
    #stim2 = h.NetCon(netStim[0],Glutneurons[1].synE,0.5,0,4./(100000))
    stim2 = h.NetCon(netStim[1],Glutneurons[1].synE,0.5,0,4./(100000))
    #stim3 = h.NetCon(netStim[0],Glutneurons[2].synE,0.5,0,4./(100000))
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

    Fs =  1/2.5e-05
    if with_V1_L4:
        #try again
        meanV1input = np.mean(Glutneurons_rec2,0)
        meanV1input = meanV1input- np.mean(meanV1input)
        
        
        print "length of time axis"
        print len(timeaxis)

        #lensignal = len(tmpmean)
        #freq = np.fft.fftfreq(lensignal,2.5e-05)

        #lenfft = len(tmpfft)

        #print "signal length"
        #print len(tmpmean)

        #fs = 1000
        #f = fs*(np.asarray(range(lenfft))-1)/lenfft

        #print "len fft"
        #print lenfft
        #print "len freq"
        #print len(freq)

        #plt.figure(5)
        #plt.plot(freq[1:lenfft/2], abs(tmpfft[1:lenfft/2]))
        #plt.xlim([0,150])
        #plt.title('PSD of V1 LFP (FFT)')


        (Pxx, freqpsd) = psd(meanV1input, 20000/2, Fs) #args are signal, nfft, Fs
        plt.figure(6)
        plt.plot(freqpsd, Pxx)
        plt.xlim([0,150])
        plt.title('PSD of V1 LFP (PSD)')


    ####################################
    #LGN LFP
    meanLGN = np.mean(Glutneurons_rec,0)
    meanLGN = meanLGN- np.mean(meanLGN)
    
    #plt.figure(7)
    #plt.plot(freq[1:lenfft/2], abs(tmpfft[1:lenfft/2]))
    #plt.xlim([0,150])
    #plt.title('PSD of LGN LFP (FFT)')
    
    
    (Pxx, freqpsd) = psd(meanLGN, 20000/2, Fs) #args are signal, nfft, Fs
    plt.figure(8)
    plt.plot(freqpsd, Pxx)
    plt.xlim([0,150])
    plt.title('PSD of LGN LFP (PSD)')

    if with_TRN:
        ####################################
        #TRN LFP 
        meanTRN = np.mean(GABAneurons_trn_rec,0)
        meanTRN = meanTRN- np.mean(meanTRN)
        
        #plt.figure(7)
        #plt.plot(freq[1:lenfft/2], abs(tmpfft[1:lenfft/2]))
        #plt.xlim([0,150])
        #plt.title('PSD of LGN LFP (FFT)')


        (Pxx, freqpsd) = psd(meanTRN, 20000/2, Fs) #args are signal, nfft, Fs
        plt.figure(9)
        plt.plot(freqpsd, Pxx)
        plt.xlim([0,150])
        plt.title('PSD of TRN LFP (PSD)')

        #plot trn neurons

        plt.figure(3)
        b = plt.subplot(211)
        
        plt.plot(timeaxis,meanTRN)
        plt.title('average membrane potential of TRN cells' )

        a = plt.subplot(212)
        for neuron_i in range(0,NGABA_trn):
            plt.plot(timeaxis,GABAneurons_trn_rec[neuron_i])
            plt.title('2-rest inputs net 2')

        plt.ylim([-100,50])
        plt.title('GABAergic TRN neurons')

        plt.xlim([0,h.tstop])

    if with_V1_L6:
        ####################################
        #L6 LFP 
        meanV1output = np.mean(Glutneurons_recL6,0)
        meanV1output = meanV1output- np.mean(meanV1output)

        (Pxx, freqpsd) = psd(meanV1output, 20000/2, Fs) #args are signal, nfft, Fs
        plt.figure(10)
        plt.plot(freqpsd, Pxx)
        plt.xlim([0,150])
        plt.title('PSD of V1 L6 LFP (PSD)')

        #plot L6 neurons

        plt.figure(4)
        a = plt.subplot(311)

        plt.plot(timeaxis,meanV1output)
        plt.title('average membrane potential of L6 cells' )

        a = plt.subplot(312)
        for neuron_i in range(2,Nneurons2):
            plt.plot(timeaxis,Glutneurons_recL6[neuron_i])
        plt.title('Glut L6')

        a = plt.subplot(313)
        for neuron_i in range(0,NGABAn2):
            plt.plot(timeaxis,GABAneurons_recL6[neuron_i])

        plt.ylim([-100,50])
        plt.title('GABAergic L6 neurons')
        plt.xlim([0,h.tstop])

    #a = plt.subplot(512)

    #tmpfft = np.fft.fft(tmpmean, 1024)

    #print "fftsize"
    #print np.shape(tmpfft)

    #signal_pow = np.abs(tmpfft)

    #fs = (1/ (timeaxis[2] - timeaxis[1]))*10
    #t = timeaxis
    #
    #print "printing time axis"
    #print timeaxis[0]
    #print timeaxis[1]
    #
    #print "printing fs"
    #print fs
    #
    #tmpmean = tmpmean - np.mean(tmpmean)
    #
    #xF = np.fft.fft(tmpmean,512)
    #N = len(xF)
    #
    #print "printing len fft"
    #print N
    #
    ##xF = xF[0:N/2]
    #fr = np.linspace(0,fs/2,N/2)
    #
    #print "printing freq axis"
    #print fr
    #print "end freq axis"ize()
    #
    #print "printing fft"
    #print xF
    #
    ##plt.plot(fr,abs(xF)**2)
    #plt.plot(abs(xF)**2)

    plt.figure(1)
    a = plt.subplot(411)

    #get fft of pyramidal neurons
    #print "printings size"
    #print np.shape(Glutneurons_rec)
    #print np.shape(tmpmean)
    #print "done printing size"
    plt.plot(timeaxis,meanLGN)
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

    if with_V1_L4:
        plt.figure(2)
        b = plt.subplot(411)

        #get fft of pyramidal neurons
        #print "printings size"
        #print np.shape(Glutneurons_rec)
        #print np.shape(tmpmean)
        #print "done printing size"
        plt.plot(timeaxis,meanV1input)
        plt.title('average membrane potential of PY cells net in V1 (L4)' )
        
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
        plt.title('GABAergic net in V1 (L4)')
        
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
    
    ofname = "../data_files/" "sim" + str(nsim+96) + ".txt"
    
    indx = np.arange(1,20001,40)
    
    w = np.array(meanLGN)
    #u = np.array(meanTRN)
    x = np.array(meanV1input)
    #y = np.array(meanV1output)
    z = np.array(timeaxis)
    np.savetxt(ofname, (w[indx], x[indx], z[indx]))
    #np.savetxt(ofname, (tmpmean_TRN, tmpmean_V1, tmpmean_LGN, timeaxis))
    
    print nsim


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

    
