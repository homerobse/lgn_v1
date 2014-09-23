import matplotlib.pyplot as plt
from matplotlib.mlab import psd
import numpy as np
from neuron import h

def plot_all(timeaxis, stim_rec, with_V1_L4, with_V1_L6, with_TRN, Glutneurons_rec2, GABAneurons_trn_rec,
             Glutneurons_recL6, Glutneurons_rec, GABAneurons_recL6, GABAneurons, GABAneurons2, GABAneurons_rec, GABAneurons_rec2,
             Nneurons, NneuronsL6, Nneurons2, NGABA_trn, NGABA_L6):
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

    Fs = 1/2.5e-05
    if with_V1_L4:
        meanV1input = np.mean(Glutneurons_rec2, 0)
        meanV1input = meanV1input - np.mean(meanV1input)

        #lensignal = len(tmpmean)
        #freq = np.fft.fftfreq(lensignal,2.5e-05)

        #lenfft = len(tmpfft)

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


        (Pxx, freqpsd) = psd(meanV1input, 20000/2, Fs)  #args are signal, nfft, Fs
        plt.figure(6)
        plt.plot(freqpsd, Pxx)
        plt.xlim([0, 150])
        plt.title('PSD of V1 LFP (PSD)')


    ####################################
    #LGN LFP
    meanLGN = np.mean(Glutneurons_rec, 0)
    meanLGN = meanLGN - np.mean(meanLGN)

    #plt.figure(7)
    #plt.plot(freq[1:lenfft/2], abs(tmpfft[1:lenfft/2]))
    #plt.xlim([0,150])
    #plt.title('PSD of LGN LFP (FFT)')


    (Pxx, freqpsd) = psd(meanLGN, 20000/2, Fs)  # args are signal, nfft, Fs
    plt.figure(8)
    plt.plot(freqpsd, Pxx)
    plt.xlim([0,150])
    plt.title('PSD of LGN LFP (PSD)')

    if with_TRN:
        ####################################
        #TRN LFP
        meanTRN = np.mean(GABAneurons_trn_rec,0)
        meanTRN = meanTRN - np.mean(meanTRN)

        #plt.figure(7)
        #plt.plot(freq[1:lenfft/2], abs(tmpfft[1:lenfft/2]))
        #plt.xlim([0,150])
        #plt.title('PSD of LGN LFP (FFT)')


        (Pxx, freqpsd) = psd(meanTRN, 20000/2, Fs)  # args are signal, nfft, Fs
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
        for neuron_i in range(0, NGABA_trn):
            plt.plot(timeaxis, GABAneurons_trn_rec[neuron_i])
            plt.title('2-rest inputs net 2')

        plt.ylim([-100,50])
        plt.title('GABAergic TRN neurons')

        plt.xlim([0, h.tstop])

    if with_V1_L6:
        ####################################
        #L6 LFP
        meanV1output = np.mean(Glutneurons_recL6,0)
        meanV1output = meanV1output- np.mean(meanV1output)

        (Pxx, freqpsd) = psd(meanV1output, 20000/2, Fs)  # args are signal, nfft, Fs
        plt.figure(10)
        plt.plot(freqpsd, Pxx)
        plt.xlim([0, 150])
        plt.title('PSD of V1 L6 LFP (PSD)')

        #plot L6 neurons

        plt.figure(4)
        a = plt.subplot(311)

        plt.plot(timeaxis,meanV1output)
        plt.title('average membrane potential of L6 cells' )

        a = plt.subplot(312)
        for neuron_i in range(2, NneuronsL6):
            plt.plot(timeaxis, Glutneurons_recL6[neuron_i])
        plt.title('Glut L6')

        a = plt.subplot(313)
        for neuron_i in range(0, NGABA_L6):
            plt.plot(timeaxis, GABAneurons_recL6[neuron_i])

        plt.ylim([-100,50])
        plt.title('GABAergic L6 neurons')
        plt.xlim([0, h.tstop])

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
    plt.plot(timeaxis, meanLGN)
    plt.title('average membrane potential of excitatory cells in LGN')

    a = plt.subplot(413)
    for neuron_i in range(2, Nneurons):
        plt.plot(timeaxis, Glutneurons_rec[neuron_i])
        plt.title('Not receiving inputs')

    a = plt.subplot(412, sharex=a)
    for neuron_i in range(0,2):
        plt.plot(timeaxis,Glutneurons_rec[neuron_i])
        plt.title('Receiving inputs')

    plt.subplot(414, sharex=a)

    for neuron_i in range(len(GABAneurons)):
        plt.plot(timeaxis, GABAneurons_rec[neuron_i])
    plt.ylim([-100, 50])
    plt.title('inhibitory cells in LGN')

    plt.xlim([0, h.tstop])

    if with_V1_L4:
        plt.figure(2)
        b = plt.subplot(411)

        #get fft of pyramidal neurons
        #print "printings size"
        #print np.shape(Glutneurons_rec)
        #print np.shape(tmpmean)
        #print "done printing size"
        plt.plot(timeaxis, meanV1input)
        plt.title('average membrane potential of PY cells net in V1 (L4)')

        a = plt.subplot(413)
        for neuron_i in range(2, Nneurons2):
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

        plt.xlim([0, h.tstop])


    #plot input
    plt.figure(5)
    plt.plot(stim_rec, np.ones([stim_rec.size(),1]), 'ob')
    plt.xlim([timeaxis[0], timeaxis[-1]])
    plt.ylim([0., 2.])

    return meanLGN, meanTRN, meanV1input, meanV1output