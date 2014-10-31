from simulation import simulate
import matplotlib.pyplot as plt
import numpy as np
from utils import exponential_connect

with_V1_L4 = True
with_V1_L6 = True
with_TRN = True

connect_E_LGN_E_L4 = True
connect_E_LGN_I_L4 = True
connect_E_L4_E_LGN = False
connect_E_LGN_E_L6 = False  # not implemented yet
connect_E_L6_E_LGN = True
connect_E_L4_TRN = False
connect_E_L6_TRN = True

nruns = 5

# number of cells should be divisible by 4, otherwise python will truncate (search simulation for "*1/4")
Nneurons = 20
NGABAn = 6
NneuronsL6 = 20
NGABA_L6 = 6
NneuronsL4 = 20
NGABAL4 = 6
NGABA_trn = 10

delay_distbtn_E_L6_LGN = (np.random.exponential(1, NneuronsL4**2)*4)+4
delay_E_L4_E_LGN = 8
delay_E_LGN_I_L4 = 8  ## check this

input = {'stimrate': 6,
         'input': 0.5,
         'nstims': 5
         }

input_glut_threshold = 0.5
input_glut_delay = 0
input_glut_weight = 4.375/100000  # PSP 5 ~ 20 mV (Wang & Hirsch 2010)
input_gaba_threshold = 0.5
input_gaba_delay = 0
input_gaba_weight = 4.375/100000  # PSP 5 ~ 20  (Wang & Hirsch 2010)

#### connections

W_E_LGN_E_LGN = np.random.exponential(1, Nneurons*Nneurons)*4/100000.
W_E_LGN_E_LGN = W_E_LGN_E_LGN.reshape((Nneurons, Nneurons))
W_E_LGN_E_LGN = W_E_LGN_E_LGN - np.diag(np.diag(W_E_LGN_E_LGN))

W_I_LGN_I_LGN = np.random.exponential(1, NGABAn*NGABAn)*4/100000.
W_I_LGN_I_LGN = W_I_LGN_I_LGN.reshape((NGABAn, NGABAn))
W_I_LGN_I_LGN = W_I_LGN_I_LGN - np.diag(np.diag(W_I_LGN_I_LGN))

W_I_LGN_E_LGN = 1.75/100000*np.ones(NGABAn*Nneurons).reshape(NGABAn, Nneurons)  #  PSP 5 mV  (Wang & Hirsch 2010)

W_I_L4_I_L4 = np.random.exponential(1, NGABAL4*NGABAL4)*4/100000.
W_I_L4_I_L4 = W_I_L4_I_L4.reshape((NGABAL4, NGABAL4))
W_I_L4_I_L4 = W_I_L4_I_L4 - np.diag(np.diag(W_I_L4_I_L4))

W_E_L4_E_L4 = np.random.exponential(1, NneuronsL4*NneuronsL4)*4/100000.
W_E_L4_E_L4 = W_E_L4_E_L4.reshape((NneuronsL4, NneuronsL4))
W_E_L4_E_L4 = W_E_L4_E_L4 - np.diag(np.diag(W_E_L4_E_L4))

W_I_L4_E_L4 = 6./10000*np.ones(NGABAL4*NneuronsL4).reshape(NGABAL4, NneuronsL4)
W_E_L4_I_L4 = 1./100000*np.ones(NGABAL4*NneuronsL4).reshape(NneuronsL4, NGABAL4)

W_I_L6_I_L6 = np.random.exponential(1, NGABA_L6*NGABA_L6)*4/100000.
W_I_L6_I_L6 = W_I_L6_I_L6.reshape((NGABA_L6, NGABA_L6))
W_I_L6_I_L6 = W_I_L6_I_L6 - np.diag(np.diag(W_I_L6_I_L6))

W_E_L6_E_L6 = np.random.exponential(1, NneuronsL6*NneuronsL6)*4/100000.
W_E_L6_E_L6 = W_E_L4_E_L4.reshape((NneuronsL6, NneuronsL6))
W_E_L6_E_L6 = W_E_L4_E_L4 - np.diag(np.diag(W_E_L6_E_L6))

W_I_L6_E_L6 = 6./10000*np.ones(NGABA_L6*NneuronsL6).reshape(NGABA_L6, NneuronsL6)
W_E_L6_I_L6 = 1./100000*np.ones(NGABA_L6*NneuronsL6).reshape(NneuronsL6, NGABA_L6)

W_TRN_TRN = np.random.exponential(1, NGABA_trn*NGABA_trn)*14/1000000.  # assuming linearity of psp with weights
W_TRN_TRN = W_TRN_TRN.reshape((NGABA_trn, NGABA_trn))
W_TRN_TRN = W_TRN_TRN - np.diag(np.diag(W_TRN_TRN))

W_E_LGN_TRN = np.random.exponential(1, Nneurons*NGABA_trn)*1/1000000.
W_E_LGN_TRN = W_E_LGN_TRN.reshape((Nneurons, NGABA_trn))

W_TRN_E_LGN = 1/1000000.*np.random.exponential(1, NGABA_trn*Nneurons)
W_TRN_E_LGN = W_TRN_E_LGN.reshape(NGABA_trn, Nneurons)

W_E_LGN_E_L4 = np.random.exponential(1, Nneurons*NneuronsL4)*4/100000.
W_E_LGN_E_L4 = W_E_LGN_E_L4.reshape((Nneurons, NneuronsL4))

W_E_LGN_I_L4 = np.random.exponential(1, Nneurons*NGABAL4)*4/1000000.
W_E_LGN_I_L4 = W_E_LGN_I_L4.reshape(Nneurons, NGABAL4)

W_E_L4_E_LGN = np.random.exponential(1, Nneurons*NneuronsL4)*4/1000000.
#W_E_L4_E_LGN = np.random.exponential(1, Nneurons*NneuronsL4)*4/10000000000.  # turn off connection
W_E_L4_E_LGN = W_E_L4_E_LGN.reshape((Nneurons, NneuronsL4))

W_E_L6_E_LGN = np.random.exponential(1, Nneurons*NneuronsL6)*4.375/1000000.  # 4 mv / 10 (Wang & Hirsch 2010) 300 pA (Granseth & Lindstrom 2002)
#W_E_L6_E_LGN = np.random.exponential(1, Nneurons*NneuronsL6)*4/10000000000.  # turn off connection
W_E_L6_E_LGN = W_E_L6_E_LGN.reshape(Nneurons, NneuronsL6)

W_E_LGN_E_L6 = exponential_connect(4.375/1000000., Nneurons, NneuronsL6)

W_E_L6_TRN = np.random.exponential(1, NneuronsL6*NGABA_trn)*4/1000000.
#W_E_L6_TRN = np.random.exponential(1,NneuronsL6*NGABA_trn)*4/100000000000.  # turn off connection
W_E_L6_TRN = W_E_L6_TRN.reshape((NneuronsL6, NGABA_trn))

W_E_L4_E_L6 = np.random.exponential(1, NneuronsL4*NneuronsL6)*4/100000.
W_E_L4_E_L6 = W_E_L4_E_L6.reshape((NneuronsL4, NneuronsL6))

W_E_L4_TRN = np.random.exponential(1, NneuronsL4*NGABA_trn)*4/1000000.
#W_E_L4_TRN = np.random.exponential(1, NneuronsL4*NGABA_trn)*4/100000000000.  #turn off connection
W_E_L4_TRN = W_E_L4_TRN.reshape((NneuronsL4, NGABA_trn))

W_E_L6_TRN = W_E_L4_TRN

simulate(nruns, with_V1_L4, with_V1_L6, with_TRN,
         input, input_glut_threshold, input_glut_delay, input_glut_weight, input_gaba_threshold, input_gaba_delay, input_gaba_weight,
         Nneurons, NGABAn, NneuronsL6, NGABA_L6, NneuronsL4, NGABAL4, NGABA_trn,
         delay_distbtn_E_L6_LGN, delay_E_L4_E_LGN, delay_E_LGN_I_L4,
         W_E_LGN_E_LGN, W_I_LGN_I_LGN, W_E_L4_E_L4, W_I_L4_I_L4, W_E_L6_E_L6, W_I_L6_I_L6, W_TRN_TRN,
         W_I_LGN_E_LGN, W_I_L4_E_L4, W_E_L4_I_L4, W_I_L6_E_L6, W_E_L6_I_L6,
         W_E_LGN_TRN, W_TRN_E_LGN, W_E_L6_TRN, W_E_L4_E_L6, W_E_LGN_E_L4, W_E_L4_E_LGN, W_E_L6_E_LGN, W_E_LGN_E_L6, W_E_LGN_I_L4, W_E_L4_TRN,
         connect_E_LGN_E_L4, connect_E_LGN_I_L4, connect_E_L4_E_LGN, connect_E_LGN_E_L6, connect_E_L6_E_LGN, connect_E_L4_TRN, connect_E_L6_TRN)

plt.show()