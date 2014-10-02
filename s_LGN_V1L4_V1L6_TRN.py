from simulation import simulate
import matplotlib.pyplot as plt
import numpy as np

with_V1_L4 = False
with_V1_L6 = False
with_TRN = True

connect_E_LGN_E_L4 = False
connect_E_LGN_I_L4 = True
connect_E_L4_E_LGN = False
connect_E_LGN_E_L6 = False
connect_E_L6_E_LGN = False
connect_E_L4_TRN = False
connect_E_L6_TRN = False

nruns = 2

Nneurons = 20
NGABAn = 6
NneuronsL6 = 20
NGABA_L6 = 6
NneuronsL4 = 20
NGABAL4 = 6
NGABA_trn = 10

delay_distbtn_E_L6_LGN = (np.random.exponential(1, NneuronsL4**2)*4)+4
delay_E_L4_E_LGN = 8
delay_E_LGN_I_L4 = 4  ## check this

input = {'stimrate': 6,
         'input': 0.5,
         'nstims': 5
         }

input_glut_threshold = 0.5
input_glut_delay = 0
input_glut_weight = 4./100000
input_gaba_threshold = 0.5
input_gaba_delay = 0
input_gaba_weight = 4./100000

#### connections

W_E_LGN_TRN = np.random.exponential(1, Nneurons*NGABA_trn)*1/1000000.
W_E_LGN_TRN = W_E_LGN_TRN.reshape((Nneurons, NGABA_trn))

W_E_L6_TRN = np.random.exponential(1, NneuronsL6*NGABA_trn)*4/1000000.
#W_E_L6_TRN = np.random.exponential(1,NneuronsL6*NGABA_trn)*4/100000000000.  #turn off connection
W_E_L6_TRN = W_E_L6_TRN.reshape((NneuronsL6, NGABA_trn))

W_E_L4_E_L6 = np.random.exponential(1, NneuronsL4*NneuronsL6)*4/100000.
W_E_L4_E_L6 = W_E_L4_E_L6.reshape((NneuronsL4, NneuronsL6))

W_E_LGN_E_L4 = np.random.exponential(1, Nneurons*NneuronsL4)*4/100000.
W_E_LGN_E_L4 = W_E_LGN_E_L4.reshape((Nneurons, NneuronsL4))

W_E_L4_E_LGN = np.random.exponential(1, Nneurons*NneuronsL4)*4/1000000.
#W_E_L4_E_LGN = np.random.exponential(1, Nneurons*NneuronsL4)*4/10000000000.  #turn off connection
W_E_L4_E_LGN = W_E_L4_E_LGN.reshape((Nneurons, NneuronsL4))

W_E_L6_E_LGN = np.random.exponential(1, Nneurons*NneuronsL6)*4/1000000.
#W_E_L6_E_LGN = np.random.exponential(1, Nneurons*NneuronsL6)*4/10000000000.  #turn off connection
W_E_L6_E_LGN = W_E_L6_E_LGN.reshape(Nneurons, NneuronsL6)

W_E_LGN_I_L4 = np.random.exponential(1, Nneurons*NGABAL4)*4/1000000.
W_E_LGN_I_L4 = W_E_LGN_I_L4.reshape(Nneurons, NGABAL4)

W_E_L4_TRN = np.random.exponential(1, NneuronsL4*NGABA_trn)*4/1000000.
#W_E_L4_TRN = np.random.exponential(1, NneuronsL4*NGABA_trn)*4/100000000000.  #turn off connection
W_E_L4_TRN = W_E_L4_TRN.reshape((NneuronsL4, NGABA_trn))

W_E_L6_TRN = W_E_L4_TRN

simulate(nruns, with_V1_L4, with_V1_L6, with_TRN,
         input, input_glut_threshold, input_glut_delay, input_glut_weight, input_gaba_threshold, input_gaba_delay, input_gaba_weight,
         Nneurons, NGABAn, NneuronsL6, NGABA_L6, NneuronsL4, NGABAL4, NGABA_trn,
         delay_distbtn_E_L6_LGN, delay_E_L4_E_LGN, delay_E_LGN_I_L4,
         W_E_LGN_TRN, W_E_L6_TRN, W_E_L4_E_L6, W_E_LGN_E_L4, W_E_L4_E_LGN, W_E_L6_E_LGN, W_E_LGN_I_L4, W_E_L4_TRN,
         connect_E_LGN_E_L4, connect_E_LGN_I_L4, connect_E_L4_E_LGN, connect_E_LGN_E_L6, connect_E_L6_E_LGN, connect_E_L4_TRN, connect_E_L6_TRN)

plt.show()