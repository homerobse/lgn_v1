from simulation import simulate
import matplotlib.pyplot as plt
import numpy as np

with_V1_L4 = True
with_V1_L6 = True
with_TRN = True

connect_E_LGN_E_L4 = True
connect_E_L4_E_LGN = False
connect_E_LGN_E_L6 = False
connect_E_L6_E_LGN = True
connect_E_L4_TRN = False
connect_E_L6_TRN = True

nruns = 2

Nneurons = 20
NGABAn = 6
NneuronsL6 = 20
NGABA_L6 = 6
NneuronsL4 = 20
NGABAL4 = 6

glut_L6_nt1_delay_distbtn = (np.random.exponential(1, NneuronsL4**2)*4)+4


input = {'stimrate': 6,
         'input': 0.5,
         'nstims': 5
         }

simulate(nruns, with_V1_L4, with_V1_L6, with_TRN,
         input, Nneurons, NGABAn, NneuronsL6, NGABA_L6, NneuronsL4, NGABAL4,
         glut_L6_nt1_delay_distbtn,
         connect_E_LGN_E_L4, connect_E_L4_E_LGN, connect_E_LGN_E_L6, connect_E_L6_E_LGN, connect_E_L4_TRN, connect_E_L6_TRN)

plt.show()