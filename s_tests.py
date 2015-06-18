from simulation import simulate
import matplotlib.pyplot as plt
import numpy as np
from utils import exponential_connect, constant_connect

with_V1_L4 = True
with_V1_L6 = True
with_TRN = True

#FF
connect_E_LGN_E_L4 = True  # connected Hirsch connectivity
connect_E_LGN_I_L4 = True  # connected Hirsch connectivity
connect_E_LGN_E_L6 = True  # connected Hirsch connectivity
connect_E_LGN_I_L6 = True  # connected Hirsch connectivity
connect_E_L4_E_L6 = True

#FB
connect_E_L6_E_LGN = True
connect_E_L6_TRN = True

#TRN loop
connect_E_LGN_TRN = True
connect_TRN_E_LGN = True

#only set to true if include V1 L4 but not V1 L6
connect_E_L4_E_LGN = False
connect_E_L4_TRN = False

nruns = 2
total_time = 500

# number of cells should be divisible by 4, otherwise python will truncate (search simulation for "*1/4")
n_e_lgn = 20
n_i_lgn = 6
n_E_L6 = 20
n_I_L6 = 6
n_E_L4 = 20
n_I_L4 = 6
n_trn = 10

delay_distbtn_E_L6_LGN = (np.random.exponential(1, n_E_L4**2)*4)+4  #  Briggs F, Usrey W M. 2009, Neuron TODO: check if delay should be different between two different net2 cells connected to the same net1 cell
# delay_E_L6_TRN = 5
delay_distbtn_E_L6_TRN = (np.random.exponential(1, n_E_L4**2)*4)+4  #  Briggs F, Usrey W M. 2009, Neuron
delay_E_LGN_TRN = 1
delay_E_LGN_E_L4 = 3
delay_E_LGN_E_L6 = 3
delay_E_LGN_I_L6 = 3
delay_E_L4_E_LGN = 8
delay_E_LGN_I_L4 = 8  ## check this
delay_E_L4_TRN = 5
delay_E_LGN_I_LGN = 1
delay_I_LGN_E_LGN = 1

# stimrate 100 = 10 Hz, 10 = 100 Hz, 1 = 1000Hz, 5 = 200 Hz, 6 = 150 Hz
#2 spikes stimrate = 400
input = {'stimrate': 6,
         'input': 0.5,
         'nstims': 5
         }

input_glut_threshold = 0.5
input_glut_delay = 0
input_glut_weight = 5./100000  # PSP 5 ~ 20 mV (Wang & Hirsch 2010)
#input_glut_weight = 1./1000000  # ONLY FOR TESTING
input_gaba_threshold = 0.5
input_gaba_delay = 0
input_gaba_weight = 5./100000  # PSP 5 ~ 20  (Wang & Hirsch 2010)
#input_gaba_weight = 1./100000  # PSP 5 ~ 20  (Wang & Hirsch 2010)

#### INTRINSIC CONNECTIONS

lgn_params = {
    'w_e_lgn_e_lgn': constant_connect(0., n_e_lgn, n_e_lgn, False),  # TODO: reference for this -> LGN doesn't have intrinsic connections
    'w_i_lgn_i_lgn': exponential_connect(0., n_i_lgn, n_i_lgn, False),  # TODO: reference for this -> LGN doesn't have intrinsic connections
    'w_i_lgn_e_lgn': exponential_connect(1.75/100000, n_i_lgn, n_e_lgn),  # PSP 5 mV  (Wang & Hirsch 2010)
    #'w_i_lgn_e_lgn': exponential_connect(5/1000000, n_i_lgn, n_e_lgn),  # PSP 5 mV  (Wang & Hirsch 2010)
    #'w_i_lgn_e_lgn': exponential_connect(1/100000, n_i_lgn, n_e_lgn),
    'w_e_lgn_i_lgn': constant_connect(0., n_e_lgn, n_i_lgn)
}

# W_E_LGN_E_LGN = constant_connect(0., n_e_lgn, n_e_lgn, False)  # TODO: reference for this -> LGN doesn't have intrinsic connections
# W_I_LGN_I_LGN = exponential_connect(0., n_i_lgn, n_i_lgn, False)  # TODO: reference for this -> LGN doesn't have intrinsic connections
# W_I_LGN_E_LGN = exponential_connect(1.75/100000, n_i_lgn, n_e_lgn)  # PSP 5 mV  (Wang & Hirsch 2010)
# #W_I_LGN_E_LGN = exponential_connect(5/1000000, n_i_lgn, n_e_lgn)  # PSP 5 mV  (Wang & Hirsch 2010)
# #W_I_LGN_E_LGN = exponential_connect(1/100000, n_i_lgn, n_e_lgn)
# w_e_lgn_i_lgn = constant_connect(0., n_e_lgn, n_i_lgn)

W_I_L4_I_L4 = exponential_connect(8/100000., n_I_L4, n_I_L4, False)
W_E_L4_E_L4 = exponential_connect(4/100000., n_E_L4, n_E_L4, False)
W_I_L4_E_L4 = exponential_connect(6./10000, n_I_L4, n_E_L4)
W_E_L4_I_L4 = exponential_connect(1./100000, n_E_L4, n_I_L4)

W_I_L6_I_L6 = exponential_connect(8/100000., n_I_L6, n_I_L6, False)
W_E_L6_E_L6 = exponential_connect(4/100000., n_E_L6, n_E_L6, False)
W_I_L6_E_L6 = exponential_connect(6./10000, n_I_L6, n_E_L6)
W_E_L6_I_L6 = exponential_connect(1./100000, n_E_L6, n_I_L6)

W_TRN_TRN = exponential_connect(14/1000000., n_trn, n_trn, False)  # assuming linearity of psp with weights  // testing effect of connection weight
#W_TRN_TRN = exponential_connect(4/1000000., n_trn, n_trn, False)  # assuming linearity of psp with weights  // testing effect of connection weight

#EXTRINSIC CONNECTIONS

# W_E_LGN_TRN = exponential_connect(4/100000., n_e_lgn, n_trn)
# W_TRN_E_LGN = exponential_connect(4/100000., n_trn, n_e_lgn)
W_E_LGN_TRN = exponential_connect(14/1000000., n_e_lgn, n_trn)
W_TRN_E_LGN = exponential_connect(14/1000000., n_trn, n_e_lgn)

W_E_LGN_E_L4 = exponential_connect(4/100000., n_e_lgn, n_E_L4)
W_E_LGN_I_L4 = exponential_connect(4/1000000., n_e_lgn, n_I_L4)
W_E_L4_E_LGN = exponential_connect(4/1000000., n_E_L4, n_e_lgn)
# TODO: E_L4 is NOT connected to I_LGN. Find reference

W_E_LGN_E_L6 = exponential_connect(2/100000., n_e_lgn, n_E_L6)
W_E_LGN_I_L6 = exponential_connect(2/1000000., n_e_lgn, n_I_L6)
W_E_L6_E_LGN = exponential_connect(4.375/1000000.,  n_E_L6, n_e_lgn)  # 4 mv / 10 (Wang & Hirsch 2010) 300 pA (Granseth & Lindstrom 2002)
# TODO: check if E_L6 should be connected to I_LGN

W_E_L4_E_L6 = exponential_connect(4/100000., n_E_L4, n_E_L6)

W_E_L6_TRN = exponential_connect(4/1000000., n_E_L6, n_trn)
W_E_L4_TRN = W_E_L6_TRN  # only if net include V1 L4 but not V1 L6

simulate(nruns, total_time, with_V1_L4, with_V1_L6, with_TRN,
         input, input_glut_threshold, input_glut_delay, input_glut_weight, input_gaba_threshold, input_gaba_delay, input_gaba_weight,
         n_e_lgn, n_i_lgn, n_E_L6, n_I_L6, n_E_L4, n_I_L4, n_trn,
         delay_distbtn_E_L6_LGN, delay_E_L4_E_LGN, delay_E_LGN_I_L4, delay_E_LGN_E_L4, delay_E_LGN_E_L6,
         delay_E_LGN_TRN, delay_E_L4_TRN, delay_distbtn_E_L6_TRN, delay_E_LGN_I_LGN, delay_I_LGN_E_LGN, delay_E_LGN_I_L6,
         lgn_params, W_E_L4_E_L4, W_I_L4_I_L4, W_E_L6_E_L6, W_I_L6_I_L6, W_TRN_TRN,
         W_I_L4_E_L4, W_E_L4_I_L4, W_I_L6_E_L6, W_E_L6_I_L6,
         W_E_LGN_TRN, W_TRN_E_LGN, W_E_L6_TRN, W_E_L4_E_L6, W_E_LGN_E_L4, W_E_L4_E_LGN,
         W_E_L6_E_LGN, W_E_LGN_E_L6, W_E_LGN_I_L6, W_E_LGN_I_L4, W_E_L4_TRN,
         connect_E_LGN_E_L4, connect_E_LGN_I_L4, connect_E_L4_E_LGN, connect_E_LGN_I_L6, connect_E_LGN_E_L6, connect_E_L6_E_LGN, connect_E_L4_TRN, connect_E_L6_TRN,
         connect_E_LGN_TRN, connect_TRN_E_LGN, connect_E_L4_E_L6)

plt.show()
