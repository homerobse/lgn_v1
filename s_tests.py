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

nruns = 3
total_time = 500

# number of cells should be divisible by 4, otherwise python will truncate (search simulation for "*1/4")
n_e_lgn = 20
n_i_lgn = 6
n_e_l6 = 20
n_i_l6 = 6
n_e_l4 = 20
n_i_l4 = 6
n_trn = 10

delay_distbtn_E_L6_LGN = (np.random.exponential(1, n_e_l4**2)*4)+4  #  Briggs F, Usrey W M. 2009, Neuron TODO: check if delay should be different between two different net2 cells connected to the same net1 cell
# delay_E_L6_TRN = 5
delay_distbtn_E_L6_TRN = (np.random.exponential(1, n_e_l4**2)*4)+4  #  Briggs F, Usrey W M. 2009, Neuron
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
input = {
    'stimrate': 6,
    'input': 0.5,  # TODO: what does this parameter mean?
    'nstims': 4
}

con_input_lgn = {
    'glut_threshold': 0,
    'glut_delay': 0,
    'glut_weight': 5./100000,  # PSP 5 ~ 20 mV (Wang & Hirsch 2010) - Using 15 mV
    'gaba_threshold':  0,
    'gaba_delay': 0,
    'gaba_weight': 5./100000  # PSP 5 ~ 20  (Wang & Hirsch 2010) - Using 15 mV
}

#### INTRINSIC CONNECTIONS

lgn_params = {
    'w_e_lgn_e_lgn': constant_connect(0., n_e_lgn, n_e_lgn, False),  # TODO: reference for this -> LGN doesn't have intrinsic connections
    'w_i_lgn_i_lgn': exponential_connect(0., n_i_lgn, n_i_lgn, False),  # TODO: reference for this -> LGN doesn't have intrinsic connections
    'w_i_lgn_e_lgn': exponential_connect(1.75/100000, n_i_lgn, n_e_lgn),  # PSP 5 mV  (Wang & Hirsch 2010)
    'w_e_lgn_i_lgn': constant_connect(0., n_e_lgn, n_i_lgn)
}

# w_i_l4_i_l4 = 0.5/100000.  # PSP = 1.55 mV Haeusler and Maass 2006
# w_e_l4_e_l4 = 0.5/100000.  # PSP = 1.1 mV Haeusler and Maass 2006
# w_i_l4_e_l4 = 0.3/100000.  # PSP = 0.85 mV Haeusler and Maass 2006
# w_e_l4_i_l4 = 1.3/100000.  # PSP = 3.7 mV Haeusler and Maass 2006

# # Backup
w_i_l4_i_l4 = 8/100000.
w_e_l4_e_l4 = 4/100000.
w_i_l4_e_l4 = 6./10000.
w_e_l4_i_l4 = 1./100000

l4_params = {
    'w_i_l4_i_l4': exponential_connect(w_i_l4_i_l4, n_i_l4, n_i_l4, False),
    'w_e_l4_e_l4': exponential_connect(w_e_l4_e_l4, n_e_l4, n_e_l4, False),
    'w_i_l4_e_l4': exponential_connect(w_i_l4_e_l4, n_i_l4, n_e_l4),
    'w_e_l4_i_l4': exponential_connect(w_e_l4_i_l4, n_e_l4, n_i_l4)
}

# w_i_l6_i_l6 = 0.5/100000  # PSP = 1.2 mV Haeusler and Maass 2006 (heuristic from L5)
# w_e_l6_e_l6 = 0.5/100000  # PSP = 1.7 mV Haeusler and Maass 2006 (heuristic from L5)
# w_e_l6_i_l6 = 0.3/100000  # PSP = 1.7 mV Haeusler and Maass 2006 (heuristic from L5)
# w_i_l6_e_l6 = 0.5/100000  # PSP = 1.2 mV Haeusler and Maass 2006 (heuristic from L5)

# Backup
w_i_l6_i_l6 = 8/100000.
w_e_l6_e_l6 = 4/100000.
w_e_l6_i_l6 = 1./100000.
w_i_l6_e_l6 = 6./10000.

l6_params = {
    'w_i_l6_i_l6': exponential_connect(w_i_l6_i_l6, n_i_l6, n_i_l6, False),
    'w_e_l6_e_l6': exponential_connect(w_e_l6_e_l6, n_e_l6, n_e_l6, False),
    'w_e_l6_i_l6': exponential_connect(w_e_l6_i_l6, n_e_l6, n_i_l6),
    'w_i_l6_e_l6': exponential_connect(w_i_l6_e_l6, n_i_l6, n_e_l6)
}

W_TRN_TRN = exponential_connect(14/1000000., n_trn, n_trn, False)  # assuming linearity of psp with weights  // testing effect of connection weight
#W_TRN_TRN = exponential_connect(4/1000000., n_trn, n_trn, False)  # assuming linearity of psp with weights  // testing effect of connection weight

#EXTRINSIC CONNECTIONS

# W_E_LGN_TRN = exponential_connect(4/100000., n_e_lgn, n_trn)
# W_TRN_E_LGN = exponential_connect(4/100000., n_trn, n_e_lgn)
W_E_LGN_TRN = exponential_connect(14/1000000., n_e_lgn, n_trn)
W_TRN_E_LGN = exponential_connect(14/1000000., n_trn, n_e_lgn)

W_E_LGN_E_L4 = exponential_connect(4/100000., n_e_lgn, n_e_l4)
W_E_LGN_I_L4 = exponential_connect(4/1000000., n_e_lgn, n_i_l4)
W_E_L4_E_LGN = exponential_connect(4/1000000., n_e_l4, n_e_lgn)
# TODO: E_L4 is NOT connected to I_LGN. Find reference

W_E_LGN_E_L6 = exponential_connect(2/100000., n_e_lgn, n_e_l6)
W_E_LGN_I_L6 = exponential_connect(2/1000000., n_e_lgn, n_i_l6)
W_E_L6_E_LGN = exponential_connect(4.375/1000000.,  n_e_l6, n_e_lgn)  # 4 mv / 10 (Wang & Hirsch 2010) 300 pA (Granseth & Lindstrom 2002)
# TODO: E_L6 is NOT connected to I_LGN. Find reference

W_E_L4_E_L6 = exponential_connect(4/100000., n_e_l4, n_e_l6)

W_E_L6_TRN = exponential_connect(4/1000000., n_e_l6, n_trn)
W_E_L4_TRN = W_E_L6_TRN  # only if net include V1 L4 but not V1 L6

simulate(nruns, total_time, with_V1_L4, with_V1_L6, with_TRN, input, con_input_lgn,
         n_e_lgn, n_i_lgn, n_e_l6, n_i_l6, n_e_l4, n_i_l4, n_trn,
         delay_distbtn_E_L6_LGN, delay_E_L4_E_LGN, delay_E_LGN_I_L4, delay_E_LGN_E_L4, delay_E_LGN_E_L6,
         delay_E_LGN_TRN, delay_E_L4_TRN, delay_distbtn_E_L6_TRN, delay_E_LGN_I_LGN, delay_I_LGN_E_LGN, delay_E_LGN_I_L6,
         lgn_params, l4_params, l6_params, W_TRN_TRN, W_E_LGN_TRN, W_TRN_E_LGN, W_E_L6_TRN, W_E_L4_E_L6,
         W_E_LGN_E_L4, W_E_L4_E_LGN, W_E_L6_E_LGN, W_E_LGN_E_L6, W_E_LGN_I_L6, W_E_LGN_I_L4, W_E_L4_TRN,
         connect_E_LGN_E_L4, connect_E_LGN_I_L4, connect_E_L4_E_LGN, connect_E_LGN_I_L6, connect_E_LGN_E_L6, connect_E_L6_E_LGN, connect_E_L4_TRN, connect_E_L6_TRN,
         connect_E_LGN_TRN, connect_TRN_E_LGN, connect_E_L4_E_L6)

plt.show()
