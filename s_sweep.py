import os
from datetime import datetime
from simulation import simulate
import matplotlib.pyplot as plt
import numpy as np
from utils import exponential_connect, constant_connect

OUTPUT_DIR = "../data_files/"
temperature = 6.3  # temperature of Hodgkin & Huxley original experiments
dt = 0.025

with_V1_L4 = False
with_V1_L6 = False
with_TRN = False

#FF
connect_E_LGN_E_L4 = False  # connected Hirsch connectivity
connect_E_LGN_I_L4 = False  # connected Hirsch connectivity
connect_E_LGN_E_L6 = False  # connected Hirsch connectivity
connect_E_LGN_I_L6 = False  # connected Hirsch connectivity
connect_E_L4_E_L6 = False

#FB
connect_E_L6_E_LGN = False
connect_E_L6_TRN = False

#TRN loop
connect_E_LGN_TRN = False
connect_TRN_E_LGN = False

# only set to true if include V1 L4 but not V1 L6
connect_E_L4_E_LGN = False
connect_E_L4_TRN = False

nruns = 10
total_time = 500

# number of cells should be divisible by 4, otherwise python will truncate (search simulation for "*1/4")
n_e_lgn = 40
n_i_lgn = 12
n_e_l6 = 40
n_i_l6 = 12
n_e_l4 = 40
n_i_l4 = 12
n_trn = 20

# threshold for spike generation
threshold = 0

# delays of connection
delay = 1
delay_E_LGN_TRN = 1
delay_E_LGN_E_L4 = 3
delay_E_LGN_I_L4 = 8  ## check this
delay_E_LGN_E_L6 = 3
delay_E_LGN_I_L6 = 3
delay_E_L4_E_LGN = 8  # only if net include V1 L4 but not V1 L6
delay_E_L4_TRN = 5  # only if net include V1 L4 but not V1 L6
delay_distbtn_E_L6_LGN = (np.random.exponential(1, n_e_l4**2)*4)+4  #  Briggs F, Usrey W M. 2009, Neuron TODO: check if delay should be different between two different net2 cells connected to the same net1 cell
delay_distbtn_E_L6_TRN = (np.random.exponential(1, n_e_l4**2)*4)+4  #  Briggs F, Usrey W M. 2009, Neuron

# stimrate 100 = 10 Hz, 10 = 100 Hz, 1 = 1000Hz, 5 = 200 Hz, 6 = 150 Hz
input = {          # TODO: find reference saying that the input goes both to excitatory and inhibitory
    # 'stimrate': 4,
    #Backup
    'stimrate': 6,  # number of milliseconds of interval between stimuli
    'position': 0.5,  # position parameter, in case the input would be located in a NEURON section, which is not the case here. This parameter is irrelevant
    'nstims': 40  # number of stimuli
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

# w_i_l4_i_l4 = 0.5/100000.  # PSP = 1.55 mV Haeusler and Maass 2006
# w_e_l4_e_l4 = 0.5/100000.  # PSP = 1.1 mV Haeusler and Maass 2006
# w_i_l4_e_l4 = 0.3/100000.  # PSP = 0.85 mV Haeusler and Maass 2006
# w_e_l4_i_l4 = 1.3/100000.  # PSP = 3.7 mV Haeusler and Maass 2006

# Backup
w_i_l4_i_l4 = 8/100000.
w_e_l4_e_l4 = 4/100000.
w_i_l4_e_l4 = 6./10000.
w_e_l4_i_l4 = 1./100000

l4_params = {
    'w_i_l4_i_l4': exponential_connect(w_i_l4_i_l4, n_i_l4, n_i_l4, False),
    'w_e_l4_e_l4': exponential_connect(w_e_l4_e_l4, n_e_l4, n_e_l4, False),
    'w_i_l4_e_l4': exponential_connect(w_i_l4_e_l4, n_i_l4, n_e_l4),
    'w_e_l4_i_l4': exponential_connect(w_e_l4_i_l4, n_e_l4, n_i_l4),
    'p_i_i': 0.7,  # 50% connectivity Hauesler and Maass 2007
    'p_e_e': 0.3,  # 17% connectivity Hauesler and Maass 2007
    'p_i_e': 0.7,  # 10% connectivity Hauesler and Maass 2007
    'p_e_i': 0.3  # 19% connectivity Hauesler and Maass 2007
}

# w_i_l6_i_l6 = 0.5/100000  # PSP = 1.2 mV Haeusler and Maass 2007 (heuristic from L5)
# w_e_l6_e_l6 = 0.5/100000  # PSP = 1.7 mV Haeusler and Maass 2007 (heuristic from L5)
# w_e_l6_i_l6 = 0.3/100000  # PSP = 0.9 mV Haeusler and Maass 2007 (heuristic from L5)
# w_i_l6_e_l6 = 0.5/100000  # PSP = 1.2 mV Haeusler and Maass 2007 (heuristic from L5)

# Backup
w_i_l6_i_l6 = 8/100000.
w_e_l6_e_l6 = 4/100000.
w_e_l6_i_l6 = 1./100000.
w_i_l6_e_l6 = 6./10000.

l6_params = {
    'w_i_l6_i_l6': exponential_connect(w_i_l6_i_l6, n_i_l6, n_i_l6, False),
    'w_e_l6_e_l6': exponential_connect(w_e_l6_e_l6, n_e_l6, n_e_l6, False),
    'w_e_l6_i_l6': exponential_connect(w_e_l6_i_l6, n_e_l6, n_i_l6),
    'w_i_l6_e_l6': exponential_connect(w_i_l6_e_l6, n_i_l6, n_e_l6),
    'p_i_i': 0.7,  # 60% connectivity Hauesler and Maass 2007 (heuristic from L5)
    'p_e_e': 0.3,  # 9% connectivity Hauesler and Maass 2007 (heuristic from L5)
    'p_i_e': 0.7,  # 12% connectivity Hauesler and Maass 2007 (heuristic from L5)
    'p_e_i': 0.3   # 10% connectivity Hauesler and Maass 2007 (heuristic from L5)
}

trn_params = {
    'w_trn_trn': exponential_connect(14/1000000., n_trn, n_trn, False),  # assuming linearity of psp with weights  // testing effect of connection weight
    #'w_trn_trn: exponential_connect(4/1000000., n_trn, n_trn, False),  # assuming linearity of psp with weights  // testing effect of connection weight
    'p_i_i': 1,
    'delay_i_i': 1
}

#EXTRINSIC CONNECTIONS

# W_E_LGN_TRN = exponential_connect(1.4/1000000., n_e_lgn, n_trn)
# W_TRN_E_LGN = exponential_connect(1.4/1000000., n_trn, n_e_lgn)
#
# W_E_LGN_E_L4 = exponential_connect(.4/100000., n_e_lgn, n_e_l4)
# W_E_LGN_I_L4 = exponential_connect(.4/1000000., n_e_lgn, n_i_l4)
# W_E_L4_E_LGN = exponential_connect(.4/1000000., n_e_l4, n_e_lgn)
# # TODO: E_L4 is NOT connected to I_LGN. Find reference

# Backup
W_E_LGN_TRN = exponential_connect(14/1000000., n_e_lgn, n_trn)
W_TRN_E_LGN = exponential_connect(14/1000000., n_trn, n_e_lgn)

W_E_LGN_E_L4 = exponential_connect(4/100000., n_e_lgn, n_e_l4)
W_E_LGN_I_L4 = exponential_connect(4/1000000., n_e_lgn, n_i_l4)
W_E_L4_E_LGN = exponential_connect(4/1000000., n_e_l4, n_e_lgn)
# TODO: E_L4 is NOT connected to I_LGN. Find reference

# W_E_LGN_E_L6 = exponential_connect(0.2/100000., n_e_lgn, n_e_l6)
# W_E_LGN_I_L6 = exponential_connect(0.2/1000000., n_e_lgn, n_i_l6)
#Backup
W_E_LGN_E_L6 = exponential_connect(2/100000., n_e_lgn, n_e_l6)
W_E_LGN_I_L6 = exponential_connect(2/1000000., n_e_lgn, n_i_l6)
W_E_L6_E_LGN = exponential_connect(0.4375/100000.,  n_e_l6, n_e_lgn)  # 4 mv / 10 (Wang & Hirsch 2010) 300 pA (Granseth & Lindstrom 2002)  # TODO: isn't it Wang & Hirsch 2011?
# TODO: E_L6 is NOT connected to I_LGN. Find reference

# Backup
W_E_L4_E_L6 = exponential_connect(4/100000., n_e_l4, n_e_l6)

# PSP = 1.4 mV Haeusler and Maass 2006 (heuristic from E L2/3 to E L5)
# W_E_L4_E_L6 = exponential_connect(0.5/100000., n_e_l4, n_e_l6)

# W_E_L6_TRN = exponential_connect(0.4/1000000., n_e_l6, n_trn)
#Backup
W_E_L6_TRN = exponential_connect(4/1000000., n_e_l6, n_trn)
W_E_L4_TRN = W_E_L6_TRN  # only if net include V1 L4 but not V1 L6

weight = 4./100000
                    # E-E, I-I, I-E, E-I
connectivity = [
                    [weight, 0, 0, 0],                  # E-E
                    [0, weight, 0, 0],                  #     I-I
                    [0, 0, weight, 0],                  #         I-E
                    [0, 0, 0, weight],                  #             E-I
                    [weight, weight, 0, 0],             # E-E I-I
                    [weight, 0, weight, 0],             # E-E     I-E
                    [weight, 0, 0, weight],             # E-E         E-I
                    [0, weight, weight, 0],             #     I-I I-E
                    [0, weight, 0, weight],             #     I-I     E-I
                    [0, 0, weight, weight],             #         I-E E-I
                    [weight, weight, weight, 0],        # E-E I-I I-E
                    [weight, weight, 0, weight],        # E-E I-I     E-I
                    [weight, 0, weight, weight],        # E-E     I-E E-I
                    [0, weight, weight, weight],        #     I-I I-E E-I
                    [weight, weight, weight, weight]    # E-E I-I I-E E-I
                ]
time_sim_start = str(datetime.now())
fname_peaks = os.path.join(OUTPUT_DIR, "sweep_" + time_sim_start + ".txt")
open(fname_peaks, 'a').close()

for conn in connectivity:

    lgn_params = {
        'w_e_lgn_e_lgn': exponential_connect(conn[0], n_e_lgn, n_e_lgn, False),  # TODO: reference for this -> LGN doesn't have intrinsic connections
        'w_i_lgn_i_lgn': exponential_connect(conn[1], n_i_lgn, n_i_lgn, False),  # TODO: reference for this -> LGN doesn't have intrinsic connections
        'w_i_lgn_e_lgn': exponential_connect(conn[2], n_i_lgn, n_e_lgn),  # PSP 5 mV  (Wang & Hirsch 2010)
        'w_e_lgn_i_lgn': exponential_connect(conn[3], n_e_lgn, n_i_lgn),
        'delay_e_i': 1,
        'delay_i_e': 1
    }
    fname_lfps_prefix = os.path.join(OUTPUT_DIR, time_sim_start + "_" + str(conn) + "_sim-")
    print conn

    simulate(conn, OUTPUT_DIR, fname_peaks, fname_lfps_prefix, dt, nruns, total_time, temperature, with_V1_L4, with_V1_L6, with_TRN, input, con_input_lgn,
             n_e_lgn, n_i_lgn, n_e_l6, n_i_l6, n_e_l4, n_i_l4, n_trn,
             threshold, delay, delay_distbtn_E_L6_LGN, delay_E_L4_E_LGN, delay_E_LGN_I_L4, delay_E_LGN_E_L4, delay_E_LGN_E_L6,
             delay_E_LGN_TRN, delay_E_L4_TRN, delay_distbtn_E_L6_TRN, delay_E_LGN_I_L6,
             lgn_params, l4_params, l6_params, trn_params, W_E_LGN_TRN, W_TRN_E_LGN, W_E_L6_TRN, W_E_L4_E_L6,
             W_E_LGN_E_L4, W_E_L4_E_LGN, W_E_L6_E_LGN, W_E_LGN_E_L6, W_E_LGN_I_L6, W_E_LGN_I_L4, W_E_L4_TRN,
             connect_E_LGN_E_L4, connect_E_LGN_I_L4, connect_E_L4_E_LGN, connect_E_LGN_I_L6, connect_E_LGN_E_L6, connect_E_L6_E_LGN, connect_E_L4_TRN, connect_E_L6_TRN,
             connect_E_LGN_TRN, connect_TRN_E_LGN, connect_E_L4_E_L6)
