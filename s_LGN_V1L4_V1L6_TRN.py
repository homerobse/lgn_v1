from simulation import simulate
import matplotlib.pyplot as plt

with_V1_L4 = True
with_V1_L6 = True
with_TRN = True

connect_E_LGN_E_L4 = True
connect_E_L4_E_LGN = False
connect_E_LGN_E_L6 = False
connect_E_L6_E_LGN = True
connect_E_L4_TRN = False
connect_E_L6_TRN = True

nruns = 1

simulate(nruns, with_V1_L4, with_V1_L6, with_TRN, connect_E_LGN_E_L4, connect_E_L4_E_LGN,
             connect_E_LGN_E_L6, connect_E_L6_E_LGN, connect_E_L4_TRN, connect_E_L6_TRN)

plt.show()