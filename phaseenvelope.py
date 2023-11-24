import Initialization
import matplotlib.pyplot as plt
import numpy as np

T_twop, P_twop, V_twop = Initialization.eos.get_envelope_twophase(1e5, Initialization.z, calc_v=True)


plt.plot(1/V_twop, T_twop)
for T in Initialization.T_values:
    plt.plot(np.array(Initialization.density[T])/Initialization.Avogadro, T*np.ones_like(Initialization.density[T]))
# Lables etc
plt.savefig("Plots/PhaseEnvelope.png")



