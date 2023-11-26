import Initialization
import matplotlib.pyplot as plt
import numpy as np

T_twop, P_twop, V_twop = Initialization.eos.get_envelope_twophase(1e5, Initialization.z, calc_v=True)


plt.plot(1/V_twop, T_twop, label='Phase Envelope')
for T in Initialization.T_values:
    label = f'Temperature {T:.2f}'
    plt.plot(np.array(Initialization.density[T])/Initialization.Avogadro, T*np.ones_like(Initialization.density[T]), label=label)
plt.xlabel('Density')
plt.ylabel('T [K]')
plt.title('Phase Envelope')
plt.legend()
plt.savefig("Plots/PhaseEnvelope.png")



