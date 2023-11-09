from pykingas.MieKinGas import MieKinGas
import numpy as np
import matplotlib.pyplot as plt

from thermopack.cubic import cubic
from thermopack.saftvrmie import saftvrmie
from scipy.constants import Boltzmann
from scipy.constants import R
from scipy.constants import Avogadro

# Create a MieKinGas instance
#eos = cubic('AR', 'PR') 
eos = saftvrmie('AR') 
kin = MieKinGas('AR', use_eos=eos)  # AR, modeled with RET-Mie

# Define temperature and pressure ranges
T_values = [200, 250, 300]  # Kelvin
pressure_start = 10e5  # Pascal -> 1 bar
pressure_end = 100e5
pressure_increment = 5e5
pressure = np.arange(pressure_start, pressure_end + 1, pressure_increment)
x = [0.5, 0.5]  # Molar composition
z = [1]
m_Ar = 39.948*1.67377e-27 # kg
molar_mass_argon = 39.948e-3 # kg/mol
kb = Boltzmann # J/K

# Store Residual entropy with thermopack and by calculating it
entropy = {T: [] for T in T_values}

# Initialize lists to store viscosity data for different temperatures
visc_data = {T: [] for T in T_values}

# Density
density = {T: [] for T in T_values}

for T in T_values:
    for p in pressure:
        visc = kin.viscosity_tp(T, p, x, N=2)
        visc_data[T].append(visc)
        spes_volume = eos.specific_volume(T, p, z, eos.VAPPH)
        density[T].append(1/spes_volume[0])
        entro = eos.entropy_tv(T, spes_volume[0], z, property_flag='R')
        entropy[T].append(entro) 


# Scaled viscosity
def dim_less_visc(T):
    if T in density:
        return np.log((np.array(density[T])*Avogadro)**(-2/3) * (kb * T* m_Ar)**(-0.5) * np.array(visc_data[T]))
    else:
        return None 


# Create a new figure
plt.figure(figsize=(8, 6))

# Plot calculated residual entropy with log dimless viscosity
for T in T_values:
    label = f'Temperature {T} K (Ideal Gas)'
    plt.plot(entropy[T], dim_less_visc(T), 'o', label=label)

# Add labels and a legend
plt.xlabel(r'S$_{res}\,$ [ J$\cdot$mol$^{-1}\cdot$K$^{-1}$]')
plt.ylabel(r'$\eta^*$')
plt.title(r'Dimensionless viscosity ,$\eta^*$, plotted against residual entropy, S$_{res}$, for Argon')
plt.legend()
#plt.ylim(0, 1e-7)

# Show the plot
plt.grid(True)
plt.savefig("Visc.png")