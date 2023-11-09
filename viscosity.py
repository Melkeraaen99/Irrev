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
kin_id = MieKinGas('AR', is_idealgas = True)

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
entropy_2 = {T: [] for T in T_values}
entro_3 = {T: [] for T in T_values}

# Initialize lists to store viscosity data for different temperatures
visc_data = {T: [] for T in T_values}
visc_data_id = {T: [] for T in T_values}

# Density
density = {T: [] for T in T_values}

for T in T_values:
    for p in pressure:
        entro, = eos.entropy(T, p, z, eos.VAPPH, residual=True) # Feil i denne beregningen da residual tror p er v
        entro_full, = eos.entropy(T, p, z, eos.VAPPH)
        entro_null, = eos.entropy(T, 0.1, z, eos.VAPPH)
        entropy_2[T].append(entro_full-entro_null) # Ikke bruk, 
        entropy[T].append(entro)
        visc = kin.viscosity_tp(T, p, x, N=2)
        visc_id = kin_id.viscosity_tp(T, 1, x, N=2)
        visc_data[T].append(visc)
        visc_data_id[T].append(visc_id)
        spes_volume = eos.specific_volume(T, p, z, eos.VAPPH)
        density[T].append(1/spes_volume[0])
        entropy_3 = eos.entropy_tv(T, spes_volume[0], z, property_flag='R')
        entro_3[T].append(entropy_3) 

print(visc_data)

# Scaled viscosity
def dim_less_visc(T): 
    return np.log(np.array(visc_data[T])/visc_data_id[T])

def dim_less_visc_2(T):
    if T in density:
        return np.log((np.array(density[T])*Avogadro)**(-2/3) * (kb * T* m_Ar)**(-0.5) * np.array(visc_data[T]))
    else:
        return None 


# Create a new figure
plt.figure(figsize=(8, 6))

# Plot Residual entropy with resiudal return value from thermopack
'''for T in T_values:
    label = f'Temperature {T} K (Non-Ideal Gas)'
    plt.plot(entropy[T], dim_less_visc, 'x', label=label)'''

# Plot calculated residual entropy with log scaled viscosity, try different scaling
for T in T_values:
    label = f'Temperature {T} K (Ideal Gas)'
    plt.plot(entro_3[T], dim_less_visc_2(T), 'o', label=label)

# Add labels and a legend
plt.xlabel(r'Residual entropy [J·mol$^{-1}$·K$^{-1}$]')
plt.ylabel(r'$\eta$/$\eta_0$')
plt.title('Argon')
plt.legend()
#plt.ylim(0, 1e-7)

# Show the plot
plt.grid(True)
plt.savefig("Visc.png")