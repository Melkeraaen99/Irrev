from pykingas.MieKinGas import MieKinGas
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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
pressure_start = 1e5  # Pascal -> 1 bar
pressure_end = 100e5
pressure_increment = 5e5
pressure = np.arange(pressure_start, pressure_end + 1, pressure_increment)
x = [0.5, 0.5]  # Molar composition
z = [1]
m_Ar = 39.948*1.67377e-27 # kg
molar_mass_argon = 39.948e-3 # kg/mol
kb = Boltzmann # J/K

# Experimental data
df_200 = pd.read_csv('nist_data_200.csv', delimiter='\t')
df_250 = pd.read_csv('nist_data_250.csv', delimiter='\t')
df_300 = pd.read_csv('nist_data_300.csv', delimiter='\t')

viscosity_200_exp = df_200['Viscosity (Pa*s)'].to_numpy()
thermal_conductivity_200_exp = df_200['Therm. Cond. (W/m*K)'].to_numpy()

viscosity_250_exp = df_250['Viscosity (Pa*s)'].to_numpy()
thermal_conductivity_250_exp = df_250['Therm. Cond. (W/m*K)'].to_numpy()

viscosity_300_exp = df_300['Viscosity (Pa*s)'].to_numpy()
thermal_conductivity_300_exp = df_300['Therm. Cond. (W/m*K)'].to_numpy()

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

percentage_difference_200 = np.abs((visc_data[200] - viscosity_200_exp) / ((visc_data[200] + viscosity_200_exp) / 2)) * 100
mean_percentage_difference_200 = np.mean(percentage_difference_200)
print('Mean Percentage Difference of viscosity calculated and viscosity from ecperimental data for Argon at 200 K:', mean_percentage_difference_200)

percentage_difference_250 = np.abs((visc_data[250] - viscosity_250_exp) / ((visc_data[250] + viscosity_250_exp) / 2)) * 100
mean_percentage_difference_250 = np.mean(percentage_difference_250)
print('Mean Percentage Difference of viscosity calculated and viscosity from ecperimental data for Argon at 250 K:', mean_percentage_difference_250)

percentage_difference_300 = np.abs((visc_data[300] - viscosity_300_exp) / ((visc_data[300] + viscosity_300_exp) / 2)) * 100
mean_percentage_difference_300 = np.mean(percentage_difference_300)
print('Mean Percentage Difference of viscosity calculated and viscosity from ecperimental data for Argon at 300 K:', mean_percentage_difference_300)

# Scaled viscosity
def dim_less_visc(T):
    if T in density:
        return np.log((np.array(density[T])*Avogadro)**(-2/3) * (kb * T* m_Ar)**(-0.5) * np.array(visc_data[T]))
    else:
        return None 

plt.figure(figsize=(8, 6))

plt.plot(pressure/1e5, visc_data[200], color='b', label='Visc thermopack 200 K')
plt.plot(pressure/1e5, viscosity_200_exp, color='r', label='Visc exp 200 K')
plt.plot(pressure/1e5, visc_data[250], color='b', label='Visc thermopack 250 K')
plt.plot(pressure/1e5, viscosity_250_exp, color='r', label='Visc exp 250 K')
plt.plot(pressure/1e5, visc_data[300], color='b', label='Visc thermopack 300 K')
plt.plot(pressure/1e5, viscosity_300_exp, color='r', label='Visc exp 300 K')

plt.text(pressure[0]/1e5, visc_data[200][0] + 0.1e-5, f' Mean Diff: {mean_percentage_difference_200:.2f}%', color='k')
plt.text(pressure[0]/1e5, visc_data[250][0] + 0.1e-5, f' Mean Diff: {mean_percentage_difference_250:.2f}%', color='k')
plt.text(pressure[0]/1e5, visc_data[300][0] + 0.1e-5, f' Mean Diff: {mean_percentage_difference_300:.2f}%', color='k')

plt.xlabel('Pressure (bar)')
plt.ylabel('Viscosity (Pa*s)')
plt.title('Viscosity (exp & calc) vs Pressure')
plt.legend()

plt.savefig("Viscosity_exp_viscosity.png")


'''# Create a new figure
plt.figure(figsize=(8, 6))

# Plot calculated residual entropy with log dimless viscosity
for T in T_values:
    label = f'Temperature {T} K (Ideal Gas)'
    plt.plot(entropy[T], dim_less_visc(T), 'o', label=label)

# Add labels and a legend
plt.xlabel(r'S$_{res}\,$ [ J$\cdot$mol$^{-1}\cdot$K$^{-1}$]')
plt.ylabel(r'$\eta^*$')
plt.title(r'Dimensionless viscosity, $\eta^*$, plotted against residual entropy, S$_{res}$, for Argon')
plt.legend()
#plt.ylim(0, 1e-7)

# Show the plot
plt.grid(True)
plt.savefig("Visc.png")'''