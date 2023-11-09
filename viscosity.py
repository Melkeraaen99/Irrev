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

# Experimental data for viscosity and conductivity
df_200 = pd.read_csv('nist_data_200.csv', delimiter='\t')
df_250 = pd.read_csv('nist_data_250.csv', delimiter='\t')
df_300 = pd.read_csv('nist_data_300.csv', delimiter='\t')

viscosity_200_exp = df_200['Viscosity (Pa*s)'].to_numpy()
thermal_conductivity_200_exp = df_200['Therm. Cond. (W/m*K)'].to_numpy()

viscosity_250_exp = df_250['Viscosity (Pa*s)'].to_numpy()
thermal_conductivity_250_exp = df_250['Therm. Cond. (W/m*K)'].to_numpy()

viscosity_300_exp = df_300['Viscosity (Pa*s)'].to_numpy()
thermal_conductivity_300_exp = df_300['Therm. Cond. (W/m*K)'].to_numpy()

visc_exp = [viscosity_200_exp, viscosity_250_exp, viscosity_300_exp]
cond_exp = [thermal_conductivity_200_exp, thermal_conductivity_250_exp, thermal_conductivity_300_exp]


# Store Residual entropy with thermopack and by calculating it
entropy = {T: [] for T in T_values}

# Initialize lists to store viscosity data for different temperatures
visc_data = {T: [] for T in T_values}
cond_data = {T: [] for T in T_values}

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
        cond = kin.thermal_conductivity(T, spes_volume[0], x, N=2)
        cond_data[T].append(cond)

def calculate_mean_percentage_difference(data_calculated, data_exp, property_name):
    percentage_difference = np.abs((data_calculated - data_exp) / ((data_calculated + data_exp) / 2)) * 100
    mean_percentage_difference = np.mean(percentage_difference)
    return mean_percentage_difference

mean_diff_visc_200 = calculate_mean_percentage_difference(visc_data[200], viscosity_200_exp, 'viscosity')
mean_diff_visc_250 = calculate_mean_percentage_difference(visc_data[250], viscosity_250_exp, 'viscosity')
mean_diff_visc_300 = calculate_mean_percentage_difference(visc_data[300], viscosity_300_exp, 'viscosity')
mean_diff_visc = [mean_diff_visc_200, mean_diff_visc_250, mean_diff_visc_300]

mean_diff_cond_200 = calculate_mean_percentage_difference(cond_data[200], thermal_conductivity_200_exp, 'conductivity')
mean_diff_cond_250 = calculate_mean_percentage_difference(cond_data[250], thermal_conductivity_250_exp, 'conductivity')
mean_diff_cond_300 = calculate_mean_percentage_difference(cond_data[300], thermal_conductivity_300_exp, 'conductivity')
mean_diff_cond = [mean_diff_cond_200, mean_diff_cond_250, mean_diff_cond_300]

# Scaled viscosity
def dim_less_visc(T):
    if T in density:
        return np.log((np.array(density[T])*Avogadro)**(-2/3) * (kb * T* m_Ar)**(-0.5) * np.array(visc_data[T]))
    else:
        return None 
    

'''# Viscosity
plt.figure(figsize=(8, 6))

plt.plot(pressure/1e5, visc_data[200], color='b', label='Visc thermopack 200 K')
plt.plot(pressure/1e5, viscosity_200_exp, color='r', label='Visc exp 200 K')
plt.plot(pressure/1e5, visc_data[250], color='b', label='Visc thermopack 250 K')
plt.plot(pressure/1e5, viscosity_250_exp, color='r', label='Visc exp 250 K')
plt.plot(pressure/1e5, visc_data[300], color='b', label='Visc thermopack 300 K')
plt.plot(pressure/1e5, viscosity_300_exp, color='r', label='Visc exp 300 K')

plt.text(pressure[0]/1e5, visc_data[200][0] + 0.1e-5, f' Mean Diff: {mean_diff_visc_200:.2f}%', color='k')
plt.text(pressure[0]/1e5, visc_data[250][0] + 0.1e-5, f' Mean Diff: {mean_diff_visc_250:.2f}%', color='k')
plt.text(pressure[0]/1e5, visc_data[300][0] + 0.1e-5, f' Mean Diff: {mean_diff_visc_300:.2f}%', color='k')

plt.xlabel('Pressure (bar)')
plt.ylabel('Viscosity (Pa*s)')
plt.title('Viscosity (exp & calc) vs Pressure')
plt.legend()

plt.savefig("Viscosity_exp_viscosity.png")

# Thermal Conductivity
plt.figure(figsize=(8, 6))

plt.plot(pressure/1e5, cond_data[200], color='b', label='Cond thermopack 200 K')
plt.plot(pressure/1e5, thermal_conductivity_200_exp, color='r', label='Cond exp 200 K')
plt.plot(pressure/1e5, cond_data[250], color='b', label='Cond  thermopack 250 K')
plt.plot(pressure/1e5, thermal_conductivity_250_exp, color='r', label='Cond  exp 250 K')
plt.plot(pressure/1e5, cond_data[300], color='b', label='Cond  thermopack 300 K')
plt.plot(pressure/1e5, thermal_conductivity_300_exp, color='r', label='Cond  exp 300 K')

plt.text(pressure[0]/1e5, cond_data[200][0] + 0.001, f' Mean Diff: {mean_diff_cond_200:.2f}%', color='k')
plt.text(pressure[0]/1e5, cond_data[250][0] + 0.001, f' Mean Diff: {mean_diff_cond_250:.2f}%', color='k')
plt.text(pressure[0]/1e5, cond_data[300][0] + 0.001, f' Mean Diff: {mean_diff_cond_300:.2f}%', color='k')

plt.xlabel('Pressure (bar)')
plt.ylabel('Conductivity [W / m K]')
plt.title('Conductivity (exp & calc) vs Pressure')
plt.legend()

plt.savefig("Conductivity_exp_viscosity.png")'''

# Plotting function
def plot_and_annotate(pressure, data_calculated, data_exp, property_name, mean_diff_calculated, temp, linestyle_calculated='-', linestyle_exp='--', marker_calculated='o', marker_exp='s', y_offset=0):
    plt.figure(figsize=(8, 6))

    colors = plt.cm.viridis(np.linspace(0, 1, len(temp)))

    for i, T in enumerate(temp):
        color = colors[i]
        plt.plot(pressure/1e5, data_calculated[T], linestyle=linestyle_calculated, marker=marker_calculated, color=color, label=f'{property_name} thermopack {property_name} at {T} K')
        plt.text(pressure[0]/1e5, data_calculated[T][0] + y_offset, f' Mean Diff: {mean_diff_calculated[i]:.2f}%', color='k')
        i += 1

    i = 0
    for data_ex in data_exp:
        color = colors[i]
        plt.plot(pressure/1e5, data_ex, linestyle=linestyle_exp, marker=marker_exp, color=color, label=f'{property_name} exp {property_name}')
        i += 1

    plt.xlabel('Pressure (bar)')
    plt.ylabel(f'{property_name}')
    plt.title(f'{property_name} (exp & calc) vs Pressure')
    plt.legend()

    plt.savefig(f"{property_name}.png")

# viscosity
plot_and_annotate(pressure, visc_data, visc_exp, 'Viscosity', mean_diff_visc, T_values, linestyle_calculated='-', linestyle_exp='--', marker_calculated='o', marker_exp='s', y_offset=0.1e-5)

# Thermal conductivity
plot_and_annotate(pressure, cond_data, cond_exp, 'Conductivity', mean_diff_cond, T_values, linestyle_calculated='-', linestyle_exp='--', marker_calculated='o', marker_exp='s', y_offset=0.001)


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