from pykingas.MieKinGas import MieKinGas
import numpy as np
import pandas as pd

from thermopack.saftvrmie import saftvrmie
from scipy.constants import Boltzmann
from scipy.constants import R
from scipy.constants import Avogadro

# Create a MieKinGas instance
eos = saftvrmie('AR') 
kin_id = MieKinGas('AR', is_idealgas=True)
kin = MieKinGas('AR', use_eos=eos)  # AR, modeled with RET-Mie

'''print(f'sigma : {kin.sigma}')
print(f'epsilon : {kin.epsilon}')
print(f'mass : {kin.m}')'''

tc, vc, pc = eos.critical([1])

# Define temperature and pressure ranges
T_values = [1.05*tc, 1.1*tc, 1.15*tc]  # Kelvin
pressure_start = 1e5  # Pascal -> 1 bar
pressure_end = 50e5
pressure = np.linspace(pressure_start, pressure_end, 20)
x = [0.5, 0.5]  # Molar composition used for conductivity and viscosity in pykingas ()
z = [1]
m_Ar = 39.948*1.67377e-27 # kg
molar_mass_argon = 39.948e-3 # kg/mol
kb = Boltzmann # J/K

# Experimental data for viscosity and conductivity
df_105 = pd.read_csv('Data/nist_data_tc_1_05.csv', delimiter='\t')
df_110 = pd.read_csv('Data/nist_data_tc_1_10.csv', delimiter='\t')
df_115 = pd.read_csv('Data/nist_data_tc_1_15.csv', delimiter='\t')

# Experimental heat capacity
cp_105_exp = df_105['Cp (J/mol*K)'].to_numpy()
cp_110_exp = df_110['Cp (J/mol*K)'].to_numpy()
cp_115_exp = df_115['Cp (J/mol*K)'].to_numpy()
cp_exp = [cp_105_exp, cp_110_exp, cp_115_exp] # [J/mol*K]

# Experimental entropy
entro_105_exp = df_105['Entropy (J/mol*K)'].to_numpy()
entro_110_exp = df_110['Entropy (J/mol*K)'].to_numpy()
entro_115_exp = df_115['Entropy (J/mol*K)'].to_numpy()
entro_exp = [entro_105_exp, entro_110_exp, entro_115_exp] # Entropy (J/mol*K)

# Experimental data for viscosity and conductivity
df_200 = pd.read_csv('Data/nist_data_200.csv', delimiter='\t')
df_250 = pd.read_csv('Data/nist_data_250.csv', delimiter='\t')
df_300 = pd.read_csv('Data/nist_data_300.csv', delimiter='\t')

# Experimental heat capacity
cp_200_exp = df_200['Cp (J/mol*K)'].to_numpy()
cp_250_exp = df_250['Cp (J/mol*K)'].to_numpy()
cp_300_exp = df_300['Cp (J/mol*K)'].to_numpy()
cp_exp = [cp_200_exp, cp_250_exp, cp_300_exp] # [J/mol*K]

# Initialize lists to store conductivity data for different temperatures
cond_data = {T: [] for T in T_values} # [W/m*K]

# Initialize lists to store viscosity data for different temperatures
visc_data = {T: [] for T in T_values} # [Pa*s]

# Store Residual entropy with thermopack in entropy list
entropy = {T: [] for T in T_values}
entropy_id = {T: [] for T in T_values}
residual_entropy_exp = {T: [] for T in T_values}

# Density
density = {T: [] for T in T_values}

# Specific volume
spes_volume = {T: [] for T in T_values}

# Heat capacity for the different pressures
cp_dict = {T: [] for T in T_values} # [J/mol*K]

# Calculate basic properties for Argon
for T in T_values:
    i = 0
    for p in pressure:
        s_volume = eos.specific_volume(T, p, z, eos.VAPPH)
        spes_volume[T].append(s_volume[0])
        density[T].append((1/s_volume[0])*Avogadro)
        entro = eos.entropy_tv(T, s_volume[0], z, property_flag='R') # Flag R to calculate residual entropy
        j = eos.getcompindex('AR')
        entro_id = eos.idealentropysingle(T, p, j)
        entropy[T].append(entro)
        entropy_id[T].append(entro_id)
        _, Cp_vap = eos.enthalpy(T, p, x, eos.VAPPH, dhdt=True) # From thermopack documentation 
        cp_dict[T].append(Cp_vap)

        # Viscosity
        visc = kin.viscosity_tp(T, p, x, N=2)
        visc_data[T].append(visc)

        # Conductivity
        cond = kin.thermal_conductivity(T, spes_volume[T][i], x, N=2)
        cond_data[T].append(cond)
        i+=1

## Her må vi beregne residual entropy for experimental data