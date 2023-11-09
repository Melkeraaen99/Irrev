from pykingas.MieKinGas import MieKinGas
import numpy as np
import pandas as pd

from thermopack.saftvrmie import saftvrmie
from scipy.constants import Boltzmann
from scipy.constants import R
from scipy.constants import Avogadro

# Create a MieKinGas instance
eos = saftvrmie('AR') 
kin = MieKinGas('AR', use_eos=eos)  # AR, modeled with RET-Mie

# Define temperature and pressure ranges
T_values = [200, 250, 300]  # Kelvin
pressure_start = 1e5  # Pascal -> 1 bar
pressure_end = 100e5
pressure_increment = 5e5
pressure = np.arange(pressure_start, pressure_end + 1, pressure_increment)
x = [0.5, 0.5]  # Molar composition used for conductivity and viscosity in pykingas ()
z = [1]
m_Ar = 39.948*1.67377e-27 # kg
molar_mass_argon = 39.948e-3 # kg/mol
kb = Boltzmann # J/K

# Experimental data for viscosity and conductivity
df_200 = pd.read_csv('Data/nist_data_200.csv', delimiter='\t')
df_250 = pd.read_csv('Data/nist_data_250.csv', delimiter='\t')
df_300 = pd.read_csv('Data/nist_data_300.csv', delimiter='\t')

# Initialize lists to store conductivity data for different temperatures
cond_data = {T: [] for T in T_values}

# Initialize lists to store viscosity data for different temperatures
visc_data = {T: [] for T in T_values}

# Store Residual entropy with thermopack in entropy list
entropy = {T: [] for T in T_values}

# Density
density = {T: [] for T in T_values}

# Specific volume
spes_volume = {T: [] for T in T_values}

# Calculate basic properties for Argon
for T in T_values:
    for p in pressure:
        s_volume = eos.specific_volume(T, p, z, eos.VAPPH)
        spes_volume[T].append(s_volume[0])
        density[T].append(1/s_volume[0])
        entro = eos.entropy_tv(T, s_volume[0], z, property_flag='R') # Flag R to calculate residual entropy
        entropy[T].append(entro) 
