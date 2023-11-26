import Initialization
import matplotlib.pyplot as plt

def calculate_mean_percentage_difference(data_calculated, data_exp):
    percentage_difference = Initialization.np.abs((data_calculated - data_exp) / ((data_calculated + data_exp) / 2)) * 100
    mean_percentage_difference = Initialization.np.mean(percentage_difference)
    return mean_percentage_difference

# Log-scaled dimensionless viscosity (Add some a or b scaling-values from article?)
def dim_less_visc(T, visc_data):
    return ((Initialization.np.array(Initialization.density[T]))**(-2/3) * (Initialization.kb * T* Initialization.m_Ar)**(-0.5) * Initialization.np.array(visc_data))

'''# dimensionless conductivity, the Prandtl number: (heat capacity * viscosity) / thermal conductivity
# cp_dict [J/(mol*K)]
# visc_data [Pa*s] -> [kg/(m*s)]
# cond_data [W/m*K] -> [J/(s*m*K)]
def dim_less_cond(T, visc_data, cond_data):
    dim_less_1 = ((Initialization.np.array(Initialization.cp_dict[T])*Initialization.np.array(visc_data)) / (Initialization.np.array(cond_data))) # kg/mol
    dim_less_2 = dim_less_1 / Initialization.molar_mass_argon # dimless
    return Initialization.np.log(dim_less_2)'''

def dim_less_cond(T, cond_data):
    dim_less = Initialization.np.log((Initialization.np.array(Initialization.density[T]))**(-2/3) * ((Initialization.kb * T)/Initialization.m_Ar)**(-0.5) * Initialization.np.array(cond_data)/Initialization.kb)
    return dim_less #/T

def dim_less_entropy(resiudal_entropy):
    return -(1/Initialization.R)*Initialization.np.array(resiudal_entropy)

# Plotting function for experimental vs calculated
def plot_and_annotate(pressure, data_calculated, data_exp, property_name, mean_diff_calculated, temp, linestyle_calculated='-', linestyle_exp='--', marker_calculated='o', marker_exp='s', y_offset=0):
    plt.figure(figsize=(8, 6))

    colors = plt.cm.viridis(Initialization.np.linspace(0, 1, len(temp)))

    for i, T in enumerate(temp):
        color = colors[i]
        plt.plot(pressure/1e5, data_calculated[T], linestyle=linestyle_calculated, marker=marker_calculated, color=color, label=f'{property_name} thermopack at {T:.2f} K')
        plt.text(pressure[15]/1e5, data_calculated[T][0] - y_offset, f' Mean Diff: ', color=color)
        plt.text(pressure[15]/1e5 + 7 , data_calculated[T][0] - y_offset, f'{mean_diff_calculated[i]:.2f}%', color='k')
        i += 1

    i = 0
    for data_ex in data_exp:
        color = colors[i]
        plt.plot(pressure/1e5, data_ex, linestyle=linestyle_exp, marker=marker_exp, color=color, label=f'{property_name} exp at {Initialization.T_values[i]:.2f} K')
        i += 1

    plt.xlabel('Pressure (bar)')
    plt.ylabel(f'{property_name}')
    plt.title(f'{property_name} vs Pressure')
    plt.legend()

    plt.savefig(f"Plots/Calc_VS_Exp_{property_name}.png")

# Plotting calculated vs residual entropy 
# Function only takes care of viscosity at the moment since dimless cond is not yet defined
def residual_plot(property_string, data_exp):
    plt.figure(figsize=(8, 6))

    i = 0
    for T in Initialization.T_values:
        label = f'Temperature {T:.2f} K'
        if property_string == "Viscosity":
            plt.plot(1/dim_less_entropy(Initialization.entropy[T]), dim_less_visc(T, Initialization.visc_data[T]), 'o', label=label)
        elif property_string == "Conductivity":
            plt.plot(1/dim_less_entropy(Initialization.entropy[T]), dim_less_cond(T, Initialization.cond_data[T]), 'o', label=label)
        elif property_string == "Viscosity_exp":
            plt.plot(1/dim_less_entropy(Initialization.entro_exp[i]), dim_less_visc(T, Initialization.np.array(data_exp[i])), 'o', label=label)
        elif property_string == "Conductivity_exp":
            plt.plot(1/dim_less_entropy(Initialization.entro_exp[i]), dim_less_cond(T, Initialization.np.array(data_exp[i])), 'o', label=label)
        else:
            print("Invalid property string! \n Set property string to either Viscosity or Conductivity")
            return None
        i += 1

    # Add labels correctly for viscosity and conductivity
    if property_string == "Viscosity":
        plt.xlabel(r'S$_{res}^*\,$') # put in this if not dimless res entro: [ J$\tcdo$mol$^{-1}\cdot$K$^{-1}$] 
        plt.ylabel(r'$\eta^*$')
        plt.title(f'Dimensionless {property_string}' + r' plotted against dimensionless residual entropy for Argon')
        plt.legend()
    elif property_string == "Conductivity":
        plt.xlabel(r'S$_{res}^*\,$')
        plt.ylabel(r'$k^*$')
        plt.title(f'Dimensionless {property_string}' + r' plotted against dimensionless residual entropy for Argon')
        plt.legend()
    elif property_string == "Viscosity_exp":
        plt.xlabel(r'S$_{res}^*\,$')
        plt.ylabel(r'$k^*$')
        plt.title(f'Dimensionless {property_string}' + r' plotted against dimensionless residual entropy for Argon')
        plt.legend()
    elif property_string == "Conductivity_exp":
        plt.xlabel(r'S$_{res}^*\,$')
        plt.ylabel(r'$k^*$')
        plt.title(f'Dimensionless {property_string}' + r' plotted against dimensionless residual entropy for Argon')
        plt.legend()

    # Show the plot
    plt.grid(True)
    plt.savefig(f"Plots/{property_string}_residual.png")