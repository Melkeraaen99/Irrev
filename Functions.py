import Initialization
import matplotlib.pyplot as plt

def calculate_mean_percentage_difference(data_calculated, data_exp):
    percentage_difference = Initialization.np.abs((data_calculated - data_exp) / ((data_calculated + data_exp) / 2)) * 100
    mean_percentage_difference = Initialization.np.mean(percentage_difference)
    return mean_percentage_difference

# Log-scaled dimensionless viscosity
def dim_less(T, property_data):
    if T in Initialization.density:
        return Initialization.np.log((Initialization.np.array(Initialization.density[T])*Initialization.Avogadro)**(-2/3) * (Initialization.kb * T* Initialization.m_Ar)**(-0.5) * Initialization.np.array(property_data))
    else:
        return None 
    
# Plotting function for experimental vs calculated
def plot_and_annotate(pressure, data_calculated, data_exp, property_name, mean_diff_calculated, temp, linestyle_calculated='-', linestyle_exp='--', marker_calculated='o', marker_exp='s', y_offset=0):
    plt.figure(figsize=(8, 6))

    colors = plt.cm.viridis(Initialization.np.linspace(0, 1, len(temp)))

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

    plt.savefig(f"Plots/{property_name}.png")

# Plotting calculated vs residual entropy 
# Function only takes care of viscosity at the moment since dimless cond is not yet defined
def residual_plot(property_string):
    plt.figure(figsize=(8, 6))

    for T in Initialization.T_values:
        label = f'Temperature {T} K (Ideal Gas)'
        plt.plot(Initialization.entropy[T], dim_less(T, Initialization.visc_data[T]), 'o', label=label) # Change visc_data here

    # Add labels and a legend
    plt.xlabel(r'S$_{res}\,$ [ J$\cdot$mol$^{-1}\cdot$K$^{-1}$]')
    plt.ylabel(r'$\eta^*$')
    plt.title(f'Dimensionless {property_string}' + r', $\eta^*$, plotted against residual entropy, S$_{res}$, for Argon')
    plt.legend()

    # Show the plot
    plt.grid(True)
    plt.savefig(f"Plots/{property_string}_resiudal.png")
