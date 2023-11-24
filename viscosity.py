import Initialization
import Functions

viscosity_105_exp = Initialization.df_105['Viscosity (Pa*s)'].to_numpy()
viscosity_110_exp = Initialization.df_110['Viscosity (Pa*s)'].to_numpy()
viscosity_115_exp = Initialization.df_115['Viscosity (Pa*s)'].to_numpy()
visc_exp = [viscosity_105_exp, viscosity_110_exp, viscosity_115_exp] 

mean_diff_visc_105 = Functions.calculate_mean_percentage_difference(Initialization.visc_data[161.48776760869768], viscosity_105_exp)
mean_diff_visc_110 = Functions.calculate_mean_percentage_difference(Initialization.visc_data[169.17766130434995], viscosity_110_exp)
mean_diff_visc_115 = Functions.calculate_mean_percentage_difference(Initialization.visc_data[176.8675550000022], viscosity_115_exp)
mean_diff_visc = [mean_diff_visc_105, mean_diff_visc_110, mean_diff_visc_115]

# Initialize plots
Functions.plot_and_annotate(Initialization.pressure, Initialization.visc_data, visc_exp, 'Viscosity', mean_diff_visc, Initialization.T_values, linestyle_calculated='-', linestyle_exp='--', marker_calculated='o', marker_exp='s', y_offset=0.1e-6)
Functions.residual_plot("Viscosity", Initialization.visc_data)
Functions.residual_plot("Viscosity_exp", visc_exp)