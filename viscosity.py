import Initialization
import Functions

viscosity_200_exp = Initialization.df_200['Viscosity (Pa*s)'].to_numpy()
viscosity_250_exp = Initialization.df_250['Viscosity (Pa*s)'].to_numpy()
viscosity_300_exp = Initialization.df_300['Viscosity (Pa*s)'].to_numpy()
visc_exp = [viscosity_200_exp, viscosity_250_exp, viscosity_300_exp] 

'''mean_diff_visc_200 = Functions.calculate_mean_percentage_difference(Initialization.visc_data[200], viscosity_200_exp)
mean_diff_visc_250 = Functions.calculate_mean_percentage_difference(Initialization.visc_data[250], viscosity_250_exp)
mean_diff_visc_300 = Functions.calculate_mean_percentage_difference(Initialization.visc_data[300], viscosity_300_exp)
mean_diff_visc = [mean_diff_visc_200, mean_diff_visc_250, mean_diff_visc_300]
'''
# Initialize plots
#Functions.plot_and_annotate(Initialization.pressure, Initialization.visc_data, visc_exp, 'Viscosity', mean_diff_visc, Initialization.T_values, linestyle_calculated='-', linestyle_exp='--', marker_calculated='o', marker_exp='s', y_offset=0.1e-5)
Functions.residual_plot("Viscosity", visc_exp)
Functions.residual_plot("Viscosity_exp", visc_exp)