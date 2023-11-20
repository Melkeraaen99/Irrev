import Initialization
import Functions

thermal_conductivity_200_exp = Initialization.df_200['Therm. Cond. (W/m*K)'].to_numpy()
thermal_conductivity_250_exp = Initialization.df_250['Therm. Cond. (W/m*K)'].to_numpy()
thermal_conductivity_300_exp = Initialization.df_300['Therm. Cond. (W/m*K)'].to_numpy()
cond_exp = [thermal_conductivity_200_exp, thermal_conductivity_250_exp, thermal_conductivity_300_exp]

mean_diff_cond_200 = Functions.calculate_mean_percentage_difference(Initialization.cond_data[200], thermal_conductivity_200_exp)
mean_diff_cond_250 = Functions.calculate_mean_percentage_difference(Initialization.cond_data[250], thermal_conductivity_250_exp)
mean_diff_cond_300 = Functions.calculate_mean_percentage_difference(Initialization.cond_data[300], thermal_conductivity_300_exp)
mean_diff_cond = [mean_diff_cond_200, mean_diff_cond_250, mean_diff_cond_300]

# Print them to see how close they are in values
'''print(Initialization.cp_dict)
print(Initialization.cp_exp)'''

# Thermal conductivity
Functions.plot_and_annotate(Initialization.pressure, Initialization.cond_data, cond_exp, 'Conductivity', mean_diff_cond, Initialization.T_values, linestyle_calculated='-', linestyle_exp='--', marker_calculated='o', marker_exp='s', y_offset=0.001)
Functions.residual_plot("Conductivity", cond_exp)
Functions.residual_plot("Conductivity_exp", cond_exp)