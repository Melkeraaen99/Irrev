import Initialization
import Functions

thermal_conductivity_105_exp = Initialization.df_105['Therm. Cond. (W/m*K)'].to_numpy()
thermal_conductivity_110_exp = Initialization.df_110['Therm. Cond. (W/m*K)'].to_numpy()
thermal_conductivity_115_exp = Initialization.df_115['Therm. Cond. (W/m*K)'].to_numpy()
cond_exp = [thermal_conductivity_105_exp, thermal_conductivity_110_exp, thermal_conductivity_115_exp] 

mean_diff_cond_105 = Functions.calculate_mean_percentage_difference(Initialization.cond_data[161.48776760869768], thermal_conductivity_105_exp)
mean_diff_cond_110 = Functions.calculate_mean_percentage_difference(Initialization.cond_data[169.17766130434995], thermal_conductivity_110_exp)
mean_diff_cond_115 = Functions.calculate_mean_percentage_difference(Initialization.cond_data[176.8675550000022], thermal_conductivity_115_exp)
mean_diff_cond = [mean_diff_cond_105, mean_diff_cond_110, mean_diff_cond_115]

# Print them to see how close they are in values
'''print(Initialization.cp_dict)
print(Initialization.cp_exp)'''

# Thermal conductivity
Functions.plot_and_annotate(Initialization.pressure, Initialization.cond_data, cond_exp, 'Conductivity', mean_diff_cond, Initialization.T_values, linestyle_calculated='-', linestyle_exp='--', marker_calculated='o', marker_exp='s', y_offset=0.0001)
Functions.residual_plot("Conductivity", cond_exp)
Functions.residual_plot("Conductivity_exp", cond_exp)