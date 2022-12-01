import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import curve_fit
from qtt.algorithms.functions import gauss_ramsey, fit_gauss_ramsey, plot_gauss_ramsey_fit

load_path_plus = '/home/sushanth/Wigner_Lab_Group_3/Group3/2022-10-06-0948_Ramsey_presto.npz' #control freq: 100MHz+200KHz
load_path_minus = '/home/sushanth/Wigner_Lab_Group_3/Group3/2022-10-06-0950_Ramsey_presto.npz' #control freq: 100MHz+200KHz

with np.load(load_path_plus) as npzfile:
    data_plus = npzfile['demod_I_p']
    x_axis = npzfile['time_LUT']
    nr = npzfile['nr']
    frequency_qubit_pulse = npzfile['control_freq'] + npzfile['drive_LO']
    control_freq=npzfile['control_freq']
    drive_LO=npzfile['drive_LO']
MW_plus=frequency_qubit_pulse

with np.load(load_path_minus) as npzfile:
    data_minus = npzfile['demod_I_p']
    x_axis = npzfile['time_LUT']
    nr = npzfile['nr']
    frequency_qubit_pulse = npzfile['control_freq'] + npzfile['drive_LO']
    control_freq=npzfile['control_freq']
    drive_LO=npzfile['drive_LO']
MW_minus=frequency_qubit_pulse

print("MW_plus:%f MW_plus:%f"%(MW_plus,MW_minus))
print("2*delta_f: %f"%(MW_plus-MW_minus))

fig, ax = plt.subplots()
ax.plot(x_axis,data_plus,ls='-', marker='', color='red', label='MW_plus data')
ax.plot(x_axis,data_minus,ls='-', marker='', color='green', label='MW_minus data')
plt.legend()
#plt.show()

#Scipy Curve Fit
def fit_func_plus(t,A,q_freq,T2,B):
    return A*np.cos(2*np.pi*(q_freq)*t)*np.exp(-1*(t/T2)**2)+B

def fit_func_minus(t,A,q_freq,T2,B):
    return A*np.cos(2*np.pi*(q_freq)*t)*np.exp(-1*(t/T2)**2)+B

popt_plus, pcov_plus=curve_fit(fit_func_plus,x_axis,data_plus)
print("\nP_opt_Plus:")
print(popt_plus)
fig, ax = plt.subplots()
ax.plot(x_axis,data_plus,ls='-', marker='', color='red', label='MW_plus data')
ax.plot(x_axis,fit_func_plus(x_axis,*popt_plus),ls='-', marker='', color='blue', label='Fitting data')
plt.legend()
plt.show()

popt_minus, pcov_minus=curve_fit(fit_func_minus,x_axis,data_minus)
print("\nP_opt_Minus:")
print(popt_minus)
fig, ax = plt.subplots()
ax.plot(x_axis,data_minus,ls='-', marker='', color='red', label='MW_minus data')
ax.plot(x_axis,fit_func_minus(x_axis,*popt_minus),ls='-', marker='', color='blue', label='Fitting data')
plt.legend()
plt.show()

#x_axis=1e5*x_axis

#TU Delft QTT
fit_parameters_plus, _ = fit_gauss_ramsey(x_axis, data_plus)
print("\nP_opt_Plus:")
print(fit_parameters_plus)
freq_fit = abs(fit_parameters_plus[2]*1e-6)
t2star_fit = abs(fit_parameters_plus[1]*1e6)

plt.figure()
plot_gauss_ramsey_fit(x_axis, data_plus, fit_parameters_plus, fig=1)
plt.show()

fit_parameters_minus, _ = fit_gauss_ramsey(x_axis, data_minus)
print("\nP_opt_Minus:")
print(fit_parameters_minus)
freq_fit = abs(fit_parameters_minus[2]*1e-6)
t2star_fit = abs(fit_parameters_minus[1]*1e6)

plt.figure()
plot_gauss_ramsey_fit(x_axis, data_minus, fit_parameters_minus, fig=1)
plt.show()