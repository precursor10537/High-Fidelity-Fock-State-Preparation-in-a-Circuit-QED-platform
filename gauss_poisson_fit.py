from statistics import mean
import matplotlib.pyplot as plt
import numpy as np
from qtt.algorithms.functions import gauss_ramsey, fit_gauss_ramsey, plot_gauss_ramsey_fit
from qtt.algorithms.fitting import fit_gaussian
from scipy.optimize import curve_fit
from scipy.special import factorial
from scipy.stats import poisson
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator

#Model for Gaussian function
#$$y = offset + amplitude * np.exp(-(1/2)*(x-mean)^2/s^2)$$
#Args:
#    x (array): data points
#    mean, std, amplitude, offset: parameters
#Returns:
#    y (array)
def gaussian(x, mean, std, amplitude, offset):
    return offset + amplitude * np.exp(-(1 / 2) * (x - mean) ** 2 / std ** 2)

def fit_function(k, lamb):
    '''poisson function, parameter lamb is the fit parameter'''
    return poisson.pmf(k, lamb)


load_path = '/home/sushanth/Wigner_Lab_Group_3/Common_groups1_2_3/2022-10-09-1142_cavity_population.npz'
#load_path='2022-10-09-1142_cavity_population.npz'
with np.load(load_path) as npzfile:
    raw_data = npzfile['raw_data']
    x_axis = npzfile['x_axis']+npzfile['drive_LO']
    y_axis = npzfile['disp_amp_LUT']
    nr = npzfile['nr']
    nr1 = npzfile['nr1']
    readout_freq = npzfile['readout_freq']
    t_template_match = npzfile['t_template_match']
    phase = npzfile['phase']
    p2_phase_shift = npzfile['p2_phase_shift']
    start = npzfile['start']
    end = npzfile['end']

data = np.zeros((nr1,nr),dtype=np.complex_)
for ii in range(nr1):
    _result = raw_data[ii]
    demod_I_p = []
    demod_Q_p = []
    for i in range(nr):
        carrier_m1_p1 = np.cos(2 * np.pi * readout_freq * t_template_match + np.pi * phase)
        carrier_m1_p2 = np.sin(2 * np.pi * readout_freq * t_template_match + np.pi * (phase + p2_phase_shift))
        carrier_m2_p1 = -np.sin(2 * np.pi * readout_freq * t_template_match + np.pi * phase)
        carrier_m2_p2 = np.cos(2 * np.pi * readout_freq * t_template_match + np.pi * (phase + p2_phase_shift))
        demod_I_p.append(
            np.sum(carrier_m1_p1 * _result[i, 0, start:end].real + carrier_m1_p2 * _result[i, 0, start:end].imag))
        demod_Q_p.append(
            np.sum(carrier_m2_p1 * _result[i, 0, start:end].real + carrier_m2_p2 * _result[i, 0, start:end].imag))

    demod_I_p = np.array(demod_I_p)
    demod_Q_p = np.array(demod_Q_p)
    data[ii] = demod_I_p + 1j*demod_Q_p

if nr1 == 1:
    plt.plot(x_axis, data.real[0])
    plt.ylabel('Readout voltage Voltage')
    plt.xlabel('Frequencies')
else:
    plt.pcolor(x_axis,y_axis,data.real,cmap='RdBu_r')
    plt.colorbar()
    plt.ylabel('Applied Voltage')
    plt.xlabel('Frequencies')

start=x_axis[len(x_axis)-1]
x_lines=[start]
chi=2.2e+6-1.5e+6
while 1>0:
    start=start-chi
    if start<=x_axis[0]:
        x_lines.append(x_axis[0])
        break
    x_lines.append(start)
#print(x_lines)

segment=[]
for i in range(0,len(x_lines),2):
    segment.append((x_lines[i]+x_lines[i+1])/2)
#print(segment)

count=0
crop_x_axis=[]
for i in x_axis:
    if i<=segment[0] and i>=segment[len(segment)-1]:
        crop_x_axis.append(i)

#print(x_axis)
#print('\n')
#print(crop_x_axis)
#print(len(crop_x_axis))
crop_data=[[],[],[],[],[]]
for j in range(5):
    for i in range(195):
        crop_data[j].append(data.real[j][i])
segment_data=[[],[],[],[],[]]
for j in range(5):
    for i in range(2,len(crop_data[j])-1,1):
        segment_data[j].append(crop_data[j][i])


#print(len(segment_data))

#for i in range(5):
#    fit_params = fit_gaussian(x_axis, data.real[i])
    #print(fit_params[0])
#    y = gaussian(x_axis, *fit_params[0])
#    plt.figure()
#    plt.scatter(x_axis, data.real[i])
#    plt.plot(x_axis, y,color='green')
    #plt.vlines(x_lines, -0.5, 1, linestyles='dashed', colors='red')
#    plt.vlines(segment, -0.5, 1, linestyles='dashed', colors='green')
#    plt.ylabel('Read out voltages')
#    plt.xlabel('Frequencies')
#plt.show()


#TRUE PLOTS def gaussian(x, mean, std, amplitude, offset):
for i in range(5):
    fit_params = fit_gaussian(crop_x_axis, segment_data[i])
    #print(fit_params[0])
    y = gaussian(crop_x_axis, *fit_params[0])
    #plt.figure()
    #plt.scatter(crop_x_axis, segment_data[i])
    #plt.plot(crop_x_axis, y,color='green')
    #plt.vlines(x_lines, -0.5, 1, linestyles='dashed', colors='red')
    #plt.vlines(segment, -0.5, 1, linestyles='dashed', colors='green')
    #plt.ylabel('Read out voltages')
    #plt.xlabel('Frequencies')
#plt.show()

new_x_axis=[]
for k in range(len(segment)-1):
    temp=[]
    start=int(segment[k])
    #print(start)
    stop=int(segment[k+1])-50000
    #print(stop)
    #print('\n')
    for z in range(start,stop,-50000):
        #print(z)
        temp.append(float(z))
    if k==0:
        new_x_axis.append(temp)
    else:
        temp.pop(0)
        new_x_axis.append(temp)

#print(len(new_x_axis))
#print(len(new_x_axis[0]))
#print(len(new_x_axis[1]))
#print(len(new_x_axis[2]))
#print(len(new_x_axis[3]))
#print(len(new_x_axis[4]))
#print(len(new_x_axis[5]))
#print(len(new_x_axis[6]))

#print('\n')
#print(crop_x_axis)
#print('\n')
#print(new_x_axis)
#length=0
#for i in range(len(new_x_axis)):
#    length=len(new_x_axis[i])+length
#print(length)

new_y_data=[]
for i in range(5):
    new_data=[]
    new_data.append(segment_data[i][0:29])
    new_data.append(segment_data[i][29:57])
    new_data.append(segment_data[i][57:85])
    new_data.append(segment_data[i][85:113])
    new_data.append(segment_data[i][113:141])
    new_data.append(segment_data[i][141:169])
    new_data.append(segment_data[i][169:192]) 
    new_y_data.append(new_data)  

#print(new_y_data[0])
#length=0
#for i in range(len(new_y_data[0])):
#    length=len(new_y_data[0][i])+length
#print(length) 

#for k in range(5):
#    for i in range(7):
#        print(len(new_x_axis[i]))
#        print(len(new_y_data[k][i]))
#        print('\n')

fitted_params=[]
y=[]
y_bin=[]
for i in range(5):
    fitted_params.append([])
    y.append([])
    y_bin.append([])
    for j in range(7):
            fit_params = fit_gaussian(new_x_axis[j], new_y_data[i][j])
            #print(fit_params[0])
            fitted_params[i].append(fit_params[0])
            g = gaussian(new_x_axis[j], *fit_params[0])
            y_bin[i].append(g)
            for c in range(len(g)):
                y[i].append(g[c])

#print(len(y[0]))
#print(y_bin[0][1])
#print(len(y_bin[0]))
for i in range(5):
    #fit_params = fit_gaussian(crop_x_axis, segment_data[i])
    #print(fit_params[0])
    #y = gaussian(crop_x_axis, *fit_params[0])
    plt.figure()
    plt.scatter(crop_x_axis, segment_data[i],label='Experimental Data')
    plt.plot(crop_x_axis, y[i],color='red',label='Fitted Data',linewidth=2.0)
    #plt.vlines(x_lines, -0.5, 1, linestyles='dashed', colors='red')
    plt.vlines(segment, -0.3, 0.15, linestyles='dashed', colors='green')
    plt.ylabel('Readout Amplitudes (a.u.)')
    plt.xlabel('Frequencies')
    plt.legend()
#plt.show()

#print(fitted_params[0])
#plt.show()

pdf=[]
for i in range(5):
    pdf.append([])
    for j in range(7):
        pdf[i].append(np.max(y_bin[i][j]))

#pdf[0][5]=pdf[0][5]+0.05
#print('\n')
#print(pdf[0])
pdf_abs=[]
for i in range(5):
    pdf_abs.append([])
    for j in range(7):
        pdf_abs[i].append(-1*(np.min(pdf[i]))+pdf[i][j])
#print(y_bin[0][0])
#print('\n')
#print(pdf_abs[0])
x_discrete=[6,5,4,3,2,1,0]
normalizing_factor=[]
for i in range(5):
    normalizing_factor.append(1/np.sum(pdf_abs[i]))
    pdf_abs[i]=np.multiply(normalizing_factor[i],pdf_abs[i])

#print(normalizing_factor)
#print(pdf_abs[0])
#for i in range(5):
#    print(np.sum(pdf_abs[i]))
#for i in range(5):
#    plt.figure()
#    plt.bar(x_discrete,pdf_abs[i])
#plt.show()

poisson_fit_param=[]
for i in range(5):
    parameters, cov_matrix = curve_fit(fit_function, x_discrete,pdf_abs[i])
    poisson_fit_param.append((parameters))
    plt.figure()
    plt.bar(x_discrete,pdf_abs[i])
    plt.plot(x_discrete,fit_function(x_discrete,*parameters),marker='o',color='red',label='Fitted Data')
    plt.ylabel('Probability')
    plt.xlabel('No. of events')
    plt.legend()
print(poisson_fit_param)
print(y_axis)
plt.show()







  