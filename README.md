# Fock State 1 Preparation using interleaved Displacement-SNAP gate sequence
## State Preparation
The pulse sequence for the state preparation runs as follow: 
1) A pulse with displacement amplitude α1 is applied the cavity to generate a superposition of Fock states. The value of α corresponds to the applied cavity voltage 
2) A sequence of pulses are applied on the qubit as follows: a qubit π pulse: Rπ y,fn --> a qubit π pulse selective on cavity Fock number: Rπ y+θ,fn. In the presented work, SNAP gate
is applied to the Fock state n = 0 with θ = π
3) We then apply another pulse with displacement amplitude α2 to get the Fock state |1⟩. Furthermore, we have derived a generalized form for the applied gate sequence
Displacement-SNAP-Displacement:

![Screenshot 2022-12-01 160600](https://user-images.githubusercontent.com/34755328/205087421-c73f5858-8788-4553-9722-f3d3561dab75.png)

## Determining lookup table for Voltage (Vin) to Displacement (α)
To generate a lookup table, we require data that correlates the applied Vin to the corresponding α in the cavity. Our CQED setup is operating in the Dispersive regime of the Jaynes-Cummings Hamiltonian. Hence, the readout
amplitude of the qubit’s frequency sweep encompasses the data regarding the Fock states and their possible superpositions (a coherent state) in the cavity. Therefore, we applied a set of input voltages to the cavity which
is equivalent to displacing it with a set of corresponding α’s. Based on the coherent state present in the cavity, the qubit’s readout data is expected to exhibit distinctive peaks spaced by the Dispersive frequency shift χ. The
measured data exhibited frequency spread peaks separated by χ = 1.4Mhz which agreed with our CQED setup’s specification. Furthermore, the amplitude peaks assumed a Gaussian-like distribution whose envelope resembled
a Poissonian-like distribution. It was an expected consequence because a coherent state displays a Poissonian-like distribution corresponding to the value of α.

Thus, we now have data that demonstrates a subtle correlation between Vin and α because our qubit is Dispersively coupled to the cavity. Furthermore, by fitting the ensemble of the Gaussian-like peaks to the discrete Poisson
distribution, we can determine the corresponding coherent state mean amplitude α. To fit the experimental data, we binned it accordingly to the known value of χ. Individual Gaussian fitting of the binned data will result in finer
confidence and prediction bounds than multiple peak Gaussian fitting of the readout amplitudes. We implemented the function fit_gaussian on the binned data from the QTT Python package. Furthermore, another advantage of binning is that it also discretizes the continuous measured data making it suitable for fitting a Poisson distribution which is a discrete probability
density function. We then normalized the fitted data based on the amplitudes of the Gaussian peaks to obtain the Poisson-like probability distribution. From the python package scipy.stats, we utilized the function poisson
to perform the Poissonian fitting. Following figure shows the results we obtained for a sample set of the measured data. 

![Figure_6](https://user-images.githubusercontent.com/34755328/205089309-015027c0-ebd2-41b7-9faa-0f1a55aa2486.png)

![Figure_11](https://user-images.githubusercontent.com/34755328/205089406-1f825c29-fe68-416f-ab4f-f057c38d4077.png)

## Tomography and Simulation of Fock state |1⟩
The Poissonian fitting followed by the Gaussian fitting of the binned data resulted in a lookup table with correspondence between Voltage (Vin) and the Displacement α
Table 1: Lookup Table
Vin(V) ---> α
 0.0   ---> 0.100
 0.1   ---> 0.223
 0.2   ---> 0.561
 0.3   ---> 1.135
 0.4   ---> 1.923
To validate, our results, we implemented the D(α2) SNAP (n,θ)D(α1) |0⟩ gate sequence with Vin1 = 0.3 V,θ = π, |n⟩ = |0⟩, Vin2 = −0.152 V in our CQED setup. The applied gate sequence is expected to yield a Fock state |1⟩. 
We performed Wigner Tomography to confirm the same. The output indicated Fock state |1⟩ characterized by its distinctive negativity in the Wigner function.
Assuming a linear relationship to hold between the parameters in Table 1, Vin,1 = 0.3 V will correspond to α1 = 1.135 and Vin,2 = −0.152 V will corresponds to α2 = −0.336. We validated the accuracy of the lookup table
by plugging these values in our (unoptimized gate sequence parameter) Qutip Simulation. The simulation results exhibit excellent conformity with the experimental data

![Wigner_tomography_reversed_cmap](https://user-images.githubusercontent.com/34755328/205091204-a6b6ba31-3629-4c2e-a989-0b1df4c54fbc.png)

![Wigner_tomography_simulation](https://user-images.githubusercontent.com/34755328/205091294-5f55cd81-5e49-43d4-a0b6-0dcb5645feb0.png)

## Optimizing Gate sequence parameters
The optimization script outputs all the possible optimal parameters that result in an implicit conditional fidelity >= 0.95 (this can be changed in the code). We have implemented least_squares() method from the SciPy library. This method was chosen as it gives the user the flexibility to adjust ftol [tolerance for change in cost function], xtol [tolerance for change independent variables/parameters], and gtol [tolerance for change in gradient norm] simultaneously. 
The gtol parameter has been disabled by default as it causes the optimization to terminate once the script finds a local minima. The other tolerances can be set by the user.

Sample:
Input(s) 𝛼1=2, 𝛼2=1, 𝜃=0.5𝜋, 𝑛=0;
Output(s): 𝛼1=−1.14, 𝛼2=0.58, 𝜃=𝜋, 𝑛=0; Fidelity: 0.981;
𝛼1=1.14, 𝛼2=−0.58, 𝜃=𝜋, 𝑛=0; Fidelity: 0.981
