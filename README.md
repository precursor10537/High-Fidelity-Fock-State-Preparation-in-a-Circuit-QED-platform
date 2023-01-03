# Preparing fock state 1 using interleaved Displacement-SNAP-Displacement gate sequence

Generation of a variety of high-fidelity Wignernegative states were experimentally demonstrated by the application of interleaved displacement and SNAP gates [1]. Following a similar approach, we prepare and validate the Fock state |1⟩ in the shown cavity quantum electrodynamics setup.

![image](https://user-images.githubusercontent.com/34755328/210422807-f341f762-0969-4e91-8ca2-568f0ec6edbf.png)

Firstly, the Ramsey Fringes experiment was carried out to probe the qubit drifting frequency and to measure its dephasing time T2*. 
After determining these parameters, the Fock state |1⟩ was prepared by applying Displacement-SNAP-Displacement gate sequence on the vacuum state |0⟩. Subsequently, Wigner tomography was performed on the prepared state. The obtained experimental data was utilized to generate a lookup table that maps the applied input voltages (Vin) to their corresponding applied displacement (α) on the cavity. Furthermore, a custom script is presented to optimize the applied gate sequence parameters.  

## Ramsey Fringes Experiment: Determining Qubit Drift Frequency (f_drift) and T2*

The Ramsey experiment was performed twice with fixed incremental delay τ + dt at two different frequencies f_Ramsey+ = f_estimate + Δf and f_Ramsey− = festimate − Δf, where Δf is known, to track the qubit resonant frequency (fqubit). 
Suppose f_estimate = f_qubit, the output signals would oscillate at the known frequency Δf.
However, if the qubit drifts (by an arbitrary f_drift), then f_estimate ̸= f_qubit, and the output signals of the Ramsey+ and Ramsey− experiments would oscillate at the unknown frequencies fout+ and fout−, respectively. 
The qubit drift frequency can then be calculated by
f_drift = (|f_out+ − Δf| + |f_out− − Δf|)/2
Furthermore, if fout+ < Δf, then f_estimate = f_qubit + f_drift, and if f_out+ > Δf, then f_estimate = f_qubit − f_drift.
The inverse of these conditions would apply to f_out−.
The constant resulting in the fitting of the exponentially decaying amplitudes of both the output signals would inherently give us T2∗.
A custom code with curve_fit function from the Python package SciPy was implemented to fit the experimental data. Additionally, another custom code was implemented utilizing the function fit_gauss_ramsey from the Python package Quantum Technology Toolbox (QTT) developed by TU Delft [3]. 
Similar results were obtained in both implementations. 
The frequency at which these signals oscillate is close to the known Δf = 200KHz, and the drift frequency was calculated to be f_drift = 0.83KHz, with an average estimate of T2* to be 41.59μs. 
As f_qubit is operating at GHz range, the calculated f_drift is inconsequential indicating that the qubit has drifted nominally from its resonant frequency (f_qubit).
Furthermore, to validate the results, Fourier Transform was applied to the Ramsey experiment output signals. 
The resulting frequency components were in agreement with the results obtained through the fitted data.

![image](https://user-images.githubusercontent.com/34755328/210425383-db9ed0c5-aebc-4664-92d6-0baa48d424c1.png)

## Jaynes-Cummings Model: Dispersive Regime

Jaynes-Cummings Hamiltonian can be derived from the Rabi Hamiltonian by using the electric field operator

![image](https://user-images.githubusercontent.com/34755328/210426210-c35168ea-468f-42f9-b896-a15487c7fc77.png)

where g is the coupling strength, σ+ = |1⟩ ⟨0| and σ− = |0⟩ ⟨1| are atomic raising and lowering operators, respectively. For a large detuning (Δ = ω0 −ωc ≫ g) between the atomic transition frequency (ω0) and the cavity frequency (ωc), the Jaynes-Cumming Hamiltonian in the dispersive regime becomes

![image](https://user-images.githubusercontent.com/34755328/210426286-ae240a26-40c9-4465-8dd2-9c2fc5efddf0.png)

where χ = g2/Δ is the dispersive frequency shift. 
Therefore, the effective frequency of the cavity (ωeff = ωc−χˆσz) is dependent on whether the qubit is in the ground state or excited state. Furthermore, the interaction between the cavity and the qubit in the dispersive Hamiltonian is shown by the term that is proportional to σza†a. 
Hence, the qubit frequency changes depending on the number of photons present in the cavity whose energy states are separated by χ (i.e. fn = nχ where n is the Fock state).

## SNAP
Selective arbitrary phase gate (SNAP) applies an arbitrary phase θ, conditioned on the |n⟩th Fock state through its projector [2]. 
The SNAP gate is experimentally realized by consecutive application of a qubit π pulse: Rπ,y,fn, and a selective qubit π pulse on the cavity Fock number: Rπ,y+θ,fn.
The SNAP gate is mathematically described by

![image](https://user-images.githubusercontent.com/34755328/210426556-90a20907-f515-4826-a83a-5d130532a36c.png)

## State Preparation
The pulse sequence for the state preparation runs as follows: 
1) A pulse with displacement amplitude α1 is applied on the cavity to generate a superposition of Fock states. The value of α corresponds to the voltage of the applied
pulse. 
2) SNAP gate is applied to the Fock state n = 0 with θ = π
3) Another pulse with displacement amplitude α2 is then applied to get the Fock state |1⟩.

Furthermore, we have derived a generalized form for the applied gate sequence Displacement-SNAP-Displacement as follows:

![Screenshot 2022-12-01 160600](https://user-images.githubusercontent.com/34755328/205087421-c73f5858-8788-4553-9722-f3d3561dab75.png)

## Wigner Tomography
The pulse sequence for the Wigner tomography runs as follows: 
1) A pulse is applied to the cavity to displace the prepared state. 
2) Two pulses separated by a time 1/2χ are applied to the qubit (similar to the Ramsey fringes experiment) to obtain the parity state of the photon in the cavity mode. It is possible because each Fock state rotates with a different frequency during the free evolution.
3) A pulse is applied to the Readout line to measure the resultant qubit states. The Wigner function can be calculated by

W(α) = (2/π)Tr[D†(α)ρD(α)P], where D(α) = exp{αa†−α∗a} is the Displacement operator and P = exp{iπa†a} is the Parity state operator. 

## Determining lookup table for Voltage (Vin) to Displacement (α)
Data that correlate the applied Vin to the corresponding α in the cavity is required to generate a lookup table.
Our CQED setup is operating in the Dispersive regime of the Jaynes-Cummings Hamiltonian. Hence, the readout
amplitude of the qubit’s frequency sweep encompasses the data regarding the Fock states and their possible
superpositions (a coherent state) in the cavity. Therefore, a set of input voltages were applied to the cavity that
is equivalent to displacing it with a set of corresponding α’s. Based on the coherent state present in the cavity,
the qubit’s readout data is expected to exhibit distinctive peaks spaced by the Dispersive frequency shift χ. The
measured data exhibited frequency spread peaks separated by χ = 1.4Mhz that agreed with our CQED setup’s
specification. Furthermore, the amplitude peaks assumed a Gaussian-like distribution whose envelope resembled
a Poissonian-like distribution. It was an expected consequence because a coherent state displays a Poissonian-like
distribution corresponding to the value of α given by

![image](https://user-images.githubusercontent.com/34755328/210428463-95ef7c94-a258-4fa4-a2b4-764a199aa9f8.png)

Thus, the collected data demonstrates a subtle correlation between Vin and α because the qubit is Dispersively
coupled to the cavity. Furthermore, by fitting the ensemble of the Gaussian-like peaks to the discrete Poisson
distribution, the corresponding coherent state mean amplitude α can be determined. The experimental data was
binned accordingly to the known value of χ. Individual Gaussian fitting of the binned data resulted in finer
confidence and prediction bounds than the multiple peak Gaussian fitting of the readout amplitudes. The function
fit_gaussian was implemented on the binned data from the QTT Python package.  
Furthermore, another advantage of binning was that it alsodiscretized the measured data making it suitable for fitting a Poisson distribution which is a discrete probability density function. The fitted data was then normalized based on the amplitudes of the Gaussian peaks to obtain the
Poisson-like probability distribution. From the python package scipy.stats, the function poisson was utilized
to perform the Poissonian fitting. The figures show the results obtained for a sample set of the measured data. 

![image](https://user-images.githubusercontent.com/34755328/210428876-6bdaa46d-5e49-49e3-807e-7f93de00fbda.png)

## Preparation, Tomography, and Simulation of the Fock state |1⟩

The Poissonian fitting followed by the Gaussian fitting of the binned data resulted in a lookup table with correspondence between Voltage (Vin) and the Displacement α
To validate the results, the D(α2)SNAP(n,θ)D(α1)|0⟩ gate sequence with
Vin1 = 0.3V,θ = π, |n⟩ = |0⟩, Vin2 = −0.152V was implemented in our CQED
setup. The applied gate sequence is expected to yield a Fock state |1⟩.
Wigner Tomography was performed to confirm the same. 
The output indicated Fock state |1⟩ characterized by its distinctive
negativity in the Wigner function. Assuming a
linear relationship to hold between the parameters in the Table, Vin1 = 0.3V
will correspond to α1 = 1.135, and Vin2 = −0.152V will corresponds to
α2 = −0.336.

![image](https://user-images.githubusercontent.com/34755328/210430281-96e0d954-05ee-4204-96b2-46152ff6578d.png) 

The accuracy of the lookup table was validated by plugging
these values in our Qutip Simulation [without optimization]. Based on the lookup table bounds
and color bar values, the simulation results exhibit good conformity with the
experimental data.

![image](https://user-images.githubusercontent.com/34755328/210431808-541fc0cb-fa23-45e4-a90e-0deb5c1ecdc8.png)

[Note that there is always an offset of ≃ 0.3 units in the +ve direction of the Re(α) axis in QuTip Wigner simulations
due to an arbitrary internal normalizing factor].

## Script to optimize the gate sequence parameters
The optimization script outputs all the possible optimal parameters that result in an implicit conditional fidelity >= 0.95 (this can be changed in the code). The least_squares() method from the SciPy library was implemented to minimize the infidelity (1-fidelity). This method was chosen as it gives the user the flexibility to adjust ftol [tolerance for change in cost function], xtol [tolerance for change in variables/parameters], and gtol [tolerance for change in gradient norm] simultaneously. The gtol parameter has been disabled by default as it causes the optimization to terminate at a local minima. The other tolerances can be set by the user.

Sample:

Input:
𝛼1=2, 𝛼2=1, 𝜃=0.5𝜋, 𝑛=0;

Output(s): 
1) 𝛼1=−1.14, 𝛼2=0.58, 𝜃=𝜋, 𝑛=0; Fidelity: 0.981
2) 𝛼1=1.14, 𝛼2=−0.58, 𝜃=𝜋, 𝑛=0; Fidelity: 0.981 

The input voltages corresponding to the computed optimized gate sequence parameters can be applied to the cQED setup to prepare an optimal fock state. 

## References
[1] Kudra M, Kervinen M, Strandberg I, Ahmed S, Scigliuzzo M, Osman A, et al. Robust Preparation of Wigner-
Negative States with Optimized SNAP-Displacement Sequences. PRX Quantum. 2022 Jul;3:030301. Available
from: https://link.aps.org/doi/10.1103/PRXQuantum.3.030301.

[2] Heeres RW, Vlastakis B, Holland E, Krastanov S, Albert VV, Frunzio L, et al. Cavity State Manipulation
Using Photon-Number Selective Phase Gates. Phys Rev Lett. 2015 Sep;115:137002. Available from: https:
//link.aps.org/doi/10.1103/PhysRevLett.115.137002.

[3] Quantum Technology Toolbox;. Available from: https://qtt.readthedocs.io/en/latest/index.html.


