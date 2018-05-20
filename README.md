# QMToolBox
Simulation ToolBox for Quadrature Measurement of Quantum Harmonic Oscillator using Qubit

## 1. Introduction

This toolbox is developed to simulate quadrature detection protocol of a harmonic
oscilator (HO) coupled to a superconducting flux qubit (quantum bit) with 
a time-dependent Hamiltonian (S15) in the supplamentary material [1]. 
It performs high fedelity measurement and generate squeezed states of the HO
which can be either electrical or mechanical resonator.

The measurement protocol, which is based on quantum nondemolishing measurement 
type protocol, consists of repeating procedure in Fig.S1 [1]. 
The qubit and HO start in a separable state where qubit is in ground and HO 
in thermal (or coherent or Fock) state respectively.

The qubit is controlled by resonant pulses (e.g spin echo or CPMG sequence)
as shown in Fig.S1. The evolution of combined system described by 
master equation in Lindblad form (S14) within control pulse and 
free precession times is numerically computed by *mesolve* in 
[QuTiP](www.qutip.org) governed by the Hamiltonian (S15). 
The influence of the qubit dephasing induced by its intrinsic environment
or 1/f noise (S13) is considered in the Hamiltonian (S15). Similarly, 
the contribution of the dissipated effects of the coupled system 
are considered in the master equation (S14).

The outcome of measurement trajectory of the quadrature and variance 
are obtained by series repetition of the protocol illustrated in Fig.S1.
A sample trajectory is obtained by this tool is shown "_qtrajectory_sample_result.png_".



## 2. Description of Simulation ToolBox (QMToolBox)

### 2.1 Classes and Modules in package directory "_qmtools/_":
	
#### (a) "_qmtools/quad.py_" consists of class hierarchy as shown in UML below: 
     
                       +-------------+
                       | QuadMeasure |
                       +------.------+
                             /_\
                              |
               --------------------------------
               |              |               |
       +-----------+    +----------+   +-------------+
       | ThermalQM |    | FockQM(*)|   |CoherentQM(*)|
       +-----------+    +----------+   +-------------+
	   
    (*) For simplicity, the code for FockQM and CoherentQM "qmtools/quad.py" does not exist.
	
   
   - Class **QuadMeasure** defines a quadrature measurement protocol for the coupled system: 
	 - Setups quantum operators, Hamiltonian, and Lindblad Master equation for qubit-HO system, 
	 - Algorithms for 
		- Quadrature measurement protocol(s), 
		- Evolution of density matrix in CPMG
		- Solving the master equation
	 - Generates the data for a single trajectory
	
   - Class **ThermalQM**, inherited from **QuadMeasure** describes thermal state of HO
    at a given temperature.
	
#### (b) "_qmtools/qsim.py_" consists of class **QuadSim** whose major component is **QuadMeasure** as shown in the UML.

                       +-------------+                +-------------+
                       | QuadMeasure |--------------<>|   QuadSim   |
                       +------.------+                +-------------+
                             /_\
                              |
               --------------------------------
               |              |               |
       +-----------+    +----------+   +-------------+
       | ThermalQM |    |  FockQM  |   | CoherentQM  |
       +-----------+    +----------+   +-------------+
	 
   - **QuadSim** 
        * Manages the simulation for various input parameters
        * Instantiates an object from **ThermalQM**, **FockQMP**, or **CoherentQMP**. 
		* Computes required parameters from the dictionary 
		* Performs the simulation and returns outcomes, 
		* Generates presentation file in LaTeX (beamer). 
		
	
#### (c) "_qmtools/qutil.py_" consists of modules for 
   - Call-back functions for Hamiltonian (S15) 
   - Generation of flux noise
   - Data processing and  visual representations using matplotlib
   - Saving data/figures in the destination directory, 
		
#### (d) "_qmtools/qpar.py_" includes system parameters in the form of dictionary (_pars_).


### 2.2 Main File "_qmtrajectory_sample_run.py_"
 When we run this file, it
 - creates output directory and data file names.
 - checks whether the data for 1/f noise exist in the current directory. 
 If not, it generates them using (S13).
 - simulates the measurement trajectory and generates data. 
 - plots the figure as shown in "_qtrajectory_sample_result.png_" and stores
 in the directory.
 
 Notice that if you run this file once more, it simulates the trajectory and plot a new figure 
 without generating the noise data.
	
	
## 3. REFERENCE: 

[1] Mehmet Canturk and Adrian Lupascu, Quadrature readout and generation 
of squeezed states of a harmonic oscillator using a qubit-based indirect 
measurement, [arXiv:1704.04533, (2017)](https://arxiv.org/abs/1704.04533)


