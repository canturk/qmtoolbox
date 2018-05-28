# QMToolBox
Simulation ToolBox for Quadrature Measurement of Quantum Harmonic Oscillator using Qubit

## 1. Introduction

This toolbox is developed to simulate quadrature detection protocol of a harmonic
oscilator (HO) coupled to a superconducting flux qubit (quantum bit) with 
a time-dependent Hamiltonian (S15) in the supplamentary material [1]. 
Repeated measurement of the qubit (based on quantum nondemolishing measurement)
leads to gradually increasing information on the quadrature I, and 
a corresponding reduction in the uncertainty
$\Delta$I corresponding to squeezing state of HO.
The detection protocol in this toolbox realistically simulates
the measurement of superconducting
electromagnetic resonators and nano-mechanical
resonators, taking into account non-idealities including
decoherence and qubit detection errors. 

The detection protocol consists of repeating procedure as shown in Fig.S1 [1]. 
The qubit and HO start in a separable state where qubit is in ground and HO 
in thermal (or coherent or Fock) state respectively.
The qubit is controlled by resonant pulses (e.g spin echo or CPMG sequence)
as shown in Fig.1(a) or Fig.S1. The evolution of combined system described by 
master equation in Lindblad form (S14) within control pulse and 
free precession times is numerically computed by mesolve in [QuTiP](www.qutip.org) 
governed by the Hamiltonian (S15). 
The influence of the qubit dephasing induced by its intrinsic environment
or 1/f noise (S13) is considered in the Hamiltonian (S15). Similarly, 
the contribution of the dissipated effects of the coupled system 
are considered in the master equation (S14).

The outcome of measurement trajectory of the quadrature and variance 
are obtained by series repetition of the protocol illustrated in Fig.S1.
A sample trajectory is obtained by this tool is shown 
"[qtrajectory_sample_result.png](https://github.com/canturk/qmtoolbox/blob/master/qtrajectory_sample_result.png)".



## 2. Description of Simulation ToolBox (QMToolBox)

### 2.1 Classes and Modules in Package directory ("[qmtools/](https://github.com/canturk/qmtoolbox/blob/master/qmtools)"):
	
#### (a) "[quad.py](https://github.com/canturk/qmtoolbox/blob/master/qmtools/quad.py)" consists of class hierarchy, shown in UML: 
     
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
	   
    (*) The class definitions for FockQM and CoherentQM were removed from "quad.py".
	
   
   - Class **QuadMeasure** defines a quadrature measurement protocol for the coupled system: 
	 - Setups quantum operators, Hamiltonian, and Lindblad Master equation for qubit-HO system, 
	 - Algorithms for 
		- Quadrature measurement protocol(s), 
		- Evolution of density matrix in CPMG
		- Solving the master equation
	 - Generates the data for a single trajectory
	
   - Class **ThermalQM**, inherited from **QuadMeasure** describes thermal state of HO
    at a given temperature.
	
#### (b) "[qsim.py](https://github.com/canturk/qmtoolbox/blob/master/qmtools/qsim.py)" consists of class QuadSim whose major component is **QuadMeasure** as shown in the UML.

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
        * Instantiates a single object from **ThermalQM**, **FockQMP**, or **CoherentQMP**. 
		* Computes required parameters from the dictionary 
		* Performs the simulation and returns outcomes, 
		* ~~Generates presentation file in LaTeX (beamer).~~ 
		
	
#### (c) "[qutil.py](https://github.com/canturk/qmtoolbox/blob/master/qmtools/qutil.py)" consists of modules for 
   - Call-back functions for Hamiltonian (S15) 
   - Generation of flux noise
   - Data processing and  visual representations using matplotlib
   - Saving data/figures in the destination directory, 
		
#### (d) "[qpar.py](https://github.com/canturk/qmtoolbox/blob/master/qmtools/qpar.py)" includes system parameters in the form of dictionary (pars).


### 2.2 Main File ("[qmtrajectory_sample_run.py](https://github.com/canturk/qmtoolbox/blob/master/qmtrajectory_sample_run.py)")
 When we run this file, 
 - First, creates output directory and data file names.
 - Second, checks whether the data for 1/f noise exist in the current directory. 
 If not, it generates them using (S13).
 - Third, simulates the measurement trajectory and generates data. 
 - Finally plots the figure as shown in "qtrajectory_sample_result.png" and stores
 in the directory.
 
 Notice that if you run this file once more, it simulates the trajectory and plot a new figure 
 without generating the noise data.
	
	
## 3. REFERENCE: 

[1] Mehmet Canturk and Adrian Lupascu, Quadrature readout and generation 
of squeezed states of a harmonic oscillator using a qubit-based indirect 
measurement, [arXiv:1704.04533, (2017)](https://arxiv.org/pdf/1704.04533.pdf)


