# QMToolBox
Simulation ToolBox for Quadrature Measurement of Quantum Harmonic Oscillator using Qubit

## 1. INTRODUCTION

This toolbox is developed to simulate quadrature detection protocol of a harmonic
oscilator (HO) coupled to a superconducting flux qubit (quantum bit) with 
a time-dependent Hamiltonian (S15) in the supplamentary material [1]. 
It performs high fedelity measurement and generate squeezed states of the HO
which can be either electrical or mechanical resonator.

The measurement protocol consists of repeating procedure in Fig.S1 [1]. 
The qubit and HO start in a separable state where qubit is in ground and HO 
in thermal (or coherent or Fock) state respectively.

The qubit is controlled by resonant pulses (e.g spin echo or CPMG sequence)
as shown in Fig.1(a) or Fig.S1. The evolution of combined system described by 
master equation in Lindblad form (S14) within control pulse and 
free precession times is numerically computed by mesolve in QuTiP 
[www.qutip.org] governed by the Hamiltonian (S15). 
The influence of the qubit dephasing induced by its intrinsic environment
(i.e. 1/f noise) is considered in the Hamiltonian (S15). Similarly, 
the contribution of the dissipated effects of the coupled system 
are considered in the master equation (S14).

The outcome of measurement trajectory of the quadrature and variance 
are obtained by series repetition of the protocol illustrated in Fig.S1.
A sample trajectory is obtained by this tool is shown "qtrajectory_sample_result.png".



2. DESCRIPTION OF SIMULATION TOOLBOX ("qmtoolbox/")

 2.1 RUN FILE ("qmtrajectory_sample_run.py")
 When we run this file, 
 - First, creates output directory and data file names.
 - Second, checks whether the data for 1/f noise exist in the current directory. 
 If not, then generates them using (S13).
 - Third, simulates the trajectory of the quadrature measurement and generates the data. 
 - Finally plots the figure as shown in "qtrajectory_sample_result.png" and stores
 in the directory.
 
 If you run this file once more, it simulates the trajectory and plot a new figure 
 without generating the noise data.


 2.2 PACKAGE FILES IN THE SUBDIRECTORY ("qmtools/"):
	
(a) "quad.py" consists of class hierarchy, shown in UML: 
     
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
	   
    (*) The code for FockQM and CoherentQM "quad.py" does not exist.
	
   
   - Class QuadMeasure defines a quadrature measurement protocol for the coupled system: 
	 - Setups quantum operators, Hamiltonian, and Lindblad Master equation for qubit-HO system, 
	 - Simulates the measurement, and 
	 - Generates a measurement trajectory
	
   - Class ThermalQM inherited from QuadMeasure describes thermal inital (or thermal displacement) state of HO
    at a given temperature.
	
(b) "qsim.py" consists of class QuadSim which has a major component class QuadMeasure.
	The relationship between them is swon in the UML form.

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
	 
   - QuadSim 
        * Manages the simulation for various input parameters
        * Instantiates a single object from ThermalQM, FockQMP, or CoherentQMP. 
		* Computes required parameters from the dictionary 
		* Obtains simulation results, 
		* Generates presentation file in LaTeX (beamer). 
		
	
(c) "qutil.py" consists of 
		- Call-back functions for Hamiltonian (S15) 
		- Modules for data processing and  visual representations 
		- Saves data/plots in the destination directory, 
		
(d) "qpar.py" contains a dictionary (pars) for System Parameters.

(e) "\_\_init\_\_.py" marks "qmtools/" as a package directory
	
	
3. REFERENCE: 

[1] Mehmet Canturk and Adrian Lupascu, Quadrature readout and generation 
of squeezed states of a harmonic oscillator using a qubit-based indirect 
measurement, arXiv:1704.04533, (2017)


