# qmtoolbox
Simulation ToolBox for Quadrature Measurement of Harmonic Oscillator using Qubit

 1. INTRODUCTION

This package is devoted to simulate quadrature detection protocol of a harmonic
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
 
 If you run this file once more it simulates the trajectory and plot a new figure 
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
	
   
   - QuadMeasure, base class, defines a quadrature measurement protocol for the coupled system. 
	It setups qubit-HO systems, and simulates quadrature measurement and generates a measurement trajectory.
	
	- ThermalQM, inherited from QuadMeasure, describes either 
	thermal inital state or initial thermal displacement state of HO
    at a given temperature.
	

	(b) "qsim.py" contains a class called QuadSim which has an object 
	instantiated from one of FockQMP, CoherentQMP or ThermalQMP. The 
	UML below shows the relationship between QuadSim and QuadMeasure.
	

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
	
	- QuadSim manages the simulation of the measurement for various input parameters. 
		* Computes required parameters from the dictionary to support QuadMeasure
		* Obtains simulation results, 
		* Generates presentation file in LaTeX (beamer) (Removed from this version for simplicity). 
		
	
	(c) "qutil.py" consists of 
		- Call-back functions for Hamiltonian (S15) 
		- Module to creates visual and data representations 
		- Saves data/plots in the destination directory, 
		
	(d) "qpar.py" contains a dictionary (pars) for System Parameters.

	(e) "\_\_init\_\_.py" used to mark "qmtools/" as Python package directory
	
	
3. REFERENCE: 

[1] Mehmet Canturk and Adrian Lupascu, Quadrature readout and generation 
of squeezed states of a harmonic oscillator using a qubit-based indirect 
measurement, arXiv:1704.04533, (2017)


