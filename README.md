# qmtoolbox
Simulation ToolBox for Quadrature Measurement of Harmonic Oscillator using Qubit

Introduction

This package is devoted to simulate quadrature detection protocol of a harmonic
oscilator (HO) coupled to a superconducting flux qubit (quantum bit) with 
a time-dependent Hamiltonian (S15) in supplamentary material in [1]. 
It performs high fedelity measurement and generate squeezed states of the HO
which can be either electrical or mechanical resonator.

The measurement protocol consists of repeating procedure in Fig.S1 [1]. 
The qubit and HO start in a separable state where qubit is in ground and HO 
in thermal (coherent and Fock state as well) states respectively.

The qubit is controlled by resonant pulses (i.g spin echo or CPMG sequence)
as shown in Fig.1(a) or Fig.S1. The evolution of combined system described by 
master equation in Lindblad form (S14) within control pulse and 
free precession times is numerically obtained by mesolve in QuTiP 
[www.qutip.org] governed by the Hamiltonian (S15). 
The influence of the qubit dephasing induced by its environment is considered 
in the Hamiltonian (S15) whereas the contribution of the dissipated effects 
of the coupled system are considered in the master equation (S14).
The measurement trajectory of the quadrature and variance are obtained by 
series repetition of the protocol illustrated in Fig.S1.
A sample trajectory is obtained by this tool is shown "qtraj.png".

Contents of "qmtoolbox/"

1) "qmtrajectory_sample_run.py" is a basic simulation file. 
When we runn this file it will creates output directory and data file names.
Next, it generates 1/f noise (S13) and outputs the files into the destination 
directory in which unless the noise data already exist. This is followed by 
the generation of the trajectory of the quadrature measurement. Finally 
the outcome as a data and associated plot files are stored in the destination. 
A sample outcome of the measurement result is given "qtrajectory_sample_result.png".
If you repeat the experiment several times, the 

2) "qmtools/" consists of required package files
	
	(a) "__init__.py" used to mark "qmtools/" as Python package directory
	
	(b) "quad.py" contains the class hierarchy in the form of UML: 
	

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
	
	- QuadMeasure, a base class, defines a quadrature measurement protocol 
	for HO coupled to the qubit. It setups qubit-HO systems, and simulates
	quadrature measurement and generates a measurement trajectory.
	
	- ThermalQMP, inherited from QuadMeasure, describes either 
	thermal inital state or initial thermal displacement state of HO
    at a given temperature.
	

	(c) "qsim.py" contains a class called QuadSim which has an object 
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
		
	
	(d) "qutil.py" consists of 
		- Call-back functions for Hamiltonian (S15) 
		- Module to creates visual and data representations 
		- Saves data/plots in the destination directory, 
		
	(e) "qpar.py" contains a dictionary (pars) for System Parameters.

	
Reference: 

[1] Mehmet Canturk and Adrian Lupascu, Quadrature readout and generation 
of squeezed states of a harmonic oscillator using a qubit-based indirect 
measurement, arXiv:1704.04533, (2017)


