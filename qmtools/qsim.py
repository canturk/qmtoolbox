#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# QND Measurement of a Nanomechanical Resonator Coupled to a Flux Qubit
###############################################################################
# Author : Mehmet Canturk 
# Date   : April 14, 2017
###############################################################################
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
from qutip import matrix_histogram, steadystate, hinton
from qutip import propagator, rotation, tensor, fidelity, tracedist, expect, mesolve, commutator
from qutip import sigmax, sigmay, sigmaz, sigmap, sigmam, qeye
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
from qmtools.qpar import pars
from qmtools.qutil import fz, makedir, save_fig
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
from qmtools.quad import ThermalQM #FockQM, CoherentQM, 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

import numpy as np
import sys # for manual exception handling
import matplotlib.pyplot as plt

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
###############################################################################
class QuadSim(object):
    """
    This class consists of a single object instantiated one of 
    FockQMP, CoherentQMP and ThermalQMP inherited from a class called QuadMeasure.
    QuadSim manages the simulation of quadrature measurement 
    for various parameters. 
    
    QuadSim is capable of obtaining simulation results, creating (and displaying) 
    figures, saving data, generating presentation file in LaTeX (beamer)
    (Removed from this version for simplicity). 
    1. Define new parameters for the dictionary (pars) and check whether 
    they are valid 
    2. Prepare initial state
    3. Setup Hamiltonian (and collapse operators if the simulation is for 
    Linblad master equation) 
    5. Check fidelity of Usim in terms of Uavg using propagator for qubit
    6. Analyze Bloch oscillation for qubit within a CPMG or real time 
    6.(b) Consider Bloch oscillation for HO+Qubit
    7. Define filename extension in the dictionary and then create a proper 
    filename in the module that generates results. After that call 
    the plot function located in qutil.py by passing filename as
    parameter.
    8. Fidelity for state between avg*rho0*Uavg.dag and rho = lindblad(rho0)
        Uavg = qs.qm.piY_half
             * rotation(qs.qmp.sx,np.pi)
             * (-1j * 2.0/np.pi * param['lambdaz'] * param['Te'] 
             * qs.qmp.sz).expm()*qs.qmp.piX_half;
        print fidelity(Uavg*rho0*Uavg.dag(),rho) 
    9. Implement ideal/actual quadrature measurement using the proper algorithms.
  """
    
    def __init__(self, state='Fock', isLind=False):
        
        #--------------------------------
        # SETUP SYSTEM PARAMETERS
        #---------------------------------
        
        # Qubit transition frequency in [GHz]
        pars['omegaGE'] = np.sqrt(pars['Delta']**2 + pars['Epsilon']**2);
        
        # The time interval between (pi/2)x and (pi/2)y pulses
        pars['Te']      = (pars['Np'] / 2.) * (2. * np.pi / pars['omegaR']);
        #
        # Free precession (evolution) time[ns]
        pars['tpr'] = (pars['Te'] / pars['Np']) / (1.0 + pars['tau2tp']);
        
        # (pi)x pulse duration [ns]
        # Redefine tau in terms of Half Period of HO
        pars['tau'] = pars['tau2tp'] * pars['tpr'];
        
        #Check whether the pulse duration (tau) must be less than half period of HO
        if pars['tau'] >= (2 * np.pi / pars['omegaR'] / 2.0):
            raise ValueError("tau exceeds  HO's half period");
        
        #-----------------------------------------
        from scipy.constants import physical_constants
        ktuple = physical_constants['Boltzmann constant']
        htuple = physical_constants['Planck constant']
        kB     = ktuple[0];#[Joule/Kelvin]
        hbar   = htuple[0] / (2*np.pi);#[Joule * sec]
        # Parameters for collapse operators
        Wr2Tho          = (hbar * pars['omegaR'] * 1e9) / (kB * pars['Tho']);
        pars['Wge2Tqb'] = (hbar * pars['omegaGE']* 1e9) / (kB * pars['Tqb']);
        Wge2Tqb         = pars['Wge2Tqb'];
        
        # The number of average photons in thermal equilibrium for the resonator
        # with the given frequency and temperature.
        pars['nHO'] = 1.0 / (np.exp(Wr2Tho) - 1.0);
        
        # Thermal constant for qubit
        pars['nqb'] = 1.0/ (1 + np.exp(- Wge2Tqb));
        
        # Decay rate for the resonator
        Kappa = pars['omegaR'] / pars['Q'];
        print("K=Wr/Q  = 2pi * %.3lg [GHz], Q=%.1lf"%(Kappa/(2*np.pi), pars['Q']))
        # Damping rates (Up and Down) for HO
        pars['kDwn'] = Kappa * (pars['nHO'] + 1.0);
        pars['kUp']  = Kappa * pars['nHO'];
        
        # Damping rate (Gamma_e2g/g2e) for qubit
        pars['Gamma_E2G'] = pars['nqb'] / (pars['T1'] / 1e-9);    #[rad/\mu sec]
        pars['Gamma_G2E'] = pars['Gamma_E2G'] * np.exp(- Wge2Tqb);#[rad/sec]
        
        #----------------------------------------------------------------------
        # Pure dephasing (Gamma_phi)[DISREGARD IT WHEN FLUX NOISE ISE CONSIDERED]
        Gamma_varphi = 0.0;#2*pi / (100e-9) / MHz; # [rad/sec]
        #----------------------------------------------------------------------
        
        print ("QuadMeasure object is created!");
        #----------------------------------------------------------------------
        self.qm    = []
        if state == 'Coherent':
            self.qm = CoherentQM(isLind);
        elif state == 'Thermal':
            self.qm = ThermalQM(isLind);
        elif state == 'Fock':
            self.qm = FockQM(isLind);
        else:
            sys.exit("Initial state of resonator was not proberly specified.");
        
        #----------------------------------------------------------------------
        print ("   Setup System Parameters in", self.__class__.__name__)
        
        #
        # Weighted coupling strentgh(s)
        pars['lambdat'] = pars['g'] * (pars['Delta']   / pars['omegaGE']);
        pars['lambdaz'] = pars['g'] * (pars['Epsilon'] / pars['omegaGE']);
        
        #
        # Coefficients for flux noise
        
        sqrtAf  = pars['sqrtAf']; #Dimensionless coefficient which is experimentally chosen from range
                #sqrt(Af)=1e-6 to sqrt(A)=2e-6.
        Ip      = 300e-9; #Persistent current [nA]
        etuple  = physical_constants['elementary charge'];
        echarge = etuple[0];#echarge = 1.60217662e-19;#[C]
        
        #Update System parameters
        pars['Traj'] = self.qm.time[-1];
        print ("Ttraj = %10.4lf [ns]"%pars['Traj'])
        
        pars['Wmin'] = 2*np.pi / pars['Traj'] / pars['v']; #Lower bound of omega in PSD
        pars['Wmax'] = 2*np.pi / pars['Te'] * pars['Np'] * pars['u'];#Upper bound of omega in PSD;
        pars['Nf']   = int(np.ceil(pars['Wmax']/pars['Wmin']));
        
        sqrtA_E     = 2 * np.pi * Ip / echarge * sqrtAf/1.0e9;#Coefficient in PSD [rad/ns]
        delEPS      = np.sqrt(2*np.pi / pars['Traj'] / pars['Wmin']) * sqrtA_E;
        print ("sqrtA_eps=%.5lg (rad/ns)"%(sqrtA_E))
        pars['deltat'] = 0.5 * delEPS * pars['Delta']   / pars['omegaGE'];#TANGANTIAL
        pars['deltaz'] = 0.5 * delEPS * pars['Epsilon'] / pars['omegaGE'];#LONGITUDIONAL
        
        #
        # Maximum amplitude of a pulse
        pars['LambdaT'] = np.pi / pars['tau'];#Amplitude of square PI pulse;
        pars['LambdaZ'] = pars['LambdaT'] * pars['Epsilon'] / pars['Delta'];#LONGIDUTIONAL
        
        #
        # Check whether the pulse amplitude, LambdaT < omegaGE such that LambdaT=pi/(2tau)
        # The accuracy is limited if the Rabi frequency of the control is comparable with the transition
        # frequency of the qubit due to the breakdown of the rotating wave approximation (RWA)
        if (pars['LambdaT'] >= pars['omegaGE'])\
        and (pars['LambdaT']==np.pi/(2*pars['tau'])):
            raise ValueError("Error in PI pulse amplitude: LambT < omegaGE such that LambT=pi/(2tau) NOT satisfied")
        #-----------------------------------------------------------------
        #self.print_System_Parameters();
        #
        return;#Construction
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def print_System_Parameters(self):
        print ("\n------------------------------------------------------")
        print ("Coupling Strength:")
        print("   g = 2pi * %.4lf [GHz]"\
              %(np.sqrt(pars['lambdat']**2 + pars['lambdaz']**2)/(2*np.pi)))
        print ("lamt = 2pi * %.4lf [GHz],\nlamz = 2pi * %.4lf [GHz]"\
               %(pars['lambdat']/(2*np.pi), pars['lambdaz']/(2*np.pi)))
        
        print("\ndelT = 2pi * %.5lf [GHz]\ndelZ = 2pi * %.5lf [GHz]"\
              %(pars['deltat']/2/np.pi, pars['deltaz']/2/np.pi));
        print ("Wmin = 2pi * %.5lf [GHz]"%(pars['Wmin']/(2*np.pi)));
        print ("Wmax = 2pi * %.5lf [GHz]"%(pars['Wmax']/(2*np.pi)));
        print ("Nf[Calc]= ",str(pars['Nf']));
        
        print ("\nMaximum Amplitude of pulses:\nLambdaT = 2pi * %5.3lf [GHz],\nLambdaZ = 2pi * %5.3lf [GHz],\nomegaGE = 2pi * %5.2lf [GHz]"\
               %(pars['LambdaT']/(2*np.pi), pars['LambdaZ']/(2*np.pi), pars['omegaGE']/(2*np.pi)));
        
        print ("\nReal-time for quadrature measurement is synchronized in module")
        print ("generate_RealTimeMeasurementSequence()")
        print ("Move required parameters  for the update")
        print ("Traj    = %8.3lf [ns]"%(pars['Traj']))
        print ("Te      = %8.3lf [ns],\ntpr     = %8.3lf [ns],\ntau     = %8.3lf [ns] => HO's half period Tr/2=pi/omegar=%8.3lf [ns]"\
               %(pars['Te'], pars['tpr'],  pars['tau'], 2*np.pi / (2 * pars['omegaR'])))
        print ("tau/2   = %8.3lf [ns]"%(pars['tau']/2))
        print ("s       = %d,\nNp      = %d,\ntau2t   = %0.4lf"%(pars['s'], pars['Np'], pars['tau2tp']))
        
        phi = 2 / np.pi * pars['g'] * pars['Te'];
        print ("varI[Asymp]=1/(1+4phi^2*s)=%.3lf where phi=%.3lf"%(1/(1+4*phi**2*pars['s']), phi))
        print ("------------------------------------------------------")
        
        k_Dwn = pars['kDwn'];
        k_Up  = pars['kUp'];
        G_e2g = pars['Gamma_E2G'];
        G_g2e = pars['Gamma_G2E'];
        nHO   = pars['nHO'];
        
        print("n_HO    = %.2lf,\nn_qb    = %.2lf"%(nHO, pars['nqb']));
        print("\nkDwn = 2pi * %.3lg [GHz],\nkUp  = 2pi * %.3lg [GHz]"%(k_Dwn/(2*np.pi), k_Up/(2*np.pi)))
        print("Ge2g =       %.3lg   [GHz],\nGg2e =       %.3lg [GHz]"%(G_e2g, G_g2e))
        print ("------------------------------------------------------")
        return;
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def plot_qubitBlochOscillations(self, tlist, exlist, eylist, ezlist, title, dirname, filename):
        #Update the matplotlib configuration parameters:
        plt.rcParams.update({'font.size': 16, 'font.family': 'serif'})
        #---------------------------------------------------------------------#
        plt.ioff() #to prevent the  multiple pop up figures for each repetition
        fig, axes = plt.subplots(1, 1, figsize=(8,6), sharey=True)
        axes.plot(tlist, np.real(ezlist), 'r-', label=r'$\left\langle\sigma_z\right\rangle$');
        axes.plot(tlist, np.real(eylist), 'b-', label=r'$\left\langle\sigma_y\right\rangle$');
        axes.plot(tlist, np.real(exlist), 'g-', label=r'$\left\langle\sigma_x\right\rangle$');
        axes.set_xlabel(r'$\tau$ [ns]', fontsize=16)
        axes.set_xlim([tlist[0], tlist[-1]]);
        axes.set_title(title)
        axes.grid(True);
        axes.set_ylabel('Bloch Oscillation', fontsize=16)
        # Now add the legend with some customizations.
        legend = axes.legend(loc='lower left', shadow=False)
        # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
        frame = legend.get_frame()
        frame.set_facecolor('0.90')
        # Set the fontsize
        for label in legend.get_texts():
            label.set_fontsize('x-large')
        # Set the line width
        for label in legend.get_lines():
            label.set_linewidth(2.0)  # the legend line width
        #-----------------------------------
        #Create directory for results
        if (dirname == [] and filename == []):
            plt.show()
        else:
            #makedir('../'+dirname)
            save_fig('../'+dirname+filename+'_bloch', ext='png', close=True, verbose=True)
        return;
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def plot_cpmg_BlochOscilations(self, tlst, exlst, eylst, ezlst, eIlst,
                                   title_list, filename, fidlst=[]):
        plt.close('all');
        
        color_list = [];
        color_list.append('r-');
        color_list.append('b-');
        color_list.append('g-');
        color_list.append('y-');
        
        label_list = [];
        label_list.append(r'$\left\langle\sigma_z\right\rangle$')
        label_list.append(r'$\left\langle\sigma_y\right\rangle$')
        label_list.append(r'$\left\langle\sigma_x\right\rangle$')
        label_list.append(r'$\left\langle I\right\rangle$')
                
        #Update the matplotlib configuration parameters:
        plt.rcParams.update({'font.size': 14, 'font.family': 'serif'})
        ##-----------------------------------------------------------------------##
        plt.ioff() #to prevent the  multiple pop up figures for each repetition
        fig, axes = plt.subplots(1, len(ezlst), figsize=(14,6), sharey=True)
        for m in range(len(axes)):
            axes[m].plot(tlst[m], np.real(ezlst[m]), color_list[0], label=label_list[0]);
            axes[m].plot(tlst[m], np.real(eylst[m]), color_list[1], label=label_list[1]);
            axes[m].plot(tlst[m], np.real(exlst[m]), color_list[2], label=label_list[2]);
            axes[m].set_xlabel(r'$\tau$ [ns]', fontsize=16)
            axes[m].set_xlim([tlst[m][0], tlst[m][-1]]);
            axes[m].set_title(title_list[m])
            axes[m].grid(True);
            
        axes[0].set_ylabel('Bloch Oscillations', fontsize=16)
        # Now add the legend with some customizations.
        legend = axes[len(axes)-1].legend(loc='lower right', shadow=False)
        # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
        frame = legend.get_frame()
        frame.set_facecolor('0.90')
        # Set the fontsize
        for label in legend.get_texts():
            label.set_fontsize('xx-large')
        # Set the line width
        for label in legend.get_lines():
            label.set_linewidth(2.0)  # the legend line width
        #-----------------------------------
        #Create directory for results
        if (filename == []):
            plt.show()
        else:
            save_fig('../'+filename+'_xyz', ext='png', close=True, verbose=True)
        
        plt.ioff() #to prevent the  multiple pop up figures for each repetition
        fig, axes = plt.subplots(1, len(eIlst), figsize=(14,6), sharey=True)
        for m in range(len(axes)):
            axes[m].plot(tlst[m], eIlst[m], color_list[0], label=label_list[3]);
            axes[m].set_xlabel(r'$\tau$ [ns]', fontsize=16)
            axes[m].set_xlim([tlst[m][0], tlst[m][-1]]);
            axes[m].set_title(title_list[m])
            axes[m].grid(True);
            
        axes[0].set_ylabel(r'$\langle I\rangle$ Oscillations', fontsize=16)
        # Now add the legend with some customizations.
        legend = axes[len(axes)-1].legend(loc='lower right', shadow=False)
        # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
        frame = legend.get_frame()
        frame.set_facecolor('0.90')
        # Set the fontsize
        for label in legend.get_texts():
            label.set_fontsize('small')
        # Set the line width
        for label in legend.get_lines():
            label.set_linewidth(2.0)  # the legend line width
        #-----------------------------------
        #Create directory for results
        if (filename == []):
            plt.show()
        else:
            save_fig('../'+filename+'_quaI', ext='png', close=True, verbose=True)
        return;
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def analyze_mesolve_halfPiX(self, rho0, ti, tf, dirname=[], filename=[]):
        
        # TEST WITHOUT COLLAPSE OPERATORS 
        self.qm.collapse = [];
        # TEST FOR SQUEARE PULSE
        pars.update({'betarf':0.0});
        
        ## Check out operator form of qb+HO Hamiltonian
        t = tf/4;# ti <= t<= tf
        HpiX = self.qm.operator_HpiX(t,1.0)
        #Visualize the Hamiltonian
        fig, ax = matrix_histogram(HpiX)
        ax.set_title(r"Hamiltonian in $(\pi/2)_x$", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        if (dirname == [] and filename == []):
            plt.show()
        else:
            #dirname  = '../' + dirname + '/'; makedir(dirname)
            save_fig('../' + dirname+filename+'_Hpix', ext='png', close=True, verbose=True)
        
        ## Solve the system for the following settings
        rho1 = self.qm.mesolve_PiX(rho0, ti, tf);
        #Visualize the Hamiltonian
        fig, ax = matrix_histogram(rho1)
        ax.set_title(r"Density matrix in $(\pi/2)_x$", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        #Create directory for results
        if (dirname == [] and filename == []):
            plt.show()
        else:
            #dirname  = '../' + dirname + '/'; makedir(dirname)
            save_fig('../' + dirname+filename+'_rho', ext='png', close=True, verbose=True)
        # Qubit rotation in (pi/2)x pulse
        Upix = rotation(sigmax(), np.pi/2)
        rho_qb0 = rho0.ptrace(0)
        rho_qbx = Upix * rho_qb0 * Upix.dag()
        # Qubit+HO rotation in (pi/2)x pulse
        U_rx = rotation(self.qm.sx, np.pi/2);
        rho_rx = U_rx * rho0 * U_rx.dag();
        print ("State Fidelity(rho_mesolve, rho_rotX) in (pi/2)x: ", fidelity(rho1, rho_rx))
        print ("Trace distance(rho_mesolve, rho_rotX) in (pi/2)x: ",  tracedist(rho1, rho_rx))
        
        ## Qubit state after evolution
        #rho_qb1 = rho1.ptrace(0) 
        #print (rho_qb1, rho_qbx)
        
        
        # Qubit Fidelity for Unitary evolution operator
        tlist = np.linspace(ti, tf, int(pars['Num']/2));#end point not included
        # Setup Hamiltonian for (pi/2)x or (pi)x pulse
        #pars.update({'betarf':0.01});
        HpiX_qb, Hargs_qb = self.qm.setup_HpiX_qb(ti, tf, 1.0);
        Us = propagator(HpiX_qb, tlist, [], Hargs_qb, opts=[]);
        Fid = 0.5 * (Us[-1].dag() * Upix).tr()
        print ("Qubit Fidelity for Unitary Evolution Operator in (pi/2)x: %.4lf"\
               %Fid.real)
        #
        # Qubit Bloch Oscillation and Expectation values of Pauli operators
        rho_qb = self.qm.prepare_qubit_state("Ideal")
        # Solve me for qubit
        out_qb = mesolve(HpiX_qb, rho_qb, tlist, [], [], Hargs_qb, self.qm.options, progress_bar=self.qm.progress_bar);
        ezlist  = expect(sigmaz(), out_qb.states);
        eylist  = expect(sigmay(), out_qb.states);
        exlist  = expect(sigmax(), out_qb.states);
        #-------------------------------------
        title = r'$\left(\frac{\pi}{2}\right)_x$,   $F=$' +str(round(Fid.real, 3));
        self.plot_qubitBlochOscillations(tlist, exlist, eylist, ezlist, title, dirname, filename+'_bloch')
        return;
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def analyze_mesolve_FreeEvolution(self, rho0, ti, tf, dirname=[], filename=[]):
        
        ## Check out operator form of qb+HO Hamiltonian
        t = tf;# ti <= t<= tf
        HwithNoise, HwoutNoise = self.qm.operator_Hfree(t)
        #Visualize the Hamiltonian with flux noise
        fig, ax = matrix_histogram(HwithNoise)
        ax.set_title(r"Hamiltonian in free evolution with flux noise", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        
        #Visualize the Hamiltonian without flux noise
        fig, ax = matrix_histogram(HwoutNoise)
        ax.set_title(r"Hamiltonian in free evolution with flux noise", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        if (dirname == [] and filename == []):
            plt.show()
        else:
            #dirname  = '../' + dirname + '/';makedir(dirname)
            save_fig('../' + dirname+filename+'_Hpix', ext='png', close=True, verbose=True)
        
        # WITHOUT COLLAPSE OPERATORS 
        self.qm.collapse = [];
        ## Solve the system for the following settings
        rho_woutNoise = self.qm.mesolve_freeEvolution(rho0, ti, tf);
        #Visualize the Hamiltonian
        fig, ax = matrix_histogram(rho_woutNoise)
        ax.set_title(r"Density matrix in half free evolution WITHOUT flux noise", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        
        # WITH COLLAPSE OPERATORS 
        self.qm.collapse = self.qm.setup_collapse_operators();
        ## Solve the system for the following settings
        rho_withNoise = self.qm.mesolve_freeEvolution(rho0, ti, tf);
        #Visualize the Hamiltonian
        fig, ax = matrix_histogram(rho_withNoise)
        ax.set_title(r"Density matrix in full free evolution WITH flux noise", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        #Create directory for results
        if (dirname == [] and filename == []):
            plt.show()
        else:
            #dirname  = '../' + dirname + '/';makedir(dirname)
            save_fig('../' + dirname+filename+'_rho', ext='png', close=True, verbose=True)
        # Qubit rotation in (pi/2)x pulse
                
        ## Qubit state after evolution
        #print (rho_withNoise.ptrace(0), rho_woutNoise.ptrace(0))
        print ("State fidelity(rho_withNoise, rho_withoutNoise) in free ev: ",fidelity(rho_withNoise.ptrace(0), rho_woutNoise.ptrace(0)))
        print ("Trace distance(rho_withNoise, rho_withoutNoise) in free ev: ",tracedist(rho_withNoise.ptrace(0), rho_woutNoise.ptrace(0)))
         
        #
        # Bloch Oscillation and Expectation values of Pauli operators
        Hfree, Hargs = self.qm.setup_Hfree(ti, tf)
        tlist = np.linspace(ti, tf, pars['Num']);
        out_free = mesolve(Hfree, rho0, tlist, [], [], Hargs, self.qm.options, progress_bar=self.qm.progress_bar);
        ezlist  = expect(self.qm.sz, out_free.states);
        eylist  = expect(self.qm.sy, out_free.states);
        exlist  = expect(self.qm.sx, out_free.states);
        title = r'Free evolution time';
        self.plot_qubitBlochOscillations(tlist, exlist, eylist, ezlist, title, dirname, filename+'_bloch')
        return;
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def analyze_mesolve_fullPiX(self, rho0, ti, tf, sign, dirname=[], filename=[]):
        
        # TEST WITHOUT COLLAPSE OPERATORS 
        self.qm.collapse = [];
        # TEST FOR SQUEARE PULSE
        pars.update({'betarf':0.0});
        
        ## Check out operator form of qb+HO Hamiltonian
        t = ti;# ti <= t<= tf
        HpiX = self.qm.operator_HpiX(t, sign)
        #Visualize the Hamiltonian
        fig, ax = matrix_histogram(HpiX)
        ax.set_title(r"Visualization of Hamiltonian in $(\pi)_x$", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        if (dirname == [] and filename == []):
            plt.show()
        else:
            save_fig('../' + dirname+filename+'_HpiX', ext='png', close=True, verbose=True)
        
        ### Solve the system for the following settings
        ##pars.update({'Num':pars['Num']})
        #self.qm.options.atol     = 1e-5;#Absolute tolerance
        #self.qm.options.rtol     = 1e-6;#Relative tolerance
        rho1 = self.qm.mesolve_PiX(rho0, ti, tf, sign);
        #Visualize the density matrix 
        fig, ax = matrix_histogram(rho1)
        ax.set_title(r"$\rho_{me}$ in $(\pi)_x$", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        #Create directory for results
        if (dirname == [] and filename == []):
            plt.show()
        else:
            save_fig('../' + dirname+filename+'_rho', ext='png', close=True, verbose=True)
        
        # Qubit rotation in (pi)x pulse
        Upix = rotation(sigmax(), sign*np.pi)
        rho_qb0 = rho0.ptrace(0)
        rho_qbx = Upix * rho_qb0 * Upix.dag()
        # Qubit+HO rotation in (pi)x pulse
        U_rx = rotation(self.qm.sx, sign*np.pi);
        rho_rx = U_rx * rho0 * U_rx.dag();
        #Visualize the density matrix rho_rx
        fig, ax = matrix_histogram(rho_rx)
        ax.set_title(r"$\rho_{\pi_x}$ in $(\pi)_x$", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        print ("State Fidelity(rho_mesolve, rho_rotX) in (pi)x: ", fidelity(rho1, rho_rx))
        print ("Trace distance(rho_mesolve, rho_rotX) in (pi)x: ", tracedist(rho1, rho_rx))
        
        # Qubit state after evolution
        #rho_qb1 = rho1.ptrace(0) 
        #print (rho_qb1, rho_qbx)       
        
        # Qubit Fidelity for Unitary evolution operator
        tlist = np.linspace(ti, tf, pars['Num']);#end point not included
        # Setup Hamiltonian for (pi/2)x or (pi)x pulse
        #pars.update({'betarf':0.01});
        HpiX_qb, Hargs_qb = self.qm.setup_HpiX_qb(ti, tf, sign);
        Us = propagator(HpiX_qb, tlist, [], Hargs_qb, opts=[]);
        Fid = 0.5 * (Us[-1].dag() * Upix).tr()
        print ("Qubit Fidelity of Unitary Evolution Operator in (pi)x pulse: %.4lf"\
               %Fid.real)
        
        #
        # Qubit Bloch Oscillation and Expectation values of Pauli operators
        rho_qb = self.qm.prepare_qubit_state("Ideal")
        # Solve me for qubit
        out_qb = mesolve(HpiX_qb, rho_qb, tlist, [], [], Hargs_qb, self.qm.options, progress_bar=self.qm.progress_bar);
        ezlist  = expect(sigmaz(), out_qb.states);
        eylist  = expect(sigmay(), out_qb.states);
        exlist  = expect(sigmax(), out_qb.states);
        #-------------------------------------
        title = r'$\left(\pi\right)_x$,   $F=$' +str(round(Fid.real, 3));
        self.plot_qubitBlochOscillations(tlist, exlist, eylist, ezlist, title,\
                                         dirname, filename+'_bloch')
        return;
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def analyze_mesolve_halfPiY(self, rho0, ti, tf, dirname=[], filename=[]):
        
        # TEST WITHOUT COLLAPSE OPERATORS 
        self.qm.collapse = [];
        # TEST FOR SQUEARE PULSE
        pars.update({'betarf':0.0});
        
        ## Check out operator form of qb+HO Hamiltonian
        t = ti;# ti <= t<= tf
        HpiY = self.qm.operator_HpiY(t)
        #Visualize the Hamiltonian
        fig, ax = matrix_histogram(HpiY)
        ax.set_title(r"Visualization of Hamiltonian in $(\frac{\pi}{2})_y$", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        if (dirname == [] and filename == []):
            plt.show()
        else:
            #dirname  = '../' + dirname + '/';makedir(dirname)
            save_fig('../' + dirname+filename+'_HpiY', ext='png', close=True, verbose=True)
        
        ## Solve the system for the following settings
        rho1 = self.qm.mesolve_PiY(rho0, ti, tf);
        #Visualize the Hamiltonian
        fig, ax = matrix_histogram(rho1)
        ax.set_title(r"Visualization of density matrix in $(\pi/2)_y$", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        #Create directory for results
        if (dirname == [] and filename == []):
            plt.show()
        else:
            #dirname  = '../' + dirname + '/';makedir(dirname)
            save_fig('../' + dirname+filename+'_rhoY', ext='png', close=True, verbose=True)
        
        # Qubit rotation in (pi/2)y pulse
        U_qby = rotation(sigmay(), np.pi/2)
        rho_qb0 = rho0.ptrace(0)
        rho_qby = U_qby * rho_qb0 * U_qby.dag()
        # Qubit+HO rotation in (pi/2)y pulse
        U_ry = rotation(self.qm.sy, np.pi/2);
        rho_ry = U_ry * rho0 * U_ry.dag();
        print ("State Fidelity(rho_mesolve, rho_rotY): ", fidelity(rho1, rho_ry))
        print ("Trace distance(rho_mesolve, rho_rotY): ", tracedist(rho1,rho_ry))
        
        ## Qubit state after evolution
        #rho_qb1 = rho1.ptrace(0) 
        #print (rho_qb1, rho_qby)
        
        
        # Qubit Fidelity for Unitary evolution operator
        tlist = np.linspace(ti, tf, int(pars['Num']/2));#end point not included
        # Setup Hamiltonian for (pi/2)y pulse
        HpiY_qb, Hargs_qb = self.qm.setup_HpiY_qb(ti, tf);
        Us = propagator(HpiY_qb, tlist, [], Hargs_qb, opts=[]);
        Fid = 0.5 * (Us[-1].dag() * U_qby).tr()
        print ("Qubit Fidelity of Unitary Evolution Operator in (pi/2)y pulse: %.4lf"\
               %Fid.real)
        #
        # Qubit Bloch Oscillation and Expectation values of Pauli operators
        rho_qb = self.qm.prepare_qubit_state("Ideal")
        # Solve me for qubit
        out_qb = mesolve(HpiY_qb, rho_qb, tlist, [], [], Hargs_qb, self.qm.options, progress_bar=self.qm.progress_bar);
        ezlist  = expect(sigmaz(), out_qb.states);
        eylist  = expect(sigmay(), out_qb.states);
        exlist  = expect(sigmax(), out_qb.states);
        #-------------------------------------
        title = r'$\left(\frac{\pi}{2}\right)_y$,   $F=$' +str(round(Fid.real, 3));
        self.plot_qubitBlochOscillations(tlist, exlist, eylist, ezlist, title, dirname, filename+'_bloch')
        return;
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def checkout_HO_synch(self, header):
        """
        Paramameters:
            tlist: subsequent time interval from 0 to pars['Ttraj']
            header={"With Sync", "Without Sync"}
        """
        # Generate subsequenent meas time interval
        if header == "With Sync":
            tlist = self.qm.generate_RealTimeSequence(sync=True)
        elif header == "Without Sync":
            tlist = self.qm.generate_RealTimeSequence(sync=False);
        else:
            raise ValueError("Invalid sync choice")
        #----------------------------------------------------------------------
        def cosine(ti,tf):
            """Internal function to generate real and imaginary components of the functions"""
            # Time list for free evolution 
            #tau_lst = np.linspace(ti, tf, args['Num'],  endpoint=False);
            #ACCORDING TO #
            tau_lst = np.linspace(ti, tf, pars['Num']);
            #
            # Analysis of fcx(t), fcy(t), fcz(t) in the (pi/2)x pulse duration
            hcos = np.zeros(len(tau_lst), dtype = float)
            for i in range(len(tau_lst)):
                hcos[i] = np.cos(pars['omegaR']*tau_lst[i])
            return (tau_lst, hcos);
        #---------------------------------------------------------------------#
        #Update the matplotlib configuration parameters:
        plt.rcParams.update({'font.size': 16, 'font.family': 'serif'})
        plt.ioff() #to prevent the  multiple pop up figures for each repetition
        fig, axes = plt.subplots(1, 1, figsize=(10,3), sharey=True)
        #---------------------------------------------------------------------#
        Np=pars['Np'];
        
        hcos_dot = []
        tcos_dot = []

        m = 0;
        while (m<pars['s']):
            mp = 2 * m * (pars['Np'] + 2);#index            
            tcos_dot.append(tlist[mp]);
            hcos_dot.append(np.cos(pars['omegaR']*tcos_dot[-1]))
            tcos_dot.append(tlist[mp + 2*Np+3]);
            hcos_dot.append(np.cos(pars['omegaR']*tcos_dot[-1]))
            tc, hzr = cosine(tlist[mp], tlist[mp+2*pars['Np']+4])
            
            axes.plot(tc, hzr, label=r'$cpmg_{'+str(m+1)+'}+M_{'+str(m+1)+'}$');
            m += 1;
        
        axes.plot(tcos_dot, hcos_dot, 'b^',label='cpmg-intvl')
        axes.set_title(header)
        axes.set_xlabel(r'$\eta(t_{cpmg}+T_M)$ [ns]', fontsize=14)
        axes.grid(True);
        # Add the legend with some customizations.
        legend = axes.legend(loc='lower right', shadow=False)
        # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
        frame = legend.get_frame()
        frame.set_facecolor('0.90')
        # Set the fontsize
        for label in legend.get_texts():
            label.set_fontsize('xx-small')
        # Set the line width
        for label in legend.get_lines():
            label.set_linewidth(2.0)  # the legend line width
        #-----------------------------------
        plt.show()
        return;
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def checkif_HO_synchronize_via_fz_DEPRECATED(self, tlist, header):
        """
        Paramameters:
            tlist: subsequent time interval from 0 to pars['Ttraj']
            header={"With Sync", "Without Sync"}
        """
        def gfz(ti,tf, args):
            """Internal function to generate real and imaginary components of the functions"""
            # Time list for free evolution 
            #tau_lst = np.linspace(ti, tf, args['Num'],  endpoint=False);
            #ACCORDING TO #
            tau_lst = np.linspace(ti, tf, args['Num']);
            #
            # Analysis of fcx(t), fcy(t), fcz(t) in the (pi/2)x pulse duration
            hZR = np.zeros(len(tau_lst), dtype = float)
            hZI = np.zeros(len(tau_lst), dtype = float)
            hZS = np.zeros(len(tau_lst), dtype = float)
            for i in range(len(tau_lst)):
                hZR[i] = args['lambdaz'] * np.real(fz(tau_lst[i], args))
                hZI[i] = args['lambdaz'] * np.imag(fz(tau_lst[i], args))
                hZS[i] = np.sqrt(hZR[i]**2 + hZI[i]**2);
            return (tau_lst, hZR, hZI, hZS);
            #return (tau_lst, hZR);
        #---------------------------------------------
        #Full Free evolution
        args = {
                'omegaR' : pars['omegaR'],
                'lambdaz': pars['lambdaz'],
                'omegaGE': pars['omegaGE'],
                'Num'    : pars['Num']};
        #---------------------------------------------------------------------#
        #Update the matplotlib configuration parameters:
        plt.rcParams.update({'font.size': 16, 'font.family': 'serif'})
        plt.ioff() #to prevent the  multiple pop up figures for each repetition
        fig, axes = plt.subplots(3, 1, figsize=(8,12), sharey=True)
        #---------------------------------------------------------------------#
        Np=pars['Np'];
        
        hzr_dot = []
        hzi_dot = []
        tim_dot = []
#        tim_dot.append(tlist[0]);
#        hzr_dot.append(np.real(args['lambdaz'] * fz(tim_dot[-1], args)))
#        hzi_dot.append(np.imag(args['lambdaz'] * fz(tim_dot[-1], args)))
#        tim_dot.append(tlist[2*Np+3]);
#        hzr_dot.append(np.real(args['lambdaz'] * fz(tim_dot[-1], args)))
#        hzi_dot.append(np.imag(args['lambdaz'] * fz(tim_dot[-1], args)))
        m = 0;
        while (m<pars['s']):
            mp = 2 * m * (pars['Np'] + 2);#index            
            tim_dot.append(tlist[mp]);
            hzr_dot.append(np.real(args['lambdaz'] * fz(tim_dot[-1], args)))
            hzi_dot.append(np.imag(args['lambdaz'] * fz(tim_dot[-1], args)))
            tim_dot.append(tlist[mp + 2*Np+3]);
            hzr_dot.append(args['lambdaz'] * np.real(fz(tim_dot[-1], args)))
            hzi_dot.append(args['lambdaz'] * np.imag(fz(tim_dot[-1], args)))
            
            tc, hzr, hzi, hzs = gfz(tlist[mp], tlist[mp+2*pars['Np']+4], args)
            
            axes[0].plot(tc, hzr, label=r'$\Re\{f_{z}\}_{'+str(m+1)+'}$');
            axes[1].plot(tc, hzi, label=r'$\Im\{f_{z}\}_{'+str(m+1)+'}$');
            axes[2].plot(tc, hzs, label=r'$\sqrt{|f_{z}|}_{'+str(m+1)+'}$');
            m += 1;
        
        axes[0].plot(tim_dot, hzr_dot, 'bo')
        axes[1].plot(tim_dot, hzi_dot, 'bo')
        axes[0].set_title(header)
        axes[2].set_xlabel(r'$\eta(t_{cpmg}+T_M)$ [ns]', fontsize=14)
        axes[0].grid(True);axes[1].grid(True);axes[2].grid(True);
        # Add the legend with some customizations.
        for n in range(0,3):
            legend = axes[n].legend(loc='lower left', shadow=False)
            # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
            frame = legend.get_frame()
            frame.set_facecolor('0.90')
            # Set the fontsize
            for label in legend.get_texts():
                label.set_fontsize('xx-small')
            # Set the line width
            for label in legend.get_lines():
                label.set_linewidth(2.0)  # the legend line width
        #-----------------------------------
        plt.show()
        return;
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def EulerMethod_vs_QuTiPsolver(self, rho0, filename=[]):
        """
        Ensure that mesolve in Qutip solve ME for square pulse
        and comparison will be done
        """
        pars.update({'betarf':0.0})
        Num = pars['Num'];
        
        #----------------------------
        # Compare for (pi/2)x pulse
        #----------------------------
        pars.update({'Num':int(Num/2)})
        rho_hPiX_euler = self.qm.euler_mesolve  (rho0, 0.0, pars['tau']/2, "X-pulse")
        rho_hPiX_qutip = self.qm.mesolve_PiX(rho0, 0.0, pars['tau']/2);
        # Fidelity
        rho_hPiX_fid = fidelity(rho_hPiX_euler, rho_hPiX_qutip);
        print("State Fidelity(rho_qutip, rho_euler) in (pi/2)x ",rho_hPiX_fid)
        fid_data = [];
        fid_head = [];
        fid_data.append(rho_hPiX_fid);
        fid_head.append("rho_half_PiX")
        # Visualization of  density matrices
        fig, ax = matrix_histogram(rho_hPiX_euler)
        ax.set_title(r"$\rho_{\left(\frac{\pi}{2}\right)_x}$ via Euler Method ", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        fig, ax = matrix_histogram(rho_hPiX_qutip)
        ax.set_title(r"$\rho_{\left(\frac{\pi}{2}\right)_x}$ via QuTiP solver", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        
        #----------------------------
        # Compare for (pi)x pulse
        #----------------------------
        pars.update({'Num':Num})
        rho_fPiX_euler = self.qm.euler_mesolve  (rho0, 0.0, pars['tau'], "X-pulse")
        rho_fPiX_qutip = self.qm.mesolve_PiX(rho0, 0.0, pars['tau'], 1.0);
        # Fidelity
        rho_fPiX_fid = fidelity(rho_fPiX_euler, rho_fPiX_qutip);
        print("State Fidelity(rho_qutip, rho_euler) in (pi)x",rho_fPiX_fid)
        fid_data.append(rho_fPiX_fid);
        fid_head.append("rho_full_PiX")
        
        #Visualization of  density matrices
        fig, ax = matrix_histogram(rho_fPiX_euler)
        ax.set_title(r"$\rho_{(\pi)_x}$ via Euler Method ", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        fig, ax = matrix_histogram(rho_fPiX_qutip)
        ax.set_title(r"$\rho_{(\pi)_x}$ via QuTiP Solver", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        
        
        #----------------------------
        # Compare for (pi/2)y pulse
        #----------------------------
        pars.update({'Num':int(Num/2)})
        rho_hPiY_euler = self.qm.euler_mesolve  (rho0, 0.0, pars['tau']/2, "Y-pulse")
        rho_hPiY_qutip = self.qm.mesolve_PiY(rho0, 0.0, pars['tau']/2);
        #Fidelity
        rho_hPiY_fid = fidelity(rho_hPiY_euler, rho_hPiY_qutip);
        print("State Fidelity(rho_qutip, rho_euler) in (pi/2)y ",rho_hPiY_fid)
        
        fid_data.append(rho_hPiY_fid);
        fid_head.append("rho_half_PiY")
        
        #Visualization of  density matrices
        fig, ax = matrix_histogram(rho_hPiY_euler)
        ax.set_title(r"$\rho_{\left(\frac{\pi}{2}\right)_y}$ via Euler Method ", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        
        fig, ax = matrix_histogram(rho_hPiY_qutip)
        ax.set_title(r"$\rho_{\left(\frac{\pi}{2}\right)_y}$ via QuTiP Solver", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        
        #----------------------------
        # Compare for Free Precessiong time
        #----------------------------
        rho_free_euler = self.qm.euler_mesolve(rho0, 0.0, pars['tpr'], "FreeEvolution")
        rho_free_qutip = self.qm.mesolve_halfFreeEvolution(rho0, 0.0, pars['tpr']);
        #Fidelity
        rho_free_fid = fidelity(rho_free_euler, rho_free_qutip)
        print("State Fidelity(rho_qutip, rho_euler) in free evolution",rho_free_fid)
        
        fid_data.append(rho_free_fid);
        fid_head.append("rho_free_Evol")
        
        #Visualization of  density matrix via Euler Method and Qutip            
        fig, ax = matrix_histogram(rho_free_euler)
        ax.set_title(r"$\rho_{t_{p}}$ via Euler Method", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        
        fig, ax = matrix_histogram(rho_free_qutip)
        ax.set_title(r"$\rho_{t_{p}}$ via QuTiP Solver", fontsize=14)
        ax.view_init(azim=-55, elev=45)
        
        #Save the result in the destination unlsee filename is []
        if (filename == []):
            plt.show()
        else:
            # Open a file for fidelity
            fidFileObj = open('../' + filename +'/stateFid.out','w') # open from scratch
            from datetime import datetime
            fidFileObj.write('Date and time: ' + str(datetime.now()) + '\n\n')
            for i in range(len(fid_data)):
                fidFileObj.write(fid_head[i] + " = " + str(fid_data[i])  + '\n')
            fidFileObj.close();
            
            save_fig('../' + filename+'/rho_hPiX_qutip', ext='png', close=True, verbose=True)
            
            save_fig('../' + filename+'/rho_hPiX_euler', ext='png', close=True, verbose=True)
            save_fig('../' + filename+'/rho_fPiX_qutip', ext='png', close=True, verbose=True)
            save_fig('../' + filename+'/rho_fPiX_euler', ext='png', close=True, verbose=True)
            save_fig('../' + filename+'/rho_hPiY_qutip', ext='png', close=True, verbose=True)
            save_fig('../' + filename+'/rho_hPiY_euler', ext='png', close=True, verbose=True)
            save_fig('../' + filename+'/rho_free_euler', ext='png', close=True, verbose=True)
            save_fig('../' + filename+'/rho_free_qutip', ext='png', close=True, verbose=True)
        return;
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def qubit_optimal_control_halfX(self, tauXmax):
        """
        #--------------------------------------#
        # Unitary Propagator for (pi/2)x pulse #
        #--------------------------------------#
        Quantum Optimal Control (Bilinear model):
            1. http://qutip.org/docs/4.1/guide/guide-control.html
            2. Dong 2011, Quantum control theory and applications: a survey, arXiv
        Qubit Hamiltonian:
            HpiX=LT*sx + LT*(cos(2Wget)*sx+sin(2Wget)*sy) - LZ*cos(Wget)*sz
            where
            LT(t)=LTmax/(1-betarf/4)*fA(t,args) and LZ(t)=2LT(t)*Epsilon/Delta
        Unitary Propagator:
            i* dUqb/dt = HpiX * Uqb(t) where Uqb(0)=I
        """
        # ENSURE THAT COLLAPSE OPERATORS ARE EMPTY
        self.qm.collapse = [];
        
        # Qubit rotation in (pi/2)x pulse
        Up_qbX = rotation(sigmax(), np.pi/2);
        
        # Qubit Unitary evolution operator
        tlistX = np.linspace(0, tauXmax, int(pars['Num']/2));#end point not included
        #TRIED:NODIFFERENCE#tlistX = np.linspace(0, tauXmax, int(pars['Num']));#end point not included
        # Setup Hamiltonian for (pi/2)x pulse
        #H_qbX, Hargs_qb = self.qm.setup_HpiX_qb(0, tauXmax, 1.0);
        H_qbX, Hargs_qb = self.qm.setup_Hamiltonian(0, tauXmax, evolType='(pi/2)x', theSys='qubit')
        Us_qbX = propagator(H_qbX, tlistX, [], Hargs_qb, opts=self.qm.options);
        # Compute unitary operator fidelity
        # http://qutip.org/docs/4.1/guide/guide-control.html
        F_UX = 0.5 * (Us_qbX[-1].dag() * Up_qbX).tr()
        #TEST#print ("Qubit Fidelity for Unitary Evolution Operator in (pi/2)x pulse: %.4lf"%F_UY.real)
        
        # Qubit initial state
        rho_qb0 = self.qm.prepare_qubit_state(protocol='Ideal');
        # Qubit state evolution in (pi/2)x pulse
        rho_qbX = Up_qbX * rho_qb0 * Up_qbX.dag()
        # State evolution
        rhoS_qbX = Us_qbX[-1] * rho_qb0 * Us_qbX[-1].dag();
        # Compute state fideltiy
        F_SX = fidelity(rhoS_qbX, rho_qbX);
        #TEST#print ("Qubit Fidelity for State Evolution in (pi/2)x pulse: %.4lf"%F_SY)
        # Output
        return (F_UX.real, F_SX);
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def qubit_optimal_control_fullX(self, tauXmax):
        """
        #--------------------------------------#
        #  Unitary Propagator for (pi)x pulse  #
        #--------------------------------------#
        
        Quantum Optimal Control (Bilinear model):
            1. http://qutip.org/docs/4.1/guide/guide-control.html
            2. Dong 2011, Quantum control theory and applications: a survey, arXiv
            
        Qubit Hamiltonian:
            HpiX=LT*sx + LT*(cos(2Wget)*sx+sin(2Wget)*sy) - LZ*cos(Wget)*sz
            where
            LT(t)=LTmax/(1-betarf/4)*fA(t,args) and LZ(t)=2LT(t)*Epsilon/Delta
        Unitary Propagator:
            i* dUqb/dt = HpiX * Uqb(t) where Uqb(0)=I
        """
        # ENSURE THAT COLLAPSE OPERATORS ARE EMPTY
        self.qm.collapse = [];
        
        # Qubit rotation in (pi)x pulse
        sign = 1.0;
        Up_qbXf = rotation(sigmax(), sign*np.pi);
        # Qubit Unitary evolution operator
        tlistX = np.linspace(0, tauXmax, int(pars['Num']));#end point not included
        # Setup Hamiltonian for (pi)x pulse
        #H_qbX, Hargs_qb = self.qm.setup_HpiX_qb(0, tauXmax, sign);
        H_qbX, Hargs_qb = self.qm.setup_Hamiltonian(0, tauXmax, 
                                                    evolType='(pi)x', 
                                                    theSys='qubit', sign=1.0)
        Us_qbXf = propagator(H_qbX, tlistX, [], Hargs_qb, opts=self.qm.options);
        # Compute unitary operator fidelity
        F_UXf = 0.5 * (Us_qbXf[-1].dag() * Up_qbXf).tr()
        #TEST#print ("Qubit Fidelity for Unitary Evolution Operator in (pi)x pulse: %.4lf"%F_UY.real)
                
        # Qubit initial state
        rho_qb0 = self.qm.prepare_qubit_state(protocol='Ideal');
        # Qubit state evolution in (pi)x pulse
        rho_qbXf = Up_qbXf * rho_qb0 * Up_qbXf.dag()
        # State evolution
        rhoS_qbXf = Us_qbXf[-1] * rho_qb0 * Us_qbXf[-1].dag();
        # Compute state fideltiy
        F_SXf = fidelity(rhoS_qbXf, rho_qbXf);
        #TEST#print ("Qubit Fidelity for State Evolution in (pi)x pulse: %.4lf"%F_SY)
        
        # Output
        return (F_UXf.real, F_SXf);
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def qubit_optimal_control_halfY(self, tauYmax):
        """
        #--------------------------------------#
        # Unitary Propagator for (pi/2)y pulse #
        #--------------------------------------#
        
        Quantum Optimal Control (Bilinear model):
            1. http://qutip.org/docs/4.1/guide/guide-control.html
            2. Dong 2011, Quantum control theory and applications: a survey, arXiv
            
        Qubit Hamiltonian:
            HpiY=LT*sy + LT*(cos(2Wget-pi/2)*sx+sin(2Wget-pi/2)*sy) - LZ*cos(Wget-pi/2)*sz
            where
            LT(t)=LTmax/(1-betarf/4)*fA(t,args) and LZ(t)=2LT(t)*Epsilon/Delta
        Unitary Propagator:
            i* dUqb/dt = HpiY * Uqb(t) where Uqb(0)=I
        """
        # ENSURE THAT COLLAPSE OPERATORS ARE EMPTY
        self.qm.collapse = [];
        
        # Qubit rotation in (pi/2)y pulse
        Up_qbY = rotation(sigmay(), np.pi/2);
        
        # Qubit Unitary evolution operator
        tlistY = np.linspace(0, tauYmax, int(pars['Num']/2));#end point not included
        #TRIED:NODIFFERENCE#tlistY = np.linspace(0, tauYmax, int(pars['Num']));#end point not included
        # Setup Hamiltonian for (pi/2)y pulse
        #H_qbY, Hargs_qb = self.qm.setup_HpiY_qb(0, tauYmax);
        H_qbY, Hargs_qb = self.qm.setup_Hamiltonian(0, tauYmax, evolType='(pi/2)y', theSys='qubit')
        Us_qbY = propagator(H_qbY, tlistY, [], Hargs_qb, opts=self.qm.options);
        # Compute unitary operator fidelity
        F_UY = 0.5 * (Us_qbY[-1].dag() * Up_qbY).tr()
        #TEST#print ("Qubit Fidelity for Unitary Evolution Operator in (pi/2)y pulse: %.4lf"%F_UY.real)
        
        # Qubit initial state
        rho_qb0 = self.qm.prepare_qubit_state(protocol='Ideal');
        # Qubit state evolution in (pi/2)y pulse
        rho_qbY = Up_qbY * rho_qb0 * Up_qbY.dag()
        # State evolution
        rhoS_qbY = Us_qbY[-1] * rho_qb0 * Us_qbY[-1].dag();
        # Compute state fideltiy
        F_SY = fidelity(rhoS_qbY, rho_qbY);
        #TEST#print ("Qubit Fidelity for State Evolution in (pi/2)y pulse: %.4lf"%F_SY)
        
        # Output
        return (F_UY.real, F_SY);
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
