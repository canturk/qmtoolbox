#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# ToolBox for Quadrature Measurement Protocol for a Harmonic Oscillator (HO)  #
# Coupled to a Flux Qubit                                                     #
###############################################################################
# Author : Mehmet Canturk 
# Date   : April 14, 2017
###############################################################################
# Class Definitions for Quadrature Measurement Protocol
###############################################################################
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
import numpy as np
import time
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
from qmtools.qpar import pars
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
from qmtools.qutil import fcx, fcy, fcz, fcAx, fcAy, fcAz, fA
from qmtools.qutil import fx, fy, fz, fx_conj, fy_conj, fz_conj
from qmtools.qutil import fEx, fEy, fEz

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

# Make qutip available in the rest of the code
from qutip import solver, mesolve, rotation, propagator, fidelity
from qutip import sigmax, sigmay, sigmaz, sigmap, sigmam, destroy
from qutip import thermal_dm, ket2dm, coherent, fock, basis, qeye, tensor
from qutip import expect, variance, commutator, qsave, qload
from qutip import matrix_histogram
from qutip.ui.progressbar import BaseProgressBar, TextProgressBar
from qutip.qobj import Qobj
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
import matplotlib.pyplot as plt
#-----------------------------------------------------------------------------#
# BinaryTree: https://pypi.python.org/pypi/binarytree/1.1.1
#from binarytree import Node, pprint, inspect
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

###############################################################################
class QuadMeasure(object):
    """
    Quadrature measurement is a type of QND detection protocol.
    QuadMeasure, a base class, defines a quadrature measurement protocol for HO 
	(resonator) coupled to a flux qubit.
    """
    
    def __init__(self, isLind=True):
        """
        Constructor initializes an instance (object)
        """
        try:
            print (self.__class__.__name__)
        except KeyError:
            print ("KeyError (arg4QM) encountered in class " + self.__class__.__name__)
        
        #---------------------------------------------------------------------#
        # Number of state levels in Hilbert Space for HO
        N = pars['N'];
        
        #
        # Setup Operators (Pauli and anhiliation operators) for Qubit-HO system
        # 
        # Tensor Product (Combining systems) describes the states of multipartite quantum systems
        # (e.g. two coupled qubits, qubit coupled to an oscillator, etc.)
        #
        # The operators acting on the state vectors in the combined 
        # Hilbert space (describing the coupled system) are formed by taking 
        #the tensor product of the individual operators.
        self.a   = tensor(qeye(2), destroy(N)) # A = I_{2x2}  \otimes a
        # Qubit Operators
        self.sx  = tensor(sigmax(), qeye(N)) # S_x = \sigma_x \otimes I_{NxN}
        self.sy  = tensor(sigmay(), qeye(N)) # S_y = \sigma_y \otimes I_{NxN}
        self.sz  = tensor(sigmaz(), qeye(N)) # S_z = \sigma_z \otimes I_{NxN}        
        
        #
        # The projection operator that projects state |q> \otimes |phi> on itself
        # for projective measurements. They are used to project qubit operators
        # without modifiying the state of HO.
        #
        self.PIg  = tensor(basis(2,0) * basis(2,0).dag(), qeye(N)); #|g> otimes |phi>
        self.PIe  = tensor(basis(2,1) * basis(2,1).dag(), qeye(N)); #|e> otimes |phi>
        
        #
        # State Preparation: Depending on the choice of Fock, or coherent or 
        # thermal state the ket psi_gnd will be set to initial state later in the associated derived
        # classes such as FockQMP, CoherentQMP, ThermalQMP.
        #

        # Operators for the conditional measurement
        self.quadI    =      (self.a.dag() + self.a); #Cosine component of Quadrature
        self.quadQ    = 1j * (self.a.dag() - self.a); #Sine component Quadrature
        
        #
        # Customize mesolver options
        self.options  = solver.Options();
        self.options.atol     = 1e-5;#Absolute tolerance
        self.options.rtol     = 1e-8;#Relative tolerance
        # When you set store_final_state=True in Options, 
        # the final state should be stored in result.final_state, 
        # not in result.states. This should probably be documented better.
        self.progress_bar = None;#OR#True
        #self.options.nsteps = 10000;
        #self.options.max_step=1e-4;
        
        print ("QuadMeasure Parameters were initialized!")
        
        #-------------------------------------------------
        # Setup collapse operators for Lindblad form
        # See http://qutip.org/docs/4.1/guide/dynamics/dynamics-time.html
        print("Quadrature measurement", end=" ")
        if isLind==True:
            print("with dissipation (collapse operators) and flux noise")
            self.collapse = self.setup_collapse_operators();
        else:
            print("with unitary evolution")
            self.collapse = [];#Default
        #--------------------------------------------------
        # Operators for expectation value in mesolve
        self.expectops = [];#DEFAULT
        #--------------------------------------------------
        # Generate synchronized real time intervals
        self.time = self.generate_RealTimeSequence(sync=True);
        # Convention used in qubit measurement
        self.qm = {'g':'+1','e':'-1'};#Qubit Measure
        return;#end
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def generate_RealTimeSequence(self, sync=True, disp=False):
        """
        date: 2017/04/14
        This module generates intervals of measurement time for a single quantum trajectory.
        Real-Time Quadrature Measurement of a Single-Photon Wavepacket with Continuous Temporal-Mode-Matching
        
        #-----------------------------+
        # Example:                    |
        #-----------------------------+
        from qmtools.qsim import QuadSim
        from qmtools.qpar import pars
        isLindblad  = True;#
        stateType  = 'Thermal';#'Fock';#
        qs  = QuadSim(stateType, isLindblad);
        t = qs.qm.generate_RealTimeSequence()
        Np = pars['Np'];
        t = qs.qm.time
        Np = pars['Np'];
        for m in range(pars['s']):
            mp = 2 * m * (Np + 2);#index
            print ("cpmg[%02d]: %.2lf"%(m+1,t[mp+1]-t[mp]),end=', ');
            print ("%.2lf"%(t[mp+2]-t[mp+1]),end=', ');
            n=1;
            while n<=Np:
                print ("%.2lf"%(t[mp+2*n+1]-t[mp+2*n]), end=', ');
                print ("%.2lf"%(t[mp+2*n+2]-t[mp+2*n+1]), end=', ');
                n += 1;
            print ("%.2lf"%(t[mp+2*Np+3]-t[mp+2*Np+2]), end='\n');
        #------------------------------------
        Output:
        #------------------------------------
                  (phx) hfree  +pfx ffree -pfx  ffree +pfx  ffree -pfx  ffree +pfx  ffree -pfx  ffree +pfx  ffree -pfx  ffree +pfx  ffree -pfx ffree  +pfx  ffree -pfx ffree  +pfx  ffree -pfx ffree  +pfx  ffree -pfx ffree  +pfx  ffree -pfx ffree  +pfx  ffree -pfx hfree  phy
        cpmg[01]: 0.11, 1.14, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 1.14, 0.11
        cpmg[02]: 0.11, 1.14, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 1.14, 0.11
                :
        cpmg[25]: 0.11, 1.14, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 2.27, 0.23, 1.14, 0.11
        #-----------------------------
        
        PULSE SYNCHRONIZATION:
        2017/03/19-1230pm
        Hi,
        Reforwarding two emails sent yesterday, cant tell what got lost.
        The measurement of one quadrature dependens on timing. 
        All the sequences have to start at multiples of T_HO, 
        otherwise at each measurement you look at a rotated quadrature.
        I can talk for a bit if this is useful.
        Best,
        Adrian
        
        2017/03/19-342pm
        Thanks for the information and the emails. 
        Please see the attached sketch that represents a single CPMG pulse 
        and a single Measurement. 

        0  1    2    3        4    5    6  7                     8
        +--+    +----+        +----+    +--+ +++++++++++         +--+    +----+
        |  |    |    |        |    |    |  | ||  <I>  ||         |  |    |    |
        +  +----+    +--------+    +----+  + +++++++++++         +  +----+    +
        | <---------------- T_cpmg ------> | |<- T_M ->|<-delay->|
        | <--------- eta * Tr = T_cpmc + T_M + delay ----------->|
        
        In the scheme,
        T_CPMG = Te + tau corresponds to total time for a single CPMG and 
        T_CPMG + T_M for the duration of CPMG and mesurement.
        Which of the quantities must be multiple of T_HO = 2p/omega_r?
        T_CPMG or T_CPMG + T_M?
        
        """
        Np   = pars['Np'];
        s    = pars['s'];       
        tau  = pars['tau'];#Duration of a (PI)x pulse is assumed to be ZERO
        TM   = pars['TM'];#Single measurement time
        tpr   = pars['tpr'];#precession time
        Tr   = 2*np.pi / pars['omegaR'];print ("Tr=",Tr," [ns]")
        
        #
        #
        #The total time trajectory: T_trj = s (Te + TM + tau + dt) = eta * 2p/omegaR
        # See Computation Notebook page 89-90
        #-------------------------------------------
        # Create an array for Real Time Measurement including the last meassurement time
        #-------------------------------------------
        t =  np.empty( 2 * (s - 1) * (Np + 2) + 2 * Np + 4 + 1);
        #--------------------------------------------
        # Fill the element of CPMG[1]
        #--------------------------------------------
        t[0] = 0.0;          # pre-(pi/2)x
        t[1] = tau / 2;      # pos-(pi/2)x
        t[2] = t[1] + tpr / 2;# pre (pi)x_1
        n=1;
        while (n<=Np-1):
            t[2*n+1] = t[2*n]   + tau;#pos (pi)x_n
            t[2*n+2] = t[2*n+1] + tpr; #
            n = n + 1;
        t[2*Np+1] = t[2*Np]   + tau;
        t[2*Np+2] = t[2*Np+1] + tpr / 2;
        t[2*Np+3] = t[2*Np+2] + tau / 2;
        
        #-----------------------------------------------
        # Fill the rest of the time interval (s-1) in the form of CPMG pulse sequence
        #-----------------------------------------------
        m = 1; # CPMG[2]
        while (m <= (s-1)):
            mp = 2 * m * (Np + 2);#index
            t[mp] = t[mp-1] + TM; # T_CPMG = Te + TM + tau
            
            # Synchronization checking
            # Check whether t[mp] is multiple of Tr=2pi/omegaR
            if sync is True:
                #TEST#print ("Pre-value t[mp]=t[%d]=%10.5lf"%(mp, t[mp]));
                eta_raw = t[mp] / Tr;
                #eta     = np.ceil(eta_raw);
                #print ("eta=t[%4d]/Tr=%.5lf/%.5lf=%.5lf, ceil(eta)=%4d"%(mp, t[mp],Tr,eta_raw, eta))
                eta     = np.floor(eta_raw);
                #print ("eta=t[%4d]/Tr=%.5lf/%.5lf=%.5lf, floor(eta)=%4d"%(mp, t[mp],Tr,eta_raw, eta))
                # Assign adjusted value
                t[mp]   = eta * Tr;      # pre-(pi/2)x
                #TEST#print("Pos-value t[mp]=t[%d]=%10.5lf"%(mp, t[mp]))
            t[mp+1] = t[mp]   + tau / 2;# pos-(pi/2)x
            t[mp+2] = t[mp+1] + tpr  / 2;
            n = 1;
            while (n<=Np-1):
                t[mp+2*n+1] = t[mp+2*n]   + tau;
                t[mp+2*n+2] = t[mp+2*n+1] + tpr;
                n = n + 1;
            t[mp+2*Np+1] = t[mp+2*Np]   + tau;
            t[mp+2*Np+2] = t[mp+2*Np+1] + tpr / 2;
            t[mp+2*Np+3] = t[mp+2*Np+2] + tau / 2;
            m = m + 1;
        #Add the last measurement duration
        t[-1] = t[len(t)-2] + pars['TM'];#t[len(t)-1]==t[-1]
        print("Synchronized real time sequence including TM is %.3lf [ns]"%t[-1])
        
        ##Update System parameters
        #pars['Traj'] = t[-1];
        #print ("Ttraj = %10.4lf [ns]"%pars['Traj'])
        
        """Display the list"""
        if disp is True:
            for n in range(len(t)):
                print ("t[%4d]=%10.6lg"%(n, t[n]));
        return t;
    #end:generate_RealTimeSequence(self,...)
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def prepare_qubit_state(self, protocol='Actual'):
        # Prepare initial qubit state 
        # This module is called at the beginning of qmProtocol()
        if protocol is 'Actual':
            pg = pars['nqb'];
            pe = pg * np.exp(-pars['Wge2Tqb']);
            print ("pg=%.5lf, pe=%.5lf"%(pg, pe));
            return pg * basis(2,0) * basis(2,0).dag() + pe * basis(2,1) * basis(2,1).dag();
        elif protocol is 'Ideal':
            return basis(2,0) * basis(2,0).dag();
    #end:prepare_qubit_state(self,...)
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def setup_collapse_operators(self):
        
        print ("Setup Collapse Operators in Lindblad form", end=" ");
        #Reference: http://qutip.org/docs/4.1/guide/dynamics/dynamics-time.html        
        collapse = [];#Empty list for collapse operators
        # Collapse coefficients from the dictionary
        k_Dwn = pars['kDwn'];
        k_Up  = pars['kUp'];
        G_e2g = pars['Gamma_E2G'];
        G_g2e = pars['Gamma_G2E'];
        print ("for the parameters")
        print("kDwn = 2pi * %.3lg [GHz],\nkUp  = 2pi * %.3lg [GHz]"\
              %(k_Dwn/(2*np.pi), k_Up/(2*np.pi)))
        print("Ge2g = %.3lg   [GHz],\nGg2e = %.3lg [GHz]\n"%(G_e2g, G_g2e))
        #------------------------------------------------
        # Collapse operators in Lindblad Master Equation
        #------------------------------------------------
        a   = self.a;
        sp  = tensor(sigmap(), qeye(pars['N'])) # sigma+
        sm  = tensor(sigmam(), qeye(pars['N'])) # sigma-
        
        #if isLindblad == True:
        collapse.append(np.sqrt(k_Dwn) * a);#[sqrt(2p*Gzh)]
        collapse.append(np.sqrt(k_Up ) * a.dag());#[sqrt(2p*Gzh)]
        collapse.append(np.sqrt(G_e2g) * sm);#[sqrt(2p*Gzh)]
        collapse.append(np.sqrt(G_g2e) * sp);#[sqrt(2p*Gzh)]
        #DISREGARDED#collapse.append(np.sqrt(Gamma_varphi) * sz);#[dimensionless]
        return collapse;
    #end:setup_collapse_operators(self)
    #
    # In[] - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def setup_Hamiltonian(self, tin, tfn, **kwargs):
        """
        Setup the Hamiltonian and required arguments for call-back functions
        in case of pulse durations e.g. (pi)x and (pi/2)x/y pulses as well as 
        free evolution time. If the pulse is (pi/2)x then sign is always +1
        else if the pulse is (pi)x then sign is either +1 or -1.
        
        The Hamiltonian:
        $\mathcal{H}(t) = \frac{\Lambda_\perp}{2}\left\{\cos\varphi\sigma_x
        - \sin\varphi\sigma_y\right\} + \frac{\Lambda_\perp}{2}\cos\left(2\omega_{ge}t
        + \varphi\right)\sigma_x + \frac{\Lambda_\perp}{2}\sin\left(2\omega_{ge}t 
        + \varphi\right)\sigma_y - \Lambda_z\cos\left(\omega_{ge}t 
        + \varphi\right)\sigma_z - \lambda_\perp\left\{a f_x^*(t) 
        + a^\dagger f_x(t)\right\}\sigma_x - \lambda_\perp\left\{a f_y^*(t) 
        + a^\dagger f_y(t)\right\}\sigma_y + \lambda_z    \left\{a f_z^*(t)
        + a^\dagger f_z(t)\right\}\sigma_z$ where $f_z(t)= e^{i\omega_r t}$, 
        $f_x(t) = f_z(t)\cos(\omega_{ge}t)$, and $f_y(t)= f_z(t) \sin(\omega_{ge}t)$.
        The pulse amplitudes are $\Lambda_\perp=A \frac{\Delta}{\omega_{ge}}$ 
        and $\Lambda_z=A \frac{\varepsilon}{\omega_{ge}}$.
        #-------------------------------------------
        # Ht = -lt * (a * conj(fx) + a.dag * fx) * sx
        #      -lt * (a * conj(fy) + a.dag * fy) * sy
        #      +lz * (a * conj(fz) + a.dag * fz) * sz
        #      +LT/2 * (cos(phi) * sx - sin(phi) * sy)
        #      +LT/2 * fcx * sx 
        #      +LT/2 * fcy * sy 
        #      -LZ   * fcz * sz
        # lt:lambda_t, lz:lambda_z, LT:Lambda_T, LZ:Lambda_Z
        # See http://qutip.org/docs/4.1/guide/dynamics/dynamics-time.html
        #print ("Setup Hamiltonian", end=" ")
        """
        if not(kwargs):
            raise KeyError("**kwargs=None")
                
        #-------------------------------------------
        # Input Parameters
        #-------------------------------------------
        
        # check its existance
        if 'evolType' in kwargs:
            evolType = kwargs.get('evolType');
        else:
            raise KeyError('evolType is not present');
        
        # modulate sign if the evolution is (pi)x pulse
        if evolType == '(pi)x':
            try:
                sign = kwargs.get('sign');
                LT   = sign * pars['LambdaT'];
                LZ   = sign * pars['LambdaZ'];
                #print ("sign=",sign, ' in ', evolType);
            except:
                raise KeyError("sign is not present for (pi)x pulse")
        else:
            #Default values
            LT = pars['LambdaT'];
            LZ = pars['LambdaZ'];
        
        #Default value
        varphi = 0.0;#Pulse phase in control pulse sinusoid
        
        # phase setting if the evolution is (pi/2)y
        if evolType == '(pi/2)y':
            varphi = -np.pi/2;
        #if evolType == '(pi/2)x' or evolType == '(pi)x' or evolType == '(pi/2)y':
        #   print ("varphi=",varphi, ' in ', evolType);
            
        lt = pars['lambdat'];
        lz = pars['lambdaz'];
        
        dt = pars['deltat'];
        dz = pars['deltaz'];
        
        #-----------------------------------------------------
        # Initialize dictionary for Hamiltonian functions
        Hargs = {}
        Hargs['varphi']   = varphi; 
        Hargs['omegaR']   = pars['omegaR'];
        Hargs['omegaGE']  = pars['omegaGE'];
        Hargs['2omegaGE'] = 2 * pars['omegaGE'];
        Hargs['Wmin']     = pars['Wmin'];
        Hargs['Nf']       = pars['Nf'];
        
        sx = self.sx;
        sy = self.sy;
        sz = self.sz;
        a  = self.a;
        
        #Set the default value for theSys label for interaction part of the Hamiltonian
        theSys = 'qubit-HO';
        
        if 'theSys' in kwargs:
            theSys = kwargs.get('theSys');
            #print ('theSys=', theSys, ' in ', evolType);
            
        if evolType=='free' and theSys == 'qubit':
            raise ValueError("Invalid Hamiltonian for qubit system");
            
        if theSys == 'qubit':
            Hop = [];
            # setting qubit operators
            sx = sigmax();
            sy = sigmay();
            sz = sigmaz();
        elif theSys == 'qubit-HO':
            Hop= [[-lt * sx * a.dag(), fx],
                  [-lt * sy * a.dag(), fy],
                  [+lz * sz * a.dag(), fz],
                  [-lt * sx * a,       fx_conj],
                  [-lt * sy * a,       fy_conj],
                  [+lz * sz * a,       fz_conj]
                  ];
        else:
            raise KeyError("theSys must be set to either 'qubit' or 'qubit-HO'")
            
        #Flux noise part of the Hamiltonian for free evolution
        if kwargs.get('evolType') == 'free' and self.collapse != []:
            ##CONSIDER REVISING ('m' in kwargs or 'n' in kwargs) FOR m & n BELOW
            if ('m' in kwargs or 'n' in kwargs):
                m = kwargs.get('m');
                n = kwargs.get('n');
            else:
                raise KeyError("Invalid m & n for flux noise data file")
                
            if ('dirname' in kwargs) and (kwargs.get('dirname') != []):
                dirname = kwargs.get('dirname');
            else:
                raise KeyError("Invalid dirname for flux noise data")
            data_xyz = np.loadtxt("../" + dirname + "/cpmg%02dn%02d.out"%(m + 1, n))
            #print ('m=',m, ', n=',n)#(data_xyz)
            
            # To turn these data points into a function we call the QuTiP qutip.interpolate.Cubic_Spline class 
            # using the first and last domain time points, t[0] and t[-1], respectively, 
            # as well as the entire array of data points:
            from qutip.interpolate import Cubic_Spline
            cSx = Cubic_Spline(data_xyz[0][0],data_xyz[0][-1], data_xyz[1])
            cSy = Cubic_Spline(data_xyz[0][0],data_xyz[0][-1], data_xyz[2])
            cSz = Cubic_Spline(data_xyz[0][0],data_xyz[0][-1], data_xyz[3])
            Hop.append([+dt * sx, cSx]);
            Hop.append([+dt * sy, cSy]);
            Hop.append([-dz * sz, cSz]);
            
        # Control part of the Hamiltonian for pulse sequence
        if (evolType=='(pi/2)x' or evolType=='(pi)x' or evolType=='(pi/2)y'):
            if (pars['betarf']==0.0):
                Hop.append( +LT/2 * (np.cos(varphi)*sx - np.sin(varphi)*sy));
                Hop.append([+LT/2 * sx, fcx]);
                Hop.append([+LT/2 * sy, fcy]);
                Hop.append([-LZ   * sz, fcz]);
            else:
                # Parameters for shaped (pi)x and (pi/2)x pulses
                Hargs['tb'] = tin;#tbegin (MUST BE SET FOR EACH HAMILTONIAN SETUP)
                Hargs['te'] = tfn;#tend
                Hargs['tr'] = (Hargs['te'] - Hargs['tb']) * pars['betarf'] / 4.0; #rising endge
                Hargs['tf'] = Hargs['tr'];#falling edge
                Hargs['tp'] = (Hargs['te'] - Hargs['tb']) * (1 - pars['betarf'] / 2.0);#Peak
                Hargs['pi2tf'] = np.pi / Hargs['tf'];#pi/tfall
                Hargs['pi2tr'] = np.pi / Hargs['tr'];#pi/trise
                #print ("(ti,tf)=","(",tin,Hargs['te'],") trise=",Hargs['tr'], ", pi/tr=", Hargs['pi2tr'])
                
                LT_prime = LT / (1 - pars['betarf'] / 4.0);
                LZ_prime = LZ / (1 - pars['betarf'] / 4.0);
                #print ("for shaped (pi/2)x or (pi)x pulse")
                Hop.append(
                        [+LT_prime / 2 * (np.cos(varphi) * sx - np.sin(varphi) * sy), fA])
                Hop.append(
                        [+LT_prime / 2 * sx,    fcAx])
                Hop.append(
                        [+LT_prime / 2 * sy,    fcAy])
                Hop.append(
                        [-LZ_prime * sz,    fcAz])
        return (Hop, Hargs);
    #end: setup_Hamiltonian
    #
    # In[] - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def solveME(self, rho0, ti, tf, **kwargs):
        """
        This function constructs the Hamiltonian and solves the Lindblad master equation 
        with or without collapse operators in the (pi/2)x, (pi)x, (pi/2)y pulse durations
        and free evolution time:
        
        +-------+               +--------------+               +-------+
        |       |               |              |               |       |
        |(pi/2)x|               |     (pi)x    |               |(pi/2)y|
        +-------+---------------+--------------+---------------+-------+---->

        Inputs: 
        rho0  : tensor(rho_qb, rho_HO)
        ti, tf: denotes initial and final time in real time duration
        kwargs: keyword arguments
        """
        #print ("Number of sub-intervals for (pi/2)y pulse duration: %d"%Num);
        
        # Time list for (pi/2)x pulse
        tlist = np.linspace(ti, tf, pars['Num']);
        
        #-------------------------------------------
        # Customize mesolver options
        options  = self.options
        
        #------------------------------
        # Retreive collapse operators
        #------------------------------
        cops = self.collapse;#Test#print("Collapse Operators",self.collapse)
        # Setup Hamiltonian for (pi/2)y pulse
        Hops, Hargs = self.setup_Hamiltonian(ti, tf, **kwargs);
        #-------------------------------------------------
        # Evolve system in (pi/2)x or free precess or (pi)x or (pi/2)y pulse duration
        # Master equation evolution of a density matrix for given Hamiltonian
        #-------------------------------------------------
        if self.expectops == []:
            options.store_states = False
            options.store_final_state = True
            out  = mesolve(Hops, rho0, tlist, cops, [], Hargs, options, progress_bar=self.progress_bar);
            rho  = out.final_state;#rho  = out.states[-1];
            return rho;
        else:
            options.store_states = True;
            out  = mesolve(Hops, rho0, tlist, cops, self.expectops, Hargs, options, progress_bar=self.progress_bar);
            return (tlist, out.states[-1], out.expect)
    #end: solveME(self,...)
    #
    # In[] - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def evolve_densitymatrix_in_CPMG(self, m, rho0):
        """
        Implementation of Lindblad master equation, based on Markovian approach
        (weak coupling) is carried out by considering the collapse operator
        including the influence of flux noise in the Hamiltonian.
        
        This module performs the numerical evolution of the total Hamiltonian
        within CPMG pulse sequence.
        
        Evolution for each pulse duration or the precession time duration between 
        two adjacent/consecutive pulses can be handled by using QuTiP.mesolve
        that calculates the evolution of the density matrix.
        
        The mesolve in QuTiP is capable of handling time-dependent Hamiltonians 
        and the terms with collapse operators.
        
        Note: This module is capable of simulating an arbitrary number of Np pulses.
        
        The first CPMG and measurement
        
        (pi/2)x       +(pi)x                -(pi)x       (pi/2)y
         +---+      +--------+            +--------+      +---+
         |   |      |        |            |        |      |   |++++++++++
         |   |      |        |            |        |      |   |    M1   +
         +---+------+--------+------------+--------+------+---+---------+
         t0  t1     t2       t3           t4       t5     t6  t7         


        The second CPMG and measurement
        
        (pi/2)x       +(pi)x                -(pi)x       (pi/2)y
         +---+      +--------+            +--------+      +---+
         |   |      |        |            |        |      |   |++++++++++
         |   |      |        |            |        |      |   |    M2   +
         +---+------+--------+------------+--------+------+---+---------+
         t8  t9     t10     t11          t12      t13    t14  t15
        """
        #--------------------------------------------------
        # Generate subsequent measurement time interval (synchronized)
        #--------------------------------------------------
        t = self.time;#The synchronized time is generated during parameter setup
        
        Np  = pars['Np'];
        # Index of m
        mp  = 2 * m * (Np + 2);
        
        # Number of subintervals
        NumPro = pars['Num'];
        #--------------------------------------------------------
        # Evolution in (pi/2)x between 0 and t0 for m=0; t[mp]-tau/2 and t[m] for m!=0
        #--------------------------------------------
        pars['Num'] = int(NumPro/2);#test#print("For (pi/2)x: Num=",pars['Num'])
        rho = self.solveME(rho0, t[mp], t[mp+1], evolType='(pi/2)x', theSys='qubit-HO');
        #rho = self.evolve_PiX(rho0, t[mp], t[mp+1], evolType='(pi/2)x', theSys='qubit-HO');
        
        #--------------------------------------------
        # Evolution in free precession time (half) between t[mp] and t[mp+1] 
        #--------------------------------------------
        if self.collapse == []:
            pars['dirname'] = [];
            
        pars['Num'] = int(NumPro);#test#print("For (free/2): Num=",pars['Num'])
        mm=m;
        rho = self.solveME(rho, t[mp+1], t[mp+2], evolType='free', dirname=pars['dirname'], theSys='qubit-HO', n=0, m=mm)
        #rho = self.evolve_free(rho, t[mp+1], t[mp+2], evolType='free', dirname=pars['dirname'], theSys='qubit-HO', n=0, m=mm);
                
        sgn = 1.0;
        n   = 1;
        while n<=Np:
            #--------------------------------------------
            # Evolution in (pi)x between t[mp+2n-1] and t[mp+2n]
            #--------------------------------------------
            pars['Num'] = int(NumPro);#test#print("For (pi)[",n,"]x: Num=",pars['Num'])
            rho = self.solveME(rho, t[mp+2*n], t[mp+2*n+1], 
                               evolType='(pi)x', theSys='qubit-HO', sign=sgn)
            #rho = self.evolve_PiX(rho, t[mp+2*n], t[mp+2*n+1], evolType='(pi)x', theSys='qubit-HO', sign=sgn)
            
            #------------------------------------------
            # Evolution in free precessing time between t[mp+2n] and t[mp+2n+1]
            #------------------------------------------
            if n == pars['Np']:
                pars['Num'] = int(NumPro);#test#print("For (free/2): Num=",pars['Num'])
            else:
                pars['Num'] = int(2*NumPro);#test#print("For (free): Num=",pars['Num'])
            
            mm = m;
            nn = n;
            rho = self.solveME(rho, t[mp+2*n+1], t[mp+2*n+2], evolType='free', 
                               dirname=pars['dirname'], theSys='qubit-HO', n=nn, m=mm);
            #rho = self.evolve_free(rho, t[mp+2*n+1], t[mp+2*n+2], evolType='free', dirname=pars['dirname'], theSys='qubit-HO', n=nn, m=mm);
            
            # Sign update for the next adjacent (pi)x pulse
            sgn = - sgn;
            #
            # Update n for next evolution
            n += 1;
        
        # Evolution in (pi/2)y pulse
        pars['Num'] = int(NumPro/2);#test#print("For (pi/2)y: Num=",pars['Num'])
        rho = self.solveME(rho, t[mp+2*Np+2], t[mp+2*Np+3], evolType='(pi/2)y', theSys='qubit-HO');
        #rho = self.evolve_PiY(rho, t[mp+2*Np+2], t[mp+2*Np+3], evolType='(pi/2)y', theSys='qubit-HO');
        
        #Invert Number of subintervals
        pars['Num'] =  int(NumPro);
        #------------------------------------------------
        # Fidelity check for rho in CPMG and rho_avg
        #------------------------------------------------
        if (m==0):
            Uavg   = rotation(self.sx, np.pi/2);
            #varphi= 2.0 / np.pi * param['lambdaz'] * param['Te'];
            varphi = 2.0 / np.pi * pars['g'] * pars['Te'];
            #Time averaged Hamiltonian
            Up = (-1j * varphi * self.sz * (self.a + self.a.dag())).expm();
            Uavg = Up * Uavg;
            if pars['Np']%2 != 0:
                Uavg = rotation(self.sx, np.pi) * Uavg;
            Uavg = rotation(self.sy, np.pi/2) * Uavg;
            
            print ("State Fidelity(rho,rho_avg)=%.3lf where rho_avg=Uavg*rho0*Uavg.dag()"\
                   %fidelity(rho, Uavg * rho0 * Uavg.dag()));
        
        return rho;
    # end: evolve_densitymatrix_in_CPMG(self,..)
    #
    # In[] 
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def qmprotocol_Ideal(self, filename):
        """
        Date: 2017/04/12
        
        Ideal quadrature measurement using Algoritm 1 
        in the latest long note (TheoryQuadDetection\notes\20170101_quadMeasure_Notes)
        """
        # Just to make sure the operators for expectation valuee in mesolve are NONE
        self.expectops = [];#DEFAULT
        #--------------------------------------------------
        # Open a file for fidelity
        #How to make python write instantly to the output file: 
        #   use flush() or set the file object to be unbuffered by using 0
        logFile = open('../' + filename +'.log','w', 0)# 0 means there is no buffer, so all output
                                                       # will be auto-flushed
                                                       
        
        from datetime import datetime
        logFile.write('Created on ' + str(datetime.now()) + '\n\n')
        #---------------------------------------------------
        # State Preparation: Initialize density matrix into ground state
        #----------------------------------------------------
        rho_qb = self.prepare_qubit_state('Ideal');
        rho0 = tensor(rho_qb, self.rho_HO);
        #--------------------------------------------------------
        # Create an empty list for expectation value and variance of quadrature I
        expI_ideal = [];
        varI_ideal = [];
        # Initial expectation value and variance.
        expI_ideal.append (expect(  self.quadI, rho0));
        varI_ideal.append (variance(self.quadI, rho0));
        
        qm = self.qm;
        # Print out and write the first data into log file
        out="Ideal[00]: r=%s  <I>=%2.5lf,  DeltaI=%.5lf"%(qm['g'], expI_ideal[0], varI_ideal[0])
        print (out);logFile.write( out + '\n')
        #print ("Ideal:r[00]=%s  <I>=%2.5lf,  DeltaI=%.5lf"%(qm['g'], expI_ideal[0], varI_ideal[0]));
        # Projection operatos
        PI_g        = self.PIg;
        PI_e        = self.PIe;
        # QND type measurement sequence
        for m in range(pars['s']):#parSim['Nrep'] was replaced by parSim['s']
            # ME evolves the density matrix governed by Htotal
            #rho = self.evolution_of_densitymatrix_in_CPMG(m, rho0)
            rho = self.evolve_densitymatrix_in_CPMG(m, rho0)
            
            #---------------------------------------------------
            # Conditional measurement: Probability of ground and excited states
            Pg = expect(PI_g, rho);
            Pe = expect(PI_e, rho);
            
            # Use random number between 0 and 1 representing the probability
            r = np.random.random();#rnd_list[m];#0<r<1
            
            #State projection
            if r < Pg:
                # Density matrices for g
                rho   = PI_g * rho * PI_g.dag() / Pg;# state is projected by PIg
                rho0  = rho;
                rr    = qm['g'];
            else:
                # Density matrices for e
                rho  = PI_e * rho * PI_e.dag() / Pe;
                # Flip the qubit
                rho0  = self.sx * rho * self.sx;
                rr = qm['e'];
            
            # Ideal measurements
            expI_ideal.append(expect  (self.quadI, rho));
            varI_ideal.append(variance(self.quadI, rho));
            
            out = "Ideal: r[%2d]=%s  <I>=%6.3lf,  varI=%6.3lf,  Pg=%0.2lf,  Pe=%0.2lf"\
            %(m+1,rr, expI_ideal[-1], varI_ideal[-1], Pg, Pe);
            print (out);logFile.write( out + '\n');
            #print ("Ideal:r[%2d]=%s  <I>=%2.4lf,  varI=%.5lf,  Pg=%0.2lf,  Pe=%0.2lf"%(m+1,rr, expI_ideal[-1], varI_ideal[-1], Pg, Pe));
            
        # Close log file at the end
        logFile.close()
        return (expI_ideal, varI_ideal);
    # end: qmprotocol_Ideal(self, filename)
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
    def qmprotocol_Algorithm3(self, gama, delta, alfa, beta, filename):
        """
        Date: 2017/04/14
        
        Actual (realistic) quadrature measurement protocol using Algoritm 3 
        in the latest long note (TheoryQuadDetection\notes\20170101_quadMeasure_Notes).
        
        # LOG FILE:
            from datetime import datetime
            datetime_now = str(datetime.now());
            logFile = open('../' + filename +'datetime_now.log','w') # open from scratch
            logFile.write('Date and time: ' + datetime_now + '\n\n')
            # In the algorithm
            logFile.write("Same as in the print line below \n")
            # At the end
            fidFileObj.close();
        
        """
        # Just to make sure the operators for expectation valuee in mesolve are NONE
        self.expectops = [];#DEFAULT
        #--------------------------------------------------
        # Open a file for fidelity
        logFile = open('../' + filename +'.log','w')
        
        from datetime import datetime
        logFile.write('Created on ' + str(datetime.now()) + '\n\n')
        logFile.flush()
        #---------------------------------------------------
        # State Preparation: Initialize density matrix into ground state
        #----------------------------------------------------
        rho_qm = self.prepare_qubit_state('Actual');
        rho0   = tensor(rho_qm, self.rho_HO);
        #
        # Create an empty list for expectation value and variance of quadrature I
        expI = [];
        varI = [];
        # Initial expectation value and variance.
        expI.append(expect(self.quadI, rho0));
        varI.append(variance(self.quadI, rho0));
        
        qm = self.qm;
        
        # Print out and write the first data into log file
        out="Alg3[00]: r=%s  <I>=%6.4lf,  DeltaI=%7.5lf"%(qm['g'], expI[0], varI[0])
        print (out);
        logFile.write( out + '\n')
        logFile.flush()
        
        # Projection operatos
        PI_g = self.PIg;
        PI_e = self.PIe;

        # QND type measurement sequence
        for m in range(pars['s']):
            
            # Evolution of density matrix
            #rho = self.evolution_of_densitymatrix_in_CPMG(m, rho0)
            rho = self.evolve_densitymatrix_in_CPMG(m, rho0)
            
            #---------------------------------------------------
            # Conditional measurement: Probability of ground and excited states
            Pg    = expect(PI_g, rho);
            Pe    = expect(PI_e, rho);
            
            # Density matrices for g or e
            rho_g = PI_g * rho * PI_g.dag() / Pg;# state is projected by PIg
            rho_e = PI_e * rho * PI_e.dag() / Pe;
            
            # Density matrix by projection errors
            rho_tilde_g = (1.0 - gama)  * rho_g + gama  * rho_e;
            rho_tilde_e = (1.0 - delta) * rho_e + delta * rho_g;
            
            # Generate a random binary number (0 for 'g' and 1 for 'e')
            # that mimics the quantum jump in qutbit and then use
            # weighted random choice
            
            print("P[%02d](%.3lf, %.3lf), Sum(%.25lf)"%(m+1,Pg,Pe,Pg+Pe))
            #TRY THIS!#from random import random choice
            r  = np.random.choice(['g','e'], p=[Pg, Pe]);
            # Imperfect state evolution
            if r == 'g':
                # Actual result observed in meter is randomly selected
                # from {g,e} using distribution {1-alpha, alpha}
                rprime = np.random.choice(['g','e'], p=[1-alfa, alfa]);
                #test#print("r[g][%02d] = %d "%(m+1, qm[rprime]))
            elif r == 'e':
                # Actual result observed in meter is randomly selected
                # from {g,e} using distribution {beta, 1-beta}
                rprime = np.random.choice(['g','e'], p=[beta, 1-beta]);
                #test#print("rp[e][%02d] = %d "%(m+1,qm[rprime]))
            else:
                raise ValueError("Ideal choice must be either 'g' or 'e'!");
                
            #Actual ReadOut: Expectation value and variance
            if rprime == 'g':
                expI.append(expect  (self.quadI, rho_tilde_g));
                varI.append(variance(self.quadI, rho_tilde_g));
                # Update projected  qubit state for the next iteration
                rho0 = rho_tilde_g;#state update 
                #--------------------------------------
            elif rprime == 'e':
                expI.append(expect  (self.quadI, rho_tilde_e));
                varI.append(variance(self.quadI, rho_tilde_e));
                # Update projected  qubit state for the next iteration
                rho0  = self.sx * rho_tilde_e * self.sx;#flip qubit state#
                #rho0 = tensor(rho_qm, rho_tilde_e.ptrace(1));
                #--------------------------------------
            else:
                raise ValueError("Actual choice must be either 'g' or 'e'!");
            #--------------------------------------
            #print(rho0.ptrace(0));
            #---------------------------------------
            # print out the instant data and write into log file
            out = "Alg3[%02d]: (r=%s, r'=%s)=%s,  Pg=%0.2lf,  Pe=%0.2lf,  expI=%2.4lf,  varI=%.4lf"\
            %(m+1, r, rprime, qm[rprime]==qm[r], Pg, Pe, expI[-1], varI[-1]);
            print (out);
            logFile.write(out + '\n');
            logFile.flush()
        
        # Close log file at the end
        logFile.close()
        return (expI, varI);
    # end: qmprotocol_Algorithm3(self,...)
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
    #
###############################################################################
#                          End of QuadMeasure                                 #
###############################################################################
class ThermalQM(QuadMeasure):
    """
    ThermalQMP, inherited from QuadMeasure, describes either 
	thermal inital state or initial thermal displacement state of HO
    at a given temperature
    """
    def __init__(self, isLindblad):
        """
        Construct an instance of class ThermalState.
        """
        super(ThermalQM, self).__init__(isLindblad)
        
        #
        #Some physical constants
        #The number of average photons in thermal equilibrium for the resonator
        #with the given frequency and temperature. 
        #
        print ("nHO=%.5g at thermal equilibrium"%pars['nHO']);
        # Density matrix of HO at thermal equlibrium
        self.rho_HO = thermal_dm(pars['N'], pars['nHO']);
        #print ('Thermal State: ' + str(self.rho_HO))
        #if isinstance(self.rho_HO, Qobj):raise TypeError("rho_HO is not a quantum object")
        
    ###########################################################################
    def __del__(self):# Constructor of ThermalQMP
        """
        This destructor prints the class name of an instance that is about to be destroyed
        """
        class_name = self.__class__.__name__;
        print (class_name, "destroyed");
###############################################################################
#                           End of ThermalQM                                  # 
###############################################################################
