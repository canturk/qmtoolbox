###############################################################################
# ToolBox Quadrature Measurement Protocol for a Harmonic Oscillator (HO)      #
# Coupled to a Flux Qubit                                                     #
###############################################################################
# Author : Mehmet Canturk 
# Date   : March 13, 2017
###############################################################################
#                     DICTIONARY OF SYSTEM PARAMETERS                         #
###############################################################################
import numpy as np
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
# Dictionary for System Parameters
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
pars = {
'Np'    : 20,#Number of Pix pulses between Pix/2 and Piy/2
's'     : 25,#10, #Alternative to 'Nrep' which is used for analytical computation for pure dephasing
'fock'  : 0, #Initial Fock state (Default is Vacuum State)
'omegaR' : 2 * np.pi * 0.2,#HO frequency [GHz]
'Delta'  : 2 * np.pi * 4.0,#Constant for qubit gap [GHz]
'Epsilon': 2 * np.pi * 10.0,#Energy bias [GHz]
'g'      : 2 * np.pi * 0.002,#[GHz], Coupling strength
'sqrtAf' : 1e-6, #Dimensionless coefficient which is experimentally chosen from range
'v'      : 2, #Required in Wmin
'u'      : 10, #Required in Wmax
'tau2tp': 0.1, #Ratio of pulse duration to free precession time. tau2t=0.25 leads to LambdaX=1/2
'betarf': 0.1,  #Rate of rising and falling times for shaped pulses
'Q'     : 1e4,  #Quality factor  for HO
'T1'    : 10e-6,#[s] Relaxation time of qubit
'Tho'   : 15e-3,#Finite temperature effect in HO (resonator) in Kelvin
'Tqb'   : 50e-3,#Finite temperature effect in the qubit in Kelvin
'TM'    : 100.0,#[ns] Single Measurement Time
                #Measurement or readout pulse ranges from 100ns to a few microseconds.
'N'     : 50,   #Number of state levels in Hilbert Space for HO
'Num'   : 2**10 #Number of subintervals for numerical evolutions
};
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
if pars['Num']<200:
    raise ValueError('Nnum must be greater than or equallto 200');
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#


