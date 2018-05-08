"""
Date: 20170323-1600
"""
import time
import numpy as np
#from qutip import *
from qmtools.qutil import plot_quadTraj, makedir, genData4FluxNoise_wholeTraj
from qmtools.qpar import pars
from qmtools.qsim import QuadSim
stime = time.time()#start timer
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
# Parameters
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
pars.update({'tau2tp' : 0.101325});
pars.update({'betarf' : 0.2})
#----------------------------------------------------------
pars.update({'Num'   : 2**9})#

isLindblad  = True;#False;#
stateType  = 'Thermal';#'Fock';#'Coherent';#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
# Create quadrature simulate object
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
qs  = QuadSim(stateType, isLindblad);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
# Display system parameters
qs.print_System_Parameters()
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
# Crate Filename and Dirname
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
maindir = "results_py3_20180112/"
dirname = maindir;makedir('../' + dirname + '/')
if stateType == 'Fock':
    dirname    += 'Fck' + str(pars['fock']).rjust(1,'0');
elif stateType == 'Coherent':
     dirname   += 'Ch' + str(round(pars['alpha'],2)).ljust(4,'0');
elif stateType == 'Thermal':
    dirname += 'Th' + str(int(pars['Tho']/1e-3))+'mK';
else:
    import sys
    sys.exit('ErrorFilename and ErrorInstance: Coherent, Thermal, or Fock');
dirname += "_" + str(round(pars['g']/2/np.pi*1e3,3))+"MHz";
makedir('../' + dirname + '/')
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
# Create subdirectory for quadrature measurement
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
dirname += "/s" + str(pars['s']) + "_Np" + str(pars['Np'])\
+ "_tau" + str(round(pars['tau'],2)) + "_Brf"+str(pars['betarf']);
makedir('../' + dirname + '/');
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
# Create subdirectory for quadrature measurement
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
dirname += "/Q"+str(int(pars['Q']));makedir('../' + dirname + '/');

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
if isLindblad is True:
    filename = 'qtr_Lind'
else: 
    filename = 'qtr_vN'
filename += "_Q"+str(int(pars['Q']))\
+ "_sqtAf" + str(pars['sqrtAf']);

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
if isLindblad is True:
    filename = 'qtr_Lind'    
    # Directory for flux noise data 
    pars['dirname'] = maindir+"noise_dat_" + "s" + str(pars['s'])\
    + "Np" + str(pars['Np']) + "_tau" + str(round(pars['tau'],2)) + 'ns';
    makedir('../' + pars['dirname'] + '/');
    import glob
    #Check if the flux noise data are generated before
    # and count number of files containing the noisy data only in the folder
    num_files4noise  = len(glob.glob("../"+ pars['dirname'] + '/cpmg*'))
    #Generate the data if the number of file does not exist
    if num_files4noise != ((pars['Np']+1)*pars['s'] * 2):
        genData4FluxNoise_wholeTraj(qs.qm.time, pars, pars['dirname']);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
filename += "_Num"+str(pars['Num']);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
# Generate data for quantum trajectory
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
measType = "Actual";#"Ideal";#"Actual-II";#
if measType is "Actual":
    print ("+----------------------------------------+")
    print ("Quantum Trajectory with Actual Measurement (Algorithm 3)")
    print ("+----------------------------------------+")
    filename += "_a3";
    #Parameters characterising projection errors (gmm, dlt) and measurement errors (alf,bta)
    gama = 0.01;delta = 0.01;alfa = 0.01;beta = 0.01;
    logFile = dirname+'/'+filename;
    eI, vI  = qs.qm.qmprotocol_Algorithm3(gama, delta, alfa, beta, logFile);#measure_trajectory_of_quadratureI
elif measType is "Ideal":
    print ("+----------------------------------------+")
    print ("Quantum Trajectory with Ideal Measurement")
    print ("+----------------------------------------+")
    filename += "_i1";
    logFile = dirname+'/'+filename;
    eI, vI = qs.qm.qmprotocol_Ideal(logFile);
else:
    raise ValueError();
print (dirname, filename);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
# Plot quadrature trajectories and variance
plot_quadTraj(np.arange(len(vI)), eI, vI, measType, dirname, filename)
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

####---------------------------------------------------------------------------#
print ('Total time elapsed: %.2lf min(s) '%((time.time() - stime) / 60.0))
###---------------------------------------------------------------------------#

