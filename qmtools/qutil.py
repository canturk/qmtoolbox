#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mehmet
"""
###############################################################################
# ToolBox Quadrature Measurement Protocol for a Harmonic Oscillator (HO)      #
# Coupled to a Flux Qubit                                                     #
###############################################################################
# @author : Mehmet Canturk 
# Created on Fri Apr 14 13:38:02 2017
###############################################################################
#                            Utility Functions                                #
###############################################################################
import os
import matplotlib.pyplot as plt
#import sys # for manual exception handling
import numpy as np

#-----------------------------------------------------------------------------#
# Setup function coefficients of total Hamiltonian  in (pi/2)x, (pi)x, (pi/2)y
# and free evolution times
# See http://qutip.org/docs/4.1/guide/dynamics/dynamics-time.html
#

# In[]:
#
#-----------------------------------------------------------------------------#
def fcx(t, args):
    return np.cos(args['2omegaGE'] * t + args['varphi']);
#
#-----------------------------------------------------------------------------#
def fcy(t, args):
    return np.sin(args['2omegaGE'] * t + args['varphi']);
#
#-----------------------------------------------------------------------------#
def fcz(t, args):
    return np.cos(args['omegaGE'] * t + args['varphi']);
#
#-----------------------------------------------------------------------------#
# Setup the shaped pulse functions for total Hamiltonian in (pi)x, (pi/2)x and
# (pi/2)y pulses
def fcAx(t, args):
    return fA(t, args) * fcx(t, args);
#
#-----------------------------------------------------------------------------#
#
def fcAy(t, args):
    return fA(t, args) * fcy(t, args);
#
#-----------------------------------------------------------------------------#
#
def fcAz(t, args):
    return fA(t, args) * fcz(t, args);

# In[]:
#
#-----------------------------------------------------------------------------#
#
# Setup the fA(t) for shaped pulses
def fA(t, args):
    if((args['tb'] <= t) and (t <= args['tb'] + args['tr'])):
        return 0.5 - 0.5 * np.cos(args['pi2tr'] * (t - args['tb']));
    elif(args['tb'] + args['tr'] < t and t < args['tb'] + args['tr'] + args['tp']):
        return 1.0;
    elif((args['tb'] + args['tr'] + args['tp'] <= t) and (t <= args['te'])):
        return 0.5 + 0.5 * np.cos(args['pi2tf'] * (t - args['tb'] - args['tr'] - args['tp']));
    else:return 0.0;#TO AVOID OUT OF RANGE
#
# In[] -----------------------------------------------------------------------#
# Setup the functions for total Hamiltonian in free evoltuion time 
# and (pi/2)x, (pi)x, and (pi/2)y pulses
#-----------------------------------------------------------------------------#
def fz(t, args):
    return np.exp(+1j * args['omegaR'] * t);
#
#-----------------------------------------------------------------------------#
def fz_conj(t, args):
    return np.exp(-1j * args['omegaR'] * t);
#
#-----------------------------------------------------------------------------#
def fy(t, args):
    return np.exp(+1j * args['omegaR'] * t) * np.sin(args['omegaGE'] * t);
#
#-----------------------------------------------------------------------------#
def fy_conj(t, args):
    return np.exp(-1j * args['omegaR'] * t) * np.sin(args['omegaGE'] * t);
#
#-----------------------------------------------------------------------------#
def fx(t, args):
    return np.exp(+1j * args['omegaR'] * t) * np.cos(args['omegaGE'] * t);
#
#-----------------------------------------------------------------------------#
def fx_conj(t, args):
    return np.exp(-1j * args['omegaR'] * t) * np.cos(args['omegaGE'] * t);
#
# In[] -----------------------------------------------------------------------#
#
#Setup the functions for flux noise
def fEz(t, args):
    """
    Wmin is already added here.
    Make sure the coefficient of this function does not include Wmin."""

    # Parameters
    Nf   = args['Nf'];
    Wmin = args['Wmin'];
    mu   = 0.0;
    var  = 0.5;
    dw   = np.complex128(0.0);
    dwsum= 0.0;
    # summing all terms
    n = 1;
    while (n<=Nf):
        #an = an' + i*an"
        an  = (np.random.normal(mu,var) + 1j*np.random.normal(mu,var))\
        * np.exp(1j * Wmin * n * t);
        dw += (an  + np.conj(an)) / np.sqrt(n);
        dwsum += dw.real;#print (n,dw.real, an)
        n  += 1;
    dwavg = dwsum / Nf;#Based
    return dw.real - dwavg;
#
#-----------------------------------------------------------------------------#
#
def fEy(t,args):
    return fEz(t, args) * np.sin(args['omegaGE']*t);
#
#-----------------------------------------------------------------------------#
#
def fEx(t,args):
    return fEz(t, args) * np.cos(args['omegaGE']*t);
#
#-----------------------------------------------------------------------------#
#
def plot_quadTraj(slist, iexpI, ivarI, measType, dirname, filename=[]):
    #-----------------------------------------------------------------#
    # Save data
    #np.savetxt(self.path4data + filename + '.out', np.column_stack((iexpI, ivarI)), delimiter='\t');
    #-----------------------------------------------------------------#
    # Plot and export eps figures for expectation value and variance.
    #s = param['s'];
    Nr = len(slist)-1;#slist = np.arange(1 + s);
    # Update the matplotlib configuration parameters:
    plt.rcParams.update({'font.size': 18, 'font.family': 'serif'});
    plt.ioff() #to prevent the  multiple pop up figures for each repetition 
    # <I> vs number of repetition
    fig, axes = plt.subplots(2, 1, figsize=(8,6))
    if measType == "Actual":
        labelType = "Actual";
    elif measType == "Ideal":
        labelType = "Ideal";
    color_list = [];
    color_list.append('r-');
    plt.xlim([0, Nr])
    axes[0].plot(slist, ivarI,  'r-', label=labelType,  linewidth=1.25);
    axes[1].plot(slist, iexpI,  'r-', linewidth=1.25);
    axes[0].grid(True);axes[1].grid(True);
    # Now add the legend with some customizations.
    legend = axes[0].legend(loc='upper right', shadow=True)
    # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    # Set the fontsize
    for label in legend.get_texts():
        label.set_fontsize('smaller')
    # Set the line width
    for label in legend.get_lines():
        label.set_linewidth(2.0)  # the legend line width
    
    #axes[0].set_xticks([]);#
    axes[0].set_xlim([0, Nr]);
    axes[1].set_xlim([0, Nr]);
    # Set the vertical and horizontal axes
    axes[1].set_xlabel(r'$s$', fontsize=20)
    axes[1].set_ylabel(r'$\left\langle I\right\rangle$', fontsize=18)
    axes[0].set_ylabel(r'$(\Delta I^2)$', fontsize=18)
    #---------------------------------
    #Create directory for results
    dirname  = '../' + dirname + '/';
    makedir(dirname)
    # save figure file
    save_fig(dirname+filename, ext='png', close=True, verbose=True);
    # save data file 
    np.savetxt(dirname+filename + '.out', np.column_stack((slist, iexpI, ivarI)))
    #plt.show();
    #plt.close();
    return;#end of module
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#
def list_LeafNodes(root,lst):
    """
    Q. How to print leaf nodes in the binary tree left to right?
    A. Perform a depth-first traversal of the tree, 
    handling the left sub-trees first and printing only the leaf nodes.
    The easiest way to implement this is with a recursive function.
    """
    if(root == None):
        return;
    list_LeafNodes(root.left, lst);
    # Print the leaf node
    if root.left == None and root.right == None:
        lst.append(root.value);#print(root.value);
    list_LeafNodes(root.right, lst);
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#
## BinaryTree: https://pypi.python.org/pypi/binarytree/1.1.1
#from binarytree import Node, pprint, inspect
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#
def print_LeafNodes(root):
    """
    Q. How to print leaf nodes in the binary tree left to right?
    A. Perform a depth-first traversal of the tree, 
    handling the left sub-trees first and printing only the leaf nodes.
    The easiest way to implement this is with a recursive function.
    """
    if(root == None):
        return;
    print_LeafNodes(root.left);
    # Print the leaf node
    if root.left == None and root.right == None:
        print(round(root.value,6), end=', ');
    print_LeafNodes(root.right);
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#
def sort_Lists(jntlst, explst, varlst, sortBy="Pjnt"):
    """Sorting three synchronized lists"""
    jntlst_arr = np.array(jntlst);
    explst_arr = np.array(explst);
    varlst_arr = np.array(varlst);
    if sortBy == "Pjnt":
        idx = np.argsort(jntlst);
    elif sortBy == "varI":
        idx = np.argsort(varlst);
    elif sortBy == "expI":
        idx = np.argsort(explst);
    else:
        raise ValueError("Invalid choice");
    return (np.array(jntlst_arr)[idx], np.array(explst_arr)[idx], np.array(varlst_arr)[idx])
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#
def plot_sortedlist(jntlst, explst, varlst, sortBy="Pjnt"):#default: "Pjnt"
    pj_sorted, eI_sorted, vI_sorted = sort_Lists(jntlst, explst, varlst, sortBy);
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    # Update the matplotlib configuration parameters:
    plt.rcParams.update({'font.size': 16, 'font.family': 'serif'});
    plt.ioff() #to prevent the  multiple pop up figures for each repetition 
    # <I> vs number of repetition
    fig, axes = plt.subplots(3, 1, figsize=(10,12))
    # Joint probability
    axes[0].plot(pj_sorted);
    axes[0].set_yscale('log');
    axes[0].grid(True);
    axes[0].set_ylabel("Joint Probabilities");
    axes[0].set_xlabel("State")
    axes[0].set_title("Sorted for " + sortBy)
    # Variance
    axes[1].plot(vI_sorted);
    axes[1].set_yscale('log');
    axes[1].grid(True);
    axes[1].set_ylabel("Variance");
    axes[1].set_xlabel("State")
    # Expectation value
    axes[2].plot(eI_sorted)
    #axes[2].set_yscale('log');
    axes[2].grid(True);
    axes[2].set_ylabel(r'$\langle I\rangle$')
    axes[2].set_xlabel('State')
    # To make sure the horizontal axis ticks are integers
    axes[0].xaxis.set_major_locator(MaxNLocator(integer=True))
    axes[1].xaxis.set_major_locator(MaxNLocator(integer=True))
    axes[2].xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.show()
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#
def makedir(filename):
    """   Creates a folder  """
    dir = os.path.dirname(filename)
    try:
        os.stat(dir)
    except:
        pass
        os.mkdir(dir);
    print ("'" + dir + "'/' is created.")
#
#-----------------------------------------------------------------------------#
#
def save_fig(path, ext='eps', close=True, verbose=True):
    """Save a figure from pyplot.
    Parameters
    ----------
    path : string
        The path (and filename, without the extension) to save the
        figure to.
    ext : string (default='png')
        The file extension. This must be supported by the active
        matplotlib backend (see matplotlib.backends module).  Most
        backends support 'png', 'pdf', 'ps', 'eps', and 'svg'.
    close : boolean (default=True)
        Whether to close the figure after saving.  If you want to save
        the figure multiple times (e.g., to multiple formats), you
        should NOT close it in between saves or you will have to
        re-plot it.
    verbose : boolean (default=True)
        Whether to print information about when and where the image
        has been saved.
    usage: save_fig("filename", ext="eps", close=True, verbose=True)
    """
    # Extract the directory and filename from the given path
    directory = os.path.split(path)[0]
    filename = "%s.%s" % (os.path.split(path)[1], ext)
    if directory == '':
        directory = '.';
    # If the directory does not exist, create it
    if not os.path.exists(directory):os.makedirs(directory);
    # The final path to save to
    savepath = os.path.join(directory, filename)
    if verbose:
        print("Saving figure to '%s'..." % savepath);
    #---------------------------------------------------------#
    import matplotlib.pyplot as plt #%matplotlib notebook
    # Actually save the figure
    plt.savefig(savepath)        
    # Close it
    if close:plt.close()
    #if verbose:print("Done")
#
#-----------------------------------------------------------------------------#
#
def josephson_energy(phi1, phi2, f, alpha):
    """
    The total josephson energy of a three-junction flux qubit [Orlando et al 1999]
    Parameters:
        phi1 : phase of junction 1
        phi2 : phase of junction 2
        f    : frustration (f=Phi/Phi0 where Phi external magnetic field)
        alpha: the ratio of third juction
    """
    return (2 + alpha - np.cos(phi1) - np.cos(phi2) - alpha * np.cos(2*np.pi*f+phi1-phi2));
# In[ ]:

#Plot expectztion values wrt times
def plot_Bloch_Oscillations(results, titles, filename=[]):
    #Update the matplotlib configuration parameters:
    plt.rcParams.update({'font.size': 14, 'font.family': 'serif'})
    #---------------------------------------------------------------------#
    plt.ioff() #to prevent the  multiple pop up figures for each repetition
    if len(titles) == 1:
        fig, axes = plt.subplots(1, 1, figsize=(8,6), sharey=True)
        axes.plot(results[-1].times, np.real(results[-1].expect[2]), 'r-', label=r'$\left\langle\sigma_z\right\rangle$');
        axes.plot(results[-1].times, np.real(results[-1].expect[1]), 'b-', label=r'$\left\langle\sigma_y\right\rangle$');
        axes.plot(results[-1].times, np.real(results[-1].expect[0]), 'g-', label=r'$\left\langle\sigma_x\right\rangle$');
        axes.set_xlabel(r'$\tau$ [ns]', fontsize=16)
        axes.set_xlim([results[-1].times[0], results[-1].times[-1]]);
        axes.set_title(titles[-1])
        axes.grid(True);
        axes.set_ylabel('Bloch Oscillations', fontsize=16)
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
        plt.show()
        return;
    #-----------------------------------
    fig, axes = plt.subplots(1, len(titles), figsize=(8,6), sharey=True)
    for i in range(len(titles)):
        axes[i].plot(results[i].times, np.real(results[i].expect[2]), 'r-', label=r'$\left\langle\sigma_z\right\rangle$');
        axes[i].plot(results[i].times, np.real(results[i].expect[1]), 'b-', label=r'$\left\langle\sigma_y\right\rangle$');
        axes[i].plot(results[i].times, np.real(results[i].expect[0]), 'g-', label=r'$\left\langle\sigma_x\right\rangle$');
        axes[i].set_xlabel(r'$\tau$ [ns]', fontsize=16)
        axes[i].set_xlim([results[i].times[0], results[i].times[-1]]);
        axes[i].set_title(titles[i])
        axes[i].grid(True);
    axes[0].set_ylabel('Bloch Oscillations', fontsize=16)
    # Now add the legend with some customizations.
    legend = axes[0].legend(loc='lower left', shadow=False)
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
    if (filename == []):
        plt.show()
    else:
        save_fig('../'+filename+'_xyz', ext='png', close=True, verbose=True)
    return;
#
# In[] - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#
def plot_qubit_cpmg_BlochOscilations(tlst, exlst, eylst, ezlst, title_list, filename):
    """REPLACED BY plot_Bloch_oscillations(...)"""
    plt.close('all');
    
    color_list = [];
    color_list.append('r-');
    color_list.append('b-');
    color_list.append('g-');
    
    label_list = [];
    label_list.append(r'$\left\langle\sigma_z\right\rangle$')
    label_list.append(r'$\left\langle\sigma_y\right\rangle$')
    label_list.append(r'$\left\langle\sigma_x\right\rangle$')
            
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
    return;
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#
def saveData4FluxNoise_singleSlot(t, noisex, noisey, noisez, fname):    
    import matplotlib.pyplot as plt    
    plt.plot(t,noisex, 'r-', label=r'$f_{\varepsilon_x}$')
    plt.plot(t,noisey, 'g-', label=r'$f_{\varepsilon_y}$')
    plt.plot(t,noisez, 'b-', label=r'$f_{\varepsilon_z}$')
    legend = plt.legend(loc='upper right', shadow=True)
    # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    # Set the fontsize
    for label in legend.get_texts():
        label.set_fontsize('smaller')
    # Set the line width
    for label in legend.get_lines():
        label.set_linewidth(2.0)  # the legend line width
    plt.grid(True)
    
    if fname==[]:
        plt.show()
    else:
        # save figure file
        save_fig(fname, ext='png', close=True, verbose=True);
        # save data file 
        np.savetxt(fname + '.out', np.row_stack((t, noisex, noisey, noisez)), delimiter='\t')
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#
def genData4FluxNoise_singleSlot(ti, tf, **kwargs):
    """
    Usage:
        genData4FluxNoise_singleSlot(ti,tf,Num=pars['Num'],Nf=pars['Nf'],Wmin=pars['Wmin'],omegaGE=pars['omegaGE'])
    """
    t = np.linspace(ti, tf, kwargs.get('Num'))
    args={'Nf':kwargs.get('Nf'),'Wmin':kwargs.get('Wmin'), 'omegaGE':kwargs.get('omegaGE')}
    #Generate noise functions
    fnoisy_x = lambda t: fEx(t, args)
    fnoisy_y = lambda t: fEy(t, args)
    fnoisy_z = lambda t: fEz(t, args)
    #Generate the data
    noisy_data_x = fnoisy_x(t)
    noisy_data_y = fnoisy_y(t)
    noisy_data_z = fnoisy_z(t)
    return (t, noisy_data_x, noisy_data_y, noisy_data_z)
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#
def genData4FluxNoise_wholeTraj(tlist, params, dirname):
    """
    Generates data for a single trajectory
    
       (pi/2)x       +(pi)x                -(pi)x       (pi/2)y
        +---+      +--------+            +--------+      +---+++++++++++
        |   |      |        |            |        |      |   |   <I>   +
        +---+------+--------+------------+--------+------+---+---------+
m=0     t0  t1     t2       t3           t4       t5     t6  t7         
m=1     t8  t9     t10     t11          t12      t13    t14  t15
m=:
  :
    """
    def filename(m,n):
        return "../" + dirname + "/cpmg%02dn%02d"%(m + 1, n);
    #--------------------------------------------------
    # Generate noise data for subsequent measurement time interval (synchronized)
    #--------------------------------------------------
    t = tlist
    for m in range(params['s']):
        mp   = 2 * m * (params['Np'] + 2);# Index of m
        _Num  = params['Num'];# Number of subintervals
        # Evolution in free precession time (half) between t[mp] and t[mp+1] 
        tslot, datx, daty, datz =\
        genData4FluxNoise_singleSlot(t[mp+1], t[mp+2],
                                     Num    = _Num,
                                     Nf     = params['Nf'],
                                     Wmin   = params['Wmin'],
                                     omegaGE= params['omegaGE']
                                     )
        #Plot the data
        saveData4FluxNoise_singleSlot(tslot, datx, daty, datz, filename(m,0))
        _Num = int(2*params['Num']);
        n   = 1;
        while n<=params['Np']:
            # Evolution in free precessing time between t[mp+2n] and t[mp+2n+1]
            if n == params['Np']:
                _Num = int(params['Num']);
            tslot, datx, daty, datz =\
            genData4FluxNoise_singleSlot(t[mp+2*n+1], t[mp+2*n+2],
                                         Num    =_Num,
                                         Nf     = params['Nf'],
                                         Wmin   = params['Wmin'],
                                         omegaGE= params['omegaGE']
                                         );
            saveData4FluxNoise_singleSlot(tslot, datx, daty, datz, filename(m,n))
            # Update n for next evolution
            n += 1;
    return;
# end: genData4FluxNoise_wholeTraj(..)
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#