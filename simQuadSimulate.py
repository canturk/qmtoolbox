###############################################################################
# QND Measurement of a Nanomechanical Resonator Coupled to a Flux Qubit
###############################################################################
# Author : Mehmet Canturk 
# Date   : May 16, 2015
# Update : 2015/05/16-12:30 AM
# Office : QND 2211,
# Address: Institute for Quantum Computing, University of Waterloo,
# 200 University Ave, West Waterloo, Ontario, Canada
# N2L 3G1         
###############################################################################
from   simFockState  import *
from   simCoherentState import *

#Other packages
import os
import sys

# %matplotlib inline 
import matplotlib.pyplot as plt
from   scipy.stats     import norm
import matplotlib.mlab as     mlab

###############################################################################
class QuadSimulation(object):
    """
    Class QuadSimulation consists of objects instantiated from class QuadFockState 
    and from class QuadCoherent. It simulates the quadrature measurement for
    various lambda values and it is capable of creating esp figures, saving data, 
    generating presentation file for LaTeX presentation and displaying them. 
    """
    def __init__(self, sim_args):
        """
        Construct Class QuadSimulation
        """
        self.args = sim_args # transfer all input args to this one 
        self.status     = sim_args['status'] #Status of simulation which can be either of {Vacuum, FockState, Coherent}
        self.Nrep       = sim_args['Nrep'] # Number of repetition of quadrature measurement
        self.Npulse     = sim_args['Npulse']
        
        self.qsim    = []
        if sim_args['status'] == 'Coherent':
            self.qsim = QuadCoherentState(sim_args)
        else:
            self.qsim = QuadFockState(sim_args)
        
        self.lmb0_list  = self.generate_lambda0List()
        
        self.teXfileObj = [] # Null File Stream object
        self.path4figs  = [] # directory path for figures
        self.path4data  = [] # directory path for data
        
        # Generate directory name for repository
#           'filename4TeX': 'qmResultTest', 
#           'dirname4Result': 'result4QM_',
        teXdir  = '../' + sim_args['dirname4Result'] + '/'  #'results4QM_' + sim_args['dir_label'] + sim_args['lego'] + '/'
        self.makedir(teXdir)
        # Generate file name for presentation LaTeX file
        teXfile = sim_args['filename4TeX'] + sim_args['status']+'_lmd' + str(sim_args['lambdai']) + 'To' + str(sim_args['lambdaf']) + '_Np' + str(sim_args['Npulse']) + '_Nr' + str(sim_args['Nrep']) + '.tex';
        teXfile = teXdir + teXfile
        # Generate directory name for figure repository
        self.path4figs  = teXdir  + 'figures/'
        self.makedir(self.path4figs)
        # Generate directory name for data repository
        self.path4data  = teXdir  + 'data/'
        self.makedir(self.path4data)
        # Create a LaTeX document to present the simulation results.
        self.open_LaTeX_Doc(teXfile, teXdir, sim_args)
    
    ###########################################################################
    def __del__(self):
        """
        Destructor of QuadSimulation to print class name of an instance 
        that is about to be destroyed
        """
        self.closeTeXfile(); print 'Closing the LaTex Presentation File...'
        del self.qsim
        class_name = self.__class__.__name__; print class_name, "destroyed ..."
    ###########################################################################
    # Generate Lambda list
    def generate_lambda0List(self):
        """
        Generate list of lambda0 values based on given parameters
        """
        lambdai = self.qsim.arg4QM['lambdai']
        lambdaf = self.qsim.arg4QM['lambdaf']
        dlambda = self.qsim.arg4QM['dLambda']
        #('dlambda', 'lambdai', 'lambdaf') Increment, 
        #Nlambda    = 1 + int((lambdaf - lambdai)/dlambda)
        #lambda0Lst = np.linspace(lambdai, lambdaf, Nlambda)
        lambda0Lst = []
        lambda0 = lambdai
        while(lambda0 <= lambdaf):
            lambda0Lst.append(lambda0)
            lambda0 = lambda0 + dlambda
        print "Lambda List: "
        print ", ".join("%.2f"%lambda0 for lambda0 in lambda0Lst)
        return lambda0Lst
        ###########################################################################
    def generate_eps_filenames(self, lambda0):
        """
        Generate eps figures file names for expectation value and variance
        """
        epsfname4expI  = 'expI_lmd' + str(lambda0).ljust(5,'0') + 'Np' + str(self.qsim.arg4QM['Npulse']) + '.eps';
        epsfname4varI  = 'varI_lmd' + str(lambda0).ljust(5,'0') + 'Np' + str(self.qsim.arg4QM['Npulse']) + '.eps';                
        return (epsfname4expI, epsfname4varI)
    ###########################################################################
    def simulate_expect_and_variance_of_quadI(self):
        """
        Simulates the system for a list of lambda0 values to make 
        quadrature measurement and resultant  <I> and DeltaI values are added 
        into their lists respectively.
        """
        for k in xrange(len(self.lmb0_list)):
            # Make the quadrature measurement and the resultant  <I> and DeltaI            
            #if sim_args['lego'] == 'Analytical': #Hbar via Analytical
            eI_list_Havg, vI_list_Havg = self.qsim.measure_expect_and_variance_of_quadI('Havg', self.lmb0_list[k])
            eI_list_Hbar, vI_list_Hbar = self.qsim.measure_expect_and_variance_of_quadI('Hbar', self.lmb0_list[k])
            # Create an empty double list 
            expI_list = []
            expI_list.append(eI_list_Havg)
            expI_list.append(eI_list_Hbar)
            varI_list = []
            varI_list.append(vI_list_Havg)
            varI_list.append(vI_list_Hbar)
                        
            # Generate eps figures file names for expectation value and variance            
            figName4I   = 'expI' + self.qsim.arg4QM['status'] + '_lmd' + str(self.lmb0_list[k]).ljust(5,'0') \
                            + 'Np' + str(self.qsim.arg4QM['Npulse']) + '.eps';            
            figName4Var = 'varI' + self.qsim.arg4QM['status'] + '_lmd' + str(self.lmb0_list[k]).ljust(5,'0') \
                            + 'Np' + str(self.qsim.arg4QM['Npulse']) + '.eps';                        
            
            self.save_expectVariance_of_epsFigs(expI_list, varI_list,\
                                            self.lmb0_list[k], figName4I, figName4Var)            
            self.insert_TeXslide4figs(self.lmb0_list[k], figName4I, figName4Var)
            print '\n' #To make an empty line before the next run starts
        
    ###########################################################################
    def save_expectVariance_of_epsFigs(self, expectI_list, varianceI_list, lambda0, figName4I, figName4Var):
        """
        Creates and export eps figures for expectation value and variance.
        """
        # Update the matplotlib configuration parameters:
        plt.rcParams.update({'font.size': 18, 'font.family': 'serif'})

        plt.ioff() #to prevent the  multiple pop up figures for each repetition 
        # <I> vs number of repetition
        fig, axes = plt.subplots(1, 1, figsize=(8,5))
        Nrep_list = np.arange(self.qsim.arg4QM['Nrep'])
        #axes.plot(nlist, quadICmbList[0], label='$H_{o}$')
        axes.plot(Nrep_list, expectI_list[0], "-b", label = r'$H_{avg}$')
        axes.plot(Nrep_list, expectI_list[1], "--r", label = r'$H_{bar}$')
        axes.grid(True)
        # Now add the legend with some customizations.
        legend = axes.legend(loc='upper center', shadow=True)
        # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
        frame = legend.get_frame()
        frame.set_facecolor('0.90')
        # Set the fontsize
        for label in legend.get_texts():
            label.set_fontsize('large')
        # Set the line width
        for label in legend.get_lines():
            label.set_linewidth(2.0)  # the legend line width
        # Set the vertical and horizontal axes
        axes.set_xlabel('Repetition',        fontsize=18)
        axes.set_ylabel(r'$\left<I\right>$', fontsize=18)
        axes.set_title(r'$\lambda_0$ = ' + str(lambda0) + r'$N_p$ = ' + str(self.qsim.arg4QM['Npulse']))
        fig.savefig(self.path4figs + figName4I)
        plt.clf()
        print '"' + self.path4figs + figName4I + ' is created. . .'
        ######################################################################
        plt.ioff()#to prevent the  multiple pop up figures for each repetition
        # Plot Variance deltaI wrt repetition
        fig, axes = plt.subplots(1, 1, figsize=(8,5))
        #axes.plot(nlist, deltaQuadICmbList[0], label='$H_{o}$')
        axes.plot(Nrep_list, varianceI_list[0], "-b", label = r'$H_{avg}$')
        axes.plot(Nrep_list, varianceI_list[1], "--r", label = r'$H_{bar}$')
        axes.grid(True)
        # Now add the legend with some customizations.
        legend = axes.legend(loc='upper center', shadow=True)
        # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
        frame = legend.get_frame()
        frame.set_facecolor('0.90')
        # Set the fontsize
        for label in legend.get_texts():
            label.set_fontsize('large')
        # Set the line width
        for label in legend.get_lines():
            label.set_linewidth(2.0)  # the legend line width
        # Set the vertical and horizontal axes
        axes.set_xlabel('Repetition',  fontsize=18)
        axes.set_ylabel(r'$\Delta I$', fontsize=18)
        axes.set_title(r'$\lambda_0$ = ' + str(lambda0) + ', ' + r'$N_p$ = ' +\
                        str(self.qsim.arg4QM['Npulse']))
        fig.savefig(self.path4figs + figName4Var)
        plt.clf()
        print '"' + self.path4figs + figName4Var + ' is created. . .'
    
    ###########################################################################
    def simulate_eigenenergies_of_Hbar_and_Havg(self):
        """
        Simulates the system for a list of lambda0 values to make 
        """
        for k in xrange(len(self.lmb0_list)):
            Hbar = self.qsim.hamiltaonian4Hbar(self.lmb0_list[k])
            Havg = self.qsim.hamiltaonian4Havg(self.lmb0_list[k])
            Ebar = Hbar.eigenenergies()
            Eavg = Havg.eigenenergies()
            psibar=Hbar.eigenstates()
            psiavg=Havg.eigenstates()
            # Generate file name eps figures of eigenenergies
            # Eigenenergy depends only on lambda0 
            figName4eig = 'eigEnergy' + '_lambda' + str(self.lmb0_list[k]).ljust(5,'0') + '.eps';
            self.save_eigenenergy_of_Hamiltonian_epsFigs(Ebar, Eavg, self.lmb0_list[k], figName4eig)
            Ubar = self.qsim.unitaryEvolution4HbarMagnus(self.lmb0_list[k])
            Uavg = self.qsim.unitaryEvolution4Havg(self.lmb0_list[k])
            Uavg.dag() * Uavg
    
    #########################################################################
    def save_eigenenergy_of_Hamiltonian_epsFigs(self, eig_bar, eig_avg, lambda0, figName4Eig):
        """
        Creates and export eps figures for expectation value and variance.
        """
        # Update the matplotlib configuration parameters:
        plt.rcParams.update({'font.size': 18, 'font.family': 'serif'})

        plt.ioff() #to prevent the  multiple pop up figures for each repetition 
        # <I> vs number of repetition
        fig, axes = plt.subplots(1, 1, figsize=(8,5))
        N_list = np.arange(len(eig_bar))
        #axes.plot(nlist, quadICmbList[0], label='$H_{o}$')
        axes.plot(N_list, eig_bar, "-b", label = r'$H_{bar}$')
        axes.plot(N_list, eig_avg, "--r", label = r'$H_{avg}$')
        axes.grid(True)
        # Now add the legend with some customizations.
        legend = axes.legend(loc='upper center', shadow=True)
        # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
        frame = legend.get_frame()
        frame.set_facecolor('0.90')
        # Set the fontsize
        for label in legend.get_texts():
            label.set_fontsize('large')
        # Set the line width
        for label in legend.get_lines():
            label.set_linewidth(2.0)  # the legend line width
        # Set the vertical and horizontal axes
        axes.set_xlabel(r'$2\otimes N$', fontsize=18)
        axes.set_ylabel('Eigenenergies', fontsize=18)
        axes.set_title(r'$\lambda_0$ = ' + str(lambda0))
        fig.savefig(self.path4figs + figName4Eig)
        plt.clf()
        print '"' + self.path4figs + figName4Eig + ' is created. . .'
    #########################################################################
    
    def simulate_probability_distributions(self, Nsample, hamiltonian):
        """
        Simulates the probability distribution of quadrature measurement. 
        The measurement process repeats Nsample times.
        """
        # Make a collection of Nsample measurement 
        for k in xrange(len(self.lmb0_list)):
            #Make the measurement for probability distribution
            pExp_list, pVar_list = self.qsim.measure_probDist_of_expectI_and_varianceI(hamiltonian, self.lmb0_list[k], Nsample)
            #Generate figure filenames
            epsfname4pExp = 'probExp4' + str(self.qsim.arg4QM['status']) + hamiltonian + '_lmb' + str(self.lmb0_list[k]).ljust(4,'0') +\
            'Np' + str(self.qsim.arg4QM['Npulse']) + 'Ns' + str(Nsample) + '.eps';
            epsfname4pVar = 'probVar4' + str(self.qsim.arg4QM['status']) + hamiltonian + '_lmb' + str(self.lmb0_list[k]).ljust(4,'0') + 'Np' +\
            str(self.qsim.arg4QM['Npulse']) + 'Ns' + str(Nsample) + '.eps';
            #create eps figures and save them
            self.save_probDist_epsFigs(hamiltonian, pExp_list, pVar_list, self.lmb0_list[k], epsfname4pExp, epsfname4pVar)
            # self.probDist_save_epsfig_varianceI(hamiltonian, pVar_list, self.lmb0_list[k], epsfname4pVar)
            # self.probDist_save_epsfig_expectI(hamiltonian, pExp_list, self.lmb0_list[k], epsfname4pExp)
            self.insert_TeXslide4figs(self.lmb0_list[k],  epsfname4pExp, epsfname4pVar)
            print '\n' #Insert a newline between previous runs and a new one (separate them)
        #for(k)
    # end of simulate_for_prob_distribution_and_savefigs()
    
    ###########################################################################
    def save_probDist_epsFigs(self, hamiltonian, elist, vlist, lambda0, figname4Exp, figname4Var):
        """
        Creates histograms in terms of eps format and same them.
        Important: The plot for DeltaI illustrates a weired values in vertical dimension. 
        """        
        # Update the matplotlib configuration parameters:
        plt.rcParams.update({'font.size': 18, 'font.family': 'serif'})
    
        plt.ioff() #to prevent the  multiple pop up figures for each repetetion 
        # Probability of <I> vs <I>
        #
        # insert figure plots here
        plt.rcParams.update({'font.size' : 18, 'font.family' : 'serif'})
        # the histogram of the data with histtype='step'
        (muExp, sigmaExp)  = norm.fit(elist)
        numBins          = 50
        nExp, binsExp, ptchsExp = plt.hist(elist, numBins, normed=True, histtype='stepfilled')
        plt.setp(ptchsExp, 'facecolor', 'g', 'alpha', 0.125)
        # add a line showing the expected distribution
        yExp = mlab.normpdf( binsExp, muExp, sigmaExp)
        lExp = plt.plot(binsExp, yExp, 'k-', linewidth = 2)
        plt.xlabel(r'$\left<I\right>$', fontsize = 20)
        plt.ylabel(r'$P\left(\left<I\right>\right)$', fontsize = 20)
        plt.grid(True)
        
        FWHM = 2 * sqrt(2 * log(2)) * sigmaExp ##FWHM, the full width at half maximum 
        info_title = r'$\Delta I_i$ = '+ str(round(self.qsim.variance4I0, 2)).ljust(4,'0') + \
        ', FWHM = ' + str(round(FWHM, 2));
        if hamiltonian == 'Havg':
            plt.title( info_title + ' for $H_{avg}$', fontsize = 16)
        elif hamiltonian == 'Ho':
            plt.title(info_title + ' for $H_{o}$', fontsize = 16)
        else:
            plt.title(info_title + ' for $H_{bar}$', fontsize = 16)        
        #
        #
        plt.savefig(self.path4figs + figname4Exp)
        plt.clf()
        print '"' + self.path4figs + figname4Exp + ' is created. . .'
        ######################################################################    
        plt.ioff()#to prevent the  multiple pop up figures for each repetetion      
        # Probability of DeltaI vs DeltaI
        #
        # insert figure plots here
        plt.figure()
        (muVar, sigmaVar) = norm.fit(vlist)
        nVar, binsVar, ptchsVar = plt.hist(vlist, numBins, normed=True, histtype='stepfilled')
        
        plt.setp(ptchsVar, 'facecolor', 'g', 'alpha', 0.125)
        # add a line showing the expected distribution
        yVar = mlab.normpdf( binsVar, muVar, sigmaVar)        
        lVar = plt.plot(binsVar, yVar, 'k-', linewidth = 2)
        plt.xlabel(r'$\Delta I$', fontsize = 20)
        plt.ylabel(r'$P\left(\Delta I \right)$', fontsize = 20)
        plt.grid(True)
        
        if hamiltonian == 'Havg':
            plt.title(r'$\Delta I_{i}$ = '+ str(self.qsim.variance4I0) + ' for $H_{avg}$', fontsize = 18)
        elif hamiltonian == 'Ho':
            plt.title(r'$\Delta I_{i}$ = '+ str(self.qsim.variance4I0) + ' for $H_{o}$', fontsize = 18)
        else:
            plt.title(r'$\Delta I_{i}$ = '+ str(self.qsim.variance4I0) + ' for $H_{bar}$', fontsize = 18)
        
        
        plt.savefig(self.path4figs + figname4Var)
        plt.clf()
        print '"' + self.path4figs + figname4Var + ' is created. . .'
        
        print 'mu=%lf, sigma=%lf for <I>'  %(muExp, sigmaExp)
        print 'mu=%lf, sigma=%lf for DelI' %(muVar, sigmaVar)
    
    ########################################################################### 
    def makedir(self, filename):
        """
        Creates a folder
        """
        dir = os.path.dirname(filename)
        try:    os.stat(dir)
        except: os.mkdir(dir)
        print '"' + dir + '/" is created.'

    ################################################################################
    # Create LaTex Presentation File
    def open_LaTeX_Doc(self, teXfile, teXdir, args):
        """
        def open_LaTeX_Doc(self, args, lambdai, lambdaf, Np):
        Creates repository directories (for presentation file and eps figure directory) 
        unless it exists and then it opens LaTeX file for the presentation. 
        Besides, it generates path for figure directory.
        """
        # Open teXfile for presentation: In later release, avoid writing on an existing file.
        self.teXfileObj = open(teXfile,'w') # open from scratch
        print teXfile + ' is opened'
        #######################################################################
        self.teXfileObj.write('\documentclass{beamer}\n')
        self.teXfileObj.write('\usepackage{graphicx, subfigure}\n')
        self.teXfileObj.write('\setbeamertemplate{caption}[numbered]\n')
        self.teXfileObj.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        self.teXfileObj.write('\n')
        self.teXfileObj.write('\mode<presentation>\n')
        self.teXfileObj.write('{\n')
        self.teXfileObj.write('\usetheme{Warsaw}\n')
        self.teXfileObj.write('\setbeamercovered{transparent}\n')
        self.teXfileObj.write('}\n')
        self.teXfileObj.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        self.teXfileObj.write('\n')
        self.teXfileObj.write('\\title{'+ args['title'] + '}\n')
        self.teXfileObj.write('\\author[M. Canturk \copyright 2015]{Written by Mehmet Canturk\inst{1}\inst{2}\\\Supervised by Adrian Lupascu\inst{1}\inst{3}}\n')
        inst1 = '\inst{1} Institute for Quantum Computing, University of Waterloo, On, Canada'
        inst2 = '\inst{2} Turgut Ozal University, Ankara, Turkey'
        inst3 = '\inst{3} Department of Physics and Astronomy, University of Waterloo, On, Canada'
        self.teXfileObj.write('\institute{' + inst1 + '\n \\\ ' + inst2 + '\n \\\ ' + inst3 + '}\n')
        self.teXfileObj.write('\n')
        self.teXfileObj.write('\n')
        self.teXfileObj.write('\usepackage{times}\n')
        self.teXfileObj.write('\n')
        self.teXfileObj.write('\n')
        self.teXfileObj.write('\\begin{document}\n')
        self.teXfileObj.write('\n')
        self.teXfileObj.write('\n')
        
        #title page
        self.teXfileObj.write('\\begin{frame}\n')
        self.teXfileObj.write('\\titlepage\n')
        self.teXfileObj.write('\end{frame}\n')
        
        msg = 'for $'+str(args['lambdai'])+' \leq \lambda_0 \leq' + str(args['lambdaf'])+'$';
        self.teXfileObj.write('\\begin{frame}{Outline of the Simulation Runs ' + msg +'}\n')
        self.teXfileObj.write('\\tableofcontents\n')
        self.teXfileObj.write('\end{frame}\n')
        self.teXfileObj.write('\n')
        self.teXfileObj.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        self.teXfileObj.write('\n')
        
        if args['status'] == 'Vacuum':
            self.teXfileObj.write('\\begin{frame}[fragile]{Quadratic Measurement for {\color{red} Vacuum} State}\n')
        elif args['status'] != 'FockState':
            self.teXfileObj.write('\\begin{frame}[fragile]{Quadratic Measurement for {\color{red} Fock} State}\n')
        else:
            self.teXfileObj.write('\\begin{frame}[fragile]{Quadratic Measurement for {\color{red} Coherent} State}\n')
        self.teXfileObj.write('\\begin{itemize}\n')
        from datetime import datetime
        self.teXfileObj.write('\item Date and time: ' + str(datetime.now()) + '\n')
        # SIMULATION REPORTd #
        codePath  = os.path.dirname(__file__)     # directory of the current python script
        codeFile  = os.path.basename(sys.argv[0]) # filename of the current python script
        codeReprt = '\item Directory: \\verb$' + codePath + '/$\n' + '\item File: \\verb$' + codeFile + '$\n'
        codeReprt = codeReprt + '\item Number of $\pi$-pulse sequence: $N_p = ' + str(args['Npulse']) + '$\n'
        codeReprt = codeReprt + '\item Number of repetitions: $N_{rep} = ' + str(args['Nrep']) + '$\n'
        codeReprt = codeReprt + '\item Report Dir: \\verb$' + teXdir + '$\n';
        self.teXfileObj.write( codeReprt)
        self.teXfileObj.write('\end{itemize}\n')
        #sim4Hamiltonian = args['enum']['Havg'] + ' and ' + args['enum']['Hbar']
        #self.teXfileObj.write('\\color{red} Simulation is performed for ' + sim4Hamiltonian + '\n')
        self.teXfileObj.write('\n')
        self.teXfileObj.write('\end{frame}\n')
        self.teXfileObj.write('\n')
        self.teXfileObj.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');

    ###########################################################################
    def closeTeXfile(self):
        """
        Close the LaTex Presentation File at the end of simulation

        """
        self.teXfileObj.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');    
        self.teXfileObj.write('\n')
        self.teXfileObj.write('\n')
        self.teXfileObj.write('\end{document}\n')
        self.teXfileObj.close(); # At the end close the file
        
    ###########################################################################
    def insert_TeXslide4figs(self, lambda0, figIName, figVarName):
        """
        write a presentation slide/frame in LaTeX document
        """
        if (self.Npulse % 2) == 0:
            whatfor = 'Run for {\color{red}EVEN} $N_p$ at $\lambda_0=' + str(lambda0) + '$'
        else:
            whatfor = 'Run for {\color{red}ODD} $N_p$ at $\lambda_0=' + str(lambda0) + '$'        
        
        # Generate LaTeX frame for a figure
        self.teXfileObj.write('\n');
        self.teXfileObj.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        self.teXfileObj.write('\n')
        self.teXfileObj.write('\\begin{frame}{' + whatfor +'}\n')
        self.teXfileObj.write('\\begin{figure}[!h]\n')
        self.teXfileObj.write('\t\centering \n')
        self.teXfileObj.write('\\begin{minipage}{0.45\linewidth}\n')
        self.teXfileObj.write('\t\includegraphics[width=1.15\linewidth]{figures/' + figIName + '}\n')
    
        figlabel = '\label{fig:quadI:lambda' + str(lambda0).ljust(5,'0') + '}'
        self.teXfileObj.write('\caption{ $\left<I\\right>=\left<a+a^\dagger\\right>$ ' + figlabel + ' }\n')
        self.teXfileObj.write('\end{minipage}\n')
        self.teXfileObj.write('\hspace{0.025\\textwidth}\n')
        self.teXfileObj.write('\\begin{minipage}{0.45\linewidth}\n')
        self.teXfileObj.write('\t\includegraphics[width=1.15\linewidth]{figures/' + figVarName + '}\n')
    
        figlabel = '\label{fig:varI:lambda' + str(lambda0).ljust(5,'0') + '}'
        self.teXfileObj.write('\caption{ $\Delta I$ ' + figlabel + ' }\n')
        self.teXfileObj.write('\end{minipage}\n')
        self.teXfileObj.write('\end{figure}\n')
        self.teXfileObj.write('\end{frame}\n')
        self.teXfileObj.write('\n');
    ###########################################################################
