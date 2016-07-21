""" Helper functions commonly used in tests. """

import os
import sys
import inspect
import pickle

def addSubFolders(currentframe):
    """Add import_modules and comparison_files to sys.path

    Attributes
    ----------
    currentframe: ?
    """
    cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( currentframe ))[0]))

    if cmd_folder not in sys.path:
        sys.path.insert(0,cmd_folder)

    cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( currentframe ))[0],
                                                                  "import_modules")))

    if cmd_subfolder not in sys.path:
        sys.path.insert(0,cmd_subfolder)

class SimulationTest():
    """Base class for simulation tests.
    
    TODO (ARB) - Is this even helpful?
    """

class NumericResults:
    """Parse and stores numerical data from a Proteus log file.

    Attributes
    ----------
    data_dictionary : dict
        A data_dictionary that stores the NumericalAnalytics data.
    data_dictionary_header : dict
        A dictionary that stores information about the simulation.
    """
    def __init__(self,file_name):
        """ Initializes the Numeric Results class 

        Parameters
        ----------
        file_name : str
            Provides the name of a proteus log file.
        """
        self.log_file = open(file_name,'r')
        self._parse_file()
   
    def _parse_file(self):
        """ Collect the NumericalAnalytics log entries """
        import re

        # TODO -
        # add unknowns
        # add fem spaces
        self.data_dictionary_header = {}
        self.data_dictionary = {}
        self.data_dictionary_header['Petsc'] = []

        NA_header     = re.compile("NAHeader(.*)$")
        NA_petsc      = re.compile("PETScOptions(.*)$")
        NA_refine     = re.compile("GridRefinements(.*)$")
        NA_timesteps  = re.compile("Num Time Steps(.*)$")
        NA_precond    = re.compile("Preconditioner(.*)$")
        NA_pattern    = re.compile("NumericalAnalytics(.*)$")
        NA_time       = re.compile("Time Step(.*)$")
        NA_level      = re.compile("Newton iteration for level(.*)$")
        NA_newton     = re.compile("NewtonNorm: (.*)$")
        NA_outersolve = re.compile("KSPOuterResidual: (.*)$")
        NA_innersolve = re.compile("KSPSchurResidual: (.*)$")

        for line in self.log_file:
            for match in re.finditer(NA_header,line):
                new_line = NA_header.search(line).groups()[0]

                for petsc_match in re.finditer(NA_petsc,new_line):
                    petsc_match = str(NA_petsc.search(new_line).groups()[0])
                    self.data_dictionary_header['Petsc'].append(petsc_match)

                for refine_match in re.finditer(NA_refine,new_line):
                    refine_match = int(NA_refine.search(new_line).groups()[0])
                    self.data_dictionary_header['Grid Refinements'] = refine_match

                for precon_match in re.finditer(NA_precond,new_line):
                    precon_match = str(NA_precond.search(new_line).groups()[0])
                    self.data_dictionary_header['Preconditioner'] = precon_match 
                
                for timestep_match in re.finditer(NA_timesteps,new_line):
                    timestep_match = int(NA_timesteps.search(new_line).groups()[0])
                    self.data_dictionary_header['Num Time Steps'] = timestep_match

            for match in re.finditer(NA_pattern,line):
                new_line = NA_pattern.search(line).groups()[0]

                for time_match in re.finditer(NA_time,new_line):
                    time_match = float(NA_time.search(new_line).groups()[0])
                    self.data_dictionary[time_match] = {}

                for level_match in re.finditer(NA_level,new_line):
                    level_match = int(NA_level.search(new_line).groups()[0])
                    self.data_dictionary[time_match][level_match] = [[],{}]

                for newton_match in re.finditer(NA_newton,new_line):
                    newton_match = float(NA_newton.search(new_line).groups()[0])
                    self.data_dictionary[time_match][level_match][0].append(newton_match)
                    newton_it_key = len(self.data_dictionary[time_match][level_match][0])-1
                    self.data_dictionary[time_match][level_match][1][newton_it_key] = [[],{}]
                
                for outerksp_match in re.finditer(NA_outersolve,new_line):
                    outerksp_match = float(NA_outersolve.search(new_line).groups()[0])
                    self.data_dictionary[time_match][level_match][1][newton_it_key][0].append(outerksp_match)
                    outerksp_match_key = len(self.data_dictionary[time_match][level_match][1][newton_it_key][0])-1
                    self.data_dictionary[time_match][level_match][1][newton_it_key][1][outerksp_match_key] = []
                
                for innerksp_match in re.finditer(NA_innersolve,new_line):
                    innerksp_match = float(NA_innersolve.search(new_line).groups()[0])
                    self.data_dictionary[time_match][level_match][1][newton_it_key][1][outerksp_match_key].append(innerksp_match)

    def pickle_data(self,filename):
        """Pickle the dictionary created after parsing the file.

        Arguments
        ---------
        filename : str
            filename for dictionary storeage
        """
        fileObject = open(filename,'wb')
        pickle.dump(self.data_dictionary,fileObject)
        fileObject.close()

    def print_info(self):
        """ Output a variety of information about the data-structure """
        print " **** HEADER INFORMATION ****"
        for key in self.data_dictionary_header.keys():
            if key == 'Petsc':
                self.__print_petsc_info()
            else:
                print `key` + '   :   ' + `self.data_dictionary_header[key]`
        print " *** VALID KEYS ***"
        for key in self.data_dictionary.keys():
            print `key`

    def __print_petsc_info(self):
        """ Prints the settings given in the PETSc command line """
        for term in self.data_dictionary_header['Petsc']:
            print term

    def print_header_latex(self):
        """ Prints the header information in a latex consistent format """
        pass


    def __init_ipython_plot(self,plot_data,legend_lst=[],title_str=' ',axis=None):
        """ Private member function that loads ipython plotting libraries 

        Parameters
        ----------
        plot_data : lst of lst
            A list of lists that stores the data series to be plotted.
        legend : lst
            A list of strings for the legend.
        """
        import matplotlib
        import numpy as np
        import matplotlib.pyplot as plt       
        for data_set in plot_data:
            plt.plot(data_set)
        plt.yscale("log")
        plt.legend(legend_lst)
        plt.title(title_str)
        if axis!=None:
            plt.xlim(axis[0],axis[1])
        plt.show() 

    def ipython_plot_newton_it(self,time_level):
        """ Plot the Newton iteration residual in a jupyter notebook.

        Parameters
        ----------
        time_level : lst of tuples
            A list of tuples with the structure (time,level) to be ploted.
        """
        plot_data = []
        legend = []
        title = 'Residual From Newton Iterations'

        for data_set in time_level:
            if data_set[0] in self.data_dictionary.keys():
                if data_set[1] in self.data_dictionary[data_set[0]].keys():
                    plot_data.append(self.data_dictionary[data_set[0]][data_set[1]][0])
                    legend.append((data_set[0],data_set[1]))
                else:
                    print 'The second key ' + `data_set[1]` + ' is not valid.'
            else:
                print 'The first key ' + `data_set[1]` + ' is not valid.'
    
        self.__init_ipython_plot(plot_data,legend,title)

    def ipython_plot_ksp_residual(self,time_level_it):
        """ Plot the outer KSP residual in a jupyter notebook.
        
        Parameters
        ----------
        time_level_it :
        """
        plot_data = []
        legend = []
        title = 'Residuals of Outer Most KSP Solve.'
        axis = [0,8]

        for data_set in time_level_it:
            if data_set[0] in self.data_dictionary.keys():
                if data_set[1] in self.data_dictionary[data_set[0]].keys():
                    if data_set[2] in self.data_dictionary[data_set[0]][data_set[1]][1].keys():
                        plot_data.append(self.data_dictionary[data_set[0]][data_set[1]][1][data_set[2]][0])
                        legend.append((data_set[0],data_set[1],data_set[2]))
                    else:
                        print 'The third key ' + `data_set[1]` + ' is not valid.'
                else:
                    print 'The second key ' + `data_set[1]` + ' is not valid.'
            else:
                print 'The first key ' + `data_set[1]` + ' is not valid.'
                
        self.__init_ipython_plot(plot_data,legend,title,axis)

    def ipython_plot_ksp_schur_residual(self,time_level_it,axis=False):
        """ Plot the inner KSP residual in a jupyter notebook.
        
        Parameters
        ----------
        time_level_it :
        """
        plot_data = []
        legend = []
        title = 'Residuals of Inner Schur KSP Solve.'
        axis_inline = axis

        for data_set in time_level_it:
            if data_set[0] in self.data_dictionary.keys():
                if data_set[1] in self.data_dictionary[data_set[0]].keys():
                    if data_set[2] in self.data_dictionary[data_set[0]][data_set[1]][1].keys():
                        if data_set[3] in self.data_dictionary[data_set[0]][data_set[1]][1][data_set[2]][1].keys():
                            plot_data.append(self.data_dictionary[data_set[0]][data_set[1]][1][data_set[2]][1][data_set[3]])
                            legend.append((data_set[0],data_set[1],data_set[2],data_set[3]))
                    else:
                        print 'The third key ' + `data_set[1]` + ' is not valid.'
                else:
                    print 'The second key ' + `data_set[1]` + ' is not valid.'
            else:
                print 'The first key ' + `data_set[1]` + ' is not valid.'
            
        if axis!=False:
            self.__init_ipython_plot(plot_data,legend,title,axis_inline)
        else:
            self.__init_ipython_plot(plot_data,legend,title)
