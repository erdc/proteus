#!/usr/bin/env python
"""
Class and script for generating a report from simulation data.

.. inheritance-diagram:: proteus.LatexReport
   :parts: 1
"""
from .Profiling import logEvent

def openLatexReport(filename,reportname):
    latexReport = open(filename,'w')
    latexReport.write(r"""\documentclass{amsart}
\usepackage{epsfig}
\setlength{\oddsidemargin}{0in}
\setlength{\evensidemargin}{0in}
\setlength{\topmargin}{-0.25in}
\setlength{\headheight}{0.1in}
\setlength{\headsep}{0.15in}
\setlength{\topskip}{0in}
\setlength{\footskip}{0.15in}
\setlength{\textwidth}{6.5in}
\setlength{\textheight}{9in}
\begin{document}
\begin{center}
{\bf"""+(' %s ' % reportname)+
r"""}\end{center}
\tableofcontents
""")
    return latexReport

def closeLatexReport(file):
    file.write('\\end{document}\n')
    file.close()


class LatexResultsSummary(object):
    """
    simple steps for taking simulation results and generating Latex
    Table of results
    """
    import pickle
    def __init__(self,resFileName,repFileName,repName=None):
        self.resFileName = resFileName
        self.repFileName = repFileName
        self.repName     = repName
        self.results     = None
        self.report      = None
        self.invalidErrorNormCombinations = ['error_velocity_TV',
                                             'error_velocity_H1',
                                             'error_velocity_H1semi',
                                             'error_velocity_W11',
                                             'error_velocity_W11semi']
    #end def
    def open(self,append=False):
        import os
        if not (os.path.exists(self.resFileName) or os.path.exists(self.resFileName + os.extsep + "dir")): #tjp hack
            logEvent("""LatexResSumm resFile= %s not found! """ % self.resFileName)
            return True
        import shelve
        self.results = shelve.open(self.resFileName)
        if self.repName is None:
            repRaw = self.results['flags']['simulationName']
            self.repName = repRaw.replace('_','-')
        self.report = openLatexReport(self.repFileName,self.repName)

        #check ok
        return False
    def close(self):
        closeLatexReport(self.report)
        return False

    def generateSpatialConvTable(self,time='Last',useRelativeError=False,relativeErrorEps=1.0e-10):
        """
        produce table that has spatial error for different levels
        at a given time. If useRelativeError is true, the error values are normalized
        by the "exact" solution norm on the finest spatial mesh with a fudge factor
        given by relativeErrorEps

        component & level & h & error & rate \\
        """
        import math #for rate calc
        nLevels = len(self.results['simulationData']['spatialMesh'])
        nErrors = len(self.results['flags']['errorNorms'])
        preamble = r"""
\begin{table}[h!]
\begin{tabular}{%s}
"""
        postamble = r"""
\hline
\end{tabular}
\end{table}
"""
        tabcols = "|l|c|c|"
        rowform0 = """%s[%d] & %s & %g """
        for enorm in self.results['flags']['errorNorms']:
            if enorm is not None:
                tabcols += "c|c|"
        #
        #treat mass conservation results differently for now
        if 'localMassBalance' in self.results['flags']['errorTypes']:
            tabcols += "c|"
        if 'globalHeavisideMassBalance' in self.results['flags']['errorTypes']:
            tabcols += "c|"
        self.report.write(preamble % tabcols)
        self.report.write("""\\hline\n""")
        rowtitle = """err. comp & level & h """
        #assumes every component and velocity gets same error norms
        for enorm in self.results['flags']['errorNorms']:
            if enorm is not None:
                rowtitle += " & "
                rowtitle += enorm.replace('_','-')
                rowtitle += """ & rate """
        if 'localMassBalance' in self.results['flags']['errorTypes']:
            rowtitle += """ & local mass err """
        if 'globalHeavisideMassBalance' in self.results['flags']['errorTypes']:
            rowtitle += """ & $L_1(H(u))$ mass err """
        self.report.write(rowtitle + """\\\\ \n""")
        self.report.write("""\\hline\n""")


        for errorQuantity in self.results['flags']['errorQuantities']:
            elabelBase = 'error '+errorQuantity
            ekeyBase   = 'error_'+errorQuantity

            #now try relative error if desired
            exKeyBase = 'exact_'+errorQuantity

            #begin section ruling out combinations that don't make sense
            #need to put this logic somewhere else once and for all
            computeLocalMassBalErr = False; computeGlobalHeavisideMassBalErr = False
            rowend = """ """
            if ('localMassBalance' in self.results['flags']['errorTypes'] and
                errorQuantity == 'u'):
                computeLocalMassBalErr = True
                rowend += """ & % g """
            elif ('localMassBalance' in self.results['flags']['errorTypes'] and
                errorQuantity != 'u'):
                rowend += """ & - """
            if ('globalHeavisideMassBalance' in self.results['flags']['errorTypes'] and
                errorQuantity == 'u'):
                computeGlobalHeavisideMassBalErr = True
                rowend += """ & % g """
            elif ('globalHeavisideMassBalance' in self.results['flags']['errorTypes'] and
                errorQuantity != 'u'):
                rowend += """ & - """

            for il in range(nLevels):
                for ci in self.results['flags']['components']:
                    h  = self.results['simulationData']['spatialMesh'][il]['h'][-1]
                    row = rowform0 % (elabelBase,ci,il,h)
                    for enorm in self.results['flags']['errorNorms']:
                        if enorm is not None:
                            elabel = elabelBase+' '+enorm
                            ekey   = ekeyBase+'_'+enorm
                            exkey  = exKeyBase+'_'+enorm
                            #need way to rule out some error-quantity combiations
                            if ekey in self.invalidErrorNormCombinations:
                                row += """ & - & - """
                                break
                            if isinstance(self.results['errorData'][ci][il][ekey],list):
                                error = self.results['errorData'][ci][il][ekey][-1]
                            else:
                                error = self.results['errorData'][ci][il][ekey]
                            if useRelativeError == True:
                                #normalize by error on finest grid?
                                if isinstance(self.results['errorData'][ci][nLevels-1][exkey],list):
                                    exact = self.results['errorData'][ci][nLevels-1][exkey][-1]
                                else:
                                    exact = self.results['errorData'][ci][nLevels-1][exkey]
                                if abs(exact) < relativeErrorEps:
                                    exact += relativeErrorEps
                                error = error/exact
                            if il == 0:
                                row += """ & %g & %s """ % (error,'-')
                            else:
                                if isinstance(self.results['errorData'][ci][il-1][ekey],list):
                                    errM=self.results['errorData'][ci][il-1][ekey][-1]
                                else:
                                    errM=self.results['errorData'][ci][il-1][ekey]
                                hM = self.results['simulationData']['spatialMesh'][il-1]['h'][-1]
                                if useRelativeError == True: #only normalizing by finest mesh val
                                    rate = math.log((errM/exact+1.0e-24)/(error+1.0e-24))/math.log(hM/h)
                                else:
                                    rate = math.log((errM+1.0e-24)/(error+1.0e-24))/math.log(hM/h)

                                row += """& %g & %g """ % (error,rate)
                    if computeLocalMassBalErr or computeGlobalHeavisideMassBalErr:
                        if computeLocalMassBalErr:
                            if isinstance(self.results['errorData'][ci][il]['localMassBalance'],list):
                                row += rowend % self.results['errorData'][ci][il]['localMassBalance'][-1]
                            else:
                                row += rowend % self.results['errorData'][ci][il]['localMassBalance']
                        if computeGlobalHeavisideMassBalErr:
                            if isinstance(self.results['errorData'][ci][il]['globalHeavisideMassF'],list):
                                row += rowend % abs(self.results['errorData'][ci][il]['globalHeavisideMassF'][-1]-
                                                    self.results['errorData'][ci][il]['globalHeavisideMass0'])
                            else:
                                row += rowend % abs(self.results['errorData'][ci][il]['globalHeavisideMassF']-
                                                    self.results['errorData'][ci][il]['globalHeavisideMass0'])
                    else:
                        row += rowend
                    row += """\\\\ \n"""
                    self.report.write(row)
                #ci
            #il
        #end error quantity
        self.report.write(postamble)
    #def
#class

if __name__ == '__main__':
    import os
    import sys
    PROTEUS_HOME = os.getenv('PROTEUS_HOME',os.getenv('HOME')+'/Public/code/proteus')
    PROTEUS_SRC  = PROTEUS_HOME+'/src'
    if PROTEUS_SRC not in sys.path:
        sys.path.insert(0,PROTEUS_SRC)
    #
    import optparse
    usage = "usage: %prog [options] resultsFile [reportFile]"
    parser = optparse.OptionParser(usage=usage)

    parser.add_option("-T","--times",
                      help="which times to calculate error: "
                      "  All, Last, tList ",
                      dest="errorTimes",
                      type="string",
                      default="Last")
    parser.add_option("--relerr",
                      help="use relative error in tables?",
                      dest="useRelativeError",
                      action="store_true",
                      default=False)
    (opts,args) = parser.parse_args()

    if len(args) < 1:
        raise RuntimeError("No results file specified")
    else:
        resfile = args[0]
    if len(args) > 1:
        repfile = args[1]
    else:
        repfile = 'dummy.tex'
    rep = LatexResultsSummary(resfile,repfile)

    rep.open()

    #mwf debug
    #print """rep.res[flags]= %s \n""" % rep.results['flags']
    #mwf debug
    #print """\n rep.res[simData]= %s \n""" % rep.results['simulationData']
    #mwf debug
    #print """\n rep.res[errData]= %s \n""" % rep.results['errorData']
    #print """nLevels = %d """ % len(rep.results['simulationData']['spatialMesh'])
    #print """nErrors = %d """ % len(rep.results['flags']['errorNorms'])


    rep.generateSpatialConvTable(useRelativeError=opts.useRelativeError)

    rep.close()

    latexcmd = """latex %s """ % repfile
    dvifile = repfile.replace('.tex','.dvi')
    os.system(latexcmd)
    dvicmd   = """xdvi %s """ % dvifile
    os.system(dvicmd)

    #raw_input('press return to exit')
