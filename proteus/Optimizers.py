"""
Class hierarchies for working with minimization problems

.. inheritance-diagram:: proteus.Optimizers
   :parts: 1
"""
#John Chrispell, Summer 07
class ObjectiveFunction_base(object):
    def __init__(self,LHS_x,RHS_x):
        self.LHS_x=LHS_x
        self.RHS_x=RHS_x
    def getResidual(self,x):
        F=0.0
        return F

class MinAlgo_base(object):
    """
    This is the base class for the minimization algorthms class,
    and can be used to find the min of a 1-D function.

    members:
    """
    def __init__(self):
        pass

class fminbound(MinAlgo_base):
    """
    This is a function that uses a golden ratio routine for finding
    the minimum of a 1-D function. It returns the pair (xmin,f(xmin))
    when its solver is called. It assumes that you have a bracketed
    already bracketed the minimum between two end points, and provide
    some good initial guess. Note there is also an endpoint hack.
    """
    def __init__(self,FuncToMinimize,Tol=1e-11):
        MinAlgo_base.__init__(self)
        self.FuncToMinimize = FuncToMinimize
        self.Tol = Tol
        self.GoldenFactor = 0.61803399
        self.maxIts = 50

    def solve(self,Guess_x=0.38196601):
        self.x0  = self.FuncToMinimize.LHS_x
        self.x3  = self.FuncToMinimize.RHS_x
        EndCheckFlag = 0
        while(EndCheckFlag <=1):

            if(abs(self.x3-Guess_x) > abs(Guess_x-self.x0)):
                self.x1 = Guess_x
                self.x2 = self.x1 + (1.0-self.GoldenFactor)*(self.x3 - self.x1)
            else:
                self.x2 = Guess_x
                self.x1 = self.x2 - (1.0-self.GoldenFactor)*(self.x2 - self.x0)

            #initial function evaluations
            F1 = self.FuncToMinimize.getResidual(self.x1)
            F2 = self.FuncToMinimize.getResidual(self.x2)
            its = 0
            while(abs(self.x3 - self.x0) > self.Tol*(abs(self.x1)+abs(self.x2)) and
                    its < self.maxIts):
                if(F2 < F1):
                    self.x0 = self.x1
                    self.x1 = self.x2
                    self.x2 = self.GoldenFactor*self.x1 + (1.0-self.GoldenFactor)*self.x3
                    F1 = F2
                    F2 = self.FuncToMinimize.getResidual(self.x2)
                else:
                    self.x3 = self.x2
                    self.x2 = self.x1
                    self.x1 = self.GoldenFactor*self.x2 + (1.0-self.GoldenFactor)*self.x0
                    F2 = F1
                    F1 = self.FuncToMinimize.getResidual(self.x1)
                its += 1
            if(F1 < F2):
                Golden = F1
                xmin   = self.x1
            else:
                Golden = F2
                xmin   = self.x2

            """ End Point Hack """
            # Look at the value on the end points
            Fleftend  = self.FuncToMinimize.getResidual(self.FuncToMinimize.LHS_x)
            Frightend = self.FuncToMinimize.getResidual(self.FuncToMinimize.RHS_x)

            if(EndCheckFlag == 1):
                # 2nd time (runs again when previous solution was an end point.
                EndCheckFlag = EndCheckFlag + 1

                Fleftend  = self.FuncToMinimize.getResidual(NewLHS_x)
                Frightend = self.FuncToMinimize.getResidual(NewRHS_x)
                if(Fleftend < Frightend):
                    Fendpt = Fleftend
                    xendpt = NewLHS_x
                else:
                    Fendpt = Frightend
                    xendpt = NewRHS_x

            if(EndCheckFlag == 0):
                # First time
                Fleftend  = self.FuncToMinimize.getResidual(self.FuncToMinimize.LHS_x)
                Frightend = self.FuncToMinimize.getResidual(self.FuncToMinimize.RHS_x)
                if(Fleftend < Frightend):
                    Fendpt = Fleftend
                    xendpt = self.FuncToMinimize.LHS_x
                    NewLHS_x = self.FuncToMinimize.LHS_x
                    # Adjust the interval properly.
                    if (self.FuncToMinimize.LHS_x > self.FuncToMinimize.RHS_x):
                        NewRHS_x = self.FuncToMinimize.LHS_x - .2*abs(abs(self.FuncToMinimize.RHS_x) - abs(self.FuncToMinimize.LHS_x))
                    else:
                        NewRHS_x = self.FuncToMinimize.LHS_x + .2*abs(abs(self.FuncToMinimize.RHS_x) - abs(self.FuncToMinimize.LHS_x))
                    self.x0  = NewLHS_x
                    self.x3  = NewRHS_x
                else:
                    Fendpt = Frightend
                    xendpt = self.FuncToMinimize.RHS_x
                    NewRHS_x = self.FuncToMinimize.RHS_x
                    # Adjust the interval properly
                    if ( self.FuncToMinimize.RHS_x < self.FuncToMinimize.LHS_x):
                        NewLHS_x = self.FuncToMinimize.RHS_x + .2*abs(abs(self.FuncToMinimize.RHS_x) - abs(self.FuncToMinimize.LHS_x))
                    else:
                        NewLHS_x = self.FuncToMinimize.RHS_x - .2*abs(abs(self.FuncToMinimize.RHS_x) - abs(self.FuncToMinimize.LHS_x))
                    self.x0  = NewLHS_x
                    self.x3  = NewRHS_x

                # Compair the ends with the value obtained.
                if(Fendpt < Golden):                          # Set to run again
                    EndCheckFlag = EndCheckFlag + 1
                    Guess_x      = 0.5*(NewLHS_x + NewRHS_x)
                else:                                         # Skip to the end
                    EndCheckFlag = EndCheckFlag + 2

            # Compair the ends with the value obtained.
            if(Fendpt < Golden):                          # Set to run again
                Golden = Fendpt
                xmin   = xendpt

        return (xmin,Golden)

if __name__ == '__main__':
     # This is a test driver for the Search Algorithm function
    from .ObjectiveFunctions import *
    import Gnuplot
    import numpy

    # Test parameters
    LHS_x   = -1.0
    RHS_x   = 2.0
    Guess_x = (RHS_x+LHS_x)/2.0
    Tol     = 1e-6
    RefnVal = 100

    """ Choose a Test Function from the following list """
    #F = SimpelFunc(LHS_x,RHS_x)
    F = SimpelFunc2(LHS_x,RHS_x)
    #F = SimpelFunc3(LHS_x,RHS_x)

    # Find the minimum on the given interval.
    solver = fminbound(FuncToMinimize=F,Tol=Tol)
    (xmin,fmin) = solver.solve(Guess_x)

    # ------ Plotting ------- #
    xLeft=LHS_x
    dx = (abs(RHS_x-LHS_x))/RefnVal
    xVec=[]
    fVec=[]
    for i in range(RefnVal+1):
        x    = xLeft+i*dx
        Fval = F.getResidual(x)
        xVec.append(x)
        fVec.append(Fval)
    #print "x",xVec
    #print "s",sVec
    g = Gnuplot.Gnuplot(debug=1)
    g.title('Plot of Function')     # (optional)
    g('set data style linespoints') # give gnuplot an arbitrary command
    d = Gnuplot.Data(numpy.array(xVec),numpy.array(fVec))
    dmin = Gnuplot.Data(xmin,fmin)
    g.xlabel('x')
    g.ylabel('f(x)')
    g.plot(d, dmin)
    input('Please press return to continue...\n')

    print("Minimum at (x,f(x)) = (",xmin,",",fmin,")")
