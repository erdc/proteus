"""
A class hierarchy for numerical integration on reference domains in 1,2, and 3D.

.. inheritance-diagram:: proteus.Quadrature
   :parts: 1
"""
from EGeometry import *
from .Profiling import logEvent
from math import *

class Q_base:
    """
    The base class for quadrature methods.
    """
    def __init__(self,order=1):
        self.points=()
        self.weights=()
        self.pointsAll=(())
        self.weightsAll=(())
        self.order=order
    def setOrder(self,k):
        self.order=k
        if self.order > len(self.pointsAll):
            logEvent("WARNING Q_base requested order=%d > max allowed=%d setting order to max" % (self.order,
                                                                                               len(self.pointsAll)))
            self.order= len(self.pointsAll)
        self.points=self.pointsAll[self.order-1]
        self.weights=self.weightsAll[self.order-1]

class GaussPoint(Q_base):
    """
    A dummy class for integrating the boundary of the unit interval (i.e. a point).
    """
    def __init__(self,order=1):
        Q_base.__init__(self,order)
        self.pointsAll=(
            (EVec(0.0),),)
        self.weightsAll=(
            (1.0,),)
        self.setOrder(order)
    def setOrder(self,k):
        self.order = k
        self.points = self.pointsAll[0]
        self.weights = self.weightsAll[0]

LobattoPoint = GaussPoint

class GaussEdge(Q_base):
    """
    Gaussian Quadrature on the unit interval.
    """
    def __init__(self,order=1):
        #mwf for convenience, see Guern and Ermond
        a1 = 0.5*sqrt((15.0+2.0*sqrt(30))/35.0)
        a2 = 0.5*sqrt((15.0-2.0*sqrt(30))/35.0)
        w1 = 0.25 - sqrt(5./6.)/12.
        w2 = 0.25 + sqrt(5./6.)/12.
        Q_base.__init__(self,order)
        self.pointsAll=(
            (EVec(0.5),),
            (EVec((sqrt(3.0)-1.0)/(2.0*sqrt(3.0))),
             EVec((sqrt(3.0)+1.0)/(2.0*sqrt(3.0)))),
            (EVec((sqrt(5.0) - sqrt(3.0))/(2.0*sqrt(5))),
             EVec(0.5),
             EVec((sqrt(5.0) + sqrt(3.0))/(2.0*sqrt(5)))),
            (EVec(0.5+a1),EVec(0.5-a1),
             EVec(0.5+a2),EVec(0.5-a2)),
            (EVec(0.5),
             EVec(0.5*(sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0) + 0.5),
             EVec(0.5*(-sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0) + 0.5),
             EVec(0.5*(sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0) + 0.5),
             EVec(0.5*(-sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0) + 0.5))
            )
        self.weightsAll=(
            (1.0,),
            (0.5,
             0.5),
            (5.0/18.0,
             8.0/18.0,
             5.0/18.0),
            (w1,w1,w2,w2),
            (0.5*(128.0/225.0),
             0.5*(322.0+13.0*sqrt(70.0))/900.0,
             0.5*(322.0+13.0*sqrt(70.0))/900.0,
             0.5*(322.0-13.0*sqrt(70.0))/900.0,
             0.5*(322.0-13.0*sqrt(70.0))/900.0)
            )
        self.setOrder(order)

    def setOrder(self,order,domain=[0.0,1.0]):
        Q_base.setOrder(self,order)

        points = self.points
        weights = self.weights
        self.points = []
        self.weights = []
        for i in range(order):
            self.points.append(EVec((domain[1]-domain[0])*points[i][0] +domain[0]))
            self.weights.append((domain[1]-domain[0])*weights[i])

class LobattoEdge(Q_base):
    """
    Gauss-Lobatto quadrature on the unit interval.
    """
    def __init__(self,order=1):
        Q_base.__init__(self,order)
        a1 = 0.2*sqrt(5.0)
        a2 = sqrt(21.0)/7.0
        self.pointsAll=(
            (EVec(0.0),EVec(1.0)),
            (EVec(0.0),EVec(0.5),EVec(1.0)),
            (EVec(0.0),EVec(0.5*(1.-a1)),EVec(0.5*(1.+a1)),EVec(1.0)),
            (EVec(0.0),EVec(0.5*(1.-a2)),EVec(0.5),EVec(0.5*(1.+a2)),EVec(1.0)),
            )
        self.weightsAll=(
            (0.5,0.5),
            (1.0/6.0, 4.0/6.0, 1.0/6.0),
            (1.0/12.0, 5.0/12.0,5.0/12.0,1.0/12.0),
            (1.0/20.0,49.0/180.0,32.0/90.0,49.0/180.0,1.0/20.0)
            )
        self.setOrder(order)

class LobattoEdgeAlt(Q_base):
    """
    Gauss-Lobatto quadrature on the [-1:1] interval.
    """
    def __init__(self,order=1):
        Q_base.__init__(self,order)
        a1 = 0.2*sqrt(5.0)
        a2 = sqrt(21)/7.0
        self.pointsAll=(
            (EVec(-1.0),EVec(1.0)),
            (EVec(-1.0),EVec(0.0),EVec(1.0)),
            (EVec(-1.0),EVec(-a1),EVec(a1) ,EVec(1.0)),
            (EVec(-1.0),EVec(-a2),EVec(0.0),EVec(a2),EVec(1.0)),)
        self.weightsAll=(
            (1.0,1.0),
            (1.0/3.0, 4.0/3.0, 1.0/3.0),
            (1.0/6.0, 5.0/6.0,5.0/6.0,1.0/6.0),
            (1.0/10.0,49.0/90.0,32.0/45.0,49.0/90.0,1.0/10.0),)
        self.setOrder(order)


class CompositeTrapezoidalEdge(Q_base):
    """
    Composite trapezoidal rule on the unit interval.
    order is number of intervals
    """
    def __init__(self,order=1,maxOrder=20):
        assert order > 0, "Composite Trapezoidal Rule requires order > 0"
        self.maxOrder = maxOrder
        Q_base.__init__(self,order)
        pointsList = []
        weightsList= []
        for nintervals in range(1,self.maxOrder+1):
            dx = 1.0/float(nintervals)
            points=numpy.arange(nintervals+1,dtype='d')*dx
            weights=numpy.zeros(nintervals+1,'d'); weights.fill(dx)
            weights[0] *= 0.5; weights[-1] *= 0.5
            pointsList.append([EVec(p) for p in points])
            weightsList.append([w for w in weights])
        self.pointsAll =tuple(tuple(pL) for pL in pointsList)
        self.weightsAll=tuple(tuple(wL) for wL in weightsList)
        self.setOrder(order)
    def setOrder(self,k):
        assert k > 0, "Composite Trapezoidal Rule requires order > 0"
        assert k <= self.maxOrder, "Composite Trapezoidal Rule k= %s maxOrder= %s need to increase in ctor " % (k,self.maxOrder)
        self.order = k
        self.points = self.pointsAll[k-1]
        self.weights = self.weightsAll[k-1]

class FaceBarycenterEdge(Q_base):
    def __init__(self,order=1):
        order=1
        Q_base.__init__(self,order)
        self.pointsAll=(
            (EVec(1.0),EVec(0.0)),)
        self.weightsAll=(
            (0.5,0.5),)
        self.setOrder(order)
    def setOrder(self,k):
        self.order = 1
        self.points = self.pointsAll[0]
        self.weights = self.weightsAll[0]

class GaussTriangle(Q_base):
    """
    Gauss quadrature on the unit triangle.
    """
    def __init__(self,order=1):
        #for convenience, see Guern and Ermond 360
        #4th order
        a1 = 0.445948490915965
        a2 = 0.091576213509771
        w1 = 0.223381589678010*0.5
        w2 = 0.109951743655322*0.5
        #5th order
        a5= ((6.-sqrt(15.0))/21.0,(6. + sqrt(15.0))/21.0)
        w5= (0.5*(155. - sqrt(15))/1200.,0.5*(155. + sqrt(15))/1200.)
        fifthOrderPoints = []; fifthOrderWeights = []
        for i in range(2):
            fifthOrderPoints.append((a5[i],a5[i])); fifthOrderWeights.append(w5[i])
            fifthOrderPoints.append((1.-2.*a5[i],a5[i])); fifthOrderWeights.append(w5[i])
            fifthOrderPoints.append((a5[i],1.-2.*a5[i])); fifthOrderWeights.append(w5[i])
        fifthOrderPoints.append((1./3.,1./3.)); fifthOrderWeights.append(0.5*9.0/40.0)

        #6th order
        a6  = (0.063089014491502,0.249286745170910)
        a6a = 0.310352451033785
        a6b = 0.053145049844816
        w6  =(0.5*0.050844906370206,0.5*0.116786275726378)
        w6ab= 0.5*0.082851075618374
        sixthOrderPoints = []; sixthOrderWeights = []
        for i in range(2):
            sixthOrderPoints.append((a6[i],a6[i])); sixthOrderWeights.append(w6[i])
            sixthOrderPoints.append((1.-2.*a6[i],a6[i])); sixthOrderWeights.append(w6[i])
            sixthOrderPoints.append((a6[i],1.-2.*a6[i])); sixthOrderWeights.append(w6[i])
        abPerms = [(a6a,a6b),(1.-a6a-a6b,a6a),(a6b,1.-a6a-a6b),
                   (a6b,a6a),(1.-a6a-a6b,a6b),(a6a,1.-a6a-a6b)]
        abWeights= [w6ab for i in range(6)]
        sixthOrderPoints.extend(abPerms); sixthOrderWeights.extend(abWeights)

        Q_base.__init__(self,order)
        self.pointsAll=(
            ( EVec(1.0/3.0,1.0/3.0),),
            ( EVec(1.0/2.0,1.0/2.0), EVec(0.0,1.0/2.0), EVec(1.0/2.0,0.0)),
            ( EVec(1.0/3.0,1.0/3.0), EVec(3.0/5.0,1.0/5.0), EVec(1.0/5.0,3.0/5.0),EVec(1.0/5.0,1.0/5.0)),
            ( EVec(a1,a1),EVec(1.0-2.0*a1,a1),EVec(a1,1.0-2.0*a1),
              EVec(a2,a2),EVec(1.0-2.0*a2,a2),EVec(a2,1.0-2.0*a2)),
            tuple(EVec(p[0],p[1]) for p in fifthOrderPoints),
            tuple(EVec(p[0],p[1]) for p in sixthOrderPoints)
            )#points All

        self.weightsAll=(
            (0.5,),
            (1.0/6.0, 1.0/6.0, 1.0/6.0),
            (-27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0),
            (w1,w1,w1,w2,w2,w2),
            tuple(fifthOrderWeights),
            tuple(sixthOrderWeights)
            )#weightsAll
        self.setOrder(order)

class LobattoTriangle(Q_base):
    """
    Gauss-Lobatto quadrature on the unit triangle.
    """
    def __init__(self,order=1):
        #mwf allow higher order now
        #order=1
        Q_base.__init__(self,order)
        #for third order Fekete rule (see Taylor etal SINUM 2000)
        #10 points below are classes of points orbit is number of permuations
        #orbit, 1st two barycentric coordinates, weight to be multiplied by area of triangle
        #1, (1/3,1/3), 0.9
        #3, (0,0),     0.0333333333
        #6, (0,a),     0.1666666667
        #where a = 0.2763932023

        a = 0.2763932023; b = 1.0-a; wv = 0.03333333333333*0.25; wab = 0.16666666666667*0.25
        self.pointsAll=(
            ( EVec(0.0,0.0),EVec(1.0,0.0),EVec(0.0,1.0)),
            ( EVec(1.0/3.0,1.0/3.0),
              EVec(0.,0.),EVec(1.0,0.0),EVec(0.0,1.0),
              EVec(0,a),EVec(b,0),EVec(a,b),EVec(a,0.),EVec(b,a),EVec(0,b)),
            )
        self.weightsAll=(
            (1.0/6.0,1.0/6.0,1.0/6.0),
            (0.5-3.0*wv-6*wab,
             wv,wv,wv,
             wab,wab,wab,wab,wab,wab),
            )
        self.setOrder(order)
    def setOrder(self,k):
        self.order = 1
        self.points = self.pointsAll[0]
        self.weights = self.weightsAll[0]
        if k > 1:
            self.order = 2
            self.points = self.pointsAll[1]
            self.weights = self.weightsAll[1]

class CompositeTrapezoidalTriangle(Q_base):
    """
    Composite trapezoidal rule on the reference triangle
    order is number of intervals
    """
    def __init__(self,order=1,maxOrder=20):
        assert order > 0, "Composite Trapezoidal Rule requires order > 0"
        self.maxOrder = maxOrder
        Q_base.__init__(self,order)
        pointsList = []
        weightsList= []
        parentArea = 0.5
        for nintervals in range(1,self.maxOrder+1):
            dx = 1.0/float(nintervals)
            #uniform subdivisions in terms of barycentric coordinates
            baryvals = numpy.arange(nintervals+1,dtype='d')*dx
            combos = []; weights = []
            for i in range(nintervals+1):
                for j in range(nintervals+1-i):
                    combos.append(EVec(baryvals[i],baryvals[j]))
                    iInBoundary = 0; jInBoundary=0;
                    if i == 0 or i == nintervals: iInBoundary = 1
                    if j == 0 or j == nintervals: jInBoundary = 1
                    if iInBoundary+jInBoundary == 2:
                        weights.append(parentArea/3.0)
                    elif iInBoundary+jInBoundary == 1:
                        weights.append(parentArea)
                    else:
                        weights.append(parentArea*2.0)
            pointsList.append([p for p in combos])
            weightsList.append([w for w in weights])
            parentArea = parentArea*0.25
        self.pointsAll =tuple(tuple(pL) for pL in pointsList)
        self.weightsAll=tuple(tuple(wL) for wL in weightsList)
        self.setOrder(order)
    def setOrder(self,k):
        assert k > 0, "Composite Trapezoidal Rule requires order > 0"
        assert k <= self.maxOrder, "Composite Trapezoidal Rule k= %s maxOrder= %s need to increase in ctor " % (k,self.maxOrder)
        self.order = k
        self.points = self.pointsAll[k-1]
        self.weights = self.weightsAll[k-1]

class FaceBarycenterTriangle(Q_base):
    def __init__(self,order=1):
        order=1
        Q_base.__init__(self,order)
        self.pointsAll=(
            ( EVec(1.0/2.0,1.0/2.0), EVec(0.0,1.0/2.0), EVec(1.0/2.0,0.0)),
            )
        self.weightsAll=(
            (1.0/6.0, 1.0/6.0, 1.0/6.0),
            )
        self.setOrder(order)
    def setOrder(self,k):
        self.order = 1
        self.points = self.pointsAll[0]
        self.weights = self.weightsAll[0]
class GaussTetrahedron(Q_base):
    """
    Gauss-Legendre quadrature on the unit tetrahedron.
    """
    def __init__(self,order=1):
        Q_base.__init__(self,order)
        #mwf for convenience, see Guern Ermond 360
        a1=(7.0-sqrt(15.))/34.0
        a2=(7.0+sqrt(15.))/34.0
        a =(10.0-2.0*sqrt(5.))/40.0
        Vf=1.0/6.0
        w1=(2665.0+14.0*sqrt(15.0))/37800.0*Vf
        w2=(2665.0-14.0*sqrt(15.0))/37800.0*Vf
        wa=10.0/189.0*Vf
        #
        self.pointsAll=(
            (EVec(0.25,0.25,0.25),),
            (EVec(0.58541020,0.13819660,0.13819660),#2nd degree
              EVec(0.13819660,0.58541020,0.13819660),
              EVec(0.13819660,0.13819660,0.58541020),
              EVec(0.13819660,0.13819660,0.13819660)),
            (EVec(0.25,0.25,0.25),#3rd degree
             EVec(1.0/2.0, 1.0/6.0, 1.0/6.0),
             EVec(1.0/6.0, 1.0/2.0, 1.0/6.0),
             EVec(1.0/6.0, 1.0/6.0, 1.0/2.0),
             EVec(1.0/6.0, 1.0/6.0, 1.0/6.0)),
            (EVec(0.50,0.50,0.0),#4th degree
             EVec(0.50,0.0,0.50),
             EVec(0.0,0.50,0.50),
             EVec(0.0,0.0,0.50),
             EVec(0.0,0.50,0.0),
             EVec(0.50,0.0,0.0),
             EVec(0.100526765225204467,0.100526765225204467,0.100526765225204467),
             EVec(0.100526765225204467,0.100526765225204467,0.698419704324386603),
             EVec(0.100526765225204467,0.698419704324386603,0.100526765225204467),
             EVec(0.698419704324386603,0.100526765225204467,0.100526765225204467),
             EVec(0.314372873493192195,0.314372873493192195,0.314372873493192195),
             EVec(0.314372873493192195,0.314372873493192195,0.568813795204234229e-1),
             EVec(0.314372873493192195,0.568813795204234229e-1,0.314372873493192195),
             EVec(0.568813795204234229e-1,0.314372873493192195,0.314372873493192195)),
            (EVec(0.333333333333333333,0.333333333333333333,0.333333333333333333),#5th degree
             EVec(0.0,0.333333333333333333,0.333333333333333333),
             EVec(0.333333333333333333,0.0,0.333333333333333333),
             EVec(0.333333333333333333,0.333333333333333333,0.0),
             EVec(0.25,0.25,0.25),
             EVec(0.909090909090909091e-1,0.909090909090909091e-1,0.909090909090909091e-1),
             EVec(0.727272727272727273,0.909090909090909091e-1,0.909090909090909091e-1),
             EVec(0.909090909090909091e-1,0.727272727272727273,0.909090909090909091e-1),
             EVec(0.909090909090909091e-1,0.909090909090909091e-1,0.727272727272727273),
             EVec(0.433449846426335728,0.665501535736642813e-1,0.665501535736642813e-1),
             EVec(0.665501535736642813e-1,0.433449846426335728,0.665501535736642813e-1),
             EVec(0.665501535736642813e-1,0.665501535736642813e-1,0.433449846426335728),
             EVec(0.433449846426335728,0.433449846426335728,0.665501535736642813e-1),
             EVec(0.433449846426335728,0.665501535736642813e-1,0.433449846426335728),
             EVec(0.665501535736642813e-1,0.433449846426335728,0.433449846426335728)),
            (EVec(0.214602871259151684,0.214602871259151684,0.214602871259151684),#6th degree
             EVec(0.214602871259151684,0.214602871259151684,0.356191386222544953),
             EVec(0.214602871259151684,0.356191386222544953,0.214602871259151684),
             EVec(0.356191386222544953,0.214602871259151684,0.214602871259151684),
             EVec(0.406739585346113397e-1,0.406739585346113397e-1,0.406739585346113397e-1),
             EVec(0.406739585346113397e-1,0.406739585346113397e-1,0.877978124396165982),
             EVec(0.406739585346113397e-1,0.877978124396165982,0.406739585346113397e-1),
             EVec(0.877978124396165982,0.406739585346113397e-1,0.406739585346113397e-1),
             EVec(0.322337890142275646,0.322337890142275646,0.322337890142275646),
             EVec(0.322337890142275646,0.322337890142275646,0.329863295731730594e-1),
             EVec(0.322337890142275646,0.329863295731730594e-1,0.322337890142275646),
             EVec(0.329863295731730594e-1,0.322337890142275646,0.322337890142275646),
             EVec(0.636610018750175299e-1,0.636610018750175299e-1,0.269672331458315867),
             EVec(0.636610018750175299e-1,0.269672331458315867,0.636610018750175299e-1),
             EVec(0.269672331458315867,0.636610018750175299e-1,0.636610018750175299e-1),
             EVec(0.636610018750175299e-1,0.636610018750175299e-1,0.603005664791649076),
             EVec(0.636610018750175299e-1,0.603005664791649076,0.636610018750175299e-1),
             EVec(0.603005664791649076,0.636610018750175299e-1,0.636610018750175299e-1),
             EVec(0.636610018750175299e-1,0.269672331458315867,0.603005664791649076),
             EVec(0.636610018750175299e-1,0.603005664791649076,0.269672331458315867),
             EVec(0.603005664791649076,0.636610018750175299e-1,0.269672331458315867),
             EVec(0.603005664791649076,0.269672331458315867,0.636610018750175299e-1),
             EVec(0.269672331458315867,0.603005664791649076,0.636610018750175299e-1),
             EVec(0.269672331458315867,0.636610018750175299e-1,0.603005664791649076)),
            (EVec(0.50,0.50,0.0),#7th degree
             EVec(0.50,0.0,0.50),
             EVec(0.0,0.50,0.50),
             EVec(0.0,0.0,0.50),
             EVec(0.0,0.50,0.0),
             EVec(0.50,0.0,0.0),
             EVec(0.25,0.25,0.25),
             EVec(0.782131923303186549e-1,0.782131923303186549e-1,0.782131923303186549e-1),
             EVec(0.782131923303186549e-1,0.782131923303186549e-1,0.765360423009044044),
             EVec(0.782131923303186549e-1,0.765360423009044044,0.782131923303186549e-1),
             EVec(0.765360423009044044,0.782131923303186549e-1,0.782131923303186549e-1),
             EVec(0.121843216663904411,0.121843216663904411,0.121843216663904411),
             EVec(0.121843216663904411,0.121843216663904411,0.634470350008286765),
             EVec(0.121843216663904411,0.634470350008286765,0.121843216663904411),
             EVec(0.634470350008286765,0.121843216663904411,0.121843216663904411),
             EVec(0.332539164446420554,0.332539164446420554,0.332539164446420554),
             EVec(0.332539164446420554,0.332539164446420554,0.238250666073834549e-2),
             EVec(0.332539164446420554,0.238250666073834549e-2,0.332539164446420554),
             EVec(0.238250666073834549e-2,0.332539164446420554,0.332539164446420554),
             EVec(0.10,0.10,0.20),
             EVec(0.10,0.20,0.10),
             EVec(0.20,0.10,0.10),
             EVec(0.10,0.10,0.60),
             EVec(0.10,0.60,0.10),
             EVec(0.60,0.10,0.10),
             EVec(0.10,0.20,0.60),
             EVec(0.10,0.60,0.20),
             EVec(0.60,0.10,0.20),
             EVec(0.60,0.20,0.10),
             EVec(0.20,0.60,0.10),
             EVec(0.20,0.10,0.60)),
            (EVec(0.25,0.25,0.25),#8th degree
             EVec(0.127470936566639015,0.127470936566639015,0.127470936566639015),
             EVec(0.127470936566639015,0.127470936566639015,0.617587190300082967),
             EVec(0.127470936566639015,0.617587190300082967,0.127470936566639015),
             EVec(0.617587190300082967,0.127470936566639015,0.127470936566639015),
             EVec(0.320788303926322960e-1,0.320788303926322960e-1,0.320788303926322960e-1),
             EVec(0.320788303926322960e-1,0.320788303926322960e-1,0.903763508822103123),
             EVec(0.320788303926322960e-1,0.903763508822103123,0.320788303926322960e-1),
             EVec(0.903763508822103123,0.320788303926322960e-1,0.320788303926322960e-1),
             EVec(0.497770956432810185e-1,0.497770956432810185e-1,0.450222904356718978),
             EVec(0.497770956432810185e-1,0.450222904356718978,0.497770956432810185e-1),
             EVec(0.450222904356718978,0.497770956432810185e-1,0.497770956432810185e-1),
             EVec(0.450222904356718978,0.450222904356718978,0.497770956432810185e-1),
             EVec(0.450222904356718978,0.497770956432810185e-1,0.450222904356718978),
             EVec(0.497770956432810185e-1,0.450222904356718978,0.450222904356718978),
             EVec(0.183730447398549945,0.183730447398549945,0.316269552601450060),
             EVec(0.183730447398549945,0.316269552601450060,0.183730447398549945),
             EVec(0.316269552601450060,0.183730447398549945,0.183730447398549945),
             EVec(0.316269552601450060,0.316269552601450060,0.183730447398549945),
             EVec(0.316269552601450060,0.183730447398549945,0.316269552601450060),
             EVec(0.183730447398549945,0.316269552601450060,0.316269552601450060),
             EVec(0.231901089397150906,0.231901089397150906,0.229177878448171174e-1),
             EVec(0.231901089397150906,0.229177878448171174e-1,0.231901089397150906),
             EVec(0.229177878448171174e-1,0.231901089397150906,0.231901089397150906),
             EVec(0.231901089397150906,0.231901089397150906,0.513280033360881072),
             EVec(0.231901089397150906,0.513280033360881072,0.231901089397150906),
             EVec(0.513280033360881072,0.231901089397150906,0.231901089397150906),
             EVec(0.231901089397150906,0.229177878448171174e-1,0.513280033360881072),
             EVec(0.231901089397150906,0.513280033360881072,0.229177878448171174e-1),
             EVec(0.513280033360881072,0.231901089397150906,0.229177878448171174e-1),
             EVec(0.513280033360881072,0.229177878448171174e-1,0.231901089397150906),
             EVec(0.229177878448171174e-1,0.513280033360881072,0.231901089397150906),
             EVec(0.229177878448171174e-1,0.231901089397150906,0.513280033360881072),
             EVec(0.379700484718286102e-1,0.379700484718286102e-1,0.730313427807538396),
             EVec(0.379700484718286102e-1,0.730313427807538396,0.379700484718286102e-1),
             EVec(0.730313427807538396,0.379700484718286102e-1,0.379700484718286102e-1),
             EVec(0.379700484718286102e-1,0.379700484718286102e-1,0.193746475248804382),
             EVec(0.379700484718286102e-1,0.193746475248804382,0.379700484718286102e-1),
             EVec(0.193746475248804382,0.379700484718286102e-1,0.379700484718286102e-1),
             EVec(0.379700484718286102e-1,0.730313427807538396,0.193746475248804382),
             EVec(0.379700484718286102e-1,0.193746475248804382,0.730313427807538396),
             EVec(0.193746475248804382,0.379700484718286102e-1,0.730313427807538396),
             EVec(0.193746475248804382,0.730313427807538396,0.379700484718286102e-1),
             EVec(0.730313427807538396,0.193746475248804382,0.379700484718286102e-1),
             EVec(0.730313427807538396,0.379700484718286102e-1,0.193746475248804382)))
#             ( EVec(0.333333333333333333,0.333333333333333333,0.333333333333333333),
#               EVec(0.0,0.333333333333333333,0.333333333333333333),
#               EVec(0.333333333333333333,0.0,0.333333333333333333),
#               EVec(0.333333333333333333,0.333333333333333333,0.0),
#               EVec(0.25,0.25,0.25),
#               EVec(0.909090909090909091e-1,0.909090909090909091e-1,0.909090909090909091e-1),
#               EVec(0.727272727272727273,0.909090909090909091e-1,0.909090909090909091e-1),
#               EVec(0.909090909090909091e-1,0.727272727272727273,0.909090909090909091e-1),
#               EVec(0.909090909090909091e-1,0.909090909090909091e-1,0.727272727272727273),
#               EVec(0.433449846426335728,0.665501535736642813e-1,0.665501535736642813e-1),
#               EVec(0.665501535736642813e-1,0.433449846426335728,0.665501535736642813e-1),
#               EVec(0.665501535736642813e-1,0.665501535736642813e-1,0.433449846426335728),
#               EVec(0.433449846426335728,0.433449846426335728,0.665501535736642813e-1),
#               EVec(0.433449846426335728,0.665501535736642813e-1,0.433449846426335728),
#               EVec(0.665501535736642813e-1,0.433449846426335728,0.433449846426335728)),
#             ( EVec(0.25,0.25,0.25),
#               EVec(a1,a1,a1),
#               EVec(1.0-2.0*a1,a1,a1),
#               EVec(a1,1.0-2.0*a1,a1),
#               EVec(a1,a1,1.0-2.0*a1),
#               EVec(a2,a2,a2),
#               EVec(1.0-2.0*a2,a2,a2),
#               EVec(a2,1.0-2.0*a2,a2),
#               EVec(a2,a2,1.0-2.0*a2),
#               EVec(0.5-a,a,a),
#               EVec(a,0.5-a,a),
#               EVec(a,a,0.5-a),
#               EVec(0.5-a,0.5-a,a),
#               EVec(a,0.5-a,0.5-a),
#               EVec(0.5-a,a,0.5-a)))
        self.weightsAll=(
            (1.0/6.0,),
            (1.0/24.0,#2nd degree
             1.0/24.0,
             1.0/24.0,
             1.0/24.0),
            (-4.0/30.0,#3rd degree
             9.0/120.0,
             9.0/120.0,
             9.0/120.0,
             9.0/120.0),
            (0.317460317460317450e-2, #4th degree
             0.317460317460317450e-2,
             0.317460317460317450e-2,
             0.317460317460317450e-2,
             0.317460317460317450e-2,
             0.317460317460317450e-2,
             0.147649707904967828e-1,
             0.147649707904967828e-1,
             0.147649707904967828e-1,
             0.147649707904967828e-1,
             0.221397911142651221e-1,
             0.221397911142651221e-1,
             0.221397911142651221e-1,
             0.221397911142651221e-1),
            (0.602678571428571597e-2,#5th degree
             0.602678571428571597e-2,
             0.602678571428571597e-2,
             0.602678571428571597e-2,
             0.302836780970891856e-1,
             0.116452490860289742e-1,
             0.116452490860289742e-1,
             0.116452490860289742e-1,
             0.116452490860289742e-1,
             0.109491415613864534e-1,
             0.109491415613864534e-1,
             0.109491415613864534e-1,
             0.109491415613864534e-1,
             0.109491415613864534e-1,
             0.109491415613864534e-1),
            (0.665379170969464506e-2,#6th degree
             0.665379170969464506e-2,
             0.665379170969464506e-2,
             0.665379170969464506e-2,
             0.167953517588677620e-2,
             0.167953517588677620e-2,
             0.167953517588677620e-2,
             0.167953517588677620e-2,
             0.922619692394239843e-2,
             0.922619692394239843e-2,
             0.922619692394239843e-2,
             0.922619692394239843e-2,
             0.803571428571428248e-2,
             0.803571428571428248e-2,
             0.803571428571428248e-2,
             0.803571428571428248e-2,
             0.803571428571428248e-2,
             0.803571428571428248e-2,
             0.803571428571428248e-2,
             0.803571428571428248e-2,
             0.803571428571428248e-2,
             0.803571428571428248e-2,
             0.803571428571428248e-2,
             0.803571428571428248e-2),
            (0.970017636684296702e-3,#7th degree
             0.970017636684296702e-3,
             0.970017636684296702e-3,
             0.970017636684296702e-3,
             0.970017636684296702e-3,
             0.970017636684296702e-3,
             0.182642234661087939e-1,
             0.105999415244141609e-1,
             0.105999415244141609e-1,
             0.105999415244141609e-1,
             0.105999415244141609e-1,
             -0.625177401143299494e-1,
             -0.625177401143299494e-1,
             -0.625177401143299494e-1,
             -0.625177401143299494e-1,
             0.489142526307353653e-2,
             0.489142526307353653e-2,
             0.489142526307353653e-2,
             0.489142526307353653e-2,
             0.275573192239850917e-1,
             0.275573192239850917e-1,
             0.275573192239850917e-1,
             0.275573192239850917e-1,
             0.275573192239850917e-1,
             0.275573192239850917e-1,
             0.275573192239850917e-1,
             0.275573192239850917e-1,
             0.275573192239850917e-1,
             0.275573192239850917e-1,
             0.275573192239850917e-1,
             0.275573192239850917e-1),
            (-0.393270066412926145e-1,#8th degree
             0.408131605934270525e-2,
             0.408131605934270525e-2,
             0.408131605934270525e-2,
             0.408131605934270525e-2,
             0.658086773304341943e-3,
             0.658086773304341943e-3,
             0.658086773304341943e-3,
             0.658086773304341943e-3,
             0.438425882512284693e-2,
             0.438425882512284693e-2,
             0.438425882512284693e-2,
             0.438425882512284693e-2,
             0.438425882512284693e-2,
             0.438425882512284693e-2,
             0.138300638425098166e-1,
             0.138300638425098166e-1,
             0.138300638425098166e-1,
             0.138300638425098166e-1,
             0.138300638425098166e-1,
             0.138300638425098166e-1,
             0.424043742468372453e-2,
             0.424043742468372453e-2,
             0.424043742468372453e-2,
             0.424043742468372453e-2,
             0.424043742468372453e-2,
             0.424043742468372453e-2,
             0.424043742468372453e-2,
             0.424043742468372453e-2,
             0.424043742468372453e-2,
             0.424043742468372453e-2,
             0.424043742468372453e-2,
             0.424043742468372453e-2,
             0.223873973961420164e-2,
             0.223873973961420164e-2,
             0.223873973961420164e-2,
             0.223873973961420164e-2,
             0.223873973961420164e-2,
             0.223873973961420164e-2,
             0.223873973961420164e-2,
             0.223873973961420164e-2,
             0.223873973961420164e-2,
             0.223873973961420164e-2,
             0.223873973961420164e-2,
             0.223873973961420164e-2))
        #             (0.602678571428571597e-2,
#              0.602678571428571597e-2,
#              0.602678571428571597e-2,
#              0.602678571428571597e-2,
#              0.302836780970891856e-1,
#              0.116452490860289742e-1,
#              0.116452490860289742e-1,
#              0.116452490860289742e-1,
#              0.116452490860289742e-1,
#              0.109491415613864534e-1,
#              0.109491415613864534e-1,
#              0.109491415613864534e-1,
#              0.109491415613864534e-1,
#              0.109491415613864534e-1,
#              0.109491415613864534e-1),
#             (16.0/135.*Vf,
#              w1,w1,w1,w1,
#              w2,w2,w2,w2,
#              wa,wa,wa,wa,wa,wa))
        self.setOrder(order)

class LobattoTetrahedron(Q_base):
    """
    Gauss-Lobatto quadrature on the unit tetrahedron.
    """
    def __init__(self,order=1):
        order=1
        Q_base.__init__(self,order)
        self.pointsAll=(
            ( EVec(0.0,0.0,0.0),
              EVec(1.0,0.0,0.0),
              EVec(0.0,1.0,0.0),
              EVec(0.0,0.0,1.0)),)
        self.weightsAll=(
            (1.0/24.0,1.0/24.0,1.0/24.0,1.0/24.0),)
        self.setOrder(order)
    def setOrder(self,k):
        self.order = 1
        self.points = self.pointsAll[0]
        self.weights = self.weightsAll[0]

class FaceBarycenterTetrahedron(Q_base):
    def __init__(self,order=1):
        order=1
        Q_base.__init__(self,order)
        self.pointsAll=(
            ( EVec(1.0/3.0,1.0/3.0,1.0/3.0), EVec(0.0,1.0/3.0,1./3.), EVec(1.0/3.0,0.0,1.0/3.0),EVec(1./3.,1./3.,0.0)),
            )
        self.weightsAll=(
            (1.0/24.0, 1.0/24.0, 1.0/24.0),
            )
        self.setOrder(order)
    def setOrder(self,k):
        self.order = 1
        self.points = self.pointsAll[0]
        self.weights = self.weightsAll[0]

class SimplexGaussQuadrature(Q_base):
    """
    A class for all quadrature on unit simplices.
    """
    def __init__(self,nd=3,order=1):
        Q_base.__init__(self,order)
        if nd == 0:
            self.quadrature = GaussPoint(order)
        if nd == 1:
            self.quadrature = GaussEdge(order)
        elif nd == 2:
            self.quadrature = GaussTriangle(order)
        elif nd == 3:
            self.quadrature = GaussTetrahedron(order)
        self.pointsAll = self.quadrature.pointsAll
        self.weightsAll = self.quadrature.weightsAll
        self.points = self.quadrature.points
        self.weights = self.quadrature.weights
    def setOrder(self,k):
        self.quadrature.setOrder(k)
        self.points = self.quadrature.points
        self.weights = self.quadrature.weights


class CubeGaussQuadrature(Q_base):
    """
    A class for all quadrature on unit simplices.
    """
    def __init__(self,nd=3,order=1):
        Q_base.__init__(self,order)
        self.nd=nd
        self.quadrature = GaussEdge(order=order)
        self.setOrder(order)


    def setOrder(self,order):
        self.quadrature.setOrder(order,[-1.0,1.0])

        if self.nd == 1:
            self.points = self.quadrature.points
            self.weights = self.quadrature.weights
        if self.nd == 2:
            self.points = []
            self.weights = []
            for i in range(order):
                for j in range(order):

                    self.points.append(EVec(self.quadrature.points[i][0],self.quadrature.points[j][0]))
                    self.weights.append(self.quadrature.weights[i]*self.quadrature.weights[j])

        if self.nd == 3:
            self.points =[]
            self.weights = []
            for i in range(order):
                for j in range(order):
                    for k in range(order):
                        self.points.append(EVec(self.quadrature.points[i][0],self.quadrature.points[j][0],self.quadrature.points[k][0]))
                        self.weights.append(self.quadrature.weights[i]*self.quadrature.weights[j]*self.quadrature.weights[k])


class SimplexLobattoQuadrature(Q_base):
    """
    A class for quadrature on unit simplices.
    """
    def __init__(self,nd=3,order=1):
        #mwf allow for higher order now
        #order=1
        Q_base.__init__(self,order)
        if nd == 0:
            self.quadrature = LobattoPoint(order)
        if nd == 1:
            self.quadrature = LobattoEdge(order)
        elif nd == 2:
            self.quadrature = LobattoTriangle(order)
        elif nd == 3:
            self.quadrature = LobattoTetrahedron(order)
        self.pointsAll = self.quadrature.pointsAll
        self.weightsAll = self.quadrature.weightsAll
        self.points = self.quadrature.points
        self.weights = self.quadrature.weights
    def setOrder(self,k):
        #mwf allow variable order now
        #order=1
        #self.quadrature.points = self.quadrature.pointsAll[0]
        #self.quadrature.weights = self.quadrature.weightsAll[0]
        self.quadrature.setOrder(k)
        self.points = self.quadrature.points
        self.weights = self.quadrature.weights

def buildUnion(quadratureDict):
        #The quadratureDict argument allows different quadrature
        #for a set of integrals (the keys). The goal of the following function is
        #to build a single array of points and an array of weights
        #for each integral (a dictionary of arrays) that matches these points--zero weight
        #if the given point isn't part of the quadrature rule
        #for that integral.
        #
        #
        #First calculate the union of all element quadrature points.
        #
    quadraturePointSet = set()
    for I,quadrature in quadratureDict.iteritems():
        quadraturePointSet |= set([(p[0],p[1],p[2]) for p in quadrature.points])
    nQuadraturePoints = len(quadraturePointSet)
    #
    #Now build a dictionary at each element quadrature point which
    #contains the weights for each integral
    #
    # e.g. quadratureWeightDict[(I,p)] is the weight at p for the
    # integral I
    #
    quadratureWeightDict={}
    #mwf check to avoid float comparison
    quadraturePointValid= {}
    for I,quadrature in quadratureDict.iteritems():
        for p in quadraturePointSet:
            quadratureWeightDict[(I,p)]=0.0
            quadraturePointValid[(I,p)]=False
        for w,p in zip(quadrature.weights,
                       quadrature.points):
            quadratureWeightDict[(I,(p[0],p[1],p[2]))]=w
            quadraturePointValid[(I,(p[0],p[1],p[2]))]=True
    #
    # Now create the desired point and weight arrays
    #
    quadraturePoints = numpy.zeros((nQuadraturePoints,3),
                                          'd')
    for k,p in enumerate(quadraturePointSet):
        quadraturePoints[k,:]=p
    quadratureWeights = {}
    #mwf add dictionary to get indeces for points corresponding to each integral type
    quadraturePointIndeces = {}
    for I in quadratureDict.keys():
        quadratureWeights[I] = numpy.zeros(
            (nQuadraturePoints,),'d')
        quadraturePointIndeces[I]= []
        for k,p in enumerate(quadraturePointSet):
            quadratureWeights[I][k] = quadratureWeightDict[(I,p)]
            if quadraturePointValid[(I,p)]:
                    #mwf is this good enough to tell me what the correct indeces are?
                assert abs(quadratureWeightDict[(I,p)]) > 1.0e-10, "valid quadrature point zero weight"
                quadraturePointIndeces[I].append(k)
    return (quadraturePoints,quadratureWeights,quadraturePointIndeces)
## @}

if __name__ == '__main__':
    from math import *
    import numpy
    from LinearAlgebraTools import *
    #define some simple functions to integrate
    a=1.1
    b=0.92
    c=1.34
    def f0(x):
        return map(lambda y: 1.0, x)
    def f1(x):
        return map(lambda y: 1.0 + a*y[X] + b*y[Y] + c*y[Z], x)
    def f2(x):
        return map(lambda y: 1.0 + a*y[X]**2 + c*y[Y]**2 + b*y[Z]**2, x)
    def f3(x):
        return map(lambda y: 1.0 + b*y[X]**3 + a*y[Y]**3 + c*y[Z]**3, x)
    def f4(x):
        return map(lambda y: 1.0 + c*y[X]**4 + b*y[Y]**4 + a*y[Z]**4, x)

    gaussEdge=GaussEdge()
    gaussTriangle=GaussTriangle()
    gaussTetrahedron=GaussTetrahedron()

    print "4th Order Polynomial"
    print "Edge"
    gaussEdge.setOrder(1)
    int0_f4 = dot(f4(gaussEdge.points),gaussEdge.weights)
    print int0_f4
    gaussEdge.setOrder(2)
    int1_f4 = dot(f4(gaussEdge.points),gaussEdge.weights)
    print int1_f4
    gaussEdge.setOrder(3)
    int2_f4 = dot(f4(gaussEdge.points),gaussEdge.weights)
    print int2_f4

    print "Triangle"
    gaussTriangle.setOrder(1)
    int0_f4 = dot(f4(gaussTriangle.points),gaussTriangle.weights)
    print int0_f4
    gaussTriangle.setOrder(2)
    int1_f4 = dot(f4(gaussTriangle.points),gaussTriangle.weights)
    print int1_f4
    gaussTriangle.setOrder(3)
    int2_f4 = dot(f4(gaussTriangle.points),gaussTriangle.weights)
    print int2_f4

    print "Tetrahedron"
    gaussTetrahedron.setOrder(1)
    int0_f4 = dot(f4(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int0_f4
    gaussTetrahedron.setOrder(2)
    int1_f4 = dot(f4(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int1_f4
    gaussTetrahedron.setOrder(3)
    int2_f4 = dot(f4(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int2_f4

    print "3rd Order Polynomial"
    print "Edge"
    gaussEdge.setOrder(1)
    int0_f3 = dot(f3(gaussEdge.points),gaussEdge.weights)
    print int0_f3
    gaussEdge.setOrder(2)
    int1_f3 = dot(f3(gaussEdge.points),gaussEdge.weights)
    print int1_f3
    gaussEdge.setOrder(3)
    int2_f3 = dot(f3(gaussEdge.points),gaussEdge.weights)
    print int2_f3

    print "Triangle"
    gaussTriangle.setOrder(1)
    int0_f3 = dot(f3(gaussTriangle.points),gaussTriangle.weights)
    print int0_f3
    gaussTriangle.setOrder(2)
    int1_f3 = dot(f3(gaussTriangle.points),gaussTriangle.weights)
    print int1_f3
    gaussTriangle.setOrder(3)
    int2_f3 = dot(f3(gaussTriangle.points),gaussTriangle.weights)
    print int2_f3

    print "Tetrahedron"
    gaussTetrahedron.setOrder(1)
    int0_f3 = dot(f3(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int0_f3
    gaussTetrahedron.setOrder(2)
    int1_f3 = dot(f3(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int1_f3
    gaussTetrahedron.setOrder(3)
    int2_f3 = dot(f3(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int2_f3

    print "2nd Order Polynomial"
    print "Edge"
    gaussEdge.setOrder(1)
    int0_f2 = dot(f2(gaussEdge.points),gaussEdge.weights)
    print int0_f2
    gaussEdge.setOrder(2)
    int1_f2 = dot(f2(gaussEdge.points),gaussEdge.weights)
    print int1_f2
    gaussEdge.setOrder(3)
    int2_f2 = dot(f2(gaussEdge.points),gaussEdge.weights)
    print int2_f2

    print "Triangle"
    gaussTriangle.setOrder(1)
    int0_f2 = dot(f2(gaussTriangle.points),gaussTriangle.weights)
    print int0_f2
    gaussTriangle.setOrder(2)
    int1_f2 = dot(f2(gaussTriangle.points),gaussTriangle.weights)
    print int1_f2
    gaussTriangle.setOrder(3)
    int2_f2 = dot(f2(gaussTriangle.points),gaussTriangle.weights)
    print int2_f2

    print "Tetrahedron"
    gaussTetrahedron.setOrder(1)
    int0_f2 = dot(f2(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int0_f2
    gaussTetrahedron.setOrder(2)
    int1_f2 = dot(f2(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int1_f2
    gaussTetrahedron.setOrder(3)
    int2_f2 = dot(f2(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int2_f2

    print "1st Order Polynomial"
    print "Edge"
    gaussEdge.setOrder(1)
    int0_f1 = dot(f1(gaussEdge.points),gaussEdge.weights)
    print int0_f1
    gaussEdge.setOrder(2)
    int1_f1 = dot(f1(gaussEdge.points),gaussEdge.weights)
    print int1_f1
    gaussEdge.setOrder(3)
    int1_f1 = dot(f1(gaussEdge.points),gaussEdge.weights)
    print int1_f1

    print "Triangle"
    gaussTriangle.setOrder(1)
    int0_f1 = dot(f1(gaussTriangle.points),gaussTriangle.weights)
    print int0_f1
    gaussTriangle.setOrder(2)
    int1_f1 = dot(f1(gaussTriangle.points),gaussTriangle.weights)
    print int1_f1
    gaussTriangle.setOrder(3)
    int1_f1 = dot(f1(gaussTriangle.points),gaussTriangle.weights)
    print int1_f1

    print "Tetrahedron"
    gaussTetrahedron.setOrder(1)
    int0_f1 = dot(f1(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int0_f1
    gaussTetrahedron.setOrder(2)
    int1_f1 = dot(f1(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int1_f1
    gaussTetrahedron.setOrder(3)
    int2_f1 = dot(f1(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int2_f1

    print "0th Order Polynomial"
    print "Edge"
    gaussEdge.setOrder(1)
    int0_f0 = dot(f0(gaussEdge.points),gaussEdge.weights)
    print int0_f0
    gaussEdge.setOrder(2)
    int1_f0 = dot(f0(gaussEdge.points),gaussEdge.weights)
    print int1_f0
    gaussEdge.setOrder(3)
    int2_f0 = dot(f0(gaussEdge.points),gaussEdge.weights)
    print int2_f0

    print "Triangle"
    gaussTriangle.setOrder(1)
    int0_f0 = dot(f0(gaussTriangle.points),gaussTriangle.weights)
    print int0_f0
    gaussTriangle.setOrder(2)
    int1_f0 = dot(f0(gaussTriangle.points),gaussTriangle.weights)
    print int1_f0
    gaussTriangle.setOrder(3)
    int2_f0 = dot(f0(gaussTriangle.points),gaussTriangle.weights)
    print int2_f0

    print "Tetrahedron"
    gaussTetrahedron.setOrder(1)
    int0_f0 = dot(f0(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int0_f0
    gaussTetrahedron.setOrder(2)
    int1_f0 = dot(f0(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int1_f0
    gaussTetrahedron.setOrder(3)
    int2_f0 = dot(f0(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int2_f0
