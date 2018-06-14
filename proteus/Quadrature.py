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
        #mql. Compute points and weigths for Gauss Quad with n=5,...,9
        #mql. TODO: we should generalize this to any order.
        #mql. Why are the points not in order? 
        [p6,w6]=numpy.polynomial.legendre.leggauss(6); p6=0.5*(p6+1); w6=0.5*w6
        [p7,w7]=numpy.polynomial.legendre.leggauss(7); p7=0.5*(p7+1); w7=0.5*w7
        [p8,w8]=numpy.polynomial.legendre.leggauss(8); p8=0.5*(p8+1); w8=0.5*w8
        [p9,w9]=numpy.polynomial.legendre.leggauss(9); p9=0.5*(p9+1); w9=0.5*w9
        [p10,w10]=numpy.polynomial.legendre.leggauss(10); p10=0.5*(p10+1); w10=0.5*w10
        [p11,w11]=numpy.polynomial.legendre.leggauss(11); p11=0.5*(p11+1); w11=0.5*w11
        [p12,w12]=numpy.polynomial.legendre.leggauss(12); p12=0.5*(p12+1); w12=0.5*w12
        [p13,w13]=numpy.polynomial.legendre.leggauss(13); p13=0.5*(p13+1); w13=0.5*w13
        [p14,w14]=numpy.polynomial.legendre.leggauss(14); p14=0.5*(p14+1); w14=0.5*w14
        [p15,w15]=numpy.polynomial.legendre.leggauss(15); p15=0.5*(p15+1); w15=0.5*w15
        [p16,w16]=numpy.polynomial.legendre.leggauss(16); p16=0.5*(p16+1); w16=0.5*w16
        [p17,w17]=numpy.polynomial.legendre.leggauss(17); p17=0.5*(p17+1); w17=0.5*w17
        [p18,w18]=numpy.polynomial.legendre.leggauss(18); p18=0.5*(p18+1); w18=0.5*w18
        [p19,w19]=numpy.polynomial.legendre.leggauss(19); p19=0.5*(p19+1); w19=0.5*w19
        [p20,w20]=numpy.polynomial.legendre.leggauss(20); p20=0.5*(p20+1); w20=0.5*w20       
        #mwf for convenience, see Ern and Guermond
        a1 = 0.5*sqrt((15.0+2.0*sqrt(30))/35.0)
        a2 = 0.5*sqrt((15.0-2.0*sqrt(30))/35.0)
        w1 = 0.25 - sqrt(5./6.)/12.
        w2 = 0.25 + sqrt(5./6.)/12.
        Q_base.__init__(self,order)
        self.pointsAll=(
            (EVec(0.5),), #n=1
            (EVec((sqrt(3.0)-1.0)/(2.0*sqrt(3.0))), #n=2
             EVec((sqrt(3.0)+1.0)/(2.0*sqrt(3.0)))),
            (EVec((sqrt(5.0) - sqrt(3.0))/(2.0*sqrt(5))), #n=3
             EVec(0.5),
             EVec((sqrt(5.0) + sqrt(3.0))/(2.0*sqrt(5)))),
            (EVec(0.5+a1),EVec(0.5-a1), #n=4
             EVec(0.5+a2),EVec(0.5-a2)),
            (EVec(0.5), #n=5
             EVec(0.5*(sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0) + 0.5),
             EVec(0.5*(-sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0) + 0.5),
             EVec(0.5*(sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0) + 0.5),
             EVec(0.5*(-sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0) + 0.5)),
            (EVec(p6[0]), #n=6
             EVec(p6[1]),
             EVec(p6[2]),
             EVec(p6[3]),
             EVec(p6[4]),
             EVec(p6[5])),            
            (EVec(p7[0]), #n=7
             EVec(p7[1]),
             EVec(p7[2]),
             EVec(p7[3]),
             EVec(p7[4]),
             EVec(p7[5]),
             EVec(p7[6])),            
            (EVec(p8[0]), #n=8
             EVec(p8[1]),
             EVec(p8[2]),
             EVec(p8[3]),
             EVec(p8[4]),
             EVec(p8[5]),
             EVec(p8[6]),             
             EVec(p8[7])),
            (EVec(p9[0]), #n=9
             EVec(p9[1]),
             EVec(p9[2]),
             EVec(p9[3]),
             EVec(p9[4]),
             EVec(p9[5]),
             EVec(p9[6]),
             EVec(p9[7]),                          
             EVec(p9[8])),
            (EVec(p10[0]), #n=10
             EVec(p10[1]),
             EVec(p10[2]),
             EVec(p10[3]),
             EVec(p10[4]),
             EVec(p10[5]),
             EVec(p10[6]),
             EVec(p10[7]),                          
             EVec(p10[8]),
             EVec(p10[9])),
            (EVec(p11[0]), #n=11
             EVec(p11[1]),
             EVec(p11[2]),
             EVec(p11[3]),
             EVec(p11[4]),
             EVec(p11[5]),
             EVec(p11[6]),
             EVec(p11[7]),                          
             EVec(p11[8]),
             EVec(p11[9]),
             EVec(p11[10])),
            (EVec(p12[0]), #n=12
             EVec(p12[1]),
             EVec(p12[2]),
             EVec(p12[3]),
             EVec(p12[4]),
             EVec(p12[5]),
             EVec(p12[6]),
             EVec(p12[7]),                          
             EVec(p12[8]),
             EVec(p12[9]),
             EVec(p12[10]),
             EVec(p12[11])),
            (EVec(p13[0]), #n=13
             EVec(p13[1]),
             EVec(p13[2]),
             EVec(p13[3]),
             EVec(p13[4]),
             EVec(p13[5]),
             EVec(p13[6]),
             EVec(p13[7]),                          
             EVec(p13[8]),
             EVec(p13[9]),
             EVec(p13[10]),
             EVec(p13[11]),
             EVec(p13[12])),            
            (EVec(p14[0]), #n=14
             EVec(p14[1]),
             EVec(p14[2]),
             EVec(p14[3]),
             EVec(p14[4]),
             EVec(p14[5]),
             EVec(p14[6]),
             EVec(p14[7]),                          
             EVec(p14[8]),
             EVec(p14[9]),
             EVec(p14[10]),
             EVec(p14[11]),
             EVec(p14[12]),
             EVec(p14[13])),
            (EVec(p15[0]), #n=15
             EVec(p15[1]),
             EVec(p15[2]),
             EVec(p15[3]),
             EVec(p15[4]),
             EVec(p15[5]),
             EVec(p15[6]),
             EVec(p15[7]),                          
             EVec(p15[8]),
             EVec(p15[9]),
             EVec(p15[10]),
             EVec(p15[11]),
             EVec(p15[12]),
             EVec(p15[13]),             
             EVec(p15[14])),            
            (EVec(p16[0]), #n=16
             EVec(p16[1]),
             EVec(p16[2]),
             EVec(p16[3]),
             EVec(p16[4]),
             EVec(p16[5]),
             EVec(p16[6]),
             EVec(p16[7]),                          
             EVec(p16[8]),
             EVec(p16[9]),
             EVec(p16[10]),
             EVec(p16[11]),
             EVec(p16[12]),
             EVec(p16[13]),
             EVec(p16[14]),
             EVec(p16[15])),
            (EVec(p17[0]), #n=17
             EVec(p17[1]),
             EVec(p17[2]),
             EVec(p17[3]),
             EVec(p17[4]),
             EVec(p17[5]),
             EVec(p17[6]),
             EVec(p17[7]),                          
             EVec(p17[8]),
             EVec(p17[9]),
             EVec(p17[10]),
             EVec(p17[11]),
             EVec(p17[12]),
             EVec(p17[13]),
             EVec(p17[14]),
             EVec(p17[15]),
             EVec(p17[16])),
            (EVec(p18[0]), #n=18
             EVec(p18[1]),
             EVec(p18[2]),
             EVec(p18[3]),
             EVec(p18[4]),
             EVec(p18[5]),
             EVec(p18[6]),
             EVec(p18[7]),                          
             EVec(p18[8]),
             EVec(p18[9]),
             EVec(p18[10]),
             EVec(p18[11]),
             EVec(p18[12]),
             EVec(p18[13]),
             EVec(p18[14]),
             EVec(p18[15]),
             EVec(p18[16]),
             EVec(p18[17])),
            (EVec(p19[0]), #n=19
             EVec(p19[1]),
             EVec(p19[2]),
             EVec(p19[3]),
             EVec(p19[4]),
             EVec(p19[5]),
             EVec(p19[6]),
             EVec(p19[7]),                          
             EVec(p19[8]),
             EVec(p19[9]),
             EVec(p19[10]),
             EVec(p19[11]),
             EVec(p19[12]),
             EVec(p19[13]),
             EVec(p19[14]),
             EVec(p19[15]),
             EVec(p19[16]),
             EVec(p19[17]),
             EVec(p19[18])),
            (EVec(p20[0]), #n=20
             EVec(p20[1]),
             EVec(p20[2]),
             EVec(p20[3]),
             EVec(p20[4]),
             EVec(p20[5]),
             EVec(p20[6]),
             EVec(p20[7]),                          
             EVec(p20[8]),
             EVec(p20[9]),
             EVec(p20[10]),
             EVec(p20[11]),
             EVec(p20[12]),
             EVec(p20[13]),
             EVec(p20[14]),
             EVec(p20[15]),
             EVec(p20[16]),
             EVec(p20[17]),
             EVec(p20[18]),
             EVec(p20[19]))            
            )
        self.weightsAll=(
            (1.0,), #n=1
            (0.5, #n=2
             0.5),
            (5.0/18.0, #n=3
             8.0/18.0,
             5.0/18.0),
            (w1,w1,w2,w2), #n=4
            (0.5*(128.0/225.0), #n=5
             0.5*(322.0+13.0*sqrt(70.0))/900.0,
             0.5*(322.0+13.0*sqrt(70.0))/900.0,
             0.5*(322.0-13.0*sqrt(70.0))/900.0,
             0.5*(322.0-13.0*sqrt(70.0))/900.0),
            (w6[0],w6[1],w6[2],w6[3],w6[4],w6[5]), #n=6
            (w7[0],w7[1],w7[2],w7[3],w7[4],w7[5],w7[6]), #n=7
            (w8[0],w8[1],w8[2],w8[3],w8[4],w8[5],w8[6],w8[7]), #n=8
            (w9[0],w9[1],w9[2],w9[3],w9[4],w9[5],w9[6],w9[7],w9[8]), #n=9
            (w10[0],w10[1],w10[2],w10[3],w10[4],w10[5],w10[6],w10[7],w10[8],w10[9]), #n=10
            (w11[0],w11[1],w11[2],w11[3],w11[4],w11[5],w11[6],w11[7],w11[8],w11[9],w11[10]), #n=11
            (w12[0],w12[1],w12[2],w12[3],w12[4],w12[5],w12[6],w12[7],w12[8],w12[9],w12[10],w12[11]), #n=12
            (w13[0],w13[1],w13[2],w13[3],w13[4],w13[5],w13[6],w13[7],w13[8],w13[9],w13[10],w13[11],w13[12]), #n=13
            (w14[0],w14[1],w14[2],w14[3],w14[4],w14[5],w14[6],w14[7],w14[8],w14[9],w14[10],w14[11],w14[12],w14[13]), #n=14
            (w15[0],w15[1],w15[2],w15[3],w15[4],w15[5],w15[6],w15[7],w15[8],w15[9],w15[10],w15[11],w15[12],w15[13],w15[14]), #n=15
            (w16[0],w16[1],w16[2],w16[3],w16[4],w16[5],w16[6],w16[7],w16[8],w16[9],w16[10],w16[11],w16[12],w16[13],w16[14],w16[15]), #n=16
            (w17[0],w17[1],w17[2],w17[3],w17[4],w17[5],w17[6],w17[7],w17[8],w17[9],w17[10],w17[11],w17[12],w17[13],w17[14],w17[15],w17[16]), #n=17
            (w18[0],w18[1],w18[2],w18[3],w18[4],w18[5],w18[6],w18[7],w18[8],w18[9],w18[10],w18[11],w18[12],w18[13],w18[14],w18[15],w18[16],w18[17]), #n=18
            (w19[0],w19[1],w19[2],w19[3],w19[4],w19[5],w19[6],w19[7],w19[8],w19[9],w19[10],w19[11],w19[12],w19[13],w19[14],w19[15],w19[16],w19[17],w19[18]), #n=19
            (w20[0],w20[1],w20[2],w20[3],w20[4],w20[5],w20[6],w20[7],w20[8],w20[9],w20[10],w20[11],w20[12],w20[13],w20[14],w20[15],w20[16],w20[17],w20[18],w20[19]) #n=20
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

########################
class GaussCompositeEdge(Q_base):
    """
    Gaussian Quadrature on the unit interval.
    """
    def __init__(self,order=1):
        [p1,w1]=numpy.polynomial.legendre.leggauss(1); #p1=0.5*(p1+1); w1=0.5*w6
        [p2,w2]=numpy.polynomial.legendre.leggauss(2); #p2=0.5*(p2+1); w2=0.5*w6
        [p3,w3]=numpy.polynomial.legendre.leggauss(3); #p3=0.5*(p3+1); w3=0.5*w6
        [p4,w4]=numpy.polynomial.legendre.leggauss(4); #p4=0.5*(p4+1); w4=0.5*w6
        [p5,w5]=numpy.polynomial.legendre.leggauss(5); #p5=0.5*(p5+1); w5=0.5*w6        
        [p6,w6]=numpy.polynomial.legendre.leggauss(6); #p6=0.5*(p6+1); w6=0.5*w6
        Q_base.__init__(self,order)
        self.pointsAll=(
            (EVec(p1[0])), #n=1
            (EVec(p2[0]),  #n=2
             EVec(p2[1])),
            (EVec(p3[0]),  #n=3
             EVec(p3[1]),
             EVec(p3[2])),
            (EVec(p4[0]),  #n=4
             EVec(p4[1]),
             EVec(p4[2]),             
             EVec(p4[3])),
            (EVec(p5[0]),  #n=5
             EVec(p5[1]),
             EVec(p5[2]),
             EVec(p5[3]),
             EVec(p5[4])),                        
            (EVec(p6[0]),  #n=6
             EVec(p6[1]),
             EVec(p6[2]),
             EVec(p6[3]),
             EVec(p6[4]),
             EVec(p6[5]))
            )
        self.weightsAll=(
            (w1[0]), #n=1
            (w2[0],w2[1]), #n=2
            (w3[0],w3[1],w3[2]), #n=3
            (w4[0],w4[1],w4[2],w4[3]), #n=4
            (w5[0],w5[1],w5[2],w5[3],w5[4]), #n=5
            (w6[0],w6[1],w6[2],w6[3],w6[4],w6[5]) #n=6
        )
        self.setOrder(order)

    def setOrder(self,order,domain=[0.0,1.0]):
        Q_base.setOrder(self,order)
        points = self.points
        weights = self.weights
        self.points = []
        self.weights = []
        for i in range(order):
            #t = 0.5*(x + 1)*(b - a) + a
            #gauss = sum(w * f(t)) * 0.5*(b - a)
            self.points.append(EVec(0.5*(domain[1]-domain[0])*(points[i][0]+1.) +domain[0]))
            self.weights.append(0.5*(domain[1]-domain[0])*weights[i])
########################

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


class CompositeTriangle(Q_base):
    """
    Composite quadrature rule on the unit triangle.
    """

    # uniform refine the reference cell until size < hk.
    def __init__(self,quad_rule,hk):

        N = int(floor(1.0/hk))
        h1= 1.0/N
        h2= 1.0/N/N
        npt = len(quad_rule.points)

        self.weights = numpy.tile(quad_rule.weights, N*N)
        self.weights *= h2

        quad_points = numpy.asarray(quad_rule.points)

        self.points = numpy.zeros((npt*N*N,3),'d')
        k = 1
        for j in range(N):
            for i in range(N-j):
                 self.points[(k-1)*npt:k*npt, :] = quad_points*h1 + numpy.array([i*h1,j*h1,0.0])
                 k += 1
        for j in range(1,N):
            for i in range(1,N-j+1):
                self.points[(k-1)*npt:k*npt, :] = -quad_points*h1 + numpy.array([i*h1,j*h1, 0.0])
                k += 1

class CompositeTetrahedron(Q_base):

    def get_detJ_and_J_from_ref(self, simplex_nodes):
        x = simplex_nodes[:, 0]
        y = simplex_nodes[:, 1]
        z = simplex_nodes[:, 2]

        J = numpy.array([[x[0] - x[2], x[1] - x[2], x[3] - x[2]],
                      [y[0] - y[2], y[1] - y[2], y[3] - y[2]],
                      [z[0] - z[2], z[1] - z[2], z[3] - z[2]]])

        detJ = numpy.linalg.det(J)

        B = numpy.array([x[2], y[2], z[2]])
        return detJ, J, B

    def get_h_of_Tetrahedron(self, simplex_nodes):
        h1 = numpy.linalg.norm(simplex_nodes[0, :] - simplex_nodes[2, :])
        h2 = numpy.linalg.norm(simplex_nodes[1, :] - simplex_nodes[2, :])
        h3 = numpy.linalg.norm(simplex_nodes[3, :] - simplex_nodes[2, :])
        h4 = numpy.linalg.norm(simplex_nodes[0, :] - simplex_nodes[1, :])
        h5 = numpy.linalg.norm(simplex_nodes[0, :] - simplex_nodes[3, :])
        h6 = numpy.linalg.norm(simplex_nodes[3, :] - simplex_nodes[2, :])
        return max([h1, h2, h3, h4, h5, h6])

    def get_max_h_of_all_tetrahedron(self, all_tet):
        n_tet = all_tet.shape[0] / 4
        h = []
        for i in range(n_tet):
            h.append(self.get_h_of_Tetrahedron(all_tet[i * 4:(i + 1) * 4, :]))

        return max(h)


    def get_8_sub_simplex(self, simplex_nodes):
        r""" Return 8 sub-simplex of the simplex givne by simplex_nodes: 4 points ordered by right hand rule.
        `<https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/22430/versions/7/previews/tetrarefine3.m/index.html?access_key=>`_
        """
        x1 = simplex_nodes[0]
        x2 = simplex_nodes[1]
        x3 = simplex_nodes[2]
        x4 = simplex_nodes[3]
        x5 = 0.5 * (x1 + x2)
        x6 = 0.5 * (x1 + x3)
        x7 = 0.5 * (x1 + x4)
        x8 = 0.5 * (x2 + x3)
        x9 = 0.5 * (x2 + x4)
        x10 = 0.5 * (x3 + x4)

        all_sub_simplex = numpy.zeros((8 * 4, 3), 'd')

        all_sub_simplex[0] = x1
        all_sub_simplex[1] = x5
        all_sub_simplex[2] = x6
        all_sub_simplex[3] = x7

        all_sub_simplex[4] = x5
        all_sub_simplex[5] = x2
        all_sub_simplex[6] = x8
        all_sub_simplex[7] = x9

        all_sub_simplex[8] = x6
        all_sub_simplex[9] = x8
        all_sub_simplex[10] = x3
        all_sub_simplex[11] = x10

        all_sub_simplex[12] = x7
        all_sub_simplex[13] = x9
        all_sub_simplex[14] = x10
        all_sub_simplex[15] = x4

        all_sub_simplex[16] = x5
        all_sub_simplex[17] = x6
        all_sub_simplex[18] = x7
        all_sub_simplex[19] = x9

        all_sub_simplex[20] = x8
        all_sub_simplex[21] = x6
        all_sub_simplex[22] = x5
        all_sub_simplex[23] = x9

        all_sub_simplex[24] = x6
        all_sub_simplex[25] = x7
        all_sub_simplex[26] = x9
        all_sub_simplex[27] = x10

        all_sub_simplex[28] = x9
        all_sub_simplex[29] = x8
        all_sub_simplex[30] = x6
        all_sub_simplex[31] = x10

        return all_sub_simplex


    def get_sub_tet_of_all_tet(self, all_tet):
        n_tet = all_tet.shape[0] / 4
        all_sub_tet = numpy.zeros((n_tet * 8 * 4, 3), 'd')
        for i in range(n_tet):
            all_sub_tet[i * 8 * 4:(i + 1) * 8 * 4,
                        :] = self.get_8_sub_simplex(all_tet[i * 4:(i + 1) * 4, :])
        return all_sub_tet

    def __init__(self, quad_rule, hk):
        """
        refine the reference tetrahedron until the largest edge < hk
        """

        all_tetrahedron = numpy.array(
            [[1, 0, 0], [0, 1, 0], [0, 0, 0], [0, 0, 1]], 'd')

        max_h = self.get_max_h_of_all_tetrahedron(all_tetrahedron)

        while max_h > hk:
            all_sub_tet = self.get_sub_tet_of_all_tet(all_tetrahedron)
            all_tetrahedron = all_sub_tet
            max_h = self.get_max_h_of_all_tetrahedron(all_tetrahedron)
#             self.plot_all_tet(all_tetrahedron)

        self.h = max_h

        n_tet = all_tetrahedron.shape[0] / 4
        quad_weights = numpy.asarray(quad_rule.weights, 'd')
        quad_points = numpy.asarray(quad_rule.points, 'd')
        n_quad_per_tet = quad_weights.shape[0]
        self.weights = numpy.zeros((n_tet * n_quad_per_tet,), 'd')
        self.points = numpy.zeros((n_tet * n_quad_per_tet, 3), 'd')

        for i in range(n_tet):
            detJ, J, B = self.get_detJ_and_J_from_ref(
                all_tetrahedron[i * 4: (i + 1) * 4, :])

            self.weights[i * n_quad_per_tet:(i + 1)
                         * n_quad_per_tet] = detJ * quad_weights
            self.points[i * n_quad_per_tet:(i + 1)
                        * n_quad_per_tet, :] = numpy.dot(quad_points, J.T) + B

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
                    if j == 0 or j == nintervals-i: jInBoundary = 1
                    if iInBoundary+jInBoundary == 2:
                        weights.append(parentArea/3.0)
                    elif iInBoundary+jInBoundary == 1:
                        weights.append(parentArea)
                    else:
                        weights.append(parentArea*2.0)
            pointsList.append([tuple(p) for p in combos])
            weightsList.append([w for w in weights])
            parentArea = 0.5*(1.0/float(nintervals+1))**2
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
    """ A class which defines quadrature on unit simplices.

    Arguments
    ---------
    nd : int
        Dimension of the finite element problem.
    order : int
        Polynomial order for which the integration is exact.
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


class CubeGaussCompositeQuadrature(Q_base):
    """
    A class for all quadrature on unit simplices.
    """
    def __init__(self,nd=3,order=1):
        Q_base.__init__(self,order)
        self.nd=nd
        self.quadratureNeg = GaussCompositeEdge(order=order)
        self.quadraturePos = GaussCompositeEdge(order=order)
        self.setOrder(order)

    def setOrder(self,order):
        self.points = []
        self.weights = []
        self.quadratureNeg.setOrder(order,[-1.0,0.0])
        self.quadraturePos.setOrder(order,[0.0,1.0])

        if self.nd == 1 or self.nd == 3:
            raise("Not implemented")
        ## Neg x Neg
        for i in range(order):
            for j in range(order):
                self.points.append(EVec(self.quadratureNeg.points[i][0],
                                        self.quadratureNeg.points[j][0]))
                self.weights.append(self.quadratureNeg.weights[i]*self.quadratureNeg.weights[j])
        ## Neg x Pos
        for i in range(order):
            for j in range(order):
                self.points.append(EVec(self.quadratureNeg.points[i][0],
                                        self.quadraturePos.points[j][0]))
                self.weights.append(self.quadratureNeg.weights[i]*self.quadraturePos.weights[j])
        ## Pos x Neg
        for i in range(order):
            for j in range(order):
                self.points.append(EVec(self.quadraturePos.points[i][0],
                                        self.quadratureNeg.points[j][0]))
                self.weights.append(self.quadraturePos.weights[i]*self.quadratureNeg.weights[j])
        ## Pos x Pos
        for i in range(order):
            for j in range(order):
                self.points.append(EVec(self.quadraturePos.points[i][0],
                                        self.quadraturePos.points[j][0]))
                self.weights.append(self.quadraturePos.weights[i]*self.quadraturePos.weights[j])

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


