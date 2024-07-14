from proteus import Domain
#         |y
#         |
#         |
#         +-------+ 
#        /|      /| 
#       +-+----+x1| 
#       | |x0  |  | 
#       | +----+--+---------x
#       |/     | /  
#       +------+   
#      /
#     /
#    / z
#===============================================================================
#  One box
#===============================================================================
def get_domain_one_box(x0=(0.0,0.0,0.0),L=(1.0,1.0,1.0),he=0.001):
    boundaries = ['left', 'right', 'bottom', 'top', 'front', 'back']
    boundaryTags = dict([(key, i + 1) for (i, key) in enumerate(boundaries)])
    x1 = (x0[0]+L[0],x0[1]+L[1],x0[2]+L[2])
    vertices = [(x0[0],x0[1],x0[2]),
                (x1[0],x0[1],x0[2]),
                (x1[0],x1[1],x0[2]),
                (x0[0],x1[1],x0[2]),
                (x0[0],x0[1],x1[2]),
                (x1[0],x0[1],x1[2]),
                (x1[0],x1[1],x1[2]),
                (x0[0],x1[1],x1[2])]
    vertexFlags=[boundaryTags['bottom'],boundaryTags['bottom'],boundaryTags['top'],boundaryTags['top'],
                boundaryTags['bottom'],boundaryTags['bottom'],boundaryTags['top'],boundaryTags['top'],
                ]
    
    planes = []
    planeFlags = []
    planeHolds = []
    planes.extend([[[0,1,2,3],],
                   [[4,5,6,7],],
                   [[0,3,7,4]],
                   [[1,2,6,5]],
                   [[0,1,5,4]],
                   [[2,3,7,6]]])
    planeFlags.extend([boundaryTags['back'],boundaryTags['front'],
                       boundaryTags['left'],boundaryTags['right'],
                       boundaryTags['bottom'],boundaryTags['top'],])
    planeHolds.extend([[],[],[],[],[],[]])
    
    regions = [[0.5*(x0[0]+x1[0]),0.5*(x0[1]+x1[1]),0.5*(x0[2]+x1[2])],
               ]
    regionFlags = [1,]
    regionConstraints=[0.5*he**2]
    
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,facets=planes,facetFlags=planeFlags,facetHoles=planeHolds,
                                       regions=regions,regionFlags=regionFlags,regionConstraints=regionConstraints)
    return domain,boundaryTags,x0,x1



#===============================================================================
# Two boxes
# the 2nd box is inside the 1st box;
# the 2nd box cannot touch the 1st box; 
#===============================================================================
def get_domain_two_box(x0=(0.0,0.0,0.0),L=(1.0,1.0,1.0),he=0.001,he2=None,x2=None,x3=None):
    boundaries = ['left', 'right', 'bottom', 'top', 'front', 'back']
    boundaryTags = dict([(key, i + 1) for (i, key) in enumerate(boundaries)])
    x1 = (x0[0]+L[0],x0[1]+L[1],x0[2]+L[2])
    if x2==None:
        x2 = (0.25*(x0[0]+x1[0]), 0.25*(x0[1]+x1[1]), 0.25*(x0[2]+x1[2]))
    if x3==None:
        x3 = (0.75*(x0[0]+x1[0]), 0.75*(x0[1]+x1[1]), 0.75*(x0[2]+x1[2]))
    if he2==None:
        he2 = he
    vertices = [(x0[0],x0[1],x0[2]),
                (x1[0],x0[1],x0[2]),
                (x1[0],x1[1],x0[2]),
                (x0[0],x1[1],x0[2]),
                (x0[0],x0[1],x1[2]),
                (x1[0],x0[1],x1[2]),
                (x1[0],x1[1],x1[2]),
                (x0[0],x1[1],x1[2]),
                (x2[0],x2[1],x2[2]),
                (x3[0],x2[1],x2[2]),
                (x3[0],x3[1],x2[2]),
                (x2[0],x3[1],x2[2]),
                (x2[0],x2[1],x3[2]),
                (x3[0],x2[1],x3[2]),
                (x3[0],x3[1],x3[2]),
                (x2[0],x3[1],x3[2])]
    vertexFlags=[boundaryTags['bottom'],boundaryTags['bottom'],boundaryTags['top'],boundaryTags['top'],
                boundaryTags['bottom'],boundaryTags['bottom'],boundaryTags['top'],boundaryTags['top'],
                0,0,0,0,
                0,0,0,0]
    
    planes = []
    planeFlags = []
    planeHolds = []
    planes.extend([[[0,1,2,3],],
                   [[4,5,6,7],],
                   [[0,3,7,4]],
                   [[1,2,6,5]],
                   [[0,1,5,4]],
                   [[2,3,7,6]]])
    planeFlags.extend([boundaryTags['back'],boundaryTags['front'],
                       boundaryTags['left'],boundaryTags['right'],
                       boundaryTags['bottom'],boundaryTags['top'],])
    planeHolds.extend([[],[],[],[],[],[]])
    
    planes.extend([[[8,9,10,11]],
                   [[12,13,14,15]],
                   [[8,11,15,12]],
                   [[9,10,14,13]],
                   [[8,9,13,12]],
                   [[10,11,15,14]]])
    planeFlags.extend([0,0,0,0,0,0])
    planeHolds.extend([[],[],[],[],[],[]])
    
    regions = [(0.5*(x2[0]+x3[0]), 0.5*(x2[1]+x3[1]), 0.5*(x2[2]+x3[2])),
               (0.5*(x3[0]+x1[0]), 0.5*(x3[1]+x1[1]), 0.5*(x3[2]+x1[2]))]
    regionFlags = [1,2]
    regionConstraints=[0.5*(he2)**2,0.5*he**2]

    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,facets=planes,facetFlags=planeFlags,facetHoles=planeHolds,
                                       regions=regions,regionFlags=regionFlags,regionConstraints=regionConstraints)
    return domain,boundaryTags,x0,x1,x2,x3

#===============================================================================
#  One box with one shelf
#===============================================================================
def get_domain_one_box_with_one_shelf(x0=(0.0,0.0,0.0),L=(1.0,1.0,1.0),he=0.001,he2=None,he3=None,L1=None,L2=None):
    boundaries = ['left', 'right', 'bottom', 'top', 'front', 'back']
    boundaryTags = dict([(key, i + 1) for (i, key) in enumerate(boundaries)])
    x1 = (x0[0]+L[0],x0[1]+L[1],x0[2]+L[2])
    if L1==None:
        L1 = 0.4*L[1]
    if L2==None:
        L2 = 0.6*L[1]
    x2 = (x0[0],x0[1]+L1,x0[2])
    x3 = (x1[0],x0[1]+L2,x1[2])
    
    if he2==None:
        he2=he
    if he3==None:
        he3=he2
    
    vertices = [(x0[0],x0[1],x0[2]),
                (x1[0],x0[1],x0[2]),
                (x1[0],x1[1],x0[2]),
                (x0[0],x1[1],x0[2]),
                (x0[0],x0[1],x1[2]),
                (x1[0],x0[1],x1[2]),
                (x1[0],x1[1],x1[2]),
                (x0[0],x1[1],x1[2]),
                
                (x2[0],x2[1],x2[2]),
                (x3[0],x2[1],x2[2]),
                (x3[0],x3[1],x2[2]),
                (x2[0],x3[1],x2[2]),
                (x2[0],x2[1],x3[2]),
                (x3[0],x2[1],x3[2]),
                (x3[0],x3[1],x3[2]),
                (x2[0],x3[1],x3[2])]
    vertexFlags=[boundaryTags['bottom'],boundaryTags['bottom'],boundaryTags['top'],boundaryTags['top'],
                boundaryTags['bottom'],boundaryTags['bottom'],boundaryTags['top'],boundaryTags['top'],
                0,0,0,0,
                0,0,0,0,]
    
    planes = []
    planeFlags = []
    planeHolds = []
    planes.extend([[[0,1,9,10,2,3,11,8],[8,9],[10,11]],
                   [[4,5,13,14,6,7,15,12],[12,13],[14,15]],
                   [[0,8,11,3,7,15,12,4],[8,12],[15,11]],
                   [[1,9,10,2,6,14,13,5],[9,13],[14,10]],
                   [[0,1,5,4]],
                   [[2,3,7,6]]])
    planeFlags.extend([boundaryTags['back'],boundaryTags['front'],
                       boundaryTags['left'],boundaryTags['right'],
                       boundaryTags['bottom'],boundaryTags['top'],])
    planeHolds.extend([[],[],[],[],[],[]])
    
    planes.extend([[[10,11,15,14]],
                   [[8,9,13,12]],
                   ])
    planeFlags.extend([0,0])
    planeHolds.extend([[],[],])
    
    regions = [(0.5*(x0[0]+x1[0]), 0.5*(x0[1]+x2[1]), 0.5*(x0[2]+x1[2])),
               (0.5*(x0[0]+x1[0]), 0.5*(x2[1]+x3[1]), 0.5*(x0[2]+x1[2])),
               (0.5*(x0[0]+x1[0]), 0.5*(x3[1]+x1[1]), 0.5*(x0[2]+x1[2]))]
    regionFlags = [1,2,3]
    regionConstraints=[0.5*he**2,0.5*(he2)**2,0.5*(he3)**2,]
    
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,facets=planes,facetFlags=planeFlags,facetHoles=planeHolds,
                                       regions=regions,regionFlags=regionFlags,regionConstraints=regionConstraints)
    return domain,boundaryTags,x0,x1,x2,x3



#===============================================================================
# For 2D cylinder problems
#===============================================================================
import math
from proteus import Domain

def circular_cross_section(center,radius,theta):
    return (radius*math.sin(theta)+center[0],
            radius*math.cos(theta)+center[1])
def get_pseudo_3D_cylinder_domain(
               x0=(0.0,0.0,0.0),L=(1.0,1.0,1.0),
               radius=0.1,
               center=(0.5,0.5),
               n_points_on_obstacle=2*21-2,
               cross_section=circular_cross_section,
               thetaOffset=0.0,
               he=1.0,
               he2=None):
    if he2==None:
        he2=he
    x1 = (x0[0]+L[0],x0[1]+L[1],x0[2]+L[2])
    boundaries=['left',
                'right',
                'bottom',
                'top',
                'front',
                'back',]
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    #work around the domain from (0.0,0.0) going counterclockwise
    vertexKeys = ['left_bottom',
                  'right_bottom',
                  'right_top',
                  'left_top']
    vertices = [[x0[0],x0[1]],
                [x1[0],x0[1]],
                [x1[0],x1[1]],
                [x0[0],x1[1]]]
    vertexFlags = [boundaryTags['bottom'],
                   boundaryTags['bottom'],
                   boundaryTags['top'],
                   boundaryTags['top']]
    nv = len(vertices)
    #cylinder
    theta=thetaOffset
    pb  = cross_section(center,radius,theta)
    vertices.append([pb[0],pb[1]])
    vertexKeys.append('obstacle_'+repr(0))
    vertexFlags.append(boundaryTags['back'])
    for gb in range(1,n_points_on_obstacle):
        theta = float(gb)/float(n_points_on_obstacle)*2.0*math.pi+thetaOffset
        pb  = cross_section(center,radius,theta)
        vertexKeys.append('obstacle_'+repr(gb))
        vertices.append([pb[0],pb[1]])
        vertexFlags.append(boundaryTags['back'])
    #convert to 3D
    vertices3dDict={}
    vertices3d=[]
    vertexFlags3d=[]
    facets3d=[]
    facetFlags3d=[]
    facetHoles3d=[]
    front_cylinder=[]
    back_cylinder=[]
    for vN,v in enumerate(vertices):
        vertices3dDict[vertexKeys[vN]+'_back'] = vN
        vertices3d.append([v[0],v[1],x0[2]])
        vertexFlags3d.append(boundaryTags['back'])
        if 'obstacle' in vertexKeys[vN]:
            back_cylinder.append(vN)
    for vN,v in enumerate(vertices):
        vertices3dDict[vertexKeys[vN]+'_front']=vN+len(vertices)
        vertices3d.append([v[0],v[1],x1[2]])
        vertexFlags3d.append(boundaryTags['front'])#note that the original tag is back
        if 'obstacle' in vertexKeys[vN]:
            front_cylinder.append(vN+len(vertices))

    #left
    facets3d.append([[vertices3dDict['left_bottom_front'],
                      vertices3dDict['left_bottom_back'],
                      vertices3dDict['left_top_back'],
                      vertices3dDict['left_top_front']]])
    facetFlags3d.append(boundaryTags['left'])
    facetHoles3d.append([])
    #right
    facets3d.append([[vertices3dDict['right_bottom_front'],
                     vertices3dDict['right_bottom_back'],
                     vertices3dDict['right_top_back'],
                     vertices3dDict['right_top_front']]])
    facetFlags3d.append(boundaryTags['right'])
    facetHoles3d.append([])
    #top
    facets3d.append([[vertices3dDict['left_top_front'],
                     vertices3dDict['right_top_front'],
                     vertices3dDict['right_top_back'],
                     vertices3dDict['left_top_back']]])
    facetFlags3d.append(boundaryTags['top'])
    facetHoles3d.append([])
    #bottom
    facets3d.append([[vertices3dDict['left_bottom_front'],
                     vertices3dDict['right_bottom_front'],
                     vertices3dDict['right_bottom_back'],
                     vertices3dDict['left_bottom_back']]])
    facetFlags3d.append(boundaryTags['bottom'])
    facetHoles3d.append([])
    #front
    facets3d.append([[vertices3dDict['left_bottom_front'],
                      vertices3dDict['right_bottom_front'],
                      vertices3dDict['right_top_front'],
                      vertices3dDict['left_top_front']],
                     front_cylinder])#add points on the front circle
    facetFlags3d.append(boundaryTags['front'])
    facetHoles3d.append([])
    #back
    facets3d.append([[vertices3dDict['left_bottom_back'],
                     vertices3dDict['right_bottom_back'],
                     vertices3dDict['right_top_back'],
                     vertices3dDict['left_top_back']],
                     back_cylinder])#add points on the back circle 
    facetFlags3d.append(boundaryTags['back'])
    facetHoles3d.append([])
    #sides of cylinder
    for fN in range(n_points_on_obstacle-1):
        facets3d.append([[front_cylinder[fN],
                         back_cylinder[fN],
                         back_cylinder[fN+1],
                         front_cylinder[fN+1]]])
        facetFlags3d.append(0)
        facetHoles3d.append([])
        
    facets3d.append([[front_cylinder[-1],
                      back_cylinder[-1],
                      back_cylinder[0],
                      front_cylinder[0]]])
    facetFlags3d.append(0)
    facetHoles3d.append([])

    
    #region 
    
    regions = [(center[0],center[1], 0.5*(x0[2]+x1[2])),
               (0.1*x0[0]+0.9*x1[0],0.1*x0[1]+0.9*x1[1],0.1*x0[2]+0.9*x1[2])]
    regionFlags = [1,2,]
    regionConstraints=[0.5*he2**2,0.5*(he)**2]
    # make domain
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices3d,
                                                 vertexFlags=vertexFlags3d,
                                                 facets=facets3d,
                                                 facetFlags=facetFlags3d,
                                                 facetHoles=facetHoles3d,
                                                 regions=regions,regionFlags=regionFlags,regionConstraints=regionConstraints)
    domain.boundaryTags = boundaryTags
    return domain,boundaryTags



#===============================================================================
# For 2D cylinder + box around this cylinder
#===============================================================================
def get_pseudo_3D_cylinder_box_domain(
               x0=(0.0,0.0,0.0),L=(1.0,1.0,1.0),
               x2=None,x3=None,
               radius=0.1,
               center=(0.5,0.5),
               n_points_on_obstacle=2*21-2,
               cross_section=circular_cross_section,
               thetaOffset=0.0,
               he=1.0,
               he2=None,
               he3=None):
    if he2==None:
        he2=he
    if he3==None:
        he3=he2
        
    x1 = (x0[0]+L[0],x0[1]+L[1],x0[2]+L[2])
    
    if x2==None:
        x2 = (0.25*(x0[0]+x1[0]), 0.25*(x0[1]+x1[1]))
    if x3==None:
        x3 = (0.75*(x0[0]+x1[0]), 0.75*(x0[1]+x1[1]))
    
    boundaries=['left',
                'right',
                'bottom',
                'top',
                'front',
                'back',]
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    #work around the domain from (0.0,0.0) going counterclockwise
    vertexKeys = ['left_bottom',
                  'right_bottom',
                  'right_top',
                  'left_top']
    vertices = [[x0[0],x0[1]],
                [x1[0],x0[1]],
                [x1[0],x1[1]],
                [x0[0],x1[1]]]
    vertexFlags = [boundaryTags['bottom'],
                   boundaryTags['bottom'],
                   boundaryTags['top'],
                   boundaryTags['top']]

    #cylinder
    theta=thetaOffset
    pb  = cross_section(center,radius,theta)
    vertices.append([pb[0],pb[1]])
    vertexKeys.append('obstacle_'+repr(0))
    vertexFlags.append(boundaryTags['back'])
    for gb in range(1,n_points_on_obstacle):
        theta = float(gb)/float(n_points_on_obstacle)*2.0*math.pi+thetaOffset
        pb  = cross_section(center,radius,theta)
        vertexKeys.append('obstacle_'+repr(gb))
        vertices.append([pb[0],pb[1]])
        vertexFlags.append(boundaryTags['back'])
    #box
    vertexKeys.extend(['box_1','box_2','box_3','box_4'])
    vertices.extend([[x2[0],x2[1]],
                     [x3[0],x2[1]],
                     [x3[0],x3[1]],
                     [x2[0],x3[1]]])
    vertexFlags.extend([boundaryTags['back'],
                        boundaryTags['back'],
                        boundaryTags['back'],
                        boundaryTags['back']])
    #convert to 3D
    vertices3dDict={}
    vertices3d=[]
    vertexFlags3d=[]
    facets3d=[]
    facetFlags3d=[]
    facetHoles3d=[]
    front_cylinder=[]
    back_cylinder=[]
    front_box=[]
    back_box=[]
    
    for vN,v in enumerate(vertices):
        vertices3dDict[vertexKeys[vN]+'_back'] = vN
        vertices3d.append([v[0],v[1],x0[2]])
        vertexFlags3d.append(boundaryTags['back'])#note that the original tag is back
        if 'obstacle' in vertexKeys[vN]:
            back_cylinder.append(vN)
        if 'box' in vertexKeys[vN]:
            back_box.append(vN)
    for vN,v in enumerate(vertices):
        vertices3dDict[vertexKeys[vN]+'_front']=vN+len(vertices)
        vertices3d.append([v[0],v[1],x1[2]])
        vertexFlags3d.append(boundaryTags['front'])#note that the original tag is back
        if 'obstacle' in vertexKeys[vN]:
            front_cylinder.append(vN+len(vertices))
        if 'box' in vertexKeys[vN]:
            front_box.append(vN+len(vertices))

    #left
    facets3d.append([[vertices3dDict['left_bottom_front'],
                      vertices3dDict['left_bottom_back'],
                      vertices3dDict['left_top_back'],
                      vertices3dDict['left_top_front']]])
    facetFlags3d.append(boundaryTags['left'])
    facetHoles3d.append([])
    #right
    facets3d.append([[vertices3dDict['right_bottom_front'],
                     vertices3dDict['right_bottom_back'],
                     vertices3dDict['right_top_back'],
                     vertices3dDict['right_top_front']]])
    facetFlags3d.append(boundaryTags['right'])
    facetHoles3d.append([])
    #top
    facets3d.append([[vertices3dDict['left_top_front'],
                     vertices3dDict['right_top_front'],
                     vertices3dDict['right_top_back'],
                     vertices3dDict['left_top_back']]])
    facetFlags3d.append(boundaryTags['top'])
    facetHoles3d.append([])
    #bottom
    facets3d.append([[vertices3dDict['left_bottom_front'],
                     vertices3dDict['right_bottom_front'],
                     vertices3dDict['right_bottom_back'],
                     vertices3dDict['left_bottom_back']]])
    facetFlags3d.append(boundaryTags['bottom'])
    facetHoles3d.append([])
    #front
    facets3d.append([[vertices3dDict['left_bottom_front'],
                      vertices3dDict['right_bottom_front'],
                      vertices3dDict['right_top_front'],
                      vertices3dDict['left_top_front']],
                      front_cylinder,
                      front_box])#add points on the front circle
    facetFlags3d.append(boundaryTags['front'])
    facetHoles3d.append([])
    #back
    facets3d.append([[vertices3dDict['left_bottom_back'],
                     vertices3dDict['right_bottom_back'],
                     vertices3dDict['right_top_back'],
                     vertices3dDict['left_top_back']],
                     back_cylinder,
                     back_box])#add points on the back circle 
    facetFlags3d.append(boundaryTags['back'])
    facetHoles3d.append([])
    #cylinder
    for fN in range(n_points_on_obstacle-1):
        facets3d.append([[front_cylinder[fN],
                         back_cylinder[fN],
                         back_cylinder[fN+1],
                         front_cylinder[fN+1]]])
        facetFlags3d.append(0)
        facetHoles3d.append([])
        
    facets3d.append([[front_cylinder[-1],
                      back_cylinder[-1],
                      back_cylinder[0],
                      front_cylinder[0]]])
    facetFlags3d.append(0)
    facetHoles3d.append([])
    
    #sides of box
    for fN in range(3):
        facets3d.append([[front_box[fN],
                         back_box[fN],
                         back_box[fN+1],
                         front_box[fN+1]]])
        facetFlags3d.append(0)
        facetHoles3d.append([])
        
    facets3d.append([[front_box[-1],
                      back_box[-1],
                      back_box[0],
                      front_box[0]]])
    facetFlags3d.append(0)
    facetHoles3d.append([])
    
    
    #region 
    
    regions = [(center[0],center[1],                    0.5*(x0[2]+x1[2])),
               (0.1*x2[0]+0.9*x3[0],0.1*x2[1]+0.9*x3[1],0.5*(x0[2]+x1[2])),
               (0.1*x0[0]+0.9*x1[0],0.1*x0[1]+0.9*x1[1],0.5*(x0[2]+x1[2]))]
    regionFlags = [1,2,3]
    regionConstraints=[0.5*he2**2,0.5*he3**2,0.5*(he)**2]
    # make domain
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices3d,
                                                 vertexFlags=vertexFlags3d,
                                                 facets=facets3d,
                                                 facetFlags=facetFlags3d,
                                                 facetHoles=facetHoles3d,
                                                 regions=regions,regionFlags=regionFlags,regionConstraints=regionConstraints)
    domain.boundaryTags = boundaryTags
    return domain,boundaryTags
