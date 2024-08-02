import math
#example 1 from COBRAS report see page 11 of appendix (Figures 1 and 2)
#original starts at x=-19m
xshift=19.0
yshift=0.0
vertices = {'bottom_left' : (-19.0,0.0),
            'bottom_inflow_right': (-18.0,0.0),#for forcing inflow
            'bottom_inflow_mid': (-18.5,0.0),#for forcing inflow
            'slope_bottom_left' : (10.0,0.0),
            'slope_shelf_left' : (23.0,0.52),
            'bottom_outflow_left': (31.0,0.52),#for forcing outflow
            'bottom_outflow_mid': (31.5,0.52), #for forcing outflow
            'slope_shelf_right' : (32.0,0.52),
            'top_right' : (32.0,1.2),
            'top_left' : (-19.0,1.2),
            'tetrapod_bottom_left' : (22.059,0.482),#tetrapod
            'tetrapod_shelf1_left' : (22.235,0.60),
            'tetrapod_shelf1_right': (22.355,0.60),
            'tetrapod_shelf2_left' : (22.88,0.95),
            'tetrapod_shelf2_right': (23,0.95),
            'tetrapod_bottom_right': (22.665,0.5062),
            'tetrapod_shelf3_left' : (22.841,0.63),
            'tetrapod_shelf3_right': (23.0,0.63),
            'rubble_mound_bottom_left': (22.713,0.5085),#rubble mound
            'rubble_mound_upper_left':  (22.85,0.60),
            'rubble_mound_upper_right': (23.433+0.15,0.60),
            'rubble_mound_bottom_right': (23.703,0.52),
            'right_filter_layer_upper_left':  (23.433,0.63),#filter layer on right side
            'right_filter_layer_upper_right': (23.433+0.159,0.63),
            'right_filter_layer_lower_right': (23.753,0.52),
            'caisson_lower_left' : (23.0,0.60),
            #duplicate of 'tetrapod_shelf2_right' 'caisson_upper_left' : (23.0,0.95),
            'caisson_upper_right': (23.433,0.95),
            'caisson_lower_right': (23.433,0.60)}

vertexId = {}
for iv,key in enumerate(vertices.keys()):
    vertexId[key] = iv
nvertices = len(vertices)

segmentLabels = {'default':0,
                 'left': 1,
                 'bottom' : 2,
                 'right'  : 3,
                 'top'    : 4,
                 'bottom_slope': 5,
                 'bottom_shelf': 6,
                 'tetrapod_exterior':7,
                 'caisson_exterior' :8,
                 'tetrapod_interior':9,
                 'caisson_interior' :9,
                 'rubble_mound': 10,
                 'filter_layer_exterior': 11,
                 'bottom_inflow' : 12,
                 'bottom_outflow' : 13}

segmentList = [('top_left','bottom_left','left'),
               ('bottom_left','bottom_inflow_mid','bottom_inflow'),
               ('bottom_inflow_mid','bottom_inflow_right','bottom_inflow'),
               ('bottom_inflow_right','slope_bottom_left','bottom'),
               ('slope_bottom_left','tetrapod_bottom_left','bottom_slope'),
               ('tetrapod_bottom_left','tetrapod_shelf1_left','tetrapod_exterior'), #loop around tetrapod
               ('tetrapod_shelf1_left' ,'tetrapod_shelf1_right','tetrapod_exterior'),
               ('tetrapod_shelf1_right' , 'tetrapod_shelf2_left','tetrapod_exterior'),
               ('tetrapod_shelf2_left'  , 'tetrapod_shelf2_right','tetrapod_exterior'),
               ('tetrapod_shelf2_right' , 'tetrapod_shelf3_right','tetrapod_interior'),
               ('tetrapod_shelf3_right' , 'tetrapod_shelf3_left','tetrapod_interior'),
               ('tetrapod_shelf3_left' , 'tetrapod_bottom_right','tetrapod_interior'),
               ('tetrapod_bottom_right' , 'tetrapod_bottom_left','tetrapod_interior' ),
               ('tetrapod_bottom_right'  , 'rubble_mound_bottom_left','bottom_slope'),
               ('rubble_mound_bottom_left' , 'slope_shelf_left','bottom_slope'),
               ('slope_shelf_left' , 'rubble_mound_bottom_right', 'bottom_shelf'),
               ('rubble_mound_bottom_left', 'rubble_mound_upper_left','rubble_mound'),
               ('rubble_mound_upper_left' , 'caisson_lower_left','rubble_mound'),
               ('caisson_lower_left' ,'tetrapod_shelf3_right','caisson_interior'),
               ('caisson_lower_left' , 'caisson_lower_right','rubble_mound'),
               ('caisson_lower_right' , 'right_filter_layer_upper_left', 'caisson_interior'),
               ('caisson_lower_right' , 'rubble_mound_upper_right', 'rubble_mound'),
               ('rubble_mound_upper_right' , 'rubble_mound_bottom_right' ,'rubble_mound'),
               ('rubble_mound_bottom_right' , 'right_filter_layer_lower_right','bottom_shelf'),
               ('right_filter_layer_lower_right' , 'right_filter_layer_upper_right','filter_layer_exterior'),
               ('right_filter_layer_upper_right' , 'right_filter_layer_upper_left','filter_layer_exterior'),
               ('right_filter_layer_upper_left' , 'caisson_upper_right','caisson_exterior'),
               ('caisson_upper_right' , 'tetrapod_shelf2_right','caisson_exterior' ),
               ('tetrapod_shelf2_right'  , 'tetrapod_shelf3_right', 'caisson_interior' ), #hopefully done with whole structure
               ('right_filter_layer_lower_right' , 'bottom_outflow_left','bottom_shelf'),
               ('bottom_outflow_left','bottom_outflow_mid','bottom_outflow'),
               ('bottom_outflow_mid','slope_shelf_right','bottom_outflow'),
               ('slope_shelf_right' , 'top_right','right'),
               ('top_right' ,'top_left','top')]

segments = set([])
for s in segmentList:
    segments.add((vertexId[s[0]],vertexId[s[1]],segmentLabels[s[2]]))

regions = {'tetrapod' : (22.6,0.6,1),             #porosity =0.5 W=370g looks like all have dm_50 0.05 [m]
           'filter_layer_left' : (22.85,0.615,2), #porosity = 0.53, W=6.5g
           'rubble_mound' : (23.5, 0.55, 3),      #porosity = 0.49, W=1.4g
           'filter_layer_right' : (23.5,0.61,2),  #porosity = 0.53, W=6.5g
           'caisson' : (23.3,0.8,4)}              #porosity=0
porous_regions = {}
solid_regions  = {}
for reg in ['tetrapod','filter_layer_left','rubble_mound','filter_layer_right']:
    porous_regions[reg]= regions[reg]
for reg in ['caisson']:
    solid_regions[reg]= regions[reg]


poly = open('cobras_wave_flume_1.poly','w')
poly.write('%d %d %d %d \n' % (nvertices,2,0,0))
#write vertices
poly.write("#vertices \n")
for key,p in vertices.items():
    poly.write('%d %12.5e %12.5e #%s \n' % (vertexId[key]+1,xshift+p[0],yshift+p[1],key))
#write segments
nSegments = len(segments)
poly.write('%d %d \n' % (nSegments,1))
poly.write("#segments \n")
for sN,s in enumerate(segments):
    poly.write('%d %d %d %d \n' % (sN+1,s[0]+1,s[1]+1,s[2]))
#if mesh just the outside of the structure insert holes here
nholes = len(solid_regions)#0
poly.write('%d \n' % (nholes,))
poly.write("#holes \n")
for i,reg in enumerate(solid_regions.keys()):
    poly.write('%d %12.5e %12.5e #%d  %s \n' % (i+1,solid_regions[reg][0]+xshift,solid_regions[reg][1]+yshift,solid_regions[reg][2],reg))
#assign region ids for porosity etc also holes if wanted
nregions = len(porous_regions)
poly.write('%d \n' % (nregions))
poly.write("#regions \n")
for i,reg in enumerate(porous_regions.keys()):
    poly.write('%d %12.5e %12.5e %d #%s \n' % (i+1,porous_regions[reg][0]+xshift,porous_regions[reg][1]+yshift,porous_regions[reg][2],reg))
poly.close()