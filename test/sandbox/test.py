import numpy
#alt
#nodes of reference hex laid out clockwise bottom to top
refNodes = [[0,0,0],
            [0,1,0],
            [1,1,0],
            [1,0,0],
            [0,0,1],
            [0,1,1],
            [1,1,1],
            [1,0,1]]
refNodes.sort()
refNodes = Numeric.array(refNodes)
#edges of reference hex laid out clockwise bottom to top
refEdges = [[[0,0,0],[0,1,0]],
            [[0,1,0],[1,1,0]],
            [[1,1,0],[1,0,0]],
            [[1,0,0],[0,0,0]],
            [[0,0,0],[0,0,1]],
            [[0,1,0],[0,1,1]],
            [[1,1,0],[1,1,1]],
            [[1,0,0],[1,0,1]],
            [[0,0,1],[0,1,1]],
            [[0,1,1],[1,1,1]],
            [[1,1,1],[1,0,1]],
            [[1,0,1],[0,0,1]]]
for e in refEdges:
    e.sort()
refEdges.sort()
refEdges = Numeric.array(refEdges)
#quadrilaterals of reference hex
refQuadrilaterals = [[[0,0,0],[0,1,0],[1,1,0],[1,0,0]],
                     [[0,0,0],[0,1,0],[0,1,1],[0,0,1]],
                     [[0,1,0],[0,1,1],[1,1,1],[1,1,0]],
                     [[1,1,0],[1,1,1],[1,0,1],[1,0,0]],
                     [[1,0,0],[0,0,0],[0,0,1],[1,0,1]],
                     [[0,0,1],[0,1,1],[1,1,1],[1,0,1]]]
for q in refQuadrilaterals:
    q.sort()
refQuadrilaterals.sort()
refQuadrilaterals = Numeric.array(refQuadrilaterals)
#now we have lexicographically sorted reference nodes, edges, and quadrilaterals
#we just need to map them to indeces
ijk=Numeric.zeros((3,),Numeric.Int)
print(ijk.shape)
for i in range(2):
    for j in range(2):
        for  k in range(2):
            ijk[0]=i
            ijk[1]=j
            ijk[2]=k
            nodeIndexes = refNodes+ijk
            edgeIndexes = refEdges+ijk
            quadrilateralIndexes = refQuadrilaterals+ijk
            print(nodeIndexes)
            print(edgeIndexes)
            print(quadrilateralIndexes)
