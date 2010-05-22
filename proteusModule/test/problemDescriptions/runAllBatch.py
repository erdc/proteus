if p.nd == 1:
    n.nLevels = min([5,n.nLevels])
elif p.nd == 2:
    n.nLevels = min([3,n.nLevels])
elif p.nd == 3:
    n.nLevels = min([2,n.nLevels])
if n.DT != None:
    n.nDTout = 5
    p.T = n.nDTout*n.DT
start
quit
