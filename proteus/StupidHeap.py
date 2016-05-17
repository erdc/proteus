"""
A stupid implementation of a heap for fast marching methods

.. inheritance-diagram:: proteus.StupidHeap
   :parts: 1
"""
class StupidHeap:
    """
    Min Heap

    heap   = [(key,val), ...] such that heap[i][1] <= heap[2*i+1][1], heap[2*i+2][1]
    heapPos= {key:i} heap[heapPos[key]]=val

    """
    #for indexing tuples
    keyIndex=0; valIndex=1;
    def __init__(self):
        self.heap   = []
        self.heapPos= {}

    def isempty(self):
        return len(self.heap) == 0

    def insert(self,key,val):
        """

        """
        last = len(self.heap)
        self.heap.append((key,val))
        self.heapPos[key]=last
        self.upHeap(last)

    def insertWithCheckForExistingKey(self,key,val):
        """

        """
        if self.heapPos.has_key(key):
            print """insertWithCheckForExistingKey found existing key=%s using updateNode!""" % key
            return self.updateNode(key,val)
        last = len(self.heap)
        self.heap.append((key,val))
        self.heapPos[key]=last
        self.upHeap(last)

    def pop(self,debugLevel=0):
        """
        return smallest value in heap, maintain heap property
        """
        assert len(self.heap) > 0, "can't pop empty heap"
        K=0; V=1;
        mnode = self.heap[0]
        bnode = self.heap.pop()
        if len(self.heap) > 0:
            self.heap[0] = bnode
            self.heapPos[bnode[K]] = 0
            self.downHeap(0,debugLevel=debugLevel)
        return mnode

    def updateNode(self,key,newval,debugLevel=0):
        """
        modify value associated with key and restore heap property
        requires that key exist in heap
        """
        K = 0; V= 1;
        pos = self.heapPos[key]
        assert self.heap[pos][K] == key, "updateVal key=%s self.heap[pos][K]=%s " % (key,self.heap[pos][K])
        self.heap[pos] = (key,newval)
        return self.downHeap(pos,debugLevel=debugLevel)

    def updateNodeWithMin(self,key,newval,debugLevel=0):
        """
        modify value associated with key and restore heap property
        requires that key exist in heap

        this version only updates value if new value is less
        """
        K = 0; V= 1;
        pos = self.heapPos[key]
        assert self.heap[pos][K] == key, "updateVal key=%s self.heap[pos][K]=%s " % (key,self.heap[pos][K])
        if newval > self.heap[pos][V]:
            return
        self.heap[pos] = (key,newval)
        return self.downHeap(pos,debugLevel=debugLevel)

    def upHeap(self,pos,debugLevel=0):
        """
        recursively move node at pos up the heap while it is less than its parent

        """
        assert pos >= 0 and pos < len(self.heap)
        K = 0; V = 1;
        Ppos = int((pos-1)/2)
        if Ppos < 0 or self.heap[Ppos][V] < self.heap[pos][V]:
            return
        #switch with parent
        tmp = (self.heap[Ppos][K],self.heap[Ppos][V])
        self.heap[Ppos] = (self.heap[pos][K],self.heap[pos][V])
        self.heapPos[tmp[K]]=pos
        self.heap[pos]  = tmp
        self.heapPos[self.heap[Ppos][K]]=Ppos
        return self.upHeap(Ppos,debugLevel)

    def downHeap(self,pos,debugLevel=0):
        """
        assume left and right children of pos are heap, but
        left or right child < pos

        promote smaller of children of node at pos, then repeat on its subheap

        then put original value at pos and work back up
        """
        assert pos >= 0 and pos < len(self.heap)
        K = 0; V=1
        startpos = pos
        endpos = len(self.heap)
        misfit = self.heap[pos]

        childpos = 2*pos+1 #start on left
        while childpos < endpos:
            #
            rightpos = childpos+1
            if rightpos < endpos and self.heap[rightpos][V] <= self.heap[childpos][V]:
                childpos = rightpos
            #promote smaller child and keep going
            tmp = (self.heap[pos][K],self.heap[pos][V])
            self.heap[pos] = (self.heap[childpos][K],self.heap[childpos][V])
            self.heapPos[self.heap[pos][K]] = pos
            self.heap[childpos] = tmp
            self.heapPos[tmp[K]]=childpos
            pos      = childpos
            childpos = 2*pos+1
        #end while

        #have pushed up one child all the way down, pos should now be empty
        self.heap[pos] = misfit
        self.heapPos[misfit[K]]=pos

        return self.upHeap(pos,debugLevel=debugLevel)

    def checkHeap(self):
        K = 0; V=1;
        for i in range(len(self.heap)/2):
            assert self.heap[i][V] <= self.heap[2*i+1][V], "failed node %d = %s left child %d =  %s " (i,self.heap[i],
                                                                                                       2*i+1,
                                                                                                       self.heap[2*i+1])
            if 2*i+2 < len(self.heap):
                assert self.heap[i][V] <= self.heap[2*i+2][V], "failed node %d = %s right child %d =  %s " (i,self.heap[i],
                                                                                                            2*i+2,
                                                                                                            self.heap[2*i+2])
        for i in range(len(self.heap)):
            assert self.heapPos[self.heap[i][K]]==i,  "failed node %d = %s heapPos[%s]=%d " % (i,self.heap[i],
                                                                                              self.heap[i][K],
                                                                                              self.heapPos[self.heap[i][K]])

    def printHeap(self):
        print "heap= %s " % self.heap
        print "heapPos= %s " % self.heapPos

if __name__ == "__main__":
    H = StupidHeap()

    #some hypothetical node-val pairs
    nn = 12
    a = [(i, 10.-i) for i in range(10)]
    a.append((11,5.)); a.append((-1,6.))
    for i in range(nn):
        print "inserting a[%s]=%s " % (i,a[i])

        H.insert(a[i][0],a[i][1])
    H.printHeap()
    H.checkHeap()


    while not H.isempty():
        h = H.pop()
        print "H.pop = (%s,%s) " % h
        H.checkHeap()

    for i in range(nn):
        print "re inserting a[%s]=%s " % (i,a[i])

        H.insert(a[i][0],a[i][1])

    #now modify some values
    H.updateNode(1,-10); H.checkHeap()
    H.updateNode(7,100.); H.checkHeap()
    while not H.isempty():
        h = H.pop()
        print "H.pop = (%s,%s) " % h
        H.checkHeap()
