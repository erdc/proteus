#ifndef STUPID_HEAP_H
#define STUPID_HEAP_H

#include <vector>
#include <map>
#include <algorithm>
#include <cassert>
#include <iostream>


class StupidHeap
{
  /***********************************************************************
    Min Heap

    heap   = [(key,val), ...] such that heap[i][1] <= heap[2*i+1][1], heap[2*i+2][1]
    heapPos= {key:i} heap[heapPos[key]]=val

   ***********************************************************************/
 public:
  typedef int KeyType; typedef double ValueType; 
  typedef std::pair<KeyType,ValueType> EntryType;
  typedef std::vector<EntryType>::size_type PositionType;
 StupidHeap(): heap(),heapPos()
    {}

  bool isEmpty() const
  { 
    return heap.empty();
  }

  bool insert(KeyType key,ValueType val, int verbose = 0)
  { 
    PositionType last = heap.size();
    heap.push_back(EntryType(key,val));
    heapPos[key] = last;
    return upHeap(last,verbose);
  }
  bool insertWithCheckForExistingKey(KeyType key, ValueType val, int verbose = 0)
  {
    std::map<KeyType,PositionType>::const_iterator pkey = heapPos.find(key);
    if (pkey != heapPos.end())
      {
	if (verbose > 0)
	  std::cout<<"insertWithCheckForExistingKey found existing key "<<key
		   <<" using updateNoe"""<<std::endl;
	return updateNode(key,val);
      }
    PositionType last = heap.size();
    heap.push_back(EntryType(key,val));
    heapPos[key] = last;
    return upHeap(last);
  }
  //return smallest value in heap, maintain heap property
  EntryType pop(int verbose = 0)
  {
    assert(!heap.empty());
    EntryType minval = heap[0];
    EntryType lastval= heap.back();
    heap.pop_back();
    heapPos.erase(minval.first);
    bool failedDown  = false;
    if (heap.size() > 0)
      {
	heap[0] = lastval;
	heapPos[lastval.first] = 0;
	failedDown = downHeap(0,verbose);
      }
    return minval;
  }

  //modify value associated with key and restore heap property
  //requires that key exist in heap
  bool updateNode(KeyType key, ValueType newval, int verbose = 0)
  {
    PositionType pos = heapPos[key];
    assert(heap[pos].first == key);
    heap[pos] = EntryType(key,newval);
    return downHeap(pos,verbose);
  }

  //modify value associated with key and restore heap property
  //requires that key exist in heap
  //this version only updates value if the new one is less than the original
  bool updateNodeWithMin(KeyType key, ValueType newval, int verbose = 0)
  {
    PositionType pos = heapPos[key];
    assert(heap[pos].first == key);
    if (newval >= heap[pos].second)
      return false;
    heap[pos] = EntryType(key,newval);
    return downHeap(pos,verbose);
  }
  //recursively move node at pos up the heap while it is less than its parent
  bool upHeap(PositionType pos, int verbose = 0);


  //reestablish heap property assuming left and right children of pos are heap,
  // but left or right child < pos
  bool downHeap(PositionType pos, int verbose = 0);

  bool checkHeap() const
  {
    //could use std::is_heap with comparison functor for pairs
    for (PositionType i = 0; i < heap.size()/2; i++)
      {
	assert(heap[i].second <= heap[2*i+1].second); //left child
	if (2*i + 2 < heap.size())
	  assert(heap[i].second <= heap[2*i+2].second);//right child
      }
    for (std::map<KeyType,PositionType>::const_iterator hp = heapPos.begin();
	 hp != heapPos.end(); hp++)
      {
	int node; PositionType pos;
	node = hp->first; pos = hp->second;
	assert(heap[pos].first == node);
      }
    return true;
  }//checkHeap
  bool printHeap() const;
 protected:
  std::vector<EntryType> heap;
  std::map<KeyType,PositionType> heapPos;
};







#endif
