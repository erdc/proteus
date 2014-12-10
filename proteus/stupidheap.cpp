#include "stupidheap.h"
#include <cassert>
#include <iostream>

bool StupidHeap::upHeap(StupidHeap::PositionType pos, int verbose)
{
  //recursively move node at pos up the heap while it is less than its parent
  assert(pos >= 0 && pos < heap.size());
  if (pos == 0)//has no parent
    return false;
  PositionType Ppos = PositionType((pos-1)/2);//parent position
  if (verbose > 0)
    {
      std::cout<<"upHeap pos= "<<pos<<" Parent position = "<<Ppos<<std::endl;
    }
  if (heap[Ppos].second < heap[pos].second)
    return false;
  //switch pos with its parent
  EntryType tmp  = heap[Ppos];//EntryType(heap[Ppos].first,heap[Ppos].second);
  heap[Ppos] = heap[pos];               //EntryType(heap[pos].first,heap[pos].second);
  heapPos[heap[Ppos].first] = Ppos;
  heap[pos]          = tmp;
  heapPos[heap[pos].first] = pos;
  return upHeap(Ppos,verbose);
}

bool StupidHeap::downHeap(StupidHeap::PositionType pos, int verbose)
{
  /**************************************************
    assume left and right children of pos are heap, but
    left or right child < pos
        
    promote smaller of children of node at pos, then repeat on its subheap

    then put original value at pos and work back up
  ***************************************************/
  assert(pos >= 0 && pos < heap.size());
  
  //PositionType startpos = pos;
  PositionType endpos   = heap.size();
  EntryType misfit      = heap[pos];

  PositionType childpos = 2*pos + 1; //start on left
  while (childpos < endpos)
    {
      //
      PositionType rightpos = childpos+1;
      if (rightpos < endpos && heap[rightpos].second <= heap[childpos].second)
	{
	  childpos = rightpos;
	}
      //promote smaller child and keep going
      EntryType tmp = heap[pos];
      heap[pos] = heap[childpos];
      heapPos[heap[pos].first] = pos;
      heap[childpos] = tmp;
      heapPos[heap[childpos].first] = childpos;
      pos      = childpos;
      childpos = 2*pos+1;
    }//while
  //have pushed up one child all the way down, pos should now be empty
  heap[pos] = misfit;
  heapPos[heap[pos].first] = pos;

  return upHeap(pos,verbose);

}

bool StupidHeap::printHeap() const
{
  std::cout<<"heap:"<<std::endl;
  for (PositionType i = 0; i < heap.size(); i++)
    {
      std::cout<<"("<<heap[i].first<<","<<heap[i].second<<")"<<std::endl;
    }
  return false;
}
