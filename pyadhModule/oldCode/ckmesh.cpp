#include <vector>
#include <map>
#include <set>
#include <list>
#include <string>
#include <fstream>
#include <iostream>
#include "ckmesh.h"
using namespace::std;

int main(int argc, char* argv[])
{
  //read file in 3dm format to extract basic element and node information
  string str;
  vector<int> elementNumbers,nodeNumbers;
  vector<vector<int> > elementNodes;
  vector<int> elementMaterials;
  vector<vector<double> > nodes;
  ifstream srcFile(argv[1]);
  srcFile>>str;
  while (srcFile)
    {
      if (str == "E4T")
        {
          int eN,mat;
          vector<int> elementNodes_element(4);
          srcFile>>eN;
          elementNumbers.push_back(eN-1);
          srcFile>>elementNodes_element[0]
                 >>elementNodes_element[1]
                 >>elementNodes_element[2]
                 >>elementNodes_element[3];
          for(int i=0;i<4;i++)
            elementNodes_element[i]-=1;
          elementNodes.push_back(elementNodes_element);
          srcFile>>mat;
          elementMaterials.push_back(mat);
        }
      else if (str== "ND")
        {
          int nN;
          vector<double> nX(3);
          srcFile>>nN;
          nodeNumbers.push_back(nN-1);
          srcFile>>nX[0]>>nX[1]>>nX[2];
          nodes.push_back(nX);
        }
      srcFile>>str;
    }
  //populate a Mesh  struct with basic  element and node information
  Mesh newMesh;
  newMesh.nNodes_global=(int)nodes.size();
  newMesh.nElements_global=(int)elementNodes.size();
  newMesh.nNodes_element=(int)elementNodes[0].size(); //assume same number  of nodes per  element for all elements in mesh
  newMesh.nodes_element = new int[newMesh.nElements_global*newMesh.nNodes_element];
  newMesh.nodeArray = new double[newMesh.nNodes_global*3];
  cout<<endl<<endl<<"elements"<<endl<<endl;
  for (int eN=0; eN<(int)elementNodes.size(); eN++)
    {
      cout<<"eN = "<<eN;
      for (int nN_element=0; nN_element<(int)elementNodes[eN].size(); nN_element++)
        {
          newMesh.nodes_element[elementNumbers[eN]*4+nN_element] = elementNodes[eN][nN_element];
          cout<<", nN["<<nN_element<<"]="<<elementNodes[eN][nN_element];
        }
      cout<<endl;
    }
  cout<<endl<<endl<<"nodes"<<endl<<endl;
  for (int nN=0; nN<(int)nodes.size(); nN++)
    {
      newMesh.nodeArray[nodeNumbers[nN]*3+0]=nodes[nN][0];
      newMesh.nodeArray[nodeNumbers[nN]*3+1]=nodes[nN][1];
      newMesh.nodeArray[nodeNumbers[nN]*3+2]=nodes[nN][2];
      cout<<"nN = "<<nN;
      cout<<", nN[0] = "<<nodes[nN][0]
          <<", nN[1] = "<<nodes[nN][1]
          <<", nN[2] = "<<nodes[nN][2]<<endl;
    }
  //print the basic node and element info from Mesh struct
  cout<<endl<<endl<<"newMesh"<<endl<<endl;
  cout<<endl<<endl<<"elements"<<endl<<endl;
  for (int eN=0;eN<newMesh.nElements_global;eN++)
    {
      cout<<"eN = "<<eN;
      for (int nN_element=0; nN_element<(int)elementNodes[eN].size(); nN_element++)
        {
          cout<<", nN["<<nN_element<<"]="<<newMesh.nodes_element[eN*4+nN_element];
        }
      cout<<endl;
    }
  cout<<endl<<endl<<"nodes"<<endl<<endl;
  for (int nN=0;nN<newMesh.nNodes_global;nN++)
    {
      cout<<"nN = "<<nN;
      cout<<", nN[0] = "<<newMesh.nodeArray[nN*3+0]
          <<", nN[1] = "<<newMesh.nodeArray[nN*3+1]
          <<", nN[2] = "<<newMesh.nodeArray[nN*3+2]<<endl;
    }
  //extract element boundaries
  newMesh.nNodes_elementBoundary = newMesh.nNodes_element-1;
  newMesh.nElementBoundaries_element = newMesh.nNodes_element;
  set<list<int> > elementBoundaryNodes;
  cout<<endl<<endl<<"Element Boundaries"<<endl;
  for  (int eN=0;eN<newMesh.nElements_global;eN++)
    for (int ebN_element=0;ebN_element<newMesh.nNodes_element;ebN_element++)
      {
		cout<<endl<<"element id = " <<eN<<" ==> ";
        list<int> nodeList;
        for (int nN_elementBoundary=0;nN_elementBoundary<newMesh.nNodes_elementBoundary;nN_elementBoundary++)
          {
            cout<<(ebN_element + nN_elementBoundary)%newMesh.nNodes_element<<" ";
            nodeList.push_back(newMesh.nodes_element[eN*newMesh.nNodes_element+(ebN_element + nN_elementBoundary)%newMesh.nNodes_element]);
          }
        nodeList.sort();
        elementBoundaryNodes.insert(nodeList);
      }
  //print element boundaries, populate nodes_elementBoundary, and extract nodeList to element boundary number mapping
  newMesh.nElementBoundaries_global = (int)elementBoundaryNodes.size();
  newMesh.nodes_elementBoundary = new int[newMesh.nElementBoundaries_global*newMesh.nNodes_elementBoundary];
  map<list<int>,int> elementBoundaryNodes_to_elementBoundaryNumber;
  int ebN=0;
  cout << endl<< endl << endl<<"Sorted Unique Element Boundaries " << endl<<endl;
  for (set<list<int> >::const_iterator nodeList=elementBoundaryNodes.begin();nodeList != elementBoundaryNodes.end();++nodeList)
    {
      cout<<ebN;
      list<int>::const_iterator node;
      elementBoundaryNodes_to_elementBoundaryNumber[*nodeList]=ebN;
      int nN_elementBoundary=0;
      for(node=nodeList->begin();node!=nodeList->end();++node)
        {
          newMesh.nodes_elementBoundary[ebN*newMesh.nNodes_elementBoundary+nN_elementBoundary] = *node;
		  nN_elementBoundary++;
          cout<<", "<<*node;
        }
      cout<<endl;
      ebN++;
    }
  // populate elementBoundaries and elements_elementBoundary
  newMesh.elementBoundaries = new  int[newMesh.nElements_global*newMesh.nElementBoundaries_element];
  newMesh.elements_elementBoundary = new int[newMesh.nElementBoundaries_global*2];

  // Initialize elements_elementBoundary array to -1 */
  for (int ebN=0 ;ebN<newMesh.nElementBoundaries_global; ebN++)
    {
      newMesh.elements_elementBoundary[ebN*2+0]=-1;
      newMesh.elements_elementBoundary[ebN*2+1]=-1;
    }
  cout<<endl<<endl<<endl<<"Populate element Boundaries and elements_elementBoundary Arrays"<<endl<<endl;
  newMesh.nInteriorElementBoundaries_global=0;
  for (int eN=0; eN<newMesh.nElements_global; eN++)
    for (int ebN_element=0; ebN_element<newMesh.nNodes_element; ebN_element++)
    {
        list<int> nodeList;
        for (int nN_elementBoundary=0; nN_elementBoundary<newMesh.nNodes_elementBoundary; nN_elementBoundary++)
        {
            nodeList.push_back(newMesh.nodes_element[eN*newMesh.nNodes_element+(ebN_element + nN_elementBoundary)%newMesh.nNodes_element]);
        }
        nodeList.sort();
        int ebN_global=elementBoundaryNodes_to_elementBoundaryNumber[nodeList];
        newMesh.elementBoundaries[eN*newMesh.nElementBoundaries_element+ebN_element] =  ebN_global;
		cout<<"elementBoundaries["<<eN*newMesh.nElementBoundaries_element+ebN_element<<"] = "<<ebN_global<<"\t,   ";
        if(newMesh.elements_elementBoundary[ebN_global*2+0] == -1)
        {
          newMesh.elements_elementBoundary[ebN_global*2+0] = eN;
		  cout << "elements_elementBoundary[" << ebN_global*2+0 << "] = " <<eN<<endl ;
        }
        else
        {
          newMesh.elements_elementBoundary[ebN_global*2+1] = eN;
          cout << "elements_elementBoundary[" << ebN_global*2+1 << "] = "<<eN<<endl ;
		  newMesh.nInteriorElementBoundaries_global++;
        }
    }
	cout<<endl<<endl<<"Elements_elementBoundaryArray"<<endl<<endl;

// populating the interior and exterior element boundaries arrays
	newMesh.nExteriorElementBoundaries_global = newMesh.nElementBoundaries_global - newMesh.nInteriorElementBoundaries_global;
	newMesh.interiorElementBoundaries = new int[newMesh.nInteriorElementBoundaries_global]();
	newMesh.exteriorElementBoundaries = new int[newMesh.nExteriorElementBoundaries_global]();

    for (int nExt=0, nInt=0, ebN=0; ebN<newMesh.nElementBoundaries_global*2; ebN++)
    { 
		if (newMesh.elements_elementBoundary[ebN*2+1] != -1)
		{
 			cout <<"elements_elementBoundary["<<ebN<<"] = "<<newMesh.elements_elementBoundary[ebN]<<endl;
			newMesh.interiorElementBoundaries[nInt] = ebN;
//			cout<<"interiorElementBoundaries[ "<<nInt<<"] = "<<ebN<<endl<<endl;
			nInt++;
		}
		else
		{
 			cout <<"elements_elementBoundary["<<ebN<<"] = "<<newMesh.elements_elementBoundary[ebN]<<endl;
			newMesh.exteriorElementBoundaries[nExt]=ebN;
//			cout<<"exteriorElementBoundaries[ "<<nExt<<"] = "<<ebN<<endl<<endl;
			nExt++;
		}
    }
//print interior and exterior element boundaries arrays
	cout<<endl<<endl<<"nEnteriorElementBoundaries_global = "<<newMesh.nExteriorElementBoundaries_global<<endl;
	cout<<endl<<endl<<"nInteriorElementBoundaries_global = "<<newMesh.nInteriorElementBoundaries_global<<endl;

	cout<<endl<<endl<<"exteriorElementBoundariesArray"<<endl<<endl;
    for (int ebN=0; ebN<newMesh.nExteriorElementBoundaries_global; ebN++)
	{
		cout<<"newMesh.exteriorElementBoundaries["<<ebN<<"] = "<<newMesh.exteriorElementBoundaries[ebN]<<endl;
	}
	cout<<endl<<endl<<"interiorElementBoundariesArray"<<endl<<endl;
    for (int ebN=0; ebN<newMesh.nInteriorElementBoundaries_global; ebN++)
	{
		cout<<"newMesh.interiorElementBoundaries["<<ebN<<"] = "<<newMesh.interiorElementBoundaries[ebN]<<endl;
	}

   return 0;
}
