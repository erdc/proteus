#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <vector>
#include <sstream>
#include "mesh.h"

using namespace std;

extern "C" int read3DM(const char* );
int write3DM(const char* );
int deleteMesh(Mesh *);

/*****************************************************************/
void error(const char* p, const char* p2="")
/*****************************************************************/
{
  cerr<<p<<' '<<p2<<'\n';
  exit(1);
}
/*****************************************************************/
/* mesh.c: */
int read3DM(const char* filename)
/*****************************************************************/
{
  ifstream meshFile(filename);
  char ch;
  string line, word;

//  string strTypes = "ND,E3T,E4T" ;
//  vector<string> v;
//  enum enumTypes = {MESH3D=1,ND,E3T,E45} ;

  bool first = true, ndfirst = true;
  int nodeid, size, newsize, oldsize=0, npos, nd=-1, ele=0;
//  int *temp;
  double gnodeid;
  Mesh *pmesh=new Mesh;

  if(!meshFile)error("cannot open inputfile ",filename);
/* Read in the header line "MESH3D") */
  getline(meshFile, line);
  cout << line << endl ;
/*  allocate  and  populate Mesh struct */
  while( !meshFile.eof())
  {
    getline(meshFile, line);
    cout << line << endl ;
    istringstream istream(line) ;
    if(istream >> word) 
      {
        if (first) 
          {
            string::size_type loc = word.find("E4T");
            if ( loc == 0 && loc != npos) 
              {
		size = 4 ;
              } 
            else 
              {
                string::size_type loc = word.find("E3T");
                if ( loc == 0 ) 
                  {
                    size = 3 ;
                  } 
                else 
                  { 
                    cout << " not E4T or E3T" << endl ; 
                  }
              }
          }
        if (ndfirst) 
          {
            string::size_type loc = word.find("ND");
            if ( loc == 0 ) 
              {
                size = 3 ;
		oldsize = 0 ;
                nd ++ ;
                ndfirst = false;
		first = true ;
              }
          }
	newsize = oldsize + size ;
	if (nd < 0)
          { 
            ele ++;
            int *temp;
            temp = pmesh->elementNodesArray;
            pmesh->elementNodesArray = new int[newsize];
    	   if (!first){
	     for (int i=0; i<size; i++){
	        pmesh->elementNodesArray[i] = temp[i] ;
             }
	   }
	}else{
	   nd ++;
	   double *gtemp;
	   gtemp = pmesh->nodeArray;
	   pmesh->nodeArray = new double[newsize];
    	   if (!first){
	     for (int i=0; i<size; i++){
	        pmesh->nodeArray[i] = gtemp[i] ;
             }
	   }
        }
    }
    int i = -1;
    if (ndfirst){
      while ( istream >> nodeid) {
        if(i>-1 && i<size){
            pmesh->elementNodesArray[size+i]=nodeid;
    	    cout << pmesh->elementNodesArray[size+i] ;
	}
	i++ ;
      }
    }else{
      while ( istream >> gnodeid) {
        if(i>-1 && i<size){
            pmesh->nodeArray[size+i]=gnodeid;
            cout << size << i << pmesh->nodeArray[size+i] << endl ;
        }
        i++ ;
      }
    }
    cout << "*** \n";
    first = false ;
    oldsize = newsize;
//    delete [] temp ;
//    delete [] gtemp ;
  }
  meshFile.close() ; 
  pmesh->nElements_global = ele ;	
  pmesh->nNodes_element = size ;
  pmesh->nNodes_global = nd ;	
  cout << "\n number of elements = " << ele << endl;
  cout << "\n number of nodes = " << nd << endl;
  return(0);
}
/*
  while (meshFile.get(ch))
    cout<<ch<<endl;
*/ 
/*****************************************************************/
int write3DM(const char* filename)
{
/* write just the 3DM info */
  ofstream meshFile(filename);
  char ch;
  string line, word;

  return(0);
}
/*****************************************************************/
int deleteMesh(Mesh* pmesh)
{
/* deallocate  Mesh structure */
   free (pmesh);
   return(0);
}
/*****************************************************************/
int main(int argc, char* argv[])
{
	int iret;
	const char *filename={"adh.3dm"};
	const char *ofilename={"adh_3dm.out"};
	Mesh *newmesh;

	/* call read, write, and delete */
	iret = read3DM(filename);
	if(iret != 0){ error("error!! reading the 3D mesh");}
	iret = write3DM(ofilename);
	if(iret != 0){ error("error!! writing the 3D mesh");}
	iret = deleteMesh(newmesh);
	if(iret != 0){ error("error!! deleting the 3D mesh");}

}
/*****************************************************************/
