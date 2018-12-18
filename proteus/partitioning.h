#ifndef PARTITIONING_H
#define PARTITIONING_H
#include <iostream>
#include <valarray>
#include "mpi.h"
#include "hdf5.h"
#include "petsc.h"
#include "petscsys.h"
#include "mesh.h"
#include "meshio.h"

namespace proteus
{
  //--memory profiling
  /*
   * Author:  David Robert Nadeau
   * Site:    http://NadeauSoftware.com/
   * License: Creative Commons Attribution 3.0 Unported License
   *          http://creativecommons.org/licenses/by/3.0/deed.en_US
   */

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

  /**
   * Returns the peak (maximum so far) resident set size (physical
   * memory use) measured in bytes, or zero if the value cannot be
   * determined on this OS.
   */
  inline size_t getPeakRSS( )
  {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
      return (size_t)0L;/* Can't open? */
    if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
      {
        close( fd );
        return (size_t)0L;/* Can't read? */
      }
    close( fd );
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;/* Unsupported. */
#endif
  }

  /**
   * Returns the current resident set size (physical memory use) measured
   * in bytes, or zero if the value cannot be determined on this OS.
   */
  inline size_t getCurrentRSS( )
  {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
                    (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
      return (size_t)0L;/* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE* fp = NULL;
    if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
      return (size_t)0L;/* Can't open? */
    if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
      {
        fclose( fp );
        return (size_t)0L;/* Can't read? */
      }
    fclose( fp );
    return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L;/* Unsupported. */
#endif
  }
  //--memory profiling
  inline int enforceMemoryLimit(const MPI_Comm& PROTEUS_COMM_WORLD, int rank, double max_rss_gb,const char* msg)
  {
    double current, current_global,gb(1.0e-9);
    PetscBarrier(NULL);
    current = double(getCurrentRSS())*gb;
    PetscBarrier(NULL);
    current_global=0.0;
    MPI_Allreduce(&current,&current_global,1,MPI_DOUBLE,MPI_MAX,PROTEUS_COMM_WORLD);
    if (current > max_rss_gb)
      {
        std::cout<<"Raising PETSC_ERR_MEM, Memory usage  on rank "<<rank<<'\t'<<current<<"GB"<<'\t'<<"limit "<<max_rss_gb<<std::endl;
        SETERRABORT(PROTEUS_COMM_WORLD,PETSC_ERR_MEM,"Exceeded Proteus memory limit");
      }
    if (rank ==  0)
      std::cout<<msg<<std::endl
               <<"Max memory usage per core "<<current_global<<"GB"<<std::endl;
    return 0;
  }

  extern int partitionElementsOriginal(const MPI_Comm& PROTEUS_COMM_WORLD, Mesh& mesh, int nElements_overlap);

  extern int partitionNodes(const MPI_Comm& PROTEUS_COMM_WORLD, Mesh& mesh, int nNodes_overlap);

  extern int partitionNodesFromTetgenFiles(const MPI_Comm& PROTEUS_COMM_WORLD, const char* filebase, int indexBase, Mesh& newMesh, int nNodes_overlap);

  extern int partitionElements(const MPI_Comm& PROTEUS_COMM_WORLD, Mesh& mesh, int nElements_overlap);

  extern int buildQuadraticSubdomain2GlobalMappings_1d(const MPI_Comm& PROTEUS_COMM_WORLD, Mesh& mesh,
                                                       const int *elementOffsets_subdomain_owned,
                                                       const int *nodeOffsets_subdomain_owned,
                                                       const int *elementNumbering_subdomain2global,
                                                       const int *nodeNumbering_subdomain2global,
                                                       int& nDOF_all_processes,//total number of dofs in whole domain
                                                       int& nDOF_subdomain,//total number of dofs in sub-domain
                                                       int& max_dof_neighbors,//maximum number of neighbors for connectivity of dofs
                                                       int *offsets_subdomain_owned, //starting point of local dofs on each processor (nProcs+1)
                                                       int * subdomain_l2g, //local to global dof mapping on subdomain
                                                       int* subdomain2global,//subdomain dof to global (parallel) numbering
                                                       double * lagrangeNodesArray);//location of nodes corresponding to dofs

  extern int buildQuadraticSubdomain2GlobalMappings_2d(const MPI_Comm& PROTEUS_COMM_WORLD, Mesh& mesh,
                                                       const int *elementBoundaryOffsets_subdomain_owned,
                                                       const int *nodeOffsets_subdomain_owned,
                                                       const int *elementBoundaryNumbering_subdomain2global,
                                                       const int *nodeNumbering_subdomain2global,
                                                       int& nDOF_all_processes,//total number of dofs in whole domain
                                                       int& nDOF_subdomain,//total number of dofs in sub-domain
                                                       int& max_dof_neighbors,//maximum number of neighbors for connectivity of dofs
                                                       int *offsets_subdomain_owned, //starting point of local dofs on each processor (nProcs+1)
                                                       int *subdomain_l2g, //local to global dof mapping on subdomain
                                                       int *subdomain2global,//subdomain dof to global (parallel) numbering
                                                       double * lagrangeNodesArray);//location of nodes corresponding to dofs

  extern int buildQuadraticSubdomain2GlobalMappings_3d(const MPI_Comm& PROTEUS_COMM_WORLD, Mesh& mesh,
                                                       const int *edgeOffsets_subdomain_owned,
                                                       const int *nodeOffsets_subdomain_owned,
                                                       const int *edgeNumbering_subdomain2global,
                                                       const int *nodeNumbering_subdomain2global,
                                                       int& nDOF_all_processes,//total number of dofs in whole domain
                                                       int& nDOF_subdomain,//total number of dofs in sub-domain
                                                       int& max_dof_neighbors,//maximum number of neighbors for connectivity of dofs
                                                       int *offsets_subdomain_owned, //starting point of local dofs on each processor (nProcs+1)
                                                       int *subdomain_l2g, //local to global dof mapping on subdomain
                                                       int *subdomain2global,//subdomain dof to global (parallel) numbering
                                                       double * lagrangeNodesArray);//location of nodes corresponding to dofs

  extern int buildQuadraticCubeSubdomain2GlobalMappings_3d(const MPI_Comm& PROTEUS_COMM_WORLD, Mesh& mesh,
                                                           const int *edgeOffsets_subdomain_owned,
                                                           const int *nodeOffsets_subdomain_owned,
                                                           const int *edgeNumbering_subdomain2global,
                                                           const int *nodeNumbering_subdomain2global,
                                                           int& nDOF_all_processes,//total number of dofs in whole domain
                                                           int& nDOF_subdomain,//total number of dofs in sub-domain
                                                           int& max_dof_neighbors,//maximum number of neighbors for connectivity of dofs
                                                           int *offsets_subdomain_owned, //starting point of local dofs on each processor (nProcs+1)
                                                           int *subdomain_l2g, //local to global dof mapping on subdomain
                                                           int *subdomain2global,//subdomain dof to global (parallel) numbering
                                                           double * lagrangeNodesArray);//location of nodes corresponding to dofs

  extern int buildDiscontinuousGalerkinSubdomain2GlobalMappings(const MPI_Comm& PROTEUS_COMM_WORLD, Mesh& mesh,
                                                                const int *elementOffsets_subdomain_owned,
                                                                const int *elementNumbering_subdomain2global,
                                                                int nDOF_element,
                                                                int& nDOF_all_processes,//total number of dofs in whole domain
                                                                int& nDOF_subdomain,//total number of dofs in sub-domain
                                                                int& max_dof_neighbors,//maximum number of neighbors for connectivity of dofs
                                                                int *offsets_subdomain_owned, //starting point of local dofs on each processor (nProcs+1)
                                                                int * subdomain_l2g, //local to global dof mapping on subdomain
                                                                int* subdomain2global);//subdomain dof to global (parallel) numbering
}//proteus
#endif
