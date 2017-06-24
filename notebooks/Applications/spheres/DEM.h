#ifndef DEM_H
#define DEM_H

#include <cmath>
#include <string>
#include <iostream>

#include "OpenMP.hpp"
#include "FileOutput.hpp"

#include "SetupGeometry.hpp"
#include "SphereCollection.hpp"
#include "Domain.hpp"
#include "ViscoElasticContact.hpp"
#include "Boundary.hpp"


class cppDEM
{
public:
  double timestep,total_time, current_time=0.0;
  int N_total,N_output,N_search,iterations = 0;

  std::shared_ptr<ContactModel> particle_contact_model, boundary_contact_model;
  std::string out_prefix;
  BoundaryConditions bcs;
  Domain spheres_search_grid;
  SphereCollection spheres;
  int get_nSolids()
  {
    return spheres.size();
  }
  cppDEM(std::string config_in)
  {
    std::cout<<"ctor"<<config_in<<std::endl;
    
    //YAML::Node config = YAML::LoadFile(config_in);
    std::cout<<"file loaded"<<config_in<<std::endl;

    SetupGeometry::RetrieveTimeParametersFromConfig(config_in, 
						    timestep, total_time, 
						    N_total, N_output, N_search);
    std::cout<<"retrieved time"<<config_in<<std::endl;

    /*
      Setup box geometry
    */
    SetupGeometry::SetupGeometry(config_in, 
				 bcs, spheres_search_grid, 
				 particle_contact_model, 
				 boundary_contact_model);
    std::cout<<"setup geometry"<<config_in<<std::endl;

    /*
    Initialize sphere collection
    */
    std::cout << "Initializing sphere collection.\n" << std::endl;
    // Create sphere object of SphereCollection class.
    SetupGeometry::SetupSpheres(config_in, spheres);
    // Add search grid and boundary conditions to collection.
    spheres.BindDomain(spheres_search_grid);
    spheres.BindBoundaryConditions(bcs);
    // Set contact parameters of particles.
    spheres.CreateDefaultStateList(particle_contact_model);
    YAML::Node config = YAML::LoadFile(config_in);
    out_prefix=config["Output prefix"].as<string>();
  }
  void step(double* force, double* torque, double dt)
  {
    double next_time = current_time+dt;
    std::cout <<"DEM step: ("<< current_time<<' ,'<<next_time<<")"<<std::endl;
    while (current_time < next_time)
    {
      // Create boundary pair lists for the spheres every N_search iterations.
      if (iterations%N_search==0)
      {
        spheres.Search(4.0*spheres[0]->GetDilationRadius());
	spheres.GenerateBoundaryPairList();
      }
      // Calculate intial force values - needed for leap-frog integration.
      if (iterations==0)
      {
	spheres.Contact();
	spheres.ContactWithBoundaries(boundary_contact_model);
      }
      // Move the particles.
      spheres.Move(timestep, boundary_contact_model);
      // Output messages, save .vtp files, calculate new timestep.
      if (iterations%N_output==0)
      {
	std::cout << "Iteration: " << iterations  << std::endl;
	WriteTimeStepVTP(out_prefix, iterations, spheres);
      }
      // Increment the iteration counter.
      current_time += timestep;
      ++iterations;
    }
  }
  double hx(double* x, double t)
  {
    return 0.0;
  }
  double hy(double* x, double t)
  {
    return 0.0;
  }
  double hz(double* x, double t)
  {
    return 0.0;
  }
  void set_center_array(double* center)
  {
    for (auto i=0; i<spheres.size(); ++i)
    {
      center[i*3+0]=spheres[i]->GetX();
      center[i*3+1]=spheres[i]->GetY();
      center[i*3+2]=spheres[i]->GetZ();
    }
  }
  void set_radius_array(double* radius)
  {
    for (auto i=0; i<spheres.size(); ++i)
    {
      radius[i]=spheres[i]->GetDilationRadius();
    }
  }
  
};

cppDEM* newDEM(const char* config_in)
{
  std::cout<<"here"<<config_in<<std::endl;
  return new cppDEM(config_in);
  std::cout<<"out"<<config_in<<std::endl;
}
#endif
