#ifndef MODELFACTORY_H
#define MODELFACTORY_H
#include <iostream>

//#define FULL_BUILD 1
#define NO_INSTANCE std::cout<<"Constructing model object from template class:"<<std::endl \
  <<"return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<" \
  <<nSpaceIn<<","                                                       \
  <<nDOF_mesh_trial_elementIn<<","                                      \
  <<nDOF_trial_elementIn<<","                                           \
  <<nDOF_test_elementIn<<">,"                                           \
  <<nSpaceIn<<","                                                       \
  <<nQuadraturePoints_elementIn<<","                                    \
  <<nDOF_mesh_trial_elementIn<<","                                      \
  <<nDOF_trial_elementIn<<","                                           \
  <<nDOF_test_elementIn<<","                                            \
  <<nQuadraturePoints_elementBoundaryIn<<">());"                        \
  <<std::endl<<std::flush

#ifdef FULL_BUILD
namespace proteus
{
  template<class Model_Base,
    template<class CompKernelType,
    int nSpace,
    int nQuadraturePoints_element,
    int nDOF_mesh_trial_element,
    int nDOF_trial_element,
    int nDOF_test_element,
    int nQuadraturePoints_elementBoundary>
    class ModelTemplate,
    template<int nSpace,
    int nDOF_mesh_trial_element,
    int nDOF_trial_element,
    int nDOF_test_element>
    class CompKernelTemplate>
    Model_Base* chooseAndAllocateDiscretization(int nSpaceIn,
                                                int nQuadraturePoints_elementIn,
                                                int nDOF_mesh_trial_elementIn,
                                                int nDOF_trial_elementIn,
                                                int nDOF_test_elementIn,
                                                int nQuadraturePoints_elementBoundaryIn,
                                                int CompKernelFlag)//0=Parametric
    {
      if (CompKernelFlag == 0)
        {
          if (nSpaceIn == 3) // 3D
            {
              if (nDOF_mesh_trial_elementIn == nDOF_trial_elementIn)//iso-parametric
                {
                  if (nDOF_mesh_trial_elementIn == 4) // P1 FE-space. Default nquad=5
                    {
                      if (nQuadraturePoints_elementIn == 5)
                        return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,4,4>,3,5,4,4,4,4>());
                      else if (nQuadraturePoints_elementIn == 4)
                        return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,4,4>,3,4,4,4,4,3>());
                      else if (nQuadraturePoints_elementIn == 15 && nQuadraturePoints_elementBoundaryIn == 7)
                        return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,4,4>,3,15,4,4,4,7>());
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }
                    }
                  else if (nDOF_mesh_trial_elementIn == 8) // Q1 FE-space. Default nquad=27
                    {
                      if (nQuadraturePoints_elementIn == 8)
                        if ( nQuadraturePoints_elementBoundaryIn == 4)
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,8,8,8>,3,8,8,8,8,4>());
                        else if ( nQuadraturePoints_elementBoundaryIn == 9)
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,8,8,8>,3,8,8,8,8,9>());
                        else
                          {
                            NO_INSTANCE;
                            abort();
                          }
                      else if (nQuadraturePoints_elementIn == 27)
                        if ( nQuadraturePoints_elementBoundaryIn == 4)
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,8,8,8>,3,27,8,8,8,4>());
                        else if ( nQuadraturePoints_elementBoundaryIn == 9)
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,8,8,8>,3,27,8,8,8,9>());
                        else
                          {
                            NO_INSTANCE;
                            abort();
                          }
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }
                    }
                  else if (nDOF_mesh_trial_elementIn == 10) // P2 FE space. Default nquad=15
                    {
                      if (nQuadraturePoints_elementIn == 15)
                        return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,10,10,10>,3,15,10,10,10,7>());
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }
                    }
                  else if (nDOF_mesh_trial_elementIn == 27) // Q2 FE space. Default nquad=27
                    {
                      if (nQuadraturePoints_elementIn == 8)
                        if ( nQuadraturePoints_elementBoundaryIn == 4)
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,27,27,27>,3,8,27,27,27,4>());
                        else if ( nQuadraturePoints_elementBoundaryIn == 9)
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,27,27,27>,3,8,27,27,27,9>());
                        else
                          {
                            NO_INSTANCE;
                            abort();
                          }
                      else if (nQuadraturePoints_elementIn == 27)
                        if ( nQuadraturePoints_elementBoundaryIn == 4)
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,27,27,27>,3,27,27,27,27,4>());
                        else if ( nQuadraturePoints_elementBoundaryIn == 9)
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,27,27,27>,3,27,27,27,27,9>());
                        else
                          {
                            NO_INSTANCE;
                            abort();
                          }
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }
                    }
                  else
                    {
                      NO_INSTANCE;
                      abort();
                    }
                }
              else if (nDOF_mesh_trial_elementIn == 4)//sub-parametric tets
                {
                  if (nDOF_trial_elementIn == 10)
                    {
                      if (nQuadraturePoints_elementIn == 4 && nQuadraturePoints_elementBoundaryIn == 3)
                        {
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,10,10>,3,4,4,10,10,3>());
                        }
                      else if (nQuadraturePoints_elementIn == 5 && nQuadraturePoints_elementBoundaryIn == 4)
                        {
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,10,10>,3,5,4,10,10,4>());
                        }
                      else if (nQuadraturePoints_elementIn == 14 && nQuadraturePoints_elementBoundaryIn == 6)
                        {
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,10,10>,3,14,4,10,10,6>());
                        }
                      else if (nQuadraturePoints_elementIn == 15 && nQuadraturePoints_elementBoundaryIn == 7)
                        {
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,10,10>,3,15,4,10,10,7>());
                        }
                      else if (nQuadraturePoints_elementIn == 24 && nQuadraturePoints_elementBoundaryIn == 12)
                        {
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,10,10>,3,24,4,10,10,12>());
                        }
                      else if (nQuadraturePoints_elementIn == 31 && nQuadraturePoints_elementBoundaryIn == 12)
                        {
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,10,10>,3,31,4,10,10,12>());
                        }
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }
                    }
                  else
                    {
                      NO_INSTANCE;
                      abort();
                    }
                }
              else if (nDOF_mesh_trial_elementIn == 8)//sub-parametric hexes
                {
                  if (nDOF_trial_elementIn == 27)
                    if (nQuadraturePoints_elementIn == 125)
                      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,8,27,27>,3,125,8,27,27,25>());
                }
              else
                {
                  NO_INSTANCE;
                  abort();
                }
            }
          else
            {
              NO_INSTANCE;
              abort();
            }
        }
      else
        {
          NO_INSTANCE;
          abort();
        }
      return NULL;
    }
  template<class Model_Base,
    template<class CompKernelType,
    int nSpace,
    int nQuadraturePoints_element,
    int nDOF_mesh_trial_element,
    int nDOF_trial_element,
    int nDOF_test_element,
    int nQuadraturePoints_elementBoundary>
    class ModelTemplate,
    template<int nSpace,
    int nDOF_mesh_trial_element,
    int nDOF_trial_element,
    int nDOF_test_element>
    class CompKernelTemplate>
    Model_Base* chooseAndAllocateDiscretization2D(int nSpaceIn,
                                                  int nQuadraturePoints_elementIn,
                                                  int nDOF_mesh_trial_elementIn,
                                                  int nDOF_trial_elementIn,
                                                  int nDOF_test_elementIn,
                                                  int nQuadraturePoints_elementBoundaryIn,
                                                  int CompKernelFlag)//0=Parametric
    {
      if (CompKernelFlag == 0)
        {
          if (nSpaceIn == 2) // 2D
            {
              if (nDOF_mesh_trial_elementIn == nDOF_trial_elementIn)//iso-parametric
                {
                  if (nDOF_mesh_trial_elementIn == 3) // P1 FE-space. Default nquad=4
                    {
                      if (nQuadraturePoints_elementIn == 1)
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 1)
                            return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,3,3>,2,1,3,3,3,1>());
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 3)
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 2)
                            return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,3,3>,2,3,3,3,3,2>());
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 4)
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 3)
                            return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,3,3>,2,4,3,3,3,3>());
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 6)
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 4)
                            return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,3,3>,2,6,3,3,3,4>());
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 7)
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 5)
                            return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,3,3>,2,7,3,3,3,5>());
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 112)//hk=0.25
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 5)
                            {
                              return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,3,3>,2,112,3,3,3,5>());
                            }
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 16)//hk=0.5
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 3)
                            {
                              return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,3,3>,2,16,3,3,3,3>());
                            }
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 36)//hk=0.3
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 3)
                            {
                              return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,3,3>,2,36,3,3,3,3>());
                            }
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 100)//hk=0.2
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 3)
                            {
                              return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,3,3>,2,100,3,3,3,3>());
                            }
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }
                    }
                  else  if(nDOF_mesh_trial_elementIn == 4)
                    {
                      if (nQuadraturePoints_elementIn == 4) // Q1 FE-space
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 2)
                            return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,4,4>,2,4,4,4,4,2>());
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }
                    }
                  else
                    {
                      NO_INSTANCE;
                      abort();
                    }
                }
              else if (nDOF_mesh_trial_elementIn == 3)
                {
                  if (nDOF_trial_elementIn == 6) // P2 FE-space. Default nquad=7
                    {
                      if (nQuadraturePoints_elementIn == 1)
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 1)
                            return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,6,6>,2,1,3,6,6,1>());
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 3)
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 2)
                            return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,6,6>,2,3,3,6,6,2>());
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 4)
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 3)
                            return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,6,6>,2,4,3,6,6,3>());
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 6)
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 4)
                            return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,6,6>,2,6,3,6,6,4>());
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 7)
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 5)
                            {
                              return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,6,6>,2,7,3,6,6,5>());
                            }
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 112)
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 5)
                            {
                              return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,6,6>,2,112,3,6,6,5>());
                            }
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }

                    }
                  else
                    {
                      NO_INSTANCE;
                      abort();
                    }
                }
              else
                {
		  //// high order hexes
		  if (nDOF_trial_elementIn == 9) //2nd order polynomials
		    {
		      if (nQuadraturePoints_elementIn == 25) //n_quad=(2*order+1)^2
			{
		          if (nQuadraturePoints_elementBoundaryIn == 5)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,9,9>,2,25,4,9,9,5>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
		        }
		      else if (nQuadraturePoints_elementIn == 36) //order=2. use n_quad=(2*order+2)^nd
			{
			  if (nQuadraturePoints_elementBoundaryIn == 6)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,9,9>,2,36,4,9,9,6>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else if (nQuadraturePoints_elementIn == 64) //order=2. use n_quad=(2*(order+1)+2)^nd
			{
			  if (nQuadraturePoints_elementBoundaryIn == 8)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,9,9>,2,64,4,9,9,8>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else if (nQuadraturePoints_elementIn == 100)
			{
			  if (nQuadraturePoints_elementBoundaryIn == 10)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,9,9>,2,100,4,9,9,10>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else if (nQuadraturePoints_elementIn == 144)
			{
			  if (nQuadraturePoints_elementBoundaryIn == 12)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,9,9>,2,144,4,9,9,12>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else if (nQuadraturePoints_elementIn == 196)
			{
			  if (nQuadraturePoints_elementBoundaryIn == 14)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,9,9>,2,196,4,9,9,14>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else if (nQuadraturePoints_elementIn == 256)
			{
			  if (nQuadraturePoints_elementBoundaryIn == 16)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,9,9>,2,256,4,9,9,16>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else if (nQuadraturePoints_elementIn == 324)
			{
			  if (nQuadraturePoints_elementBoundaryIn == 18)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,9,9>,2,324,4,9,9,18>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else if (nQuadraturePoints_elementIn == 400)
			{
			  if (nQuadraturePoints_elementBoundaryIn == 20)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,9,9>,2,400,4,9,9,20>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else //other quad rules
			{
			  NO_INSTANCE;
			  abort();
			}
		    }
		  else if (nDOF_trial_elementIn == 16) //3rd order polynomials
		    {
		      if (nQuadraturePoints_elementIn == 64) //n_quad=64 
			{
		          if (nQuadraturePoints_elementBoundaryIn == 8)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,16,16>,2,64,4,16,16,8>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else
			{
			  NO_INSTANCE;
			  abort();
			}
		    }
		  else //higher-order polynomials
		    {
		      NO_INSTANCE;
		      abort();
		    }
                }
            }
        }
      else
        {
          NO_INSTANCE;
          abort();
        }
      return NULL;
    }
}
#else
// QUICK BUILD
namespace proteus
{
  template<class Model_Base,
    template<class CompKernelType,
    int nSpace,
    int nQuadraturePoints_element,
    int nDOF_mesh_trial_element,
    int nDOF_trial_element,
    int nDOF_test_element,
    int nQuadraturePoints_elementBoundary>
    class ModelTemplate,
    template<int nSpace,
    int nDOF_mesh_trial_element,
    int nDOF_trial_element,
    int nDOF_test_element>
    class CompKernelTemplate>
    Model_Base* chooseAndAllocateDiscretization(int nSpaceIn,
                                                int nQuadraturePoints_elementIn,
                                                int nDOF_mesh_trial_elementIn,
                                                int nDOF_trial_elementIn,
                                                int nDOF_test_elementIn,
                                                int nQuadraturePoints_elementBoundaryIn,
                                                int CompKernelFlag)//0=Parametric
    {
      if (CompKernelFlag == 0)
        {
          if (nSpaceIn == 3) // 3D
            {
              if (nDOF_mesh_trial_elementIn == nDOF_trial_elementIn)//iso-parametric
                {
                  if (nDOF_mesh_trial_elementIn == 4) // P1 FE-space.
                    {
		      if (nQuadraturePoints_elementIn == 15 && nQuadraturePoints_elementBoundaryIn == 7)
			return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,4,4>,3,15,4,4,4,7>());//3D, for pressure on p1 while vel on p2
                      else if (nQuadraturePoints_elementIn == 5)
                        return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,4,4>,3,5,4,4,4,4>());
                      else if (nQuadraturePoints_elementIn == 4)
                        return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,4,4>,3,4,4,4,4,3>());//added to pass: tests/griffiths_lane_6/test_griffiths_lane6.py
		      else if (nQuadraturePoints_elementIn == 24 && nQuadraturePoints_elementBoundaryIn == 12)
			return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,4,4>,3,24,4,4,4,12>());
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }
                    }
                  else if (nDOF_mesh_trial_elementIn == 8) // Q1 FE-space
                    {
                      if (nQuadraturePoints_elementIn == 27)
                        if ( nQuadraturePoints_elementBoundaryIn == 9)
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,8,8,8>,3,27,8,8,8,9>());
                        else
                          {
                            NO_INSTANCE;
                            abort();
                          }
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }
                    }
                  else if (nDOF_mesh_trial_elementIn == 10) // P2 FE space
                    {
                      if (nQuadraturePoints_elementIn == 15)
                        return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,10,10,10>,3,15,10,10,10,7>());
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }
                    }
                  else if (nDOF_mesh_trial_elementIn == 27) // Q2 FE space
                    {
                      if (nQuadraturePoints_elementIn == 27)
                        if ( nQuadraturePoints_elementBoundaryIn == 9)
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,27,27,27>,3,27,27,27,27,9>());
                        else
                          {
                            NO_INSTANCE;
                            abort();
                          }
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }
                    }
                  else
                    {
                      NO_INSTANCE;
                      abort();
                    }
                }
              else if (nDOF_mesh_trial_elementIn == 4)//sub-parametric tets
                {
                  if (nDOF_trial_elementIn == 10)
                    {
		      if (nQuadraturePoints_elementIn == 15 && nQuadraturePoints_elementBoundaryIn == 7) //3D, for velocity on p2 while pressure on p1
			{
			  return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,10,10>,3,15,4,10,10,7>());
			}
                      else if (nQuadraturePoints_elementIn == 24 && nQuadraturePoints_elementBoundaryIn == 12)
                        {
                          return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,10,10>,3,24,4,10,10,12>());//added to pass: tests/griffiths_lane_6/test_griffiths_lane6.py
                        }
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }
                    }
                  else
                    {
                      NO_INSTANCE;
                      abort();
                    }
                }
              else
                {
                  NO_INSTANCE;
                  abort();
                }
            }
          else
            {
              NO_INSTANCE;
              abort();
            }
        }
      else
        {
          NO_INSTANCE;
          abort();
        }
      return NULL;
    }
  template<class Model_Base,
    template<class CompKernelType,
    int nSpace,
    int nQuadraturePoints_element,
    int nDOF_mesh_trial_element,
    int nDOF_trial_element,
    int nDOF_test_element,
    int nQuadraturePoints_elementBoundary>
    class ModelTemplate,
    template<int nSpace,
    int nDOF_mesh_trial_element,
    int nDOF_trial_element,
    int nDOF_test_element>
    class CompKernelTemplate>
    Model_Base* chooseAndAllocateDiscretization2D(int nSpaceIn,
                                                  int nQuadraturePoints_elementIn,
                                                  int nDOF_mesh_trial_elementIn,
                                                  int nDOF_trial_elementIn,
                                                  int nDOF_test_elementIn,
                                                  int nQuadraturePoints_elementBoundaryIn,
                                                  int CompKernelFlag)//0=Parametric
    {
      if (CompKernelFlag == 0)
        {
          if (nSpaceIn == 2) // 2D
            {
              if (nDOF_mesh_trial_elementIn == nDOF_trial_elementIn)//iso-parametric
                {
                  if (nDOF_mesh_trial_elementIn == 3) // P1 FE-space
                    {
                      if (nQuadraturePoints_elementIn == 4)
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 3)
                            return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,3,3>,2,4,3,3,3,3>());
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 6)
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 4)
                            return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,3,3>,2,6,3,3,3,4>());//added to pass: tests/solver_tests_mprans/test_bochev_pressure_stabilization.py
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 7)
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 5)
                            return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,3,3>,2,7,3,3,3,5>());//added to pass: tests/ProjScheme_with_EV
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else if (nQuadraturePoints_elementIn == 100)//hk=0.2
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 3)
                            {
                              return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,3,3>,2,100,3,3,3,3>());
                            }
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }
                    }
                  else  if(nDOF_mesh_trial_elementIn == 4) //Hexes
                    {
                      if (nQuadraturePoints_elementIn == 4) // Q1 FE-space
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 2)
                            return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,4,4>,2,4,4,4,4,2>());
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
		      else if (nQuadraturePoints_elementIn == 9) // Q1 FE-space
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 3)
                            return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,4,4>,2,9,4,4,4,3>());
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }
                    }
                  else
                    {
                      NO_INSTANCE;
                      abort();
                    }
                }
              else if (nDOF_mesh_trial_elementIn == 3)
                {
                  if (nDOF_trial_elementIn == 6) // P2 FE-space
                    {
                      if (nQuadraturePoints_elementIn == 7)
                        {
                          if (nQuadraturePoints_elementBoundaryIn == 5)
                            {
                              return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,3,6,6>,2,7,3,6,6,5>());
                            }
                          else
                            {
                              NO_INSTANCE;
                              abort();
                            }
                        }
                      else
                        {
                          NO_INSTANCE;
                          abort();
                        }
                    }
                  else
                    {
                      NO_INSTANCE;
                      abort();
                    }
                }
              else
                {
		  //// high order hexes
		  if (nDOF_trial_elementIn == 9) //2nd order polynomials
		    {
		      if (nQuadraturePoints_elementIn == 16) //order=2 with 16 quad points
			{
			  if (nQuadraturePoints_elementBoundaryIn == 4)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,9,9>,2,16,4,9,9,4>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else if (nQuadraturePoints_elementIn == 36) //order=2. use n_quad=(2*order+2)^nd
			{
			  if (nQuadraturePoints_elementBoundaryIn == 6)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,9,9>,2,36,4,9,9,6>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else if (nQuadraturePoints_elementIn == 64) //order=2. use n_quad=(2*(order+1)+2)^nd
			{
			  if (nQuadraturePoints_elementBoundaryIn == 8)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,9,9>,2,64,4,9,9,8>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else //other quad rules
			{
			  NO_INSTANCE;
			  abort();
			}
		    }
		  else if (nDOF_trial_elementIn == 16) //3rd order polynomials
		    {
		      if (nQuadraturePoints_elementIn == 36) //n_quad=36
			{
		          if (nQuadraturePoints_elementBoundaryIn ==6)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,16,16>,2,36,4,16,16,6>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else if (nQuadraturePoints_elementIn == 64) //n_quad=64 
			{
		          if (nQuadraturePoints_elementBoundaryIn == 8)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,16,16>,2,64,4,16,16,8>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else if (nQuadraturePoints_elementIn == 81) //n_quad=81
			{
		          if (nQuadraturePoints_elementBoundaryIn == 9)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,16,16>,2,81,4,16,16,9>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}		      
		      else
			{
			  NO_INSTANCE;
			  abort();
			}
		    }
		  else if (nDOF_trial_elementIn == 25) //4th order polynomials
		    {
		      if (nQuadraturePoints_elementIn == 64) //n_quad=64
			{
		          if (nQuadraturePoints_elementBoundaryIn ==8)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,25,25>,2,64,4,25,25,8>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else if (nQuadraturePoints_elementIn == 144) //n_quad=144 
			{
		          if (nQuadraturePoints_elementBoundaryIn == 12)
			    {
			      return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<2,4,25,25>,2,144,4,25,25,12>());
			    }
			  else
			    {
			      NO_INSTANCE;
			      abort();
			    }
			}
		      else
			{
			  NO_INSTANCE;
			  abort();
			}
		    }
		  else //higher-order polynomials
		    {
		      NO_INSTANCE;
		      abort();
		    }
                }
            }
        }
      else
        {
          NO_INSTANCE;
          abort();
        }
      return NULL;
    }
}
#endif
#endif
