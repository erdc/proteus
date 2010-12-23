#ifndef MODELFACTORY_H
#define MODELFACTORY_H

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
	     std::cout<<"Constructing model object from template class:"<<std::endl
		      <<"return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<"
		      <<nSpaceIn<<","
		      <<nDOF_mesh_trial_elementIn<<","
		      <<nDOF_trial_elementIn<<","
		      <<nDOF_test_elementIn<<">,"
		      <<nSpaceIn<<","
		      <<nQuadraturePoints_elementIn<<","
		      <<nDOF_mesh_trial_elementIn<<","
		      <<nDOF_trial_elementIn<<","
		      <<nDOF_test_elementIn<<","
		      <<nQuadraturePoints_elementBoundaryIn<<">());"
		      <<std::endl<<std::flush;
	     if (CompKernelFlag == 0)
	       {
		 if (nSpaceIn == 3)
		   {
		     if (nDOF_mesh_trial_elementIn == nDOF_trial_elementIn)//iso-parametric
		       {
			 if (nDOF_mesh_trial_elementIn == 4)
			   {
			     if (nQuadraturePoints_elementIn == 5)
			       return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,4,4>,3,5,4,4,4,4>());
			     else if (nQuadraturePoints_elementIn == 4)
			       return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,4,4>,3,4,4,4,4,3>());
			     else
			       abort();
			   }
			 else if (nDOF_mesh_trial_elementIn == 8)
			   {
			     if (nQuadraturePoints_elementIn == 8)
			       if ( nQuadraturePoints_elementBoundaryIn == 4)
			         return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,8,8,8>,3,8,8,8,8,4>());
			       else if ( nQuadraturePoints_elementBoundaryIn == 9)
			         return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,8,8,8>,3,8,8,8,8,9>());
			       else  
				 abort();
			     else if (nQuadraturePoints_elementIn == 27)
			       if ( nQuadraturePoints_elementBoundaryIn == 4)
			         return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,8,8,8>,3,27,8,8,8,4>());
			       else if ( nQuadraturePoints_elementBoundaryIn == 9)			     			     
			         return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,8,8,8>,3,27,8,8,8,9>());
			       else
			         abort();	 
			     else			   
			   
			     abort();
			   }
			 else if (nDOF_mesh_trial_elementIn == 10)
			   {
			     if (nQuadraturePoints_elementIn == 15)
			       return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,10,10,10>,3,15,10,10,10,7>());
			     else
			       abort();
			   }
			 else if (nDOF_mesh_trial_elementIn == 27)
			   {
			     if (nQuadraturePoints_elementIn == 8)
			       if ( nQuadraturePoints_elementBoundaryIn == 4)
			         return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,27,27,27>,3,8,27,27,27,4>());
			       else if ( nQuadraturePoints_elementBoundaryIn == 9)
			         return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,27,27,27>,3,8,27,27,27,9>());
			       else  
				 abort();
			     else if (nQuadraturePoints_elementIn == 27)
			       if ( nQuadraturePoints_elementBoundaryIn == 4)
			         return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,27,27,27>,3,27,27,27,27,4>());
			       else if ( nQuadraturePoints_elementBoundaryIn == 9)			     			     
			         return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,27,27,27>,3,27,27,27,27,9>());
			       else
			         abort();	 
			     else
			       abort();
			   }
			 else
			   abort();
		       }
		     else if (nDOF_mesh_trial_elementIn == 4)//sub-parametric tets
		       {
			 if (nDOF_trial_elementIn == 10)
			   return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,10,10>,3,15,4,10,10,7>());
			 else
			   abort();
		       }
		     else if (nDOF_mesh_trial_elementIn == 8)//sub-parametric hexes
		       {
			 if (nDOF_trial_elementIn == 27)
			   if (nQuadraturePoints_elementIn == 125)
			     return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,8,27,27>,3,125,8,27,27,25>());
		       }
		     else
		       abort();
		   }
		 else
		   abort();
	       }
	     else
	       abort();	       
	     return NULL;   
	   }	  		      
}
#endif
