# A type of -*- python -*- file
import numpy as np
cimport numpy as np
cimport postprocessing as pp
    
def postProcessRT0velocityFromP1nc(&nFreeDOF_element,
                                   &freeLocal_element,
                                   &detJ,
                                   &sqrt_det_g,
                                   &n,
                                   &elementBarycenters,
                                   &quad_a,
                                   &quad_f,
                                   &w_dV_r,
                                   &u,
                                   &gradu,
                                   &a,
                                   &f,
                                   &r,
                                   &rt0vdofs,
                                   &w_dV_m=None,
                                   &mt=None):
    if mt is not None and w_dV_m is not None:
        pp.postProcessRT0velocityFromP1nc(gradu.shape[0],#/*nElements_global*/
                                          gradu.shape[1],#/*nQuadraturePoints_element*/
                                          w_dV_r.shape[2],#/*nDOF_test_element*/
                                          n.shape[1],    #/*nElementBoundaries_element*/
                                          n.shape[2],    #/*nQuadraturePoints_elementBoundary*/
                                          gradu.shape[2], #/*spaceDim*/
                                          <int*>(nFreeDOF_element.data),
                                          <int*>(freeLocal_element.data),
                                          <double*>(detJ.data),
                                          <double*>(sqrt_det_g.data),
                                          <double*>(n.data),
                                          <double*>(elementBarycenters.data),
                                          <double*>(quad_a.data),
                                          <double*>(quad_f.data),
                                          <double*>(w_dV_r.data),
                                          <double*>(w_dV_m.data),
                                          <double*>(u.data),
                                          <double*>(gradu.data),
                                          <double*>(a.data),
                                          <double*>(f.data),
                                          <double*>(r.data),
                                          <double*>(mt.data),
                                          <double*>(rt0vdofs.data))
    else:
        pp.postProcessRT0velocityFromP1ncNoMass(gradu.shape[0],#/*nElements_global*/
                                                gradu.shape[1],#/*nQuadraturePoints_element*/
                                                w_dV_r.shape[2],#/*nDOF_test_element*/
                                                n.shape[1],    #/*nElementBoundaries_element*/
                                                n.shape[2],    #/*nQuadraturePoints_elementBoundary*/
                                                gradu.shape[2], #/*spaceDim*/
                                                <int*>(nFreeDOF_element.data),
                                                <int*>(freeLocal_element.data),
                                                <double*>(detJ.data),
                                                <double*>(sqrt_det_g.data),
                                                <double*>(n.data),
                                                <double*>(elementBarycenters.data),
                                                <double*>(quad_a.data),
                                                <double*>(quad_f.data),
                                                <double*>(w_dV_r.data),
                                                <double*>(u.data),
                                                <double*>(gradu.data),
                                                <double*>(a.data),
                                                <double*>(f.data),
                                                <double*>(r.data),
                                                <double*>(rt0vdofs.data))

def postProcessRT0velocityFromP1nc_sd(&rowptr,
                                      &colind,
                                      &nFreeDOF_element,
                                      &freeLocal_element,
                                      &detJ,
                                      &sqrt_det_g,
                                      &n,
                                      &elementBarycenters,
                                      &quad_a,
                                      &quad_f,
                                      &w_dV_r,
                                      &u,
                                      &gradu,
                                      &a,
                                      &f,
                                      &r,
                                      &rt0vdofs,
                                      &w_dV_m=None,
                                      &mt=None):
    if mt is None None and w_dV_m is not None:
        pp.postProcessRT0velocityFromP1nc_sd(gradu.shape[0],#/*nElements_global*/
                                             gradu.shape[1],#/*nQuadraturePoints_element*/
                                             w_dV_r.shape[2],#/*nDOF_test_element*/
                                             n.shape[1],    #/*nElementBoundaries_element*/
                                             n.shape[2],    #/*nQuadraturePoints_elementBoundary*/
                                             gradu.shape[2], #/*spaceDim*/
                                             <int*>(rowptr.data),
                                             <int*>(colind.data),
                                             <int*>(nFreeDOF_element.data),
                                             <int*>(freeLocal_element.data),
                                             <double*>(detJ.data),
                                             <double*>(sqrt_det_g.data),
                                             <double*>(n.data),
                                             <double*>(elementBarycenters.data),
                                             <double*>(quad_a.data),
                                             <double*>(quad_f.data),
                                             <double*>(w_dV_r.data),
                                             <double*>(w_dV_m.data),
                                             <double*>(u.data),
                                             <double*>(gradu.data),
                                             <double*>(a.data),
                                             <double*>(f.data),
                                             <double*>(r.data),
                                             <double*>(mt.data),
                                             <double*>(rt0vdofs.data))
    else:
        pp.postProcessRT0velocityFromP1ncNoMass_sd(gradu.shape[0],#/*nElements_global*/
                                                   gradu.shape[1],#/*nQuadraturePoints_element*/
                                                   w_dV_r.shape[2],#/*nDOF_test_element*/
                                                   n.shape[1],    #/*nElementBoundaries_element*/
                                                   n.shape[2],    #/*nQuadraturePoints_elementBoundary*/
                                                   gradu.shape[2], #/*spaceDim*/
                                                   <int*>(rowptr.data),
                                                   <int*>(colind.data),
                                                   <int*>(nFreeDOF_element.data),
                                                   <int*>(freeLocal_element.data),
                                                   <double*>(detJ.data),
                                                   <double*>(sqrt_det_g.data),
                                                   <double*>(n.data),
                                                   <double*>(elementBarycenters.data),
                                                   <double*>(quad_a.data),
                                                   <double*>(quad_f.data),
                                                   <double*>(w_dV_r.data),
                                                   <double*>(u.data),
                                                   <double*>(gradu.data),
                                                   <double*>(a.data),
                                                   <double*>(f.data),
                                                   <double*>(r.data),
                                                   <double*>(rt0vdofs.data))

def updateRT0velocityWithAveragedPotentialP1nc(&detJ,
                                               &quad_a,
                                               &phi,
                                               &gradphi,
                                               &a,
                                               &rt0vdofs_element):
    pp.updateRT0velocityWithAveragedPotentialP1nc(a.shape[0],
                                             a.shape[1],
                                             gradphi.shape[-1],
                                             <double*>(detJ.data),
                                             <double*>(quad_a.data),
                                             <double*>(phi.data),
                                             <double*>(gradphi.data),
                                             <double*>(a.data), 
                                             <double*>(rt0vdofs_element.data))

def updateRT0velocityWithAveragedPotentialP1nc_sd(&rowptr,
                                                  &colind,
                                                  &detJ,
                                                  &quad_a,
                                                  &phi,
                                                  &gradphi,
                                                  &a,
                                                  &rt0vdofs_element):
    pp.updateRT0velocityWithAveragedPotentialP1nc_sd(a.shape[0],
                                                     a.shape[1],
                                                     gradphi.shape[-1],
                                                     <int*>(rowptr.data),
                                                     <int*>(colind.data),
                                                     <double*>(detJ.data),
                                                     <double*>(quad_a.data),
                                                     <double*>(phi.data),
                                                     <double*>(gradphi.data),
                                                     <double*>(a.data), 
                                                     <double*>(rt0vdofs_element.data))

def getElementRT0velocityValues(&x_element,
                                &rt0vdofs_element,
                                &v_element):
    pp.getElementRT0velocityValues(v_element.shape[0],
                                   v_element[1],
                                   v_element.shape[2],
                                   <double*>(x_element.data),
                                   <double*>(rt0vdofs_element.data),
                                   <double*><v_element.data))

def getElementBoundaryRT0velocityValues(&x_elementBoundary,
                                        &rt0vdofs_element,
                                        &v_elementBoundary):
    pp.getElementBoundaryRT0velocityValues(v_elementBoundary.shape[0],
                                           v_elementBoundary.shape[1],
                                           v_elementBoundary.shape[2],
                                           v_elementBoundary.shape[3],
                                           <double*>(x_elementBoundary.data),
                                           <double*>(rt0vdofs_element.data),
                                           <double*>(v_elementBoundary.data))

def getGlobalElementBoundaryRT0velocityValues(&elementBoundaryElementsArray,
                                              &x_elementBoundary_global,
                                              &rt0vdofs_element,
                                              &v_elementBoundary_global):
    pp.getGlobalElementBoundaryRT0velocityValues(v_elementBoundary_global.shape[0],
                                                 v_elementBoundary_global.shape[1],
                                                 v_elementBoundary_global.shape[2],
                                                 <int*>(elementBoundaryElementsArray.data),
                                                 <double*>(x_elementBoundary_global.data),
                                                 <double*>(rt0vdofs_element.data),
                                                 <double*>(v_elementBoundary_global.data))
    
def getGlobalExteriorElementBoundaryRT0velocityValues(&exteriorElementBoundariesArray,
                                                      &elementBoundaryElementsArray,
                                                      &x_elementBoundary_global,
                                                      &rt0vdofs_element,
                                                      &v_elementBoundary_global):
    pp.getGlobalExteriorElementBoundaryRT0velocityValues(exteriorElementBoundariesArray.shape[0],
                                                         v_elementBoundary_global.shape[1],
                                                         v_elementBoundary_global.shape[2],
                                                         <int*>(elementBoundaryElementsArray.data),
                                                         <int*>(exteriorElementBoundariesArray.data),
                                                         <double*>(x_elementBoundary_global.data),
                                                         <double*>(rt0vdofs_element.data),
                                                         <double*>(v_elementBoundary_global.data))

def postProcessRT0potentialFromP1nc(&uQuadratureWeights_element,
                                    &elementBarycenters,
                                    &aElementQuadWeights,
                                    &detJ,
                                    &uQuadratureWeights_elementBoundary,
                                    &x,
                                    &u,
                                    &gradu,
                                    &x_elementBoundary,
                                    &u_elementBoundary,
                                    &n,
                                    &a,
                                    &f,
                                    &r,
                                    &rt0vdofs,
                                    &rt0potential):
    pp.postProcessRT0potentialFromP1nc(u.shape[0],
                                       u.shape[1],
                                       n.shape[1],
                                       n.shape[2],
                                       n.shape[3],
                                       <double*>(uQuadratureWeights_element),
                                       <double*>(elementBarycenters.data),
                                       <double*>(aElementQuadWeights.data),
                                       <double*>(detJ.data),
                                       <double*>(uQuadratureWeights_element.data),
                                       <double*>(x.data),
                                       <double*>(u.data),
                                       <double*>(gradu.data),
                                       <double*>(x_elementBoundary.data),
                                       <double*>(u_elementBoundary.data),
                                       <double*>(n.data),
                                       <double*>(a.data),
                                       <double*>(f.data),
                                       <double*>(r.data),
                                       <double*>(rt0vdofs.data),
                                       <double*>(rt0potential.data))

def postProcessRT0potentialFromP1nc_sd(&rowptr,
                                       &colind,
                                       &uQuadratureWeights_element,
                                       &elementBarycenters,
                                       &aElementQuadWeights,
                                       &detJ,
                                       &uQuadratureWeights_elementBoundary,
                                       &x,
                                       &u,
                                       &gradu,
                                       &x_elementBoundary,
                                       &u_elementBoundary,
                                       &n,
                                       &a,
                                       &f,
                                       &r,
                                       &rt0vdofs,
                                       &rt0potential):
    pp.postProcessRT0potentialFromP1nc_sd(u.shape[0],
                                          u.shape[1],
                                          n.shape[1],
                                          n.shape[2],
                                          n.shape[3],
                                          <int*>(rowptr.data),
                                          <int*>(colind.data),
                                          <double*>(uQuadratureWeights_element),
                                          <double*>(elementBarycenters.data),
                                          <double*>(aElementQuadWeights.data),
                                          <double*>(detJ.data),
                                          <double*>(uQuadratureWeights_element.data),
                                          <double*>(x.data),
                                          <double*>(u.data),
                                          <double*>(gradu.data),
                                          <double*>(x_elementBoundary.data),
                                          <double*>(u_elementBoundary.data),
                                          <double*>(n.data),
                                          <double*>(a.data),
                                          <double*>(f.data),
                                          <double*>(r.data),
                                          <double*>(rt0vdofs.data),
                                          <double*>(rt0potential.data))

def buildLocalBDM1projectionMatrices(&w_dS_f,
                                     &ebq_n,
                                     &ebq_v,
                                     &BDMmat_element):
    pp.buildLocalBDM1projectionMatrices(ebq_n.shape[0],
                                        ebq_n.shape[1],
                                        ebq_n.shape[2],
                                        ebq_n.shape[3],
                                        w_dS_f.shape[3],
                                        ebq_v.shape[3],
                                        BDMmat_element.shape[1],
                                        <double*>(w_dS_f.data),
                                        <double*>(ebq_n.data),
                                        <double*>(ebq_v.data),
                                        <double*>(BDMmat_element.data))

def buildLocalBDM2projectionMatrices(&degree,
                                     &w_dS_f,
                                     &ebq_n,
                                     &ebq_v,
                                     &q_basis_vals,
                                     &w_int_test_grads,
                                     &w_int_div_free,
                                     &piola_trial_fun,
                                     &edgeFlags,
                                     &BDMmat_element):
    pp.buildLocalBDM2projectionMatrices(<double*>(degree.data),
                                        ebq_n.shape[0],
                                        ebq_n.shape[1],
                                        ebq_n.shape[2],
                                        q_basis_vals.shape[1],
                                        ebq_n.shape[3],
                                        w_dS_f.shape[3],
                                        ebq_v.shape[3],
                                        w_int_test_grads.shape[2],
                                        BDMmat_element.shape[1],
                                        <double*>(edgeFlags.data),
                                        <double*>(w_dS_f.data),
                                        <double*>(ebq_n.data),
                                        <double*>(ebq_v.data),
                                        <double*>(BDMmat_element.data),
                                        <double*>(q_basis_vals.data),
                                        <double*>(w_int_test_grads.data),
                                        <double*>(w_int_div_free.data),
                                        <double*>(piola_trial_fun.data))

def factorLocalBDM1projectionMatrices(&BDMmat_element,
                                      &BDMmatPivots_element):
    pp.factorLocalBDM1projectionMatrices(BDMmat_element.shape[0],
                                         BDMmat_element.shape[1],
                                         <double*>(BDMmat_element.data),
                                         <int*>(BDMmatPivots_element.data))

def factorLocalBDM2projectionMatrices(&BDMmat_element,
                                      &BDMmatPivots_element):
    pp.factorLocalBDM2projectionMatrices(BDMmat_element.shape[0],
                                         BDMmat_element.shape[1],
                                         <double*>(BDMmat_element.data),
                                         <int*>(BDMmatPivots_element.data))

def solveLocalBDM1projection(&BDMmat_element,
                             &BDMmatPivots_element,
                             &w_dS_f,
                             &ebq_n,
                             &ebq_velocity,
                             &q_vdofs):
    pp.solveLocalBDM1projection(ebq_n.shape[0],
                                ebq_n.shape[1],
                                ebq_n.shape[2],
                                ebq_n.shape[3],
                                w_dS_f.shape[3],
                                BDMmat_element.shape[1],
                                <double*>(BDMmat_element.data),
                                <int*>(BDMmatPivots_element.data),
                                <double*>(w_dS_f.data),
                                <double*>(ebq_n.data),
                                <double*>(ebq_velocity.data),
                                <double*>(q_vdofs.data))

def solveLocalBDM2projection(&BDMmat_element,
                             &BDMmatPivots_element,
                             &w_dS_f,
                             &ebq_n,
                             &w_interior_gradients,
                             &ebq_velocity,
                             &q_velocity,
                             &q_vdofs):
    pp.solveLocalBDM2projection(ebq_n.shape[0],
                                ebq_n.shape[1],
                                ebq_n.shape[2],
                                ebq_n.shape[3],
                                w_dS_f.shape[3],
                                BDMmat_element.shape[1],
                                <double*>(BDMmat_element.data),
                                <int*>(BDMmatPivots_element.data),
                                <double*>(w_dS_f.data),
                                <double*>(ebq_n.data),
                                <double*>(w_interior_gradients.data),
                                <double*>(q_velocity.data),
                                <double*>(ebq_velocity.data),
                                <double*>(q_vdofs.data))

def buildBDM2rhs(&BDMmat_element,
                 &BDMmatPivots_element,
                 &w_dS_f,
                 &ebq_n,
                 &w_interior_gradients,
                 &w_interior_divfree,
                 &ebq_velocity,
                 &q_velocity,
                 &q_vdofs,
                 &edgeFlags):
    pp.buildBDM2rhs(ebq_n.shape[0],
                    ebq_n.shape[1],
                    ebq_n.shape[2],
                    q_velocity.shape[1],
                    ebq_n.shape[3],
                    w_dS_f.shape[3],
                    BDMmat_element.shape[1],
                    w_interior_gradients.shape[2],
                    <double*>(BDMmat_element.data),
                    <int*>(BDMmatPivots_element.data),
                    <double*>(edgeFlags.data),
                    <double*>(w_dS_f.data),
                    <double*>(ebq_n.data),
                    <double*>(w_interior_gradients.data),
                    <double*>(w_interior_divfree.data),
                    <double*>(ebq_velocity.data),
                    <double*>(q_velocity.data),               
                    <double*>(q_vdofs.data))

def solveLocalBDM1projectionFromFlux(&BDMmat_element,
                                     &BDMmatPivots_element,
                                     &elementBoundaryElementsArray,
                                     &elementBoundariesArray,
                                     &w_dS_f,
                                     &ebq_global_flux,
                                     &q_vdofs):
    pp.solveLocalBDM1projectionFromFlux(w_dS_f.shape[0],
                                        w_dS_f.shape[1],
                                        w_dS_f.shape[2],
                                        w_dS_f.shape[3],
                                        BDMmat_element.shape[1],
                                        <double*>(BDMmat_element.data),
                                        <int*>(BDMmatPivots_element.data),
                                        <int*>(elementBoundaryElementsArray.data),
                                        <int*>(elementBoundariesArray.data),
                                        <double*>(w_dS_f.data),
                                        <double*>(ebq_global_flux.data),
                                        <double*>(q_vdofs.data))

def getElementBDM1velocityValuesLagrangeRep(&q_v,
                                            &p1_vdofs,
                                            &q_velocity):
    cdef int nVDOFs_element = q_velocity.shape[2]*(q_velocity.shape[2]+1)
    
    if p1_vdofs.ndim > 1:
        assert nVDOFs_element == p1_vdofs.shape[1]
    pp.getElementBDM1velocityValuesLagrangeRep(q_velocity.shape[0],
                                            q_velocity.shape[1],
                                            q_velocity.shape[2],
                                            q_v.shape[2],
                                            nVDOFs_element,
                                            <double*>(q_v.data),
                                            <double*>(p1_vdofs.data),
                                            <double*>(q_velocity.data))

def getElementBDM2velocityValuesLagrangeRep(&q_v,
                                            &p1_vdofs,
                                            &q_velocity):
    if q_velocity.shape[2] == 2:
        #dimension of bdm2 in 2d
        nVDOFs_element = 12
    elif q_velocity.shape[2] == 3:
        #dimension of bdm2 in 3d
        nVDOFs_element = 30
        
    if p1_vdofs.ndim > 1:
        assert nVDOFs_element == p1_vdofs.shape[1]
        
    pp.getElementBDM2velocityValuesLagrangeRep(q_velocity.shape[0],
                                          q_velocity.shape[1],
                                          q_velocity.shape[2],
                                          q_v.shape[2],
                                          nVDOFs_element,
                                          <double*>(q_v.data),
                                          <double*>(p1_vdofs.data),
                                          <double*>(q_velocity.data))
    
def getElementLDGvelocityValuesLagrangeRep(&q_v,
                                           &vdofs,
                                           &q_velocity):
    cdef int nVDOF_element =vdofs.shape[1]
    cdef int nDOF_trial_element=q_v.shape[2]
    if  vdofs.ndim > 1:
        assert nVDOF_element == q_v.shape[2]*q_velocity.shape[2]
    pp.getElementLDGvelocityValuesLagrangeRep(q_velocity.shape[0],
                                              q_velocity.shape[1],
                                              q_velocity.shape[2],
                                              nDOF_trial_element,
                                              nVDOF_element,
                                              <double*>(q_v.data),
                                              <double*>(vdofs.data),
                                              <double*>(q_velocity.data))

def getGlobalExteriorElementBoundaryBDM1velocityValuesLagrangeRep(&elementBoundaryElementsArray,
                                                                  &exteriorElementBoundariesArray,
                                                                  &ebqe_v,
                                                                  &p1_vdofs,
                                                                  &ebqe_velocity):
    cdef int nVDOFs_element = ebqe_velocity.shape[2]*(ebqe_velocity.shape[2]+1)
    if p1_vdofs.ndim > 1:
        assert nVDOFs_element == p1_vdofs.shape[1]
    pp.getGlobalExteriorElementBoundaryBDM1velocityValuesLagrangeRep(ebqe_velocity.shape[0],
                                                                     ebqe_velocity.shape[1],
                                                                     ebqe_velocity.shape[2],
                                                                     ebqe_v.shape[2],
                                                                     nVDOFs_element,
                                                                     <int*>(elementBoundaryElementsArray.data),
                                                                     <int*>(exteriorElementBoundariesArray.data),
                                                                     <double*>(ebqe_v.data),
                                                                     <double*>(p1_vdofs.data),
                                                                     <double*>(ebqe_velocity.data))

def getGlobalElementBoundaryBDM1velocityValuesLagrangeRep(&elementBoundaryElementsArray,
                                                          &exteriorElementBoundariesArray,
                                                          &ebqe_v,
                                                          &p1_vdofs,
                                                          &ebq_global_velocity):
    cdef int nVDOFs_element = ebq_global_velocity.shape[2]*(ebq_global_velocity.shape[2]+1)
    if p1_vdofs.ndim > 1:
        assert nVDOFs_element == p1_vdofs.shape[1]
    pp.getGlobalExteriorElementBoundaryBDM1velocityValuesLagrangeRep(ebqe_v.shape[0],
                                                                     ebq_global_velocity.shape[1],
                                                                     ebq_global_velocity.shape[2],
                                                                     ebqe_v.shape[2],
                                                                     nVDOFs_element,
                                                                     <int*>(elementBoundaryElementsArray.data),
                                                                     <int*>(exteriorElementBoundariesArray.data),
                                                                     <double*>(ebqe_v.data),
                                                                     <double*>(p1_vdofs.data),
                                                                     <double*>(ebq_global_velocity.data))

def getElementBoundaryBDM1velocityValuesLagrangeRep(&elementBoundaryElementsArray,
                                                    &exteriorElementBoundariesArray,
                                                    &ebq_v,
                                                    &p1_vdofs,
                                                    &ebq_velocity):
    cdef int nVDOFs_element = ebq_velocity.shape[3]*(ebq_velocity.shape[3]+1)
    if p1_vdofs.ndim > 1:
        assert nVDOFs_element == p1_vdofs.shape[1]
        pp.getElementBoundaryBDM1velocityValuesLagrangeRep(ebq_velocity.shape[0],
                                                           ebq_velocity.shape[1],
                                                           ebq_velocity.shape[2],
                                                           ebq_velocity.shape[3],
                                                           ebq_v.shape[2],
                                                           nVDOFs_element,
                                                           <int*>(elementBoundaryElementsArray.data),
                                                           <int*>(exteriorElementBoundariesArray.data), #/*need the correct mapping here */
                                                           <double*>(ebq_v.data),
                                                           <double*>(p1_vdofs.data),
                                                           <double*>(ebq_velocity.data))

def projectElementBoundaryVelocityToRT0fluxRep(&elementBoundaryQuadratureWeights,
                                               &n,
                                               &v_elementBoundary,
                                               &rt0vdofs_element):
    pp.projectElementBoundaryVelocityToRT0fluxRep(n.shape[0],
                                                  n.shape[1],
                                                  n.shape[2],
                                                  n.shape[3],
                                                  <double*>(elementBoundaryQuadratureWeights.data),
                                                  <double*>(n.data),
                                                  <double*>(v_elementBoundary.data),
                                                  <double*>(rt0vdofs_element.data))

def projectElementBoundaryFluxToRT0fluxRep(&elementBoundaryElementsArray,
                                           &elementBoundariesArray,
                                           &elementBoundaryQuadratureWeights,
                                           &flux_elementBoundary,
                                           &rt0vdofs_element):
    pp.projectElementBoundaryFluxToRT0fluxRep(elementBoundaryQuadratureWeights.shape[0],
                                              elementBoundaryQuadratureWeights.shape[1],
                                              elementBoundaryQuadratureWeights.shape[2],
                                              rt0vdofs_element.shape[1],
                                              <int*>(elementBoundaryElementsArray.data),
                                              <int*>(elementBoundariesArray.data),
                                              <double*>(elementBoundaryQuadratureWeights.data),
                                              <double*>(flux_elementBoundary.data),
                                              <double*>(rt0vdofs_element.data))

def getElementRT0velocityValuesFluxRep(&nodeArray,
                                       &elementNodesArray,
                                       &abs_det_J,
                                       &x_element,
                                       &rt0vdofs_element,
                                       &v_element):
    pp.getElementRT0velocityValuesFluxRep(v_element.shape[0],
                                          elementNodesArray.shape[1],
                                          v_element[1],
                                          v_element.shape[2],
                                          abs_det_J.shape[1],
                                          <double*>(nodeArray.data),
                                          <int*>(elementNodesArray.data),
                                          <double*>(abs_det_J.data),
                                          <double*>(x_element.data),
                                          <double*>(rt0vdofs_element.data),
                                          <double*><v_element.data))

def getElementBoundaryRT0velocityValuesFluxRep(&nodeArray,
                                               &elementNodesArray,
                                               &abs_det_J,
                                               &x_elementBoundary,
                                               &rt0vdofs_element,
                                               &v_elementBoundary):
    pp.getElementBoundaryRT0velocityValuesFluxRep(v_elementBoundary.shape[0],
                                                  v_elementBoundary.shape[1],
                                                  v_elementBoundary.shape[2],
                                                  v_elementBoundary.shape[3],
                                                  abs_det_J.shape[1],
                                                  <double*>(nodeArray.data),
                                                  <int*>(elementNodesArray.data),
                                                  <double*>(abs_det_J.data),
                                                  <double*>(x_elementBoundary.data),
                                                  <double*>(rt0vdofs_element.data),
                                                  <double*>(v_elementBoundary.data))
    
def getGlobalElementBoundaryRT0velocityValuesFluxRep(&nodeArray,
                                                     &elementNodesArray,
                                                     &elementBoundaryElementsArray,
                                                     &abs_det_J,
                                                     &x_elementBoundary_global,
                                                     &rt0vdofs_element,
                                                     &v_elementBoundary_global):
    pp.getGlobalElementBoundaryRT0velocityValuesFluxRep(v_elementBoundary_global.shape[0],
                                                        v_elementBoundary_global.shape[1],
                                                        v_elementBoundary_global.shape[2],
                                                        abs_det_J.shape[1],
                                                        <double*>(nodeArray.data),
                                                        <int*>(elementNodesArray.data),
                                                        <int*>(elementBoundaryElementsArray.data),
                                                        <double*>(abs_det_J.data),
                                                        <double*>(x_elementBoundary_global.data),
                                                        <double*>(rt0vdofs_element.data),
                                                        <double*>(v_elementBoundary_global.data))

def getGlobalExteriorElementBoundaryRT0velocityValuesFluxRep(&nodeArray,
                                                             &elementNodesArray,
                                                             &elementBoundaryElementsArray,
                                                             &exteriorElementBoundariesArray,
                                                             &abs_det_J,
                                                             &x_ebqe,
                                                             &rt0vdofs_element,
                                                             &v_ebqe):
    pp.getGlobalExteriorElementBoundaryRT0velocityValuesFluxRep(v_ebqe.shape[0],
                                                                v_ebqe.shape[1],
                                                                v_ebqe.shape[2],
                                                                abs_det_J.shape[1],
                                                                <double*>(nodeArray.data),
                                                                <int*>(elementNodesArray.data),
                                                                <int*>(elementBoundaryElementsArray.data),
                                                                <int*>(exteriorElementBoundariesArray.data),
                                                                <double*>(abs_det_J.data),
                                                                <double*>(x_ebqe.data),
                                                                <double*>(rt0vdofs_element.data),
                                                                <double*>(v_ebqe.data))

def getRT0velocityValuesFluxRep_arbitraryElementMembership(&nodeArray,
                                                           &elementNodesArray,
                                                           &abs_det_J,
                                                           &x,
                                                           &element_locations,
                                                           &rt0vdofs_element,
                                                           &v_element):
    cdef int nPoints = x.size/x.shape[-1]
    pp.getRT0velocityValuesFluxRep_arbitraryElementMembership(elementNodesArray.shape[0],
                                                              elementNodesArray.shape[1],
                                                              nPoints,
                                                              v_element[1],
                                                              abs_det_J.shape[1],
                                                              <double*>(nodeArray.data),
                                                              <int*>(elementNodesArray.data),
                                                              <double*>(abs_det_J.data),
                                                              <double*>(x.data),
                                                              <int*>(element_locations.data),
                                                              <double*>(rt0vdofs_element.data),
                                                              <double*><v_element.data))

def fluxCorrectionVelocityUpdate(&interiorElementBoundaries,
                                 &exteriorElementBoundaries,
                                 &elementBoundaryElements,
                                 &elementBoundaryLocalElementBoundaries,
                                 &dS,
                                 &n,
                                 &fluxCorrection,
                                 &velocity,
                                 &velocity_element):
    pp.fluxCorrectionVelocityUpdate(dS.shape[0],
                                    velocity.shape[0],
                                    interiorElementBoundaries.shape[0],
                                    exteriorElementBoundaries.shape[0],
                                    dS.shape[1],
                                    dS.shape[2],
                                    velocity.shape[2],
                                    <int*>(interiorElementBoundaries.data),
                                    <int*>(exteriorElementBoundaries.data),
                                    <int*>(elementBoundaryElements.data),
                                    <int*>(elementBoundaryLocalElementBoundaries.data),
                                    <double*>(dS.data),
                                    <double*>(n.data),
                                    <double*>(fluxCorrection.data),
                                    <double*>(velocity.data),
                                    <double*>(velocity_element.data))

def computeFluxCorrectionPWC(&interiorElementBoundaries,
                             &exteriorElementBoundaries,
                             &elementBoundaryElements,
                             &pwcW,
                             &pwcV,
                             &fluxCorrection):
    pp.computeFluxCorrectionPWC(elementBoundaryElements.shape[0],
                                interiorElementBoundaries.shape[0],
                                exteriorElementBoundaries.shape[0],
                                <int*>(interiorElementBoundaries.data),
                                <int*>(exteriorElementBoundaries.data),
                                <int*>(elementBoundaryElements.data),
                                <double*>(pwcW.data),
                                <double*>(pwcV.data),
                                <double*>(fluxCorrection.data))

def sunWheelerGSsweep(&interiorElementBoundaries,
                      &exteriorElementBoundaries,
                      &elementBoundaryElements,
                      &elementBoundaryLocalElementBoundaries,
                      &dS,
                      &n,
                      &sqrt_det_g,
                      &alphaFactor,
                      &fluxCorrection,
                      &conservationResidual):
    pp.sunWheelerGSsweep(conservationResidual.shape[0],
                         elementBoundaryElements.shape[0],
                         interiorElementBoundaries.shape[0],
                         exteriorElementBoundaries.shape[0],
                         dS.shape[1],
                         dS.shape[2],
                         n.shape[2],
                         <int*>(interiorElementBoundaries.data),
                         <int*>(exteriorElementBoundaries.data),
                         <int*>(elementBoundaryElements.data),
                         <int*>(elementBoundaryLocalElementBoundaries.data),
                         <double*>(dS.data),
                         <double*>(n.data),
                         <double*>(sqrt_det_g.data),
                         <double*>(alphaFactor.data),
                         <double*>(fluxCorrection.data),
                         <double*>(conservationResidual.data))

#  put in data structure for node star solves now
cdef class NodeStarFactor:
    cdef pp.NodeStarFactor nsf
    def __cinit__(self,
                  &nElements_node,
                  &nodeStarElementsArray,
                  &nodeStarElementsNeighborsArray):
        cdef int rval=0
        rval = pp.nodeStar_init(nodeStarElementsArray.shape[0],
                                nodeStarElementsArray.shape[1],
                                nElements_node.shape[0],
                                <int*>(nElements_node.data),
                                <int*>(nodeStarElementsArray.data),
                                <int*>(nodeStarElementsNeighborsArray.data),
                                self.nsf.N,
                                self.nsf.subdomain_dim,
                                self.nsf.subdomain_L,
                                self.nsf.subdomain_R,
                                self.nsf.subdomain_U,
                                self.nsf.subdomain_pivots,
                                self.nsf.subdomain_column_pivots)
        assert rval == 0

    def __dealloc__(self):
        pp.nodeStar_free(self.nsf.N, 
                         self.nsf.subdomain_dim,
                         self.nsf.subdomain_L,
                         self.nsf.subdomain_R,
                         self.nsf.subdomain_U,
                         self.nsf.subdomain_psivots,
                         self.nsf.subdomain_column_pivots)
    
    def setU(self,
             double val):
        cdef int I
        cdef int i
        for I in range(nodeStarFactor.nsf.N):
            for i in range(nodeStarFactor.nsf.subdomain_dim[I]):
                nodeStarFactor.nsf.subdomain_U[I][i] = val
    
    def copyData(self,
                 NodeStarFactor other):
        pp.nodeStar_copy(other.nsf.N, 
                         other.nsf.subdomain_dim,
                         other.nsf.subdomain_L,
                         other.nsf.subdomain_R,
                         other.nsf.subdomain_U,
                         other.nsf.subdomain_pivots,
                         other.nsf.subdomain_column_pivots,
                         self.nsf.N, 
                         self.nsf.subdomain_dim,
                         self.nsf.subdomain_L,
                         self.nsf.subdomain_R,
                         self.nsf.subdomain_U,
                         self.nsf.subdomain_pivots,
                         self.nsf.subdomain_column_pivots)

def calculateConservationResidualPWL(&interiorElementBoundaries,
                                     &exteriorElementBoundaries,
                                     &elementBoundaryElements,
                                     &elementBoundaryLocalElementBoundaries,
                                     &elementNodes,
                                     &dofMapl2g,
                                     &nodeStarElements,
                                     &nodeStarElementNeighbors,
                                     &nElements_node,
                                     &fluxElementBoundaries,
                                     &elementResidual,
                                     &vAverage,
                                     &dx,
                                     &w,
                                     &n,
                                     NodeStarFactor nodeStarFactor,
                                     &conservationResidual,
                                     &vConservative,
                                     &vConservative_element):
    pp.calculateConservationResidualPWL(w.shape[0] 
                                        #/* nElements_global*/, 
                                        interiorElementBoundaries.shape[0] 
                                        #/*nInteriorElementBoundaries_global*/,
                                        exteriorElementBoundaries.shape[0]
                                        #/*nExteriorElementBoundaries_global*/,
                                        w.shape[1]
                                        #/*nElementBoundaries_element*/,
                                        w.shape[2]
                                        #/*nQuadraturePoints_elementBoundary*/,
                                        w.shape[3]
                                        #/*nNodes_element*/,
                                        dofMapl2g.shape[1]
                                        #/*nDOF_element*/,
                                        n.shape[2]
                                        #/*nSpace */,
                                        <int*>(interiorElementBoundaries.data),
                                        <int*>(exteriorElementBoundaries.data),
                                        <int*>(elementBoundaryElements.data),
                                        <int*>(elementBoundaryLocalElementBoundaries.data),
                                        <int*>(elementNodes.data),
                                        <int*>(dofMapl2g.data),
                                        <int*>(nodeStarElements.data),
                                        <int*>(nodeStarElementNeighbors.data),
                                        <int*>(nElements_node.data),
                                        <int*>(fluxElementBoundaries.data),
                                        <double*>(elementResidual.data),
                                        <double*>(vAverage.data),
                                        <double*>(dx.data),
                                        <double*>(w.data),
                                        <double*>(n.data),
                                        nodeStarFactor.nsf,
                                        <double*>(conservationResidual.data),
                                        <double*>(vConservative.data),
                                        <double*>(vConservative_element.data))

def calculateConservationJacobianPWL(&interiorElementBoundaries,
                                     &exteriorElementBoundaries,
                                     &elementBoundaryElements,
                                     &elementBoundaryLocalElementBoundaries,
                                     &elementNodes,
                                     &dofMapl2g,
                                     &dofStarElements,
                                     &dofStarElementNeighbors,
                                     &nElements_dof,
                                     &internalNodes,
                                     &fluxElementBoundaries,
                                     &fluxBoundaryNodes,
                                     &w,
                                     &n,
                                     &nodeStarFactor):
    pp.calculateConservationJacobianPWL(nElements_dof.shape[0],
                                        internalNodes.shape[0],
                                        w.shape[0],
                                        interiorElementBoundaries.shape[0],
                                        exteriorElementBoundaries.shape[0],
                                        w.shape[1],
                                        w.shape[2],
                                        w.shape[3],
                                        dofMapl2g.shape[1],
                                        n.shape[2],
                                        <int*>(interiorElementBoundaries.data),
                                        <int*>(exteriorElementBoundaries.data),
                                        <int*>(elementBoundaryElements.data),
                                        <int*>(elementBoundaryLocalElementBoundaries.data),
                                        <int*>(elementNodes.data),
                                        <int*>(dofMapl2g.data),
                                        <int*>(dofStarElements.data),
                                        <int*>(dofStarElementNeighbors.data),
                                        <int*>(nElements_dof.data),
                                        <int*>(internalNodes.data),
                                        <int*>(fluxElementBoundaries.data),
                                        <int*>(fluxBoundaryNodes.data),
                                        <double*>(w.data),
                                        <double*>(n.data),
                                        nodeStarFactor.nsf)

def calculateConservationFluxPWL(&nElements_node,
                                 &internalNodes,
                                 &fluxBoundaryNodes,
                                 NodeStarFactor nodeStarFactor):
    pp.calculateConservationFluxPWL(nElements_node.shape[0],
                                    internalNodes.shape[0],
                                    <int*>(nElements_node.data),
                                    <int*>(internalNodes),
                                    <int*>(fluxBoundaryNodes),
                                    nodeStarFactor.nsf)

def calculateConservationFluxPWL_noNeumannFix(&nElements_node,
                                              NodeStarFactor nodeStarFactor):
    pp.calculateConservationFluxPWL_noNeumannFix(nElements_node.shape[0],
                                                 <int*>(nElements_node.data),
                                                 nodeStarFactor.nsf)

def calculateConservationResidualPWL_opt(&nNodes_owned,
                                         &interiorElementBoundaries,
                                         &exteriorElementBoundaries,
                                         &elementBoundaryElements,
                                         &elementBoundaryLocalElementBoundaries,
                                         &elementNodes,
                                         &nodeStarElements,
                                         &nodeStarElementNeighbors,
                                         &nElements_node,
                                         &fluxElementBoundaries,
                                         &elementResidual,
                                         &vAverage,
                                         &dx,
                                         &w,
                                         &n,
                                         NodeStarFactor nodeStarFactor,
                                         &conservationResidual,
                                         &vConservative,
                                         &vConservative_element):
    pp.calculateConservationResidualPWL_opt(nNodes_owned,
                                            w.shape[0],
                                            interiorElementBoundaries.shape[0],
                                            exteriorElementBoundaries.shape[0],
                                            w.shape[1],
                                            w.shape[2],
                                            w.shape[3],
                                            n.shape[2],
                                            <int*>(interiorElementBoundaries.data),
                                            <int*>(exteriorElementBoundaries.data),
                                            <int*>(elementBoundaryElements.data),
                                            <int*>(elementBoundaryLocalElementBoundaries.data),
                                            <int*>(elementNodes.data),
                                            <int*>(nodeStarElements.data),
                                            <int*>(nodeStarElementNeighbors.data),
                                            <int*>(nElements_node.data),
                                            <int*>(fluxElementBoundaries.data),
                                            <double*>(elementResidual.data),
                                            <double*>(vAverage.data),
                                            <double*>(dx.data),
                                            <double*>(w.data),
                                            <double*>(n.data),
                                            nodeStarFactor.nsf,
                                            <double*>(conservationResidual.data),
                                            <double*>(vConservative.data),
                                            <double*>(vConservative_element.data))

def calculateConservationResidualPWL_primative(&interiorElementBoundaries,
                                               &exteriorElementBoundaries,
                                               &elementBoundaryElements,
                                               &elementBoundaryLocalElementBoundaries,
                                               &skipflag_elementBoundaries,
                                               &elementResidual,
                                               &dx,
                                               &n,
                                               &conservationResidual,
                                               &vConservative):
    pp.calculateConservationResidualPWL_primative(dx.shape[0],
                                                  interiorElementBoundaries.shape[0],
                                                  exteriorElementBoundaries.shape[0],
                                                  dx.shape[1],
                                                  dx.shape[2],
                                                  elementResidual.shape[1],
                                                  n.shape[2],
                                                  <int*>(interiorElementBoundaries.data),
                                                  <int*>(exteriorElementBoundaries.data),
                                                  <int*>(elementBoundaryElements.data),
                                                  <int*>(elementBoundaryLocalElementBoundaries.data),
                                                  <int*>(skipflag_elementBoundaries.data),
                                                  <double*>(elementResidual.data),
                                                  <double*>(dx.data),
                                                  <double*>(n.data),
                                                  <double*>(conservationResidual.data),
                                                  <double*>(vConservative.data))

def calculateConservationJacobianPWL_opt(int nNodes_owned,
                                         &interiorElementBoundaries,
                                         &exteriorElementBoundaries,
                                         &elementBoundaryElements,
                                         &elementBoundaryLocalElementBoundaries,
                                         &elementNodes,
                                         &nodeStarElements,
                                         &nodeStarElementNeighbors,
                                         &nElements_node,
                                         &internalNodes,
                                         &fluxElementBoundaries,
                                         &fluxBoundaryNodes,
                                         &w,
                                         &n,
                                         NodeStarFactor nodeStarFactor):
    pp.calculateConservationJacobianPWL_opt(nNodes_owned,
                                            nElements_node.shape[0],
                                            internalNodes.shape[0],
                                            w.shape[0],
                                            interiorElementBoundaries.shape[0],
                                            exteriorElementBoundaries.shape[0],
                                            w.shape[1],
                                            w.shape[2],
                                            w.shape[3],
                                            n.shape[2],
                                            <int*>(interiorElementBoundaries.data),
                                            <int*>(exteriorElementBoundaries.data),
                                            <int*>(elementBoundaryElements.data),
                                            <int*>(elementBoundaryLocalElementBoundaries.data),
                                            <int*>(elementNodes.data),
                                            <int*>(nodeStarElements.data),
                                            <int*>(nodeStarElementNeighbors.data),
                                            <int*>(nElements_node.data),
                                            <int*>(internalNodes.data),
                                            <int*>(fluxElementBoundaries.data),
                                            <int*>(fluxBoundaryNodes.data),
                                            <double*>(w.data),
                                            <double*>(n.data),
                                            nodeStarFactor.nsf)

def calculateConservationFluxPWL_opt(int nNodes_owned,
                                     &nElements_node,
                                     &internalNodes,
                                     &fluxBoundaryNodes,
                                     NodeStarFactor nodeStarFactor):
    pp.calculateConservationFluxPWL_opt(nNodes_owned,
                                        nElements_node.shape[0],
                                        internalNodes.shape[0],
                                        <int*>(nElements_node.data),
                                        <int*>(internalNodes.data),
                                        <int*>(fluxBoundaryNodes.data),
                                        nodeStarFactor.nsf)

def calculateConservationResidualPWL_interiorBoundaries(&interiorElementBoundaries,
                                                        &exteriorElementBoundaries,
                                                        &elementBoundaryElements,
                                                        &elementBoundaryLocalElementBoundaries,
                                                        &elementNodes,
                                                        &nodeStarElements,
                                                        &nodeStarElementNeighbors,
                                                        &nElements_node,
                                                        &fluxElementBoundaries,
                                                        &elementResidual,
                                                        &vAverage,
                                                        &dx,
                                                        &w,
                                                        &n,
                                                        NodeStarFactor nodeStarFactor,
                                                        &conservationResidual,
                                                        &vConservative,
                                                        &vConservative_element):
    assert fluxElementBoundaries.shape[0] == elementBoundaryLocalElementBoundaries.shape[0]
    pp.calculateConservationResidualPWL_interiorBoundaries(w.shape[0],
                                                           interiorElementBoundaries.shape[0],
                                                           exteriorElementBoundaries.shape[0],
                                                           w.shape[1],
                                                           w.shape[2],
                                                           w.shape[3],
                                                           n.shape[2],
                                                           <int*>(interiorElementBoundaries.data),
                                                           <int*>(exteriorElementBoundaries.data),
                                                           <int*>(elementBoundaryElements.data),
                                                           <int*>(elementBoundaryLocalElementBoundaries.data),
                                                           <int*>(elementNodes.data),
                                                           <int*>(nodeStarElements.data),
                                                           <int*>(nodeStarElementNeighbors.data),
                                                           <int*>(nElements_node.data),
                                                           <int*>(fluxElementBoundaries.data),
                                                           <double*>(elementResidual.data),
                                                           <double*>(vAverage.data),
                                                           <double*>(dx.data),
                                                           <double*>(w.data),
                                                           <double*>(n.data),
                                                           nodeStarFactor.nsf,
                                                           <double*>(conservationResidual.data),
                                                           <double*>(vConservative.data),
                                                           <double*>(vConservative_element.data))

def calculateConservationJacobianPWL_interiorBoundaries(&interiorElementBoundaries,
                                                        &exteriorElementBoundaries,
                                                        &elementBoundaryElements,
                                                        &elementBoundaryLocalElementBoundaries,
                                                        &elementNodes,
                                                        &nodeStarElements,
                                                        &nodeStarElementNeighbors,
                                                        &nElements_node,
                                                        &fluxElementBoundaries,
                                                        &w,
                                                        &n,
                                                        NodeStarFactor nodeStarFactor):
    assert fluxElementBoundaries.shape[0] == elementBoundaryLocalElementBoundaries.shape[0]
    pp.calculateConservationJacobianPWL_interiorBoundaries(nElements_node.shape[0],
                                                           w.shape[0],
                                                           interiorElementBoundaries.shape[0],
                                                           exteriorElementBoundaries.shape[0],
                                                           w.shape[1],
                                                           w.shape[2],
                                                           w.shape[3],
                                                           n.shape[2],
                                                           <int*>(interiorElementBoundaries.data),
                                                           <int*>(exteriorElementBoundaries.data),
                                                           <int*>(elementBoundaryElements.data),
                                                           <int*>(elementBoundaryLocalElementBoundaries.data),
                                                           <int*>(elementNodes.data),
                                                           <int*>(nodeStarElements.data),
                                                           <int*>(nodeStarElementNeighbors.data),
                                                           <int*>(nElements_node.data),
                                                           <int*>(fluxElementBoundaries.data),
                                                           <double*>(w.data),
                                                           <double*>(n.data),
                                                           nodeStarFactor.nsf)

def _subdomain_U_copy_global2local(int max_nN_owned,
                                   &elementNodes,
                                   &nodeStarElements,
                                   NodeStarFactor nodeStarFactor,
                                   &subdomain_U):
    pp.subdomain_U_copy_global2local(max_nN_owned,
                                     subdomain_U.shape[0],
                                     subdomain_U.shape[1],
                                     <int*>(elementNodes.data),
                                     <int*>(nodeStarElements.data),
                                     nodeStarFactor.nsf,
                                     <double*>(subdomain_U.data))

def _subdomain_U_copy_local2global(int max_nN_owned,
                                   &elementNodes,
                                   &nodeStarElements,
                                   NodeStarFactor nodeStarFactor,
                                   &subdomain_U):
    pp.subdomain_U_copy_local2global(max_nN_owned,
                                     subdomain_U.shape[0],
                                     subdomain_U.shape[1],
                                     <int*>(elementNodes.data),
                                     <int*>(nodeStarElements.data),
                                     nodeStarFactor.nsf,
                                     <double*>(subdomain_U.data))

def calculateConservationResidualGlobalBoundaries(&interiorElementBoundaries,
                                                  &exteriorElementBoundaries,
                                                  &elementBoundaryElements,
                                                  &elementBoundaryLocalElementBoundaries,
                                                  &exteriorElementBoundariesToSkip,
                                                  &dS,
                                                  &n,
                                                  &elementResidual,
                                                  &velocity,
                                                  &conservationResidual):
    pp.calculateConservationResidualGlobalBoundaries(conservationResidual.shape[0],
                                                     interiorElementBoundaries.shape[0],
                                                     exteriorElementBoundaries.shape[0],
                                                     dS.shape[1],
                                                     dS.shape[2],
                                                     elementResidual.shape[1],
                                                     velocity.shape[2],
                                                     <int*>(interiorElementBoundaries.data),
                                                     <int*>(exteriorElementBoundaries.data),
                                                     <int*>(elementBoundaryElements.data),
                                                     <int*>(elementBoundaryLocalElementBoundaries.data),
                                                     <int*>(exteriorElementBoundariesToSkip.data),
                                                     <double*>(dS.data),
                                                     <double*>(n.data),
                                                     <double*>(elementResidual.data),
                                                     <double*>(velocity.data),
                                                     <double*>(conservationResidual.data))

def updateSelectedExteriorElementBoundaryFlux(&exteriorElementBoundaries,
                                              &elementBoundaryElements,
                                              &elementBoundaryLocalElementBoundaries,
                                              &skipflag_elementBoundaries,
                                              &flux,
                                              &w,
                                              &residual):
    pp.updateSelectedExteriorElementBoundaryFlux(exteriorElementBoundaries.shape[0],
                                                 w.shape[1],
                                                 w.shape[2],
                                                 w.shape[3],
                                                 <int*>(exteriorElementBoundaries.data),
                                                 <int*>(elementBoundaryElements.data),
                                                 <int*>(elementBoundaryLocalElementBoundaries.data),
                                                 <int*>(skipflag_elementBoundaries.data),
                                                 <double*>(flux.data),
                                                 <double*>(w.data),
                                                 <double*>(residual.data))

def updateAdvectiveVelocityPointEval(double updateCoef,
                                     &advectiveVelocity,
                                     &velocity):
    cdef int nPoints = velocity.size/velocity.shape[-1]
    pp.postprocessAdvectiveVelocityPointEval(nPoints,
                                             velocity.shape[-1],
                                             updateCoef,
                                             <double*>(advectiveVelocity.data),
                                             <double*>(velocity.data))

def updateDiffusiveVelocityPointEval(double updateCoef,
                                     &diffusionTensor,
                                     &grad_phi,
                                     &velocity):
    cdef int nPoints = velocity.size/velocity.shape[-1]
    pp.postprocessDiffusiveVelocityPointEval(nPoints,
                                             velocity.shape[-1],
                                             updateCoef,
                                             <double*>(diffusionTensor.data),
                                             <double*>(grad_phi.data),
                                             <double*>(velocity.data))

def updateDiffusiveVelocityPointEval_sd(double updateCoef,
                                        &rowptr,
                                        &colind,
                                        &diffusionTensor,
                                        &grad_phi,
                                        &velocity):
    cdef int nPoints = velocity.size/velocity.shape[-1]
    pp.postprocessDiffusiveVelocityPointEval_sd(nPoints,
                                                velocity.shape[-1],
                                                updateCoef,
                                                <int*>(rowptr.data),
                                                <int*>(colind.data),
                                                <double*>(diffusionTensor.data),
                                                <double*>(grad_phi.data),
                                                <double*>(velocity.data))

def calculateElementResidualPWL(&alpha,
                                &elementResidual,
                                &elementResidualPWL):
    pp.calculateElementResidualPWL(elementResidual.shape[0],
                                   elementResidual.shape[1],
                                   elementResidualPWL.shape[1],
                                   <double*>(alpha.data),
                                   <double*>(elementResidual.data),
                                   <double*>(elementResidualPWL.data))

def copyElementBoundaryVelocityToParVec(&ebq_velocity, 
                                        &permutations, 
                                        &ebq_v_par_local):
    cdef int ebN
    cdef int k
    cdef int I
    cdef int nElementBoundaries=ebq_v_par_local.shape[0]
    cdef int nQuadraturePoints=ebq_v_par_local.shape[1]
    cdef int nSpace=ebq_v_par_local.shape[2]
    for ebN in range(nElementBoundaries):
        for k in range(nQuadraturePoints):
            for I=0 in range(nSpace):
                ebq_v_par_local[ebN,permutations[ebN,k],I] = ebq_velocity[ebN,k,I]

def addAverageToParVec(&ebq_velocityAverage, 
                       &permutations, 
                       &ebq_v_par_local):
    cdef int ebN
    cdef int k
    cdef int I
    cdef int nElementBoundaries=ebq_v_par_local.shape[0]
    cdef int nQuadraturePoints=ebq_v_par_local.shape[1]
    cdef int nSpace=ebq_v_par_local.shape[2]
    for  ebN in range(nElementBoundaries):
        for k in range(nQuadraturePoints):
            for I in range(nSpace):
                ebq_v_par_local[ebN, permutations[ebN,k], I] += ebq_velocityAverage[ebN, k, I]

def copyParVecToElementBoundaryVelocity(&ebq_velocity, 
                                        &permutations, 
                                        &ebq_v_par_local):
    cdef int ebN
    cdef int k
    cdef int I
    cdef int nElementBoundaries=ebq_v_par_local.shape[0]
    cdef int nQuadraturePoints=ebq_v_par_local.shape[1]
    cdef int nSpace=ebq_v_par_local.shape[2]
    for ebN in range(nElementBoundaries):
        for k in range(nQuadraturePoints):
            for I in range(nSpace):
                ebq_velocity[ebN, k, I] = ebq_v_par_local[ebN, permutations[ebN, k], I]
