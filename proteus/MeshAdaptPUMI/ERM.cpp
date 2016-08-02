#include "MeshAdaptPUMI.h"
#include <PCU.h>
#include <petscksp.h>

#include <apf.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <apfDynamicMatrix.h>
#include <apfNumbering.h>

#include <iostream>
#include <fstream>

//proxy variables used to make it easier to pass these variables from MeshAdaptPUMIDrvr
int approx_order;
int int_order;
int norm_order;
double nu_0,nu_1,rho_0,rho_1;
double a_kl = 0.5; //flux term weight
int testcount = 0;

int eID = -1;
//int eID = 5081;
//int eID = 3865; // water element
//int eID = 7418; // air element  
//int eID = 4482; //mixed element

void getProps(double*rho,double*nu)
{
  rho_0 = rho[0];
  nu_0 = nu[0];      
  rho_1 = rho[1];
  nu_1 = nu[1];
  return;
}

double MeshAdaptPUMIDrvr::getMPvalue(double field_val,double val_0, double val_1)
{
  return val_0*(1-field_val)+val_1*field_val;
}

apf::Vector3 getFaceNormal(apf::Mesh* mesh, apf::MeshEntity* face){ //get the normal vector
  apf::Vector3 normal;
  apf::Adjacent verts;
  mesh->getAdjacent(face,0,verts);
  apf::Vector3 vtxs[verts.getSize()]; //4 points
  for(int i=0;i<verts.getSize();i++){
    mesh->getPoint(verts[i],0,vtxs[i]); 
  } 
  apf::Vector3 a,b;
  a = vtxs[1]-vtxs[0];
  b = vtxs[2]-vtxs[0];
  normal = apf::cross(a,b);

  return normal.normalize();
}

double getDotProduct(apf::Vector3 a, apf::Vector3 b){
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

double getDotProduct(apf::Matrix3x3 a, apf::Matrix3x3 b){
  double temp =0;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      temp = temp + a[i][j]*b[i][j];
    }
  }
  return temp;
}


bool isInTet(apf::Mesh* mesh, apf::MeshEntity* ent, apf::Vector3 pt){
  bool isin=0;

  apf::Adjacent verts;
  mesh->getAdjacent(ent,0,verts);
  apf::Vector3 vtxs[4]; //4 points
  for(int i=0;i<4;i++){
    mesh->getPoint(verts[i],0,vtxs[i]); 
  } 
  apf::Vector3 c[4];
  c[0] = vtxs[1]-vtxs[0];
  c[1] = vtxs[2]-vtxs[0];
  c[2] = vtxs[3]-vtxs[0];
  c[3] = pt-vtxs[0];

  apf::Matrix3x3 K,Kinv;
  apf::Vector3 F;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      K[i][j] = getDotProduct(c[i],c[j]);
    }
    F[i] = getDotProduct(c[3],c[i]);
  }
  Kinv = invert(K);
  apf::DynamicMatrix Kinv_dyn = apf::fromMatrix(Kinv);
  apf::DynamicVector F_dyn = apf::fromVector(F);
  apf::DynamicVector uvw; //result
  apf::multiply(Kinv_dyn,F_dyn,uvw);
  if(uvw[0] >= 0 && uvw[1] >=0 && uvw[2] >=0 && (uvw[0]+uvw[1]+uvw[2])<=1) isin = 1;
  return isin;
}

double a_k(apf::Matrix3x3 u, apf::Matrix3x3 v,double nu){
  //u and v are gradients of a vector
  apf::Matrix3x3 temp_u = u+apf::transpose(u);
  apf::Matrix3x3 temp_v = v+apf::transpose(v);
  return nu*getDotProduct(temp_u,temp_v);
} 

double b_k(double a, apf::Matrix3x3 b){
  //b is a gradient of a vector
  return a*(b[0][0]+b[1][1]+b[1][1]);
}

double c_k(apf::Vector3 a, apf::Matrix3x3 b, apf::Vector3 c){
  //b is a gradient of a vector
  return getDotProduct(b*a,c);
}

void getLHS(Mat &K,apf::NewArray <apf::DynamicVector> &shdrv,int nsd,double weight, double visc_val,int nshl)
{
      PetscScalar term1[nshl][nshl], term2[nshl][nshl];
      //Calculate LHS Diagonal Block Term
      for(int s=0; s<nshl;s++){
        for(int t=0; t<nshl;t++){
          double temp=0;
          for(int j=0;j<nsd;j++){
            temp+=shdrv[s][j]*shdrv[t][j];
          }
          term1[s][t] = temp*weight*visc_val;
        }
      } 
      int idx[nshl]; //indices for PETSc Mat insertion
      for(int i=0; i<nsd;i++){
        for(int j=0;j<nshl;j++){
          idx[j] = i*nshl+j;
        }
        MatSetValues(K,nshl,idx,nshl,idx,term1[0],ADD_VALUES);
      }
      int idxr[nshl],idxc[nshl]; //indices for PETSc rows and columns
      for(int i = 0; i< nsd;i++){
        for(int j=0; j< nsd;j++){
          for(int s=0;s<nshl;s++){
            for(int t=0;t<nshl;t++){
              term2[s][t] = shdrv[s][j]*shdrv[t][i]*weight*visc_val;
            }
          }
          for(int count=0;count<nshl;count++){
            idxr[count] = i*nshl+count;
            idxc[count] = j*nshl+count;
          }
          MatSetValues(K,nshl,idxr,nshl,idxc,term2[0],ADD_VALUES);
        }
      } //end 2nd term loop
}

void getRHS(Vec &F,apf::NewArray <double> &shpval,apf::NewArray <apf::DynamicVector> &shdrv,apf::Vector3 vel_vect,apf::Matrix3x3 grad_vel,
            int nsd,double weight,int nshl,
            double visc_val,double density,apf::Vector3 grad_density,double pressure,
            double g[3])
{
      int idx[nshl];
      for( int i = 0; i<nsd; i++){
        double temp_vect[nshl];
        for( int s=0;s<nshl;s++){
          idx[s] = i*nshl+s;

          //forcing term
          //temp_vect[s] = (g[i]+0.0)*shpval[s];
          //temp_vect[s] += pressure/density*shdrv[s][i]; //pressure term
          double force = (g[i]+0.0)*shpval[s];
          double pressure_force = pressure/density*shdrv[s][i];
          double a_rho_term = 0;
          double b_rho_term = -pressure/(density*density)*grad_density[i]*shpval[s];
          double a_term = 0; 
          double c_term = 0;
          //a(u,v) and c(u,u,v) term
          for(int j=0;j<nsd;j++){
            a_term += -visc_val*shdrv[s][j]*(grad_vel[i][j]+grad_vel[j][i]);
            a_rho_term += visc_val*(grad_vel[i][j]+grad_vel[j][i])*shpval[s]*grad_density[j]/(density);
            c_term += -shpval[s]*grad_vel[i][j]*vel_vect[j];///density;
          }
/*
if(testcount==eID){
  std::cout<<"RHS i "<<i<<" s "<<s<<std::endl;
  std::cout<<"force "<<force<<" gravity? "<<g[i]<<std::endl;
  std::cout<<"pressure "<<pressure_force<<" pressure? "<<pressure<<std::endl;
  std::cout<<"a_rho_term "<<a_rho_term<<" grad density "<< grad_density<<" density "<<density<<std::endl;
  std::cout<<"b_rho_term "<<b_rho_term<<" grad density "<< grad_density<<" density "<<density<<std::endl;
  std::cout<<"a_term "<<a_term<<std::endl;
  std::cout<<"c_term "<<c_term<<" shpval? "<<shpval[s]<<" gradvel "<<grad_vel<<" vel_vect "<<vel_vect<<" density "<<density<<std::endl;
  std::cout<<"combo " <<force+pressure_force+a_term+c_term<<std::endl;
  std::cout<<"temp_vect "<< temp_vect[s]<<std::endl;
}
*/
          temp_vect[s] = force+pressure_force+a_term+c_term;
          //temp_vect[s] = force+pressure_force+a_rho_term+b_rho_term+a_term+c_term;
          temp_vect[s] = temp_vect[s]*weight;
        } //end loop over number of shape functions
        VecSetValues(F,nshl,idx,temp_vect,ADD_VALUES);
      } //end loop over spatial dimensions
}

double getL2error(apf::Mesh* m, apf::MeshEntity* ent, apf::Field* voff, apf::Field* visc,apf::Field* pref, apf::Field* velf){

    int norm_order = int_order + 1;

    int nsd = m->getDimension();
    int nshl; //assuming linear solution
    int numqpt;
    int elem_type;

    apf::FieldShape* err_shape = apf::getLagrange(approx_order);
    apf::EntityShape* elem_shape;
    elem_type = m->getType(ent);

    nshl=apf::countElementNodes(err_shape,elem_type);

    double Lz = 0.05;
    double Ly = 0.2;
    double u_exact, u_h,p_exact,p_h, dpdy;
    double L2_err=0.0; 

    apf::MeshElement* element;
    apf::Element *pres_elem,*velo_elem;
    double weight, Jdet;
    apf::Matrix3x3 J;
    apf::Vector3 qpt;

    element = apf::createMeshElement(m,ent);
    pres_elem = apf::createElement(pref,element);
    velo_elem = apf::createElement(velf,element);
  
    numqpt=apf::countIntPoints(element,norm_order); //generally p*p maximum for shape functions
    apf::Vector3 xyz;
    apf::Vector3 vel_vect;

    for(int k=0;k<numqpt;k++){
      apf::getIntPoint(element,norm_order,k,qpt);
      apf::getJacobian(element,qpt,J); 
      Jdet=fabs(apf::getJacobianDeterminant(J,nsd)); 
      weight = apf::getIntWeight(element,norm_order,k);
   
      apf::mapLocalToGlobal(element,qpt,xyz);
      apf::getVector(velo_elem,qpt,vel_vect);

      //Hardcoded Exact Solution    
//int casenum = 1;
/*
      if(casenumber==0){ 
      //Poiseuille Flow
        dpdy = -1/Ly;
        u_exact= 0.5/(nu_0*rho_0)*(dpdy)*(xyz[2]*xyz[2]-Lz*xyz[2]);
        p_exact = 1+xyz[1]*dpdy;
      }
      else if(casenumber ==1){
      //Couette Flow
        //u_exact = 1.0*xyz[2]/Lz;
        u_exact = 2e-3*xyz[2]/Lz;
        p_exact = 0;
      }
*/

      u_h = vel_vect[1]; 
      p_h = apf::getScalar(pres_elem,qpt);

      double temp=0.0;
      temp = temp + (u_exact-u_h)*(u_exact-u_h);
      temp = temp + (p_exact-p_h)*(p_exact-p_h)/rho_0/rho_0;
      L2_err = L2_err+temp *weight*Jdet;
   }

  apf::destroyMeshElement(element);apf::destroyElement(velo_elem);apf::destroyElement(pres_elem);
  return L2_err;
}

double getStarerror(apf::Mesh* m, apf::MeshEntity* ent, apf::Field* voff, apf::Field* visc,apf::Field* pref, apf::Field* velf, apf::Field* estimate){

    int norm_order = int_order + 1;

    int nsd = m->getDimension();
    int nshl_err,nshl_est;  //exact error vs estimate
    int numqpt;
    int elem_type;

    apf::FieldShape* err_shape = apf::getLagrange(approx_order);
    apf::FieldShape* est_shape = apf::getHierarchic(2);
    apf::EntityShape* elem_shape;
    elem_type = m->getType(ent);

    nshl_err=apf::countElementNodes(err_shape,elem_type);
    nshl_est=apf::countElementNodes(est_shape,elem_type);

    double Lz = 0.05;
    double Ly = 0.2;
    //double u_exact, u_h,
    double p_exact,p_h, dpdy,div_u_h;
    double star_err=0.0; 

    apf::MeshElement* element;
    apf::Element *pres_elem,*velo_elem, *est_elem;
    double weight, Jdet;
    apf::Matrix3x3 J,grad_u_exact,grad_u_h,grad_est;
    apf::Vector3 qpt;

    element = apf::createMeshElement(m,ent);
    pres_elem = apf::createElement(pref,element);
    velo_elem = apf::createElement(velf,element);
    est_elem = apf::createElement(estimate,element);
  
    numqpt=apf::countIntPoints(element,norm_order); //generally p*p maximum for shape functions
    apf::Vector3 xyz;
    apf::Vector3 vel_vect;
    apf::Vector3 est_vect;
    apf::Vector3 u_exact, u_h;

    for(int k=0;k<numqpt;k++){
      apf::getIntPoint(element,norm_order,k,qpt);
      apf::getJacobian(element,qpt,J); 
      Jdet=fabs(apf::getJacobianDeterminant(J,nsd)); 
      weight = apf::getIntWeight(element,norm_order,k);
   
      apf::mapLocalToGlobal(element,qpt,xyz);
      apf::getVector(velo_elem,qpt,vel_vect);

      //Hardcoded Exact Solution    
//int casenum = 1;
/*
      if(casenumber==0){ 
      //Poiseuille Flow
        dpdy = -1/Ly;
        u_exact[0]=0;
        u_exact[1]= 0.5/(nu_0*rho_0)*(dpdy)*(xyz[2]*xyz[2]-Lz*xyz[2]);
        u_exact[2]=0;
        p_exact = 1+xyz[1]*dpdy;
        grad_u_exact[1][2] = 0.5/(nu_0*rho_0)*(dpdy)*(2*xyz[2]-Lz);
      }
      else if(casenumber ==1){
      //Couette Flow
        u_exact[0]=0;
        u_exact[1] = 1.0*xyz[2]/Lz;
        u_exact[2]=0;
        p_exact = 0;
        grad_u_exact[1][2] = 1.0/Lz;
      }
*/

      u_h = vel_vect; 
      p_h = apf::getScalar(pres_elem,qpt);

      apf::getVector(est_elem,qpt,est_vect); 
      apf::getVectorGrad(est_elem,qpt,grad_est);
      apf::getVectorGrad(velo_elem,qpt,grad_u_h);
      grad_u_h=apf::transpose(grad_u_h);
      div_u_h = grad_u_h[0][0]+grad_u_h[1][1]+grad_u_h[2][2];            

      double temp=0.0;
      temp = temp+a_k(grad_u_exact-grad_u_h,grad_est,nu_0); //nu is hardcoded because it's not necessary to generalize yet
      temp = temp-b_k(p_exact-p_h,grad_est); 
      temp = temp-b_k(div_u_h,grad_u_exact-grad_u_h); 
      temp = temp + c_k(u_exact,grad_u_exact,est_vect);
      temp = temp - c_k(u_h,grad_u_h,est_vect);
if(k==0)
//std::cout<<"C contribution "<<c_k(u_exact,grad_u_exact,est_vect)<<" "<<c_k(u_h,grad_u_h,est_vect)<<" "<<temp<<std::endl;
//std::cout<<"Grad "<<grad_u_h<<std::endl;

      star_err = star_err+temp *weight*Jdet;
   }

  apf::destroyMeshElement(element);apf::destroyElement(velo_elem);apf::destroyElement(pres_elem); apf::destroyElement(est_elem);
  return star_err;
}

void MeshAdaptPUMIDrvr::computeDiffusiveFlux(apf::Mesh*m,apf::Field* voff, apf::Field* visc,apf::Field* pref, apf::Field* velf){
  std::cerr<<"Begin computeDiffusiveFlux()"<<std::endl;
  int numbqpt, nshl;
  int hier_off = 4;
  apf::MeshEntity* bent,*ent;
  apf::MeshIterator* iter = m->begin(nsd-1); //loop over faces

  //Need to get number of bqpt
  while(bent = m->iterate(iter)){
    apf::MeshElement* b_elem;
    b_elem = apf::createMeshElement(m,bent);
    numbqpt = apf::countIntPoints(b_elem,int_order);
    apf::destroyMeshElement(b_elem);
    break;
  }
  m->end(iter);
    
  diffFlux = m->createDoubleTag("diffFlux",numbqpt*nsd*2);
  apf::MeshElement* tempelem; apf::Element * tempvelo,*temppres,*tempvoff;
  apf::MeshElement* b_elem;
  apf::Adjacent adjFaces;
  apf::Vector3 normal,centerdir;
  int orientation;
  double tempflux[numbqpt*nsd];
  double *flux; flux = (double*) calloc(numbqpt*nsd*2,sizeof(double));
  apf::NewArray <double> shpval;
  apf::NewArray <double> shpval_temp;

  apf::FieldShape* err_shape = apf::getHierarchic(2);
  apf::EntityShape* elem_shape;
  apf::Vector3 bqpt,bqptl,bqptshp;
  double weight, Jdet;
  apf::Matrix3x3 J;
  apf::Matrix3x3 tempgrad_velo;
  apf::Matrix3x3 identity(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
  
  iter=m->begin(nsd-1);
  while(ent=m->iterate(iter))
  {
    m->setDoubleTag(ent,diffFlux,flux);
  }
  m->end(iter);
  free(flux);

std::cerr<<"Initialized flux"<<std::endl;
  //loop over regions
  PCU_Comm_Begin();
  iter = m->begin(nsd);
  int ent_count=0;
  while(ent = m->iterate(iter))
  { 
    //Shape functions of the region and not the boundaries
    if(ent_count==0){
      nshl=apf::countElementNodes(err_shape,m->getType(ent));
      shpval_temp.allocate(nshl);
      nshl= nshl-hier_off;
      shpval.allocate(nshl);
      elem_shape = err_shape->getEntityShape(m->getType(ent));
    }

    m->getAdjacent(ent,nsd-1,adjFaces);
    for(int adjcount =0;adjcount<adjFaces.getSize();adjcount++){
      bent = adjFaces[adjcount];
      normal=getFaceNormal(m,bent);
      centerdir=apf::getLinearCentroid(m,ent)-apf::getLinearCentroid(m,bent);
      if(isInTet(m,ent,apf::project(normal,centerdir)*centerdir.getLength()+apf::getLinearCentroid(m,bent)))
        orientation = 1;
      else
        orientation = 0;
//      apf::getOrientation(m,ent,bent,adjcount,orientation,0);

      //begin calculation of flux
      b_elem = apf::createMeshElement(m,bent);
      tempelem = apf::createMeshElement(m,ent);
      temppres = apf::createElement(pref,tempelem);
      tempvelo = apf::createElement(velf,tempelem);
      tempvoff = apf::createElement(voff,tempelem);

      for(int l = 0;l<numbqpt;l++)
      {
        apf::Vector3 bflux(0.0,0.0,0.0); 
        apf::Matrix3x3 tempbflux(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
        apf::getIntPoint(b_elem,int_order,l,bqpt);
        weight = apf::getIntWeight(b_elem,int_order,l);
        apf::getJacobian(b_elem,bqpt,J); //evaluate the Jacobian at the quadrature point
        Jdet=fabs(apf::getJacobianDeterminant(J,nsd-1));
        bqptl=apf::boundaryToElementXi(m,bent,ent,bqpt); 
        apf::getVectorGrad(tempvelo,bqptl,tempgrad_velo);

        apf::ModelEntity* me=m->toModel(bent);
        int tag = m->getModelTag(me);
        apf::ModelEntity* boundary_face = m->findModelEntity(nsd-1,tag);

        if(me==boundary_face && has_gBC){
          int BCtype[4];
          double fluxdata[4][numbqpt];
          //for(int i=1;i<nsd+1;i++){ //ignores 0th index because that's pressure
          //  m->getIntTag(bent,BCtag[i],&(BCtype[i]));                 
          //}
          m->getIntTag(bent,BCtag,&(BCtype[0]));                 
          if((BCtype[1]+BCtype[2]+BCtype[3] != 3) && BCtype[1] == 1 ){
            std::cerr << "diffusive flux not fully specified on face " << localNumber(bent) << '\n';
            std::cerr << "BCtype "<<BCtype[1]<<" "<<BCtype[2]<<" "<<BCtype[3]<<std::endl;
            abort();
          }        
          if(BCtype[1]+BCtype[2]+BCtype[3] == 3){
            for(int i=1;i<nsd+1;i++)
              m->getDoubleTag(bent,fluxtag[i],&(fluxdata[i][0]));
            bflux = {fluxdata[1][l],fluxdata[2][l],fluxdata[3][l]};
            bflux = bflux-identity*apf::getScalar(temppres,bqptl)/getMPvalue(apf::getScalar(tempvoff,bqptl),rho_0,rho_1)*normal;
          }
          else{
            tempbflux = (tempgrad_velo+apf::transpose(tempgrad_velo))*getMPvalue(apf::getScalar(tempvoff,bqptl),nu_0,nu_1)
                   -identity*apf::getScalar(temppres,bqptl)/getMPvalue(apf::getScalar(tempvoff,bqptl),rho_0,rho_1);
            bflux = tempbflux*normal;
          }
        }
        else{
          tempbflux = (tempgrad_velo+apf::transpose(tempgrad_velo))*getMPvalue(apf::getScalar(tempvoff,bqptl),nu_0,nu_1)
              -identity*apf::getScalar(temppres,bqptl)/getMPvalue(apf::getScalar(tempvoff,bqptl),rho_0,rho_1);
          bflux = tempbflux*normal;
//*
if(comm_rank==0 && localNumber(ent)==eID){
  std::cout<<"quadrature point "<<l<<" value "<<bqpt<<" weight "<<weight<<" Jdet "<<Jdet<<std::endl;
  std::cout<<"velocity gradient "<<std::endl;
  std::cout<<tempgrad_velo<<std::endl;
  std::cout<<"Viscosity value "<<getMPvalue(apf::getScalar(tempvoff,bqptl),nu_0,nu_1)<<std::endl;
  std::cout<<"pressure "<<apf::getScalar(temppres,bqptl)<<std::endl;
  std::cout<<"density "<<getMPvalue(apf::getScalar(tempvoff,bqptl),rho_0,rho_1)<<std::endl;
  std::cout<<"tempbflux "<<std::endl;
  std::cout<< (tempgrad_velo+apf::transpose(tempgrad_velo))*getMPvalue(apf::getScalar(tempvoff,bqptl),nu_0,nu_1)
              -identity*apf::getScalar(temppres,bqptl)/getMPvalue(apf::getScalar(tempvoff,bqptl),rho_0,rho_1)<<std::endl;
  std::cout<<"normal? "<<normal<<std::endl;
  std::cout<<"bflux "<<tempbflux*normal<<std::endl;

}
//*/ 
        } //end if boundary
        bflux = bflux*weight*Jdet;
        bflux.toArray(&(tempflux[l*nsd]));
      }
      flux = (double*) calloc(numbqpt*nsd*2,sizeof(double));
      m->getDoubleTag(bent,diffFlux,flux);
      for (int i=0;i<numbqpt*nsd;i++){
        flux[orientation*numbqpt*nsd+i] = tempflux[i];
      }
      m->setDoubleTag(bent,diffFlux,flux);
      free(flux);
      apf::destroyMeshElement(tempelem);apf::destroyElement(tempvelo);apf::destroyElement(temppres); apf::destroyElement(tempvoff);
      
      //Parallel Communications
      apf::ModelEntity* me=m->toModel(bent);
      apf::ModelEntity* boundary_face = m->findModelEntity(nsd-1,m->getModelTag(me));
      apf::Copies remotes;
      if(m->isShared(bent))
      {
        m->getRemotes(bent,remotes);
        for(apf::Copies::iterator it=remotes.begin(); it!=remotes.end();++it)
        {
          PCU_COMM_PACK(it->first, it->second);
          PCU_COMM_PACK(it->first, orientation);
          PCU_COMM_PACK(it->first, tempflux);
        }
      } //end if
      ent_count++;
    } //end loop over faces
  } //end loop over regions
  m->end(iter);
std::cerr<<"Sending flux"<<std::endl;
  PCU_Comm_Send(); 
  flux = (double*) calloc(numbqpt*nsd*2,sizeof(double));
  while(PCU_Comm_Receive())
  {
    PCU_COMM_UNPACK(bent);
    PCU_COMM_UNPACK(orientation);
    PCU_COMM_UNPACK(tempflux);
    m->getDoubleTag(bent,diffFlux,flux);
    for (int i=0;i<numbqpt*nsd;i++){
//if(localNumber(ent)==eID)
//std::cout<<"Before numbqpt i "<<i<<" flux "<<flux[orientation*numbqpt*nsd+i]<<" "<<tempflux[i]<<std::endl;
      flux[orientation*numbqpt*nsd+i] = flux[orientation*numbqpt*nsd+i]+tempflux[i];
//if(localNumber(ent)==eID)
//std::cout<<"After numbqpt i "<<i<<" flux "<<flux[orientation*numbqpt*nsd+i]<<" "<<tempflux[i]<<std::endl;
    }
    m->setDoubleTag(bent,diffFlux,flux);
  }
  PCU_Barrier();
  free(flux);
  std::cerr<<"End computeDiffusiveFlux()"<<std::endl;
}
/*
void MeshAdaptPUMIDrvr::getBoundaryFlux(apf::Mesh* m, apf::MeshEntity* ent, apf::Field* voff, apf::Field* visc,apf::Field* pref, apf::Field* velf, double * endflux){

    int nshl;
    apf::NewArray <double> shpval;
    apf::NewArray <double> shpval_temp;

    //apf::FieldShape* err_shape = apf::getLagrange(approx_order);
    apf::FieldShape* err_shape = apf::getHierarchic(2);
    apf::EntityShape* elem_shape;

    //loop over element faces
    apf::Adjacent boundaries;
    apf::Adjacent neighbors;
    apf::MeshEntity* bent;
    apf::MeshElement* b_elem;
    apf::Vector3 bqpt,bqptl,bqptshp;

    double weight, Jdet;
    apf::Matrix3x3 J;
    apf::Vector3 normal;
    apf::Vector3 centerdir;

    //Shape functions of the region and not the boundaries
    nshl=apf::countElementNodes(err_shape,m->getType(ent));
    shpval_temp.allocate(nshl);
    int hier_off = 4;
    nshl= nshl-hier_off;
    shpval.allocate(nshl);
    elem_shape = err_shape->getEntityShape(m->getType(ent));
  
    m->getAdjacent(ent,nsd-1,boundaries);
    for(int adjcount =0;adjcount<boundaries.getSize();adjcount++){

      apf::Vector3 bflux(0.0,0.0,0.0); 
      apf::Matrix3x3 tempbflux(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
      bent = boundaries[adjcount];

      apf::ModelEntity* me=m->toModel(bent);
      int tag = m->getModelTag(me);
      apf::ModelEntity* boundary_face = m->findModelEntity(nsd-1,tag);
        if(m->isShared(bent)){//is shared by a parallel entity
          std::cout<<"PARALLEL "<<std::endl;
        }
        else{
          m->getAdjacent(bent,nsd,neighbors);
          b_elem = apf::createMeshElement(m,bent);
          apf::Matrix3x3 tempgrad_velo[2];
          apf::Matrix3x3 identity(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
          apf::MeshElement* tempelem; apf::Element * tempvelo,*temppres,*tempvoff;
          
          normal=getFaceNormal(m,bent);
          centerdir=apf::getLinearCentroid(m,ent)-apf::getLinearCentroid(m,bent);
          if(isInTet(m,ent,apf::project(normal,centerdir)*centerdir.getLength()+apf::getLinearCentroid(m,bent)))
            normal = normal*-1.0; //normal needs to face the other direction

          //shape functions are from weighting function and independent of neighbors

          double flux_weight;         
          for(int idx_neigh=0; idx_neigh<neighbors.getSize();idx_neigh++){ //at most two neighboring elements
            int numqpt = apf::countIntPoints(b_elem,int_order); 
            int BCtype[nsd];
            double fluxdata[4][numqpt];
            if(me==boundary_face && has_gBC){
              m->getIntTag(bent,BCtag,&(BCtype[0]));                 
              for(int i=1;i<nsd+1;i++){ //ignores 0th index because that's pressure
                //m->getIntTag(bent,BCtag[i],&BCtype[i]);                 
                m->getDoubleTag(bent,fluxtag[i],&(fluxdata[i][0]));
              }
            }

            tempelem = apf::createMeshElement(m,neighbors[idx_neigh]);
            temppres = apf::createElement(pref,tempelem);
            tempvelo = apf::createElement(velf,tempelem);
            tempvoff = apf::createElement(voff,tempelem);
            for(int l=0; l<apf::countIntPoints(b_elem,int_order);l++){
              apf::getIntPoint(b_elem,int_order,l,bqpt);
              weight = apf::getIntWeight(b_elem,int_order,l);
              apf::getJacobian(b_elem,bqpt,J); //evaluate the Jacobian at the quadrature point
              Jdet=fabs(apf::getJacobianDeterminant(J,nsd-1));
              bqptl=apf::boundaryToElementXi(m,bent,neighbors[idx_neigh],bqpt); 
              bqptshp=apf::boundaryToElementXi(m,bent,ent,bqpt); 
              elem_shape->getValues(NULL,NULL,bqptshp,shpval_temp);
              for(int j=0;j<nshl;j++){ shpval[j] = shpval_temp[hier_off+j];}

              if(me==boundary_face){
                if((BCtype[1]+BCtype[2]+BCtype[3] != 3) && BCtype[1] == 1 ){
                  std::cerr << "diffusive flux not fully specified on face " << localNumber(bent) << '\n';
                  std::cerr << "BCtype "<<BCtype[1]<<" "<<BCtype[2]<<" "<<BCtype[3]<<std::endl;
                  abort();
                }        
                if(BCtype[1]+BCtype[2]+BCtype[3] == 3){
                  bflux = {fluxdata[1][l],fluxdata[2][l],fluxdata[3][l]};
                  bflux = bflux-identity*apf::getScalar(temppres,bqptl)/getMPvalue(apf::getScalar(tempvoff,bqptl),rho_0,rho_1)*normal;
                  bflux = bflux*weight*Jdet;
                }
                else{
                  flux_weight = 1;
                  apf::getVectorGrad(tempvelo,bqptl,tempgrad_velo[idx_neigh]);
                  tempbflux = ((tempgrad_velo[idx_neigh]+apf::transpose(tempgrad_velo[idx_neigh]))*getMPvalue(apf::getScalar(tempvoff,bqptl),nu_0,nu_1)
                    -identity*apf::getScalar(temppres,bqptl)/getMPvalue(apf::getScalar(tempvoff,bqptl),rho_0,rho_1))*weight*Jdet*flux_weight;
                  bflux = tempbflux*normal;
                }
              }
              else{
                if(neighbors[idx_neigh]==ent) flux_weight = 1-a_kl;
                else flux_weight = a_kl;
                apf::getVectorGrad(tempvelo,bqptl,tempgrad_velo[idx_neigh]);
                tempbflux = ((tempgrad_velo[idx_neigh]+apf::transpose(tempgrad_velo[idx_neigh]))*getMPvalue(apf::getScalar(tempvoff,bqptl),nu_0,nu_1)
                  -identity*apf::getScalar(temppres,bqptl)/getMPvalue(apf::getScalar(tempvoff,bqptl),rho_0,rho_1))*weight*Jdet*flux_weight;
                bflux = tempbflux*normal;
              }

              for(int i=0;i<nsd;i++){ 
                for(int s=0;s<nshl;s++){
                  endflux[i*nshl+s] = endflux[i*nshl+s]+bflux[i]*shpval[s];
                }
              } 
  
            } //end boundary integration loop
          } //end for loop of neighbors

          apf::destroyMeshElement(tempelem);apf::destroyElement(tempvelo);apf::destroyElement(temppres); apf::destroyElement(tempvoff);
        }
    } //end loop over adjacent faces
}//end function
*/

void MeshAdaptPUMIDrvr::getBoundaryFlux(apf::Mesh* m, apf::MeshEntity* ent, double * endflux){

    int nshl;
    apf::NewArray <double> shpval;
    apf::NewArray <double> shpval_temp;
    double* flux;

    apf::FieldShape* err_shape = apf::getHierarchic(2);
    apf::EntityShape* elem_shape;

    //loop over element faces
    apf::Adjacent boundaries;
    apf::MeshEntity* bent;
    apf::MeshElement* b_elem;
    apf::Vector3 bqpt,bqptl,bqptshp;

    apf::Vector3 normal;
    apf::Vector3 centerdir;

    //Shape functions of the region and not the boundaries
    nshl=apf::countElementNodes(err_shape,m->getType(ent));
    shpval_temp.allocate(nshl);
    int hier_off = 4;
    nshl= nshl-hier_off;
    shpval.allocate(nshl);
    elem_shape = err_shape->getEntityShape(m->getType(ent));

    m->getAdjacent(ent,nsd-1,boundaries);
    for(int adjcount =0;adjcount<boundaries.getSize();adjcount++){

      apf::Vector3 bflux(0.0,0.0,0.0); 
      bent = boundaries[adjcount];
          
      b_elem = apf::createMeshElement(m,bent);
      normal=getFaceNormal(m,bent);
      centerdir=apf::getLinearCentroid(m,ent)-apf::getLinearCentroid(m,bent);
      int orientation = 0;
      if(isInTet(m,ent,apf::project(normal,centerdir)*centerdir.getLength()+apf::getLinearCentroid(m,bent))){
            normal = normal*-1.0; //normal needs to face the other direction
            orientation=1;
      }
      apf::ModelEntity* me=m->toModel(bent);
      int tag = m->getModelTag(me);
      apf::ModelEntity* boundary_face = m->findModelEntity(nsd-1,tag);

      double flux_weight[2];         
      if(me==boundary_face){
        if(orientation==0){flux_weight[0]=1; flux_weight[1]=0;}
        else{flux_weight[0]=0; flux_weight[1]=-1;}
      }
      else{
        if(orientation==0){flux_weight[0]=(1-a_kl); flux_weight[1]=a_kl;}
        else{ flux_weight[0]=-a_kl; flux_weight[1]=-1*(1-a_kl);}
      }

      int numbqpt = apf::countIntPoints(b_elem,int_order); 
      flux = (double*) calloc(numbqpt*nsd*2,sizeof(double));
      m->getDoubleTag(bent,diffFlux,flux);

      for(int l=0; l<numbqpt;l++){
        apf::getIntPoint(b_elem,int_order,l,bqpt);
        bqptshp=apf::boundaryToElementXi(m,bent,ent,bqpt); 
        elem_shape->getValues(NULL,NULL,bqptshp,shpval_temp);
        for(int j=0;j<nshl;j++){shpval[j] = shpval_temp[hier_off+j];}
        for(int i=0;i<nsd;i++){
          for(int s=0;s<nshl;s++){
            endflux[i*nshl+s] = endflux[i*nshl+s]+(flux_weight[0]*flux[l*nsd+i]+flux_weight[1]*flux[numbqpt*nsd+l*nsd+i])*shpval[s];
          }
//if(localNumber(ent)==3865)
//std::cout<<"Components "<<flux_weight[0]*flux[l*nsd+i]<<" "<<flux_weight[1]*flux[numbqpt*nsd+l*nsd+i]<<std::endl;
        }
      }//end of boundary integration loop
      free(flux);
   } //end for adjacent faces
}

apf::Field* MeshAdaptPUMIDrvr::getViscosityField(apf::Field* voff)
{
  apf::Field* visc = apf::createLagrangeField(m,"viscosity",apf::SCALAR,1);
  apf::MeshEntity* ent;
  apf::MeshIterator* iter = m->begin(0);
  double vof_val, visc_val;
  while(ent = m->iterate(iter)){ //loop through all vertices
    vof_val=apf::getScalar(voff,ent,0);
    visc_val = getMPvalue(vof_val,nu_0, nu_1);
    apf::setScalar(visc, ent, 0,visc_val);
  }
  m->end(iter);
  return visc;
}


void setErrorField(apf::Field* estimate,Vec coef,apf::MeshEntity* ent,int nsd,int nshl)
{
    apf::Mesh* m = apf::getMesh(estimate);
    double coef_ez[nshl*nsd];
    int ez_idx[nshl*nsd];
    for(int ez=0;ez<nshl*nsd;ez++){ez_idx[ez]=ez;}
    VecGetValues(coef,nshl*nsd,ez_idx,coef_ez);

    //Copy coefficients onto field
    apf::Adjacent adjvert;
    m->getAdjacent(ent,0,adjvert);
    for(int idx=0;idx<4;idx++){
      double coef_sub[3]={0,0,0};
      apf::setVector(estimate,adjvert[idx],0,&coef_sub[0]);
    }
    
    apf::Adjacent adjedg;
    m->getAdjacent(ent,1,adjedg);
    for(int idx=0;idx<nshl;idx++){
      double coef_sub[3] ={coef_ez[idx],coef_ez[nshl+idx],coef_ez[nshl*2+idx]};
      apf::setVector(estimate,adjedg[idx],0,&coef_sub[0]);
    }
}

void MeshAdaptPUMIDrvr::removeBCData()
{
  if(comm_rank==0) std::cout<<"Start removing BC tags/data"<<std::endl;
  apf::MeshEntity* ent;   
  apf::MeshIterator* fIter = m->begin(2);
  while(ent=m->iterate(fIter))
  {
    if(has_gBC && m->hasTag(ent,BCtag)){
      m->removeTag(ent,BCtag);
      for(int i=0;i<4;i++)
      {
        if(i>0 && m->hasTag(ent,fluxtag[i]))
          m->removeTag(ent,fluxtag[i]);
      }
/*
    for(int i=0;i<4;i++)
    {
      if(i>0 && m->hasTag(ent,fluxtag[i]))
        m->removeTag(ent,fluxtag[i]);
      if(m->hasTag(ent,diffFlux))
        m->removeTag(ent,diffFlux);
>>>>>>> remotes/origin/local_adapt
*/
    }
    if(m->hasTag(ent,diffFlux))
      m->removeTag(ent,diffFlux);

  }
  std::cout<<"Flag 1"<<std::endl;
  m->end(fIter);
  if(has_gBC){
    m->destroyTag(BCtag);
    for(int i=0;i<4;i++)
    {
      if(i>0)
        m->destroyTag(fluxtag[i]);
    }
  }
  m->destroyTag(diffFlux);
  if(comm_rank==0) std::cout<<"Destroyed BC and flux tags"<<std::endl;
}

void MeshAdaptPUMIDrvr::get_local_error() 
//This function aims to compute error at each element via ERM.
//First get the mesh and impose a 2nd order field
//Then get the desired quadrature points
{
  getProps(rho,nu);
  approx_order = approximation_order; 
  int_order = integration_order;
  nsd = m->getDimension();

  //***** Get Solution Fields First *****//
  apf::Field* voff = m->findField("vof");
  assert(voff);
  apf::Field* velf = m->findField("velocity");
  assert(velf);
  apf::Field* pref = m->findField("p");
  assert(pref);
  //*****               *****//
  
  //***** Compute the viscosity field *****//
  apf::Field* visc = getViscosityField(voff);

  //***** Compute diffusive flux *****//
  computeDiffusiveFlux(m,voff,visc,pref,velf);

  //Initialize the Error Fields
  freeField(err_reg);
  err_reg = apf::createField(m,"ErrorRegion",apf::VECTOR,apf::getConstant(nsd));

  //Start computing element quantities
  int numqpt; //number of quadrature points
  int nshl; //number of local shape functions
  int elem_type; //what type of topology
  double weight; //value container for the weight at each qpt
  double Jdet;
  apf::FieldShape* err_shape = apf::getHierarchic(approx_order);
  apf::Field * estimate = apf::createField(m, "err_est", apf::VECTOR, err_shape);
  apf::EntityShape* elem_shape;
  apf::Vector3 qpt; //container for quadrature points
  apf::MeshElement* element;
  apf::Element* visc_elem, *pres_elem,*velo_elem,*vof_elem;
  apf::Element* est_elem;
  apf::Matrix3x3 J; //actual Jacobian matrix
  apf::Matrix3x3 invJ; //inverse of Jacobian
  apf::NewArray <double> shpval; //array to store shape function values at quadrature points
  apf::NewArray <double> shpval_temp; //array to store shape function values at quadrature points temporarily
  apf::NewArray <apf::Vector3> shgval; //array to store shape function values at quadrature points

  apf::DynamicMatrix invJ_copy;
  apf::NewArray <apf::DynamicVector> shdrv;
  apf::NewArray <apf::DynamicVector> shgval_copy;
  
  apf::MeshIterator* iter = m->begin(nsd); //loop over elements
  apf::MeshEntity* ent;
  

double L2_total=0;
double star_total=0;
double err_est = 0;
double err_est_total=0;
double u_norm_total=0;
  while(ent = m->iterate(iter)){ //loop through all elements
    
    elem_type = m->getType(ent);
    if(elem_type != m->TET){
      std::cout<<"Not a Tet present"<<std::endl;
      exit(0); 
    }
    element = apf::createMeshElement(m,ent);
    pres_elem = apf::createElement(pref,element);
    velo_elem = apf::createElement(velf,element);
    visc_elem = apf::createElement(visc,element); //at vof currently
    vof_elem = apf::createElement(voff,element);
  
    numqpt=apf::countIntPoints(element,int_order); //generally p*p maximum for shape functions
    nshl=apf::countElementNodes(err_shape,elem_type);
    shgval.allocate(nshl);
    shpval_temp.allocate(nshl);

    int hier_off = 4; //there is an offset that needs to be made to isolate the hierarchic edge modes
    nshl = nshl - hier_off;

    shpval.allocate(nshl);   shgval_copy.allocate(nshl); shdrv.allocate(nshl);

    //LHS Matrix Initialization
    int ndofs = nshl*nsd;
    Mat K; //matrix size depends on nshl, which may vary from element to element
    MatCreate(PETSC_COMM_SELF,&K);
    MatSetSizes(K,ndofs,ndofs,ndofs,ndofs);
    MatSetFromOptions(K);
    MatSetUp(K); //is this inefficient? check later

    //RHS Vector Initialization
    Vec F;
    VecCreate(PETSC_COMM_SELF,&F);
    VecSetSizes(F,ndofs,ndofs);
    VecSetUp(F);

    //loop through all qpts
    for(int k=0;k<numqpt;k++){
      apf::getIntPoint(element,int_order,k,qpt); //get a quadrature point and store in qpt
      apf::getJacobian(element,qpt,J); //evaluate the Jacobian at the quadrature point
      J = apf::transpose(J); //Is PUMI still defined in this way?
      invJ = invert(J);
      Jdet=fabs(apf::getJacobianDeterminant(J,nsd)); 
      weight = apf::getIntWeight(element,int_order,k);
      invJ_copy = apf::fromMatrix(invJ);

      //first get the shape function values for error shape functions
      elem_shape = err_shape->getEntityShape(elem_type);
      elem_shape->getValues(NULL,NULL,qpt,shpval_temp);
      elem_shape->getLocalGradients(NULL,NULL,qpt,shgval); 

      for(int i =0;i<nshl;i++){ //get the true derivative and copy only the edge modes for use
        shgval_copy[i] = apf::fromVector(shgval[i+hier_off]);
        shpval[i] = shpval_temp[i+hier_off];
        apf::multiply(shgval_copy[i],invJ_copy,shdrv[i]); 
      }

      //obtain needed values

      apf::Vector3 vel_vect;
      apf::Matrix3x3 grad_vel;
      apf::getVector(velo_elem,qpt,vel_vect);
      apf::getVectorGrad(velo_elem,qpt,grad_vel);
      grad_vel = apf::transpose(grad_vel);
      apf::Vector3 grad_vof;
      apf::getGrad(vof_elem,qpt,grad_vof);
    
      double density = getMPvalue(apf::getScalar(vof_elem,qpt),rho_0,rho_1);
      double pressure = apf::getScalar(pres_elem,qpt);
      double visc_val = apf::getScalar(visc_elem,qpt);
      apf::Vector3 grad_rho = grad_vof*(rho_1-rho_0);

/*
      if(testcount==eID && k==0 && comm_rank==0){
      
        std::cout<<std::setprecision(15);
        apf::Adjacent dbg_vadj;
        m->getAdjacent(ent,0,dbg_vadj);
        std::cout<<"adjacent verts ";
        apf::Vector3 testpt;
        for(int test_count=0;test_count<4;test_count++){
          m->getPoint(dbg_vadj[test_count],0,testpt);
          std::cout<<testpt<<" ";
        }
        std::cout<<std::endl;
  
        apf::Vector3 vectpt;
        std::cout<<"Velocity at Nodes\n";
        for(int test_count=0;test_count<4;test_count++){
          apf::getVector(velf,dbg_vadj[test_count],0,vectpt);
          std::cout<<vectpt<<" ";
        }
        std::cout<<std::endl;
        
        double pressurept;
        std::cout<<"Pressure at Nodes\n";
        for(int test_count=0;test_count<4;test_count++){
          pressurept = apf::getScalar(pref,dbg_vadj[test_count],0);
          std::cout<<pressurept<<" ";
        }
        std::cout<<std::endl;


        std::cout<<"nshl & numqpt "<<nshl<<" "<<numqpt<<std::endl;
        std::cout<<"Quadrature point "<<k<<std::endl;
        std::cout<<"quad pt" <<qpt<<std::endl;
        apf::Vector3 xyz;
        apf::mapLocalToGlobal(element,qpt,xyz);
        std::cout<<"Global qpt "<<xyz<<std::endl;
        std::cout<<"Jacobian "<<J<<std::endl;
        std::cout<<"Jdet "<<Jdet<<std::endl;
        std::cout<<"invJ "<<invJ<<std::endl;
        std::cout<<"Shape function"<<shpval[0]<<" "<<shpval[1]<<" "<<shpval[2]<<" "<<shpval[3]<<" "<<shpval[4]<<" "<<shpval[5]<<std::endl;
        std::cout<<"weight "<<weight<<std::endl;
        std::cout<<"Shape global derivatives\n";
        for(int i=0;i<nshl;i++){
          std::cout<<i<<" "<<shdrv[i][0]<<" "<<shdrv[i][1]<<" "<<shdrv[i][2]<<std::endl;
        }
        std::cout<<"Density "<<density<<std::endl;
        std::cout<<"pressure "<<pressure<<std::endl;
        std::cout<<"viscosity "<<visc_val<<std::endl;
        std::cout<<"Velocity "<<std::endl;
        std::cout<<vel_vect<<std::endl;
        std::cout<<"Velocity Gradient "<<std::endl;
        std::cout<<grad_vel<<std::endl;
      }
*/
      //Left-Hand Side
      getLHS(K,shdrv,nsd,weight,visc_val,nshl);

      //Get RHS
      getRHS(F,shpval,shdrv,vel_vect,grad_vel,nsd,weight,nshl,visc_val,density,grad_rho,pressure,g);

    } // end quadrature loop

    //to complete integration, scale by the determinant of the Jacobian

    MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);
    MatScale(K,Jdet); //must be done after assembly
    VecAssemblyBegin(F);
    VecAssemblyEnd(F); 
    VecScale(F,Jdet); //must be done after assembly
    double* bflux;
    int F_idx[ndofs];
    bflux = (double*) calloc(ndofs,sizeof(double));

    getBoundaryFlux(m, ent,bflux);
    for(int s=0;s<ndofs;s++){
      F_idx[s]=s;
    }
    VecSetValues(F,ndofs,F_idx,bflux,ADD_VALUES);
if(comm_rank==0 && testcount==eID){
  std::cout<<"What is bflux?\n";
  for(int s=0;s<ndofs;s++){
    std::cout<<bflux[s]<<std::endl;
  }
  std::cout<<"End of bflux "<<std::endl;
}
    VecAssemblyBegin(F); VecAssemblyEnd(F);
    free(bflux);
    Vec coef;
    VecCreate(PETSC_COMM_SELF,&coef);
    VecSetSizes(coef,ndofs,ndofs);
    VecSetUp(coef);

/*
if(testcount==eID && comm_rank==0){

//Save Temporarily for Debugging
      std::ofstream myfile ("stiffness.csv");
      std::ofstream myfile2 ("force.csv");
      myfile<<std::scientific<<std::setprecision(15);
      myfile2<<std::scientific<<std::setprecision(15);
      PetscScalar matstor;  
      PetscScalar vecstor;  
      int idxr[ndofs], idxc[ndofs];
      for(int ii=0;ii<ndofs;ii++){
        idxr[ii]=ii;
        for(int jj=0;jj<ndofs;jj++){
          idxc[jj]=jj;
          MatGetValues(K,1,&idxr[ii],1,&idxc[jj],&matstor);
          myfile<<matstor<<","; 
        }
        myfile<<std::endl;
        VecGetValues(F,1,&idxr[ii],&vecstor);
        myfile2<<vecstor<<std::endl;
      }
      myfile.close();
      myfile2.close();
    //MatView(K,PETSC_VIEWER_STDOUT_SELF);
    //VecView(F,PETSC_VIEWER_STDOUT_SELF);
}
*/

    KSP ksp; //initialize solver context
    KSPCreate(PETSC_COMM_SELF,&ksp);
    KSPSetOperators(ksp,K,K);
    KSPSetType(ksp,KSPPREONLY);
    PC pc;
    //PCSetOperators(pc,K,K);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
    KSPSetFromOptions(ksp);

    KSPSolve(ksp,F,coef);
    
    KSPDestroy(&ksp); //destroy ksp
    //PCDestroy(&pc);

    setErrorField(estimate,coef,ent,nsd,nshl);

    //compute the local error  
    double Acomp=0;
    double Bcomp=0;
    double visc_avg=0;
    double u_norm = 0;
    apf::Matrix3x3 phi_ij;
    apf::Matrix3x3 vel_ij;
    apf::Vector3 vel_vect;
    apf::Vector3 grad_vof;

    est_elem= apf::createElement(estimate,element);   
    for(int k=0; k<numqpt;k++){ 
      apf::getIntPoint(element,int_order,k,qpt); //get a quadrature point and store in qpt
      apf::getJacobian(element,qpt,J); //evaluate the Jacobian at the quadrature point

      invJ = invert(J);
      invJ = apf::transpose(invJ);
      Jdet=fabs(apf::getJacobianDeterminant(J,nsd)); 
      weight = apf::getIntWeight(element,int_order,k);
      invJ_copy = apf::fromMatrix(invJ);

      //first get the shape function values for error shape functions
      elem_shape = err_shape->getEntityShape(elem_type);
      elem_shape->getValues(NULL,NULL,qpt,shpval_temp);
      elem_shape->getLocalGradients(NULL,NULL,qpt,shgval); 

      for(int i =0;i<nshl;i++){ //get the true derivative and copy only the edge modes for use
        shgval_copy[i] = apf::fromVector(shgval[i+hier_off]);
        shpval[i] = shpval_temp[i+hier_off];
        apf::multiply(shgval_copy[i],invJ_copy,shdrv[i]); 
      }
      double visc_val = apf::getScalar(visc_elem,qpt);
      double pres_val = apf::getScalar(pres_elem,qpt);
      double density = getMPvalue(apf::getScalar(vof_elem,qpt),rho_0,rho_1);
      apf::getVectorGrad(est_elem,qpt,phi_ij);
      apf::getGrad(vof_elem,qpt,grad_vof);
      apf::getVector(velo_elem,qpt,vel_vect);
      apf::getVectorGrad(velo_elem,qpt,vel_ij);
      vel_ij = apf::transpose(vel_ij);
      phi_ij = apf::transpose(phi_ij);

      Acomp = Acomp + visc_val*getDotProduct(phi_ij,phi_ij+apf::transpose(phi_ij))*weight;
      Bcomp = Bcomp + apf::getDiv(velo_elem,qpt)*apf::getDiv(velo_elem,qpt)*weight;
      visc_avg = visc_avg + visc_val*weight;
      u_norm = u_norm + visc_val*getDotProduct(vel_ij,vel_ij+apf::transpose(vel_ij))*weight;

    } //end compute local error
    visc_avg = visc_avg*Jdet/apf::measure(element);
    Acomp = Acomp*Jdet/visc_avg; //nondimensionalize with average viscosity, Jacobians can cancel out, but this is done for clarity
    Bcomp = Bcomp*Jdet;
    u_norm = u_norm/visc_avg*Jdet;
    err_est = sqrt(Acomp+Bcomp); 

    apf::Vector3 err_in(err_est,Acomp,Bcomp);
    apf::setVector(err_reg,ent,0,err_in);
    err_est_total = err_est_total+(Acomp+Bcomp); //for tracking the upper bound
    u_norm_total = u_norm_total + u_norm;
   
    MatDestroy(&K); //destroy the matrix
    VecDestroy(&F); //destroy vector
    VecDestroy(&coef); //destroy vector

    apf::destroyElement(visc_elem);apf::destroyElement(pres_elem);apf::destroyElement(velo_elem);apf::destroyElement(est_elem);apf::destroyElement(vof_elem);
testcount++;

  } //end element loop

  PCU_Add_Doubles(&err_est_total,1);
  PCU_Add_Doubles(&u_norm_total,1);

  total_error = sqrt(err_est_total);
  u_norm_total = sqrt(u_norm_total);

  if(comm_rank==0){
    std::cout<<std::setprecision(10)<<std::endl;
    std::cout<<"Error estimate "<<total_error<<std::endl;
  }

  if(logging_config=="on"){
    std::cout<<"outputting error field\n";
    char namebuffer[20];
    sprintf(namebuffer,"err_reg_%i",nEstimate);
    apf::writeVtkFiles(namebuffer, m);
    nEstimate++;
  }

  m->end(iter);
  apf::destroyField(visc);
  apf::destroyField(estimate);

  printf("It cleared the function.\n");
}


