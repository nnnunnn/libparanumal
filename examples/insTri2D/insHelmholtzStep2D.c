#include "ins2D.h"

// complete a time step using LSERK4
void insHelmholtzStep2D(ins_t *ins, iint tstep,  iint haloBytes,
			dfloat * sendBuffer, dfloat * recvBuffer, 
			char   * options){
  
  mesh2D *mesh = ins->mesh; 
  solver_t *solver = ins->vSolver; 
 // dfloat t = tstep*ins->dt;
  dfloat t = tstep*ins->dt + ins->dt;
  
  iint offset = mesh->Nelements+mesh->totalHaloPairs;

  iint rhsPackingMode = (strstr(options, "VECTORHELMHOLTZ")) ? 1:0;

  iint subcycling = (strstr(options,"SUBCYCLING")) ? 1:0;

   
   occaTimerTic(mesh->device,"HelmholtzRhsForcing"); 
   // compute all forcing i.e. f^(n+1) - grad(Pr)
  ins->helmholtzRhsForcingKernel(mesh->Nelements,
                                 subcycling,
			                           rhsPackingMode,
                                 mesh->o_vgeo,
                                 mesh->o_MM,
                                 ins->idt,
                                 ins->inu,
                                 ins->a0,
                                 ins->a1,
                                 ins->a2,
                                 ins->b0,
                                 ins->b1,
                                 ins->b2,
                                 ins->c0,
                                 ins->c1,
                                 ins->c2,
                                 ins->index,
                                 offset,
                                 ins->o_U,
                                 ins->o_V,
                                 ins->o_NU,
                                 ins->o_NV,
                                 ins->o_Px,
                                 ins->o_Py,
                                 ins->o_rhsU,
                                 ins->o_rhsV);
 occaTimerToc(mesh->device,"HelmholtzRhsForcing"); 


   occaTimerTic(mesh->device,"HelmholtzRhsIpdg");   
  ins->helmholtzRhsIpdgBCKernel(mesh->Nelements,
				                        rhsPackingMode,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                ins->tau,
                                t,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_vgeo,
                                mesh->o_sgeo,
                                mesh->o_EToB,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                ins->o_rhsU,
                                ins->o_rhsV);
    occaTimerToc(mesh->device,"HelmholtzRhsIpdg");   

  //use intermediate buffer for solve storage TODO: fix this later. Should be able to pull out proper buffer in elliptic solve
  if(rhsPackingMode==0){
    iint Ntotal = offset*mesh->Np;
    ins->o_UH.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
    ins->o_VH.copyFrom(ins->o_V,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));

    // printf("Solving for Ux ... ");
    mesh->device.finish();
    occa::tic("Ux-Solve");
    
    ins->NiterU = ellipticSolveTri2D( solver, ins->lambda, ins->velTOL, ins->o_rhsU, ins->o_UH, ins->vSolverOptions);
    mesh->device.finish();
    occa::toc("Ux-Solve"); 

    // printf("%d iteration(s)\n", ins->NiterU);
   


     // printf("Solving for Uy ... ");
    mesh->device.finish();
    occa::tic("Uy-Solve");
    ins->NiterV = ellipticSolveTri2D(solver, ins->lambda, ins->velTOL, ins->o_rhsV, ins->o_VH, ins->vSolverOptions);
    mesh->device.finish();
    occa::toc("Uy-Solve");
     // printf("%d iteration(s)\n", ins->NiterV);


    //copy into next stage's storage
    int index1 = (ins->index+1)%3; //hard coded for 3 stages
    ins->o_UH.copyTo(ins->o_U,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);
    ins->o_VH.copyTo(ins->o_V,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);  
  }else{

    printf("Solving for Ux and Uy \n");
    parAlmondPrecon(ins->precon->parAlmond, ins->o_rhsV, ins->o_rhsU); // rhs in rhsU, solution in rhsV

    iint Ntotal = mesh->Np*offset;
    dfloat *tmp  = (dfloat*) calloc(2*Ntotal, sizeof(dfloat));
    dfloat *tmpU = (dfloat*) calloc(Ntotal, sizeof(dfloat));
    dfloat *tmpV = (dfloat*) calloc(Ntotal, sizeof(dfloat));
    
    ins->o_rhsV.copyTo(tmp);

    for(iint n=0;n<Ntotal;++n){
      tmpU[n] = tmp[2*n+0];
      tmpV[n] = tmp[2*n+1];
    }

    //copy into next stage's storage
    int index1 = (ins->index+1)%3; //hard coded for 3 stages
    
    ins->o_U.copyFrom(tmpU,Ntotal*sizeof(dfloat),index1*offset*mesh->Np*sizeof(dfloat));
    ins->o_V.copyFrom(tmpV,Ntotal*sizeof(dfloat),index1*offset*mesh->Np*sizeof(dfloat));

    free(tmp);
    free(tmpU);
    free(tmpV);
  }
  
}