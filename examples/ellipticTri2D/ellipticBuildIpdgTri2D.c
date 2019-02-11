
#include "ellipticTri2D.h"

int parallelCompareRowColumn(const void *a, const void *b);

void ellipticBuildIpdgTri2D(mesh2D *mesh, int basisNp, dfloat *basis,
                            dfloat tau, dfloat lambda, iint *BCType, nonZero_t **A,
                            iint *nnzA, iint *globalStarts, const char *options){

  iint size, rankM;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankM);

  int Np = mesh->Np;
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;
  iint Nelements = mesh->Nelements;

  if(!basis) { // default to degree N Lagrange basis
    basisNp = Np;
    basis = (dfloat*) calloc(basisNp*basisNp, sizeof(dfloat));
    for(iint n=0;n<basisNp;++n){
      basis[n+n*basisNp] = 1;
    }
  }

  // number of degrees of freedom on this rank
  iint Nnum = Np*Nelements;

  // create a global numbering system
  iint *globalIds = (iint *) calloc((Nelements+mesh->totalHaloPairs)*Np,sizeof(iint));

  // every degree of freedom has its own global id
  MPI_Allgather(&(Nelements), 1, MPI_IINT, globalStarts+1, 1, MPI_IINT, MPI_COMM_WORLD);
  for(iint r=0;r<size;++r)
    globalStarts[r+1] = globalStarts[r]+globalStarts[r+1]*Np;

  /* so find number of elements on each rank */
  iint *rankNelements = (iint*) calloc(size, sizeof(iint));
  iint *rankStarts = (iint*) calloc(size+1, sizeof(iint));
  MPI_Allgather(&(Nelements), 1, MPI_IINT,
    rankNelements, 1, MPI_IINT, MPI_COMM_WORLD);
  //find offsets
  for(iint r=0;r<size;++r){
    rankStarts[r+1] = rankStarts[r]+rankNelements[r];
  }
  //use the offsets to set a global id
  for (iint e =0;e<Nelements;e++) {
    for (int n=0;n<Np;n++) {
      globalIds[e*Np +n] = n + (e + rankStarts[rankM])*Np;
    }
  }

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    iint *idSendBuffer = (iint *) calloc(Np*mesh->totalHaloPairs,sizeof(iint));
    meshHaloExchange(mesh, Np*sizeof(iint), globalIds, idSendBuffer, globalIds + Nelements*Np);
    free(idSendBuffer);
  }

  iint nnzLocalBound = basisNp*basisNp*(1+Nfaces)*Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(Nfaces*Nfp*Nfp,sizeof(dfloat));
  for (iint f=0;f<Nfaces;f++) {
    for (iint n=0;n<Nfp;n++) {
      iint fn = mesh->faceNodes[f*Nfp+n];

      for (iint m=0;m<Nfp;m++) {
        dfloat MSnm = 0;

        for (iint i=0;i<Np;i++){
          MSnm += mesh->MM[fn+i*Np]*mesh->LIFT[i*Nfp*Nfaces+f*Nfp+m];
        }
        MS[m+n*Nfp + f*Nfp*Nfp]  = MSnm;
      }
    }
  }


  // reset non-zero counter
  int nnz = 0;

  *A = (nonZero_t*) calloc(nnzLocalBound, sizeof(nonZero_t));

  dfloat *SM = (dfloat*) calloc(Np*Np,sizeof(dfloat));
  dfloat *SP = (dfloat*) calloc(Np*Np,sizeof(dfloat));

  // loop over all elements
  for(iint eM=0;eM<Nelements;++eM){

    iint vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat J = mesh->vgeo[vbase+JID];

    /* start with stiffness matrix  */
    for(iint n=0;n<Np;++n){
      for(iint m=0;m<Np;++m){
        SM[n*Np+m]  = J*lambda*mesh->MM[n*Np+m];
        SM[n*Np+m] += J*drdx*drdx*mesh->Srr[n*Np+m];
        SM[n*Np+m] += J*drdx*dsdx*mesh->Srs[n*Np+m];
        SM[n*Np+m] += J*dsdx*drdx*mesh->Ssr[n*Np+m];
        SM[n*Np+m] += J*dsdx*dsdx*mesh->Sss[n*Np+m];

        SM[n*Np+m] += J*drdy*drdy*mesh->Srr[n*Np+m];
        SM[n*Np+m] += J*drdy*dsdy*mesh->Srs[n*Np+m];
        SM[n*Np+m] += J*dsdy*drdy*mesh->Ssr[n*Np+m];
        SM[n*Np+m] += J*dsdy*dsdy*mesh->Sss[n*Np+m];
      }
    }

    for (iint fM=0;fM<Nfaces;fM++) {

      for (int n=0;n<Np*Np;n++) SP[n] =0;

      // load surface geofactors for this face
      iint sid = mesh->Nsgeo*(eM*Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat hinv = mesh->sgeo[sid+IHID];
      dfloat penalty = tau*hinv;

      iint eP = mesh->EToE[eM*Nfaces+fM];
      if (eP < 0) eP = eM;

      iint vbaseP = eP*mesh->Nvgeo;
      dfloat drdxP = mesh->vgeo[vbaseP+RXID];
      dfloat drdyP = mesh->vgeo[vbaseP+RYID];
      dfloat dsdxP = mesh->vgeo[vbaseP+SXID];
      dfloat dsdyP = mesh->vgeo[vbaseP+SYID];

      int bcD = 0, bcN =0;
      int bc = mesh->EToB[fM+Nfaces*eM]; //raw boundary flag
      iint bcType = 0;

      if(bc>0) bcType = BCType[bc];          //find its type (Dirichlet/Neumann)

      // this needs to be double checked (and the code where these are used)
      if(bcType==1){ // Dirichlet
        bcD = 1;
        bcN = 0;
      } else if(bcType==2){ // Neumann
        bcD = 0;
        bcN = 1;
      }

      // reset eP
      eP = mesh->EToE[eM*Nfaces+fM];

      // mass matrix for this face
      dfloat *MSf = MS+fM*Nfp*Nfp;

      // penalty term just involves face nodes
      for(iint n=0;n<Nfp;++n){
        for(iint m=0;m<Nfp;++m){
          iint idM = eM*Nfp*Nfaces+fM*Nfp+m;
          iint nM = mesh->faceNodes[fM*Nfp+n];
          iint mM = mesh->faceNodes[fM*Nfp+m];
          iint mP  = mesh->vmapP[idM]%Np;

          dfloat MSfnm = sJ*MSf[n*Nfp+m];

          SM[nM*Np+mM] +=  0.5*(1.-bcN)*(1.+bcD)*penalty*MSfnm;
          SP[nM*Np+mP] += -0.5*(1.-bcN)*(1.-bcD)*penalty*MSfnm;
        }
      }

      // now add differential surface terms
      for(iint n=0;n<Nfp;++n){
        for(iint m=0;m<Np;++m){
          iint nM = mesh->faceNodes[fM*Nfp+n];

          for(iint i=0;i<Nfp;++i){
            iint iM = mesh->faceNodes[fM*Nfp+i];
            iint iP = mesh->vmapP[i + fM*Nfp+eM*Nfp*Nfaces]%Np;

            dfloat MSfni = sJ*MSf[n*Nfp+i]; // surface Jacobian built in

            dfloat DxMim = drdx*mesh->Dr[iM*Np+m] + dsdx*mesh->Ds[iM*Np+m];
            dfloat DyMim = drdy*mesh->Dr[iM*Np+m] + dsdy*mesh->Ds[iM*Np+m];

            dfloat DxPim = drdxP*mesh->Dr[iP*Np+m] + dsdxP*mesh->Ds[iP*Np+m];
            dfloat DyPim = drdyP*mesh->Dr[iP*Np+m] + dsdyP*mesh->Ds[iP*Np+m];

            // OP11 = OP11 + 0.5*( - mmE*Dn1)
            SM[nM*Np+m] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
            SM[nM*Np+m] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;

            SP[nM*Np+m] += -0.5*nx*(1-bcD)*(1-bcN)*MSfni*DxPim;
            SP[nM*Np+m] += -0.5*ny*(1-bcD)*(1-bcN)*MSfni*DyPim;
          }
        }
      }

      for(iint n=0;n<Np;++n){
        for(iint m=0;m<Nfp;++m){
          iint mM = mesh->faceNodes[fM*Nfp+m];
          iint mP = mesh->vmapP[m + fM*Nfp+eM*Nfp*Nfaces]%Np;

          for(iint i=0;i<Nfp;++i){
            iint iM = mesh->faceNodes[fM*Nfp+i];

            dfloat MSfim = sJ*MSf[i*Nfp+m];

            dfloat DxMin = drdx*mesh->Dr[iM*Np+n] + dsdx*mesh->Ds[iM*Np+n];
            dfloat DyMin = drdy*mesh->Dr[iM*Np+n] + dsdy*mesh->Ds[iM*Np+n];

            SM[n*Np+mM] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
            SM[n*Np+mM] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;

            SP[n*Np+mP] +=  +0.5*nx*(1-bcD)*(1-bcN)*DxMin*MSfim;
            SP[n*Np+mP] +=  +0.5*ny*(1-bcD)*(1-bcN)*DyMin*MSfim;
          }
        }
      }

      // store non-zeros for off diagonal block
      if (strstr(options,"NONSYM")) { //multiply by inverse mass matrix
        dfloat *tmpSP = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  
        for(iint j=0;j<Np;++j){
          for(iint i=0;i<Np;++i){
            dfloat val = 0;
            for(iint m=0;m<Np;++m){
              tmpSP[i*Np+j] += mesh->invMM[i*Np+m]*SP[m*Np+j];
            }
          }
        }

        for(iint j=0;j<Np;++j){
          for(iint i=0;i<Np;++i){
            SP[i*Np+j] = tmpSP[i*Np+j];
          }
        }
        free(tmpSP);
      }
      for(iint j=0;j<basisNp;++j){
        for(iint i=0;i<basisNp;++i){
          dfloat val = 0;
          for(iint n=0;n<Np;++n){
            for(iint m=0;m<Np;++m){
              val += basis[n*Np+j]*SP[n*Np+m]*basis[m*Np+i];
            }
          }

          if(fabs(val)>tol){
            (*A)[nnz].row = globalIds[eM*Np + j];
            (*A)[nnz].col = globalIds[eP*Np + i];
            (*A)[nnz].val = val;
            (*A)[nnz].ownerRank = rankM;
            ++nnz;
          }
        }
      }
    }
    // store non-zeros for diagonal block
    if (strstr(options,"NONSYM")) { //multiply by inverse mass matrix
      dfloat *tmpSM = (dfloat *) calloc(Np*Np,sizeof(dfloat));
      for(iint j=0;j<Np;++j){
        for(iint i=0;i<Np;++i){
          dfloat val = 0;
          for(iint m=0;m<Np;++m){
            val += mesh->invMM[i*Np+m]*SM[m*Np+j];
          }
          tmpSM[i*Np+j] = val;
        }
      }
      for(iint j=0;j<Np;++j){
        for(iint i=0;i<Np;++i){
          SM[i*Np+j] = tmpSM[i*Np+j];
        }
      }
      free(tmpSM);
    }
    for(iint j=0;j<basisNp;++j){
      for(iint i=0;i<basisNp;++i){
        dfloat val = 0;
        for(iint n=0;n<Np;++n){
          for(iint m=0;m<Np;++m){
            val += basis[n*Np+j]*SM[n*Np+m]*basis[m*Np+i];
          }
        }

        if(fabs(val)>tol){
          (*A)[nnz].row = globalIds[eM*Np + j];
          (*A)[nnz].col = globalIds[eM*Np + i];
          (*A)[nnz].val = val;
          (*A)[nnz].ownerRank = rankM;
          ++nnz;
        }
      }
    }
  }

  printf("nnz = %d\n", nnz);

  qsort((*A), nnz, sizeof(nonZero_t), parallelCompareRowColumn);
  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  *nnzA = nnz;

#if 0
  dfloat* Ap = (dfloat *) calloc(Np*Np*Nelements*Nelements,sizeof(dfloat));
  for (int n=0;n<nnz;n++) {
    int row = (*A)[n].row;
    int col = (*A)[n].col;

    Ap[col+row*Np*Nelements] = (*A)[n].val;
  }

  for (int i=0;i<Np*Nelements;i++) {
    for (int j =0;j<Nelements*Np;j++) {
      printf("%4.2f \t", Ap[j+i*Np*Nelements]);
    }
    printf("\n");
  }
#endif

  free(globalIds);

  free(SM); free(SP);
  free(MS);
}