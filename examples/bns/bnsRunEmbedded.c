#include "bns.h"

void bnsRunEmbedded(bns_t *bns, int haloBytes, dfloat * sendBuffer,
	                                     dfloat *recvBuffer, setupAide &options){

	mesh_t *mesh = bns->mesh;

	dfloat safe  = 0.9;   //safety factor
	//error control parameters
	dfloat beta       = 0.05; 
	dfloat factor1    = 0.2;
	dfloat factor2    = 10.0;
	dfloat exp1       = 1.0/bns->rkp -0.75*beta; 
	dfloat invfactor1 = 1.0/factor1;
	dfloat invfactor2 = 1.0/factor2;
	dfloat facold     = 1E-4;

	dfloat hmin = 1e9;
	for(dlong e=0;e<mesh->Nelements;++e){  
		for(int f=0;f<mesh->Nfaces;++f){
			dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
			dfloat sJ   = mesh->sgeo[sid + SJID];
			dfloat invJ = mesh->sgeo[sid + IJID];
			dfloat hest = .5/(sJ*invJ);
			hmin = mymin(hmin, hest);
			}
	}


	// hard code this for the moment
	dfloat outputInterval;
	options.getArgs("OUTPUT INTERVAL", outputInterval);

	dfloat nextOutputTime = outputInterval;
	dfloat outputNumber = 0;

	//initial time
	bns->time = 0.0;
	bns->tstep = 0;
	bns->atstep = 0; 
	bns->rtstep = 0;  

	if(bns->reportFlag)
		bnsReport(bns, bns->time, options);

	int done =0;
	while (!done) {

		if (bns->dt<bns->dtMIN){
			printf("ERROR: Time step became too small at time step=%d\n", bns->tstep);
			exit (-1);
		}
		if (isnan(bns->dt)) {
			printf("ERROR: Solution became unstable at time step=%d\n", bns->tstep);
			exit (-1);
		}


		//check for final timestep
		if (bns->time+bns->dt > bns->finalTime){
			bns->dt = bns->finalTime-bns->time;
			done = 1;
		}

    occaTimerTic(mesh->device, "SARK_STEP"); 
  	bnsSARKStep(bns, bns->time, haloBytes, sendBuffer, recvBuffer, options);
  	occaTimerToc(mesh->device, "SARK_STEP"); 
		
    

    occaTimerTic(mesh->device, "SARK_ERROR"); 
		//Error estimation 
		dlong Ntotal = mesh->Nelements*mesh->Np*bns->Nfields;
		// printf("Ntotal: %d \t %d\n", Ntotal, bns->Nblock);
		bns->errorEstimateKernel(Ntotal, 
				                    bns->ATOL,
				                    bns->RTOL,
				                    bns->o_q,
				                    bns->o_rkq,
				                    bns->o_rkerr,
				                    bns->o_errtmp);

		bns->o_errtmp.copyTo(bns->errtmp);
		bns->o_rkerr.copyTo(bns->rkerr);


		dfloat maxerr = 0.f;
		for(dlong e=0; e<mesh->Nelements; e++){
			for(int n =0;n<mesh->Np;n++){
				const dlong id = e*mesh->Np + n; 
          maxerr= mymax(maxerr,bns->rkerr[id]);
       }
		}

		// printf("Max err:\t%.5e\n",maxerr);
		dfloat localerr = 0, err = 0;

		for(int n=0;n<bns->Nblock;++n){
			localerr += bns->errtmp[n];
			// printf("Error\t:%02d\t%.5e\n", n, localerr);
		}

		MPI_Allreduce(&localerr, &err, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
   
    // printf("Error\t:\t%.5e\t%.5e\t\n", err, bns->dt);
		err = sqrt(err/(bns->totalElements*mesh->Np));
		
    occaTimerToc(mesh->device, "SARK_ERROR"); 




		dfloat fac1 = pow(err,exp1);
		dfloat fac  = fac1/pow(facold,beta);

		fac = mymax(invfactor2, mymin(invfactor1,fac/safe));
		dfloat dtnew = bns->dt/fac;

		if(err<1.0){
      

			if(bns->reportFlag){
				occaTimerTic(mesh->device, "SARK_OUTPUT"); 
				// check for output during this step and do a mini-step
				if(bns->time<nextOutputTime && bns->time+bns->dt>nextOutputTime){
					dfloat savedt = bns->dt;			  
					// save rkq
					bns->o_saveq.copyFrom(bns->o_rkq);
					if(bns->pmlFlag){
						bns->o_saveqx.copyFrom(bns->o_rkqx);
						bns->o_saveqy.copyFrom(bns->o_rkqy);
						if(bns->dim==3)
							bns->o_saveqz.copyFrom(bns->o_rkqz);
					}

					// change dt to match output
					bns->dt = nextOutputTime-bns->time;
					// print
					printf("Taking output mini step: %g\n", bns->dt);

					// if(options.compareArgs("TIME INTEGRATOR","SARK"))  // SA Adaptive RK 
					bnsSARKStep(bns, bns->time, haloBytes, sendBuffer, recvBuffer, options);


					// shift for output
					bns->o_rkq.copyTo(bns->o_q);
					// output  (print from rkq)
					bnsReport(bns, nextOutputTime, options);
					// restore time step
					bns->dt = savedt;
					// increment next output time
					nextOutputTime += outputInterval;

					// accept saved rkq
					bns->o_q.copyFrom(bns->o_saveq);
					if(bns->pmlFlag){
						bns->o_pmlqx.copyFrom(bns->o_saveqx);
						bns->o_pmlqy.copyFrom(bns->o_saveqy);
						if(bns->dim==3)
							bns->o_pmlqz.copyFrom(bns->o_saveqz);
					}
				} 
				else{
					// Accept rkq
					bns->o_q.copyFrom(bns->o_rkq);
					if(bns->pmlFlag){
						bns->o_pmlqx.copyFrom(bns->o_rkqx);
						bns->o_pmlqy.copyFrom(bns->o_rkqy);	
						if(bns->dim==3)
							bns->o_pmlqz.copyFrom(bns->o_rkqz);
					}	
				}
				occaTimerToc(mesh->device, "SARK_OUTPUT"); 
			}
			else{
				// Accept rkq
				bns->o_q.copyFrom(bns->o_rkq);
				if(bns->pmlFlag){
					bns->o_pmlqx.copyFrom(bns->o_rkqx);
					bns->o_pmlqy.copyFrom(bns->o_rkqy);	
					if(bns->dim==3)
						bns->o_pmlqz.copyFrom(bns->o_rkqz);
				}
			}	

			facold = mymax(err,1E-4);
			bns->time += bns->dt;

			printf("\r time = %g (%d), dt = %g accepted (ratio dt/hmin = %g)               ", bns->time, bns->atstep, bns->dt, bns->dt/hmin);
			bns->tstep++;
		}
		else{
			bns->rtstep++; 
			dtnew = bns->dt/(mymax(invfactor1,fac1/safe));
			printf("\r time = %g (%d), dt = %g rejected (ratio dt/min = %g), trying %g", bns->time,bns->atstep, bns->dt, bns->dt/hmin, dtnew);
			done =0;
		}

		bns->dt = dtnew;
		// bns->dt = bns->dt;
		bns->atstep++;

		bnsSAADRKCoefficients(bns, options);

		char fname[BUFSIZ]; sprintf(fname, "boltzmannAddaptiveDt.dat");
		FILE *fp; fp = fopen(fname, "a");
		fprintf(fp, "%.5e %.5e\n", bns->time, bns->dt); 
		fclose(fp);
	}


}


