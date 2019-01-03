#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"
#include "parAlmond.h"
#include "ellipticPreconQuad2D.h"

// block size for reduction (hard coded)
#define blockSize 256

typedef struct {

  mesh_t *mesh;

  precon_t *precon;

  ogs_t *ogs;

  ogs_t *ogsDg;

  // C0 halo gather-scatter info
  ogs_t *halo;

  // C0 nonhalo gather-scatter info
  ogs_t *nonHalo;

  char *type;

  iint Nblock;

  dfloat tau;

  // HOST shadow copies
  dfloat *Ax, *p, *r, *z, *Ap, *tmp, *grad;

  dfloat *sendBuffer, *recvBuffer;

  occa::stream defaultStream;
  occa::stream dataStream;

  occa::memory o_p; // search direction
  occa::memory o_z; // preconditioner solution
  occa::memory o_res;
  occa::memory o_Sres;
  occa::memory o_Ax; // A*initial guess
  occa::memory o_Ap; // A*search direction
  occa::memory o_tmp; // temporary
  occa::memory o_grad; // temporary gradient storage (part of A*)
  occa::memory o_rtmp;
  occa::memory o_invDegree;
  occa::memory o_EToB;

  // list of elements that are needed for global gather-scatter
  iint NglobalGatherElements;
  occa::memory o_globalGatherElementList;

  // list of elements that are not needed for global gather-scatter
  iint NlocalGatherElements;
  occa::memory o_localGatherElementList;

  occa::kernel AxKernel;
  occa::kernel partialAxKernel;
  occa::kernel rhsBCKernel;
  occa::kernel innerProductKernel;
  occa::kernel weightedInnerProduct1Kernel;
  occa::kernel weightedInnerProduct2Kernel;
  occa::kernel scaledAddKernel;
  occa::kernel dotMultiplyKernel;
  occa::kernel dotDivideKernel;

  occa::kernel gradientKernel;
  occa::kernel ipdgKernel;
  occa::kernel partialGradientKernel;
  occa::kernel partialIpdgKernel;
  occa::kernel rhsBCIpdgKernel;

}solver_t;


void ellipticSetupQuad2D(mesh2D *mesh, occa::kernelInfo &kernelInfo);

void ellipticParallelGatherScatterQuad2D(mesh2D *mesh, ogs_t *ogs, occa::memory &o_v, occa::memory &o_gsv,
          const char *type, const char *op);

void ellipticPreconditionerSetupQuad2D(solver_t *solver, ogs_t *ogs, dfloat tau, dfloat lambda, iint *BCType, const char *options, const char *parAlmondOptions);

void diagnostic(int N, occa::memory &o_x, const char *message);

void ellipticPreconditioner2D(solver_t *solver, dfloat lambda, occa::memory &o_r,occa::memory &o_z,const char *options);

int ellipticSolveQuad2D(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const char *options);

solver_t *ellipticSolveSetupQuad2D(mesh_t *mesh, dfloat tau, dfloat lambda, iint *BCType, occa::kernelInfo &kernelInfo, const char *options, const char *parAlmondOptions);

solver_t *ellipticBuildMultigridLevelQuad2D(solver_t *baseSolver, int* levelDegrees, int n, const char *options);

void ellipticStartHaloExchange2D(solver_t *solver, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer);

void ellipticInterimHaloExchange2D(solver_t *solver, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer);

void ellipticEndHaloExchange2D(solver_t *solver, occa::memory &o_q, dfloat *recvBuffer);

void ellipticParallelGatherScatterSetup(mesh_t *mesh,    // provides DEVICE
          iint Nlocal,     // number of local nodes
          iint Nbytes,     // number of bytes per node
          iint *gatherLocalIds,  // local index of nodes
          iint *gatherBaseIds,   // global index of their base nodes
          iint *gatherHaloFlags,
          ogs_t **halo,
          ogs_t **nonHalo);   // 1 for halo node, 0 for not

//smoother setups
void ellipticSetupSmootherOverlappingPatchIpdg(solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int *BCType, const char *options);
void ellipticSetupSmootherDampedJacobiIpdg    (solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int* BCType, const char *options);
void ellipticSetupSmootherFullPatchIpdg (solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int *BCType, dfloat rateTolerance, const char *options);
void ellipticSetupSmootherFacePatchIpdg (solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int *BCType, dfloat rateTolerance, const char *options);
void ellipticSetupSmootherLocalPatchIpdg(solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int *BCType, dfloat rateTolerance, const char *options);


void ellipticMultiGridSetupQuad2D(solver_t *solver, precon_t* precon, dfloat tau, dfloat lambda, iint *BCType, const char *options, const char *parAlmondOptions);
void ellipticSetupSmootherQuad2D(solver_t *solver, precon_t *precon,
                                dfloat tau, dfloat lambda, int* BCType,
                                const char *options);
dfloat maxEigSmoothAx(solver_t* solver, agmgLevel *level);