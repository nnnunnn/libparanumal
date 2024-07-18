/*

The MIT License (MIT)

Copyright (c) 2017-2024 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "cns.hpp"

// cns_t member functions defined in this file:
// 
//  cns_t::saveCheckpoint(...)
//  cns_t::loadCheckpoint(...)
//  cns_t::changeOrder(...)


// When calculating Ndofs, how to count NhaloElems?
// 
// NhaloElems = mesh.totalHaloPairs;  // #faces shared with halo elements?
// NhaloElems = mesh.NhaloElements;   // actual number of halo elements?

#define USE_TOTAL_HALO_PAIRS 1        // Note: totalHaloPairs >= NhaloElements

//---------------------------------------------------------
// Re: changeOrder(...)
//---------------------------------------------------------
// Intended use case - raise/lower the order of solution 
// fields when restarting a run from a checkpoint.
// 
// Note: 
// Current version assumes that the restart solution was 
// calculated on the same mesh and partition to be used 
// in the new run, and that same "dfloat" type is used.
// 
// TODO: 
// - enable changing the number of ranks
// - handle conversion between dfloat types


// Use struct restart_t to save/load items of interest
// as (binary) header data in checkpont file.
// 
// Note: "sim_time" is used to set start-time for new run
// Note: "sim_N" is used to raise/lower order of solution 

typedef struct {
  int     sim_prec;   // dfloat type of solution (fp32/fp64)
  int     sim_N;      // polynomial order of solution data
  double  sim_time;   // sim time of solution data
  size_t  sim_tstep;  // time step of solution data
} restart_t;


void cns_t::saveCheckpoint(const memory<dfloat>&outQ, dfloat t, int tstep) {

  if (! settings.compareSetting("CHECKPOINT SAVE", "TRUE")) { return; }

  restart_t rstData;
  rstData.sim_prec = (sizeof(dfloat) == sizeof(float)) ? 32 : 64;
  rstData.sim_N = mesh.N;
  rstData.sim_time = t;
  rstData.sim_tstep = g_tstep;

  static int s_count=0;
  stdS name, fname;
  settings.getSetting("CHECKPOINT SAVE NAME", name);
  fname = nnSTR("%s/%s_%d_%06d_t_%0.2lf.rst", g_restartDir.c_str(), name.c_str(), mesh.rank, s_count++, g_time);

  // write header info, then field data
  FILE *fp = fopen(fname.c_str(), "wb");  // binary
  LIBP_ABORT("failed to open checkpoint file: [" << fname << "].",
              fp==NULL);

  // write header data
  size_t nwrite = fwrite(&rstData, sizeof(restart_t), 1, fp);
  LIBP_ABORT("failed to write header to checkpoint file: [" << fname << "].", 
              nwrite!=1);

#if (USE_TOTAL_HALO_PAIRS)
  dlong NhaloElems = mesh.totalHaloPairs; // #faces shared with halo elements?
#else
  dlong NhaloElems = mesh.NhaloElements;  // actual number of halo elements?
#endif

  // if (mesh.size>1), data includes halo fields
  dlong NlocalFields = mesh.Nelements * (mesh.Np * Nfields);
  dlong NhaloFields = NhaloElems * (mesh.Np * Nfields);
  size_t writeLen = NlocalFields + NhaloFields;

  // write solution field data
  nwrite = fwrite(outQ.ptr(), sizeof(dfloat), writeLen, fp);
  fclose(fp);
  LIBP_ABORT("error writing " << writeLen << " data to checkpoint file: [" << fname << "].", 
              nwrite!=writeLen);

  nnMSG(1, "[%d:%d] saved checkpoint %s for time %0.3e (order %d, step %d)\n",
            mesh.rank, mesh.size, fname.c_str(), (double)t, mesh.N, g_tstep);
}


void cns_t::loadCheckpoint(dfloat& t_old) {

  if (!settings.compareSetting("CHECKPOINT LOAD", "TRUE")) { return; }

  // NBN: globals used by GUI driver:
  // g_restartDir  : directory for restart/checkpoint files
  // g_restartPath : full path/name of restart file
  // g_restartName : name of restart file

  // TODO: add menu option
#ifdef _MSC_VER
  g_restartDir = "D:/TW/bin/restart";
#else
  g_restartDir = "./restart";
#endif

  settings.getSetting("CHECKPOINT LOAD NAME", g_restartName);

  if (mesh.size > 1) {
    // build filename for this rank,
    // rst_T6 --> rst_T6_p#.rst
    g_restartPath = nnSTR("%s/%s_p%d.rst", g_restartDir.c_str(), g_restartName.c_str(), mesh.rank);
  } else {
    // rst_T6 --> rst_T6.rst
    g_restartPath = nnSTR("%s/%s.rst", g_restartDir.c_str(), g_restartName.c_str());
  }
  nnTRC(1, "[%d:%d] g_restartPath: %s\n", mesh.rank, mesh.size, g_restartPath.c_str());

  // read header info, then field data
  FILE *fp = fopen(g_restartPath.c_str(), "rb");  // binary
  LIBP_ABORT("failed to open checkpoint file: [" << g_restartPath << "].",
              fp==NULL);

  // read checkpoint info
  restart_t rstData;
  size_t nread = fread(&rstData, sizeof(restart_t), 1, fp);
  LIBP_ABORT("failed to read restart header from file: [" << g_restartPath << "].",
              nread!=1);

  // caller uses restart time as new start time
  t_old = (dfloat)rstData.sim_time;
  int N_old = rstData.sim_N;
  int tstep_old = rstData.sim_tstep;

  // To calculate the size of restart solution, we need
  // Np_old correspondng to order N of the restart data.

  int Np_old = 0;
  switch (mesh.elementType) {
  case Mesh::TRIANGLES:      Np_old = ((N_old+1)*(N_old+2)) / 2; break;
  case Mesh::QUADRILATERALS: Np_old =  (N_old+1)*(N_old+1); break;
  case Mesh::TETRAHEDRA:     Np_old = ((N_old+1)*(N_old+2)*(N_old+3)) / 6; break;
  case Mesh::HEXAHEDRA:      Np_old =  (N_old+1)*(N_old+1)*(N_old+1); break;
  default: 
    LIBP_FORCE_ABORT("unexpected element type: [" << mesh.elementType << "]");
  }

#if (USE_TOTAL_HALO_PAIRS)
  dlong NhaloElems = mesh.totalHaloPairs; // #faces shared with halo elements?
#else
  dlong NhaloElems = mesh.NhaloElements;  // actual number of halo elements?
#endif

  // if (mesh.size>1), data includes halo fields
  dlong NlocalFields = mesh.Nelements * (Np_old * Nfields);
  dlong NhaloFields  = NhaloElems * (Np_old * Nfields);
  size_t readLen     = NlocalFields + NhaloFields;

  // resize host buffer
  this->q.malloc(readLen);

  // read solution field data
  nread = fread(q.ptr(), sizeof(dfloat), readLen, fp);
  fclose(fp);

  LIBP_ABORT("error reading " << readLen << " data from checkpoint file: [" << g_restartPath << "].",
              nread!=readLen);

  if (N_old != mesh.N) {
    // interpolate o_q to new order in device memory
    changeOrder(N_old, Np_old, mesh.N, mesh.Np);
  } 
  else {
    // copy old solution directly to device memory
    o_q = platform.malloc<dfloat>(q);
  }

  nnMSG(1, "[%d:%d] loaded checkpoint %s for time %0.3e (order %d, step %d)\n",
            mesh.rank, mesh.size, g_restartPath.c_str(), (double)t_old, N_old, tstep_old);

  if (N_old != mesh.N) {
    nnMSG(1, "[%d:%d] interpolated checkpoint solution from order %d to %d\n",
              mesh.rank, mesh.size, N_old, mesh.N);
  }
}


void cns_t::changeOrder(int N_old, int Np_old, int N_new, int Np_new){

  // change order of restart solution from N_old to N_new
  // 
  // {N_old, Np_old}: {N, Np} from previous run
  // {N_new, Np_new}: {N, Np} for new run

  if (N_old == N_new) {
    nnMSG(1, "cns_t::changeOrder(%d,%d) called with same coarse/fine order?\n", N_old, N_new);
    return;
  }

  // Sequence:
  //-------------------------------------------------------
  // Old solution is currently in cns::q
  // 
  // 1. select "coarse-to-fine" or "fine-to-coarse"
  // 2. build interpolation operator
  // 3. load old q into device mem o_qOld
  // 4. build and run kernel, returning new o_q
  // 5. resize host buffer q to match size of new o_q
  //-------------------------------------------------------

  // 1. select "coarse-to-fine" or "fine-to-coarse"
  int N_coarse, Np_coarse, N_fine, Np_fine;
  if (N_old < N_new) { 
     N_coarse =  N_old;  N_fine =  N_new;   // "coarse-to-fine"
    Np_coarse = Np_old; Np_fine = Np_new;
  }
  else {
    N_fine  =  N_old;  N_coarse =  N_new;   // "fine-to-coarse"
    Np_fine = Np_old; Np_coarse = Np_new;
  }

  // 2. build interpolation operator
  memory<dfloat> qInterp;
  if (N_old < N_new) {
    switch (mesh.elementType) {
    case Mesh::TRIANGLES:  mesh.DegreeRaiseMatrixTri2D(N_coarse, N_fine, qInterp); break;
    case Mesh::TETRAHEDRA: mesh.DegreeRaiseMatrixTet3D(N_coarse, N_fine, qInterp); break;
    default:
      LIBP_FORCE_ABORT("cns_t::changeOrder - TODO: add element type: [" << mesh.elementType << "]");
    }
  }
  else {
    switch (mesh.elementType) {
    case Mesh::TRIANGLES:  mesh.DegreeLowerMatrixTri2D(N_fine, N_coarse, qInterp); break;
    case Mesh::TETRAHEDRA: mesh.DegreeLowerMatrixTet3D(N_fine, N_coarse, qInterp); break;
    default:
      LIBP_FORCE_ABORT("cns_t::changeOrder - TODO: add element type: [" << mesh.elementType << "]");
    }
  }

  deviceMemory<dfloat> o_qInterp = platform.malloc<dfloat>(qInterp);

#if (USE_TOTAL_HALO_PAIRS)
  dlong NhaloElems = mesh.totalHaloPairs; // #faces shared with halo elements?
#else
  dlong NhaloElems = mesh.NhaloElements;  // actual number of halo elements?
#endif

  // if (mesh.size>1), data includes halo fields
  dlong NlocalFields = mesh.Nelements * (Np_new * Nfields);
  dlong NhaloFields  = NhaloElems     * (Np_new * Nfields);
  dlong totalNdofs   = NlocalFields + NhaloFields;

  // number of elements handled by this rank
  dlong totalNel = mesh.Nelements + NhaloElems;

  // 3. load old q into device mem o_qOld
  deviceMemory<dfloat> o_qOld;
  o_qOld = platform.malloc<dfloat>(q);  // copy old solution to device mem

  std::string suffix = mesh.elementSuffix();
  std::string oklFilePrefix = DCNS "okl/";
  std::string oklFileSuffix = ".okl";

  // add properties for COARSE to FINE kernel
  properties_t& kernelInfo = mesh.props;
  kernelInfo["defines/" "p_NpCoarse"] = Np_coarse;
  kernelInfo["defines/" "p_NpFine"]   = Np_fine;

  // 4. build and run kernel, returning new o_q
  if (N_old < N_new) {
    //-----------------------------
    // interpolate "coarse-to-fine"
    //-----------------------------
    std::string fileName = oklFilePrefix + "cnsInterpCoarseToFine" + oklFileSuffix;
    std::string kernelName = "cnsInterpCoarseToFine";
    coarseToFineKernel = platform.buildKernel(fileName, kernelName, kernelInfo);
    coarseToFineKernel(totalNel, o_qInterp, o_qOld, o_q);
  }
  else {
    //-----------------------------
    // interpolate "fine-to-coarse"
    //-----------------------------
    std::string fileName = oklFilePrefix + "cnsInterpFineToCoarse" + oklFileSuffix;
    std::string kernelName = "cnsInterpFineToCoarse";
    fineToCoarseKernel = platform.buildKernel(fileName, kernelName, kernelInfo);
    fineToCoarseKernel(totalNel, o_qInterp, o_qOld, o_q);
  }

  // 5. resize host buffer q to match size of new o_q
  q.realloc(totalNdofs);
}
