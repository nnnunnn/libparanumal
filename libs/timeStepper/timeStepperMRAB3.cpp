/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "core.hpp"
#include "timeStepper.hpp"

namespace libp {

namespace TimeStepper {

mrab3::mrab3(dlong Nelements, dlong NhaloElements,
             int _Np, int _Nfields,
             platform_t& _platform, mesh_t& _mesh):
  mrab3(Nelements, 0, NhaloElements,
        _Np, _Nfields, 0,
        _platform, _mesh) {}

mrab3::mrab3(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
             int _Np, int _Nfields, int _Npmlfields,
             platform_t& _platform, mesh_t& _mesh):
  timeStepperBase_t(Nelements, NpmlElements, NhaloElements,
                    _Np, _Nfields, _Npmlfields,
                    _platform, _mesh.comm),
  mesh(_mesh),
  Nlevels(mesh.mrNlevels),
  Nfields(_Nfields),
  Npmlfields(_Npmlfields) {

  //Nstages = 3;

  properties_t kernelInfo = platform.props(); //copy base occa properties from solver

  const int blocksize=256;

  kernelInfo["defines/" "p_blockSize"] = blocksize;
  kernelInfo["defines/" "p_Nstages"] = Nstages;
  kernelInfo["defines/" "p_Np"] = mesh.Np;
  kernelInfo["defines/" "p_Nfp"] = mesh.Nfp;
  kernelInfo["defines/" "p_Nfaces"] = mesh.Nfaces;
  kernelInfo["defines/" "p_Nfields"] = Nfields;
  int maxNodes = std::max(mesh.Np, mesh.Nfp*mesh.Nfaces);
  kernelInfo["defines/" "p_maxNodes"] = maxNodes;

  updateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperMRAB.okl",
                                    "mrabUpdate",
                                    kernelInfo);
  pmlUpdateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                      "timeStepperMRAB.okl",
                                      "mrabPmlUpdate",
                                      kernelInfo);
  traceUpdateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperMRAB.okl",
                                    "mrabTraceUpdate",
                                    kernelInfo);

  // initialize AB time stepping coefficients
  dfloat _ab_a[Nstages*Nstages] = {
                           1.0,      0.0,    0.0,
                         3./2.,   -1./2.,    0.0,
                       23./12., -16./12., 5./12.};
  dfloat _ab_b[(Nstages+1)*Nstages] = {
                         1./2.,      0.0,    0.0,
                         5./8.,   -1./8.,    0.0,
                       17./24.,  -7./24., 2./24.};

  ab_a.malloc(Nstages*Nstages);
  ab_b.malloc(Nstages*Nstages);
  ab_a.copyFrom(_ab_a);
  ab_b.copyFrom(_ab_b);

  h_shiftIndex = platform.hostMalloc<int>(Nlevels);
  o_shiftIndex = platform.malloc<int>(Nlevels);

  mrdt.malloc(Nlevels, 0.0);
  o_mrdt = platform.malloc<dfloat>(mrdt);

  o_ab_a = platform.malloc<dfloat>(ab_a);
  o_ab_b = platform.malloc<dfloat>(ab_b);

  memory<dfloat> zeros(Nstages, 0.0);
  o_zeros = platform.malloc<dfloat>(zeros);
}

void mrab3::Run(solver_t& solver,
                deviceMemory<dfloat> o_q,
                std::optional<deviceMemory<dfloat>> o_pmlq,
                dfloat start, dfloat end) {

  dlong Ntraces = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Nfp
                  *mesh.Nfaces*Nfields;

  /*Pre-reserve memory pool space to avoid some unnecessary re-sizing*/
  platform.reserve<dfloat>(Nstages * N + Nstages * Npml
                           + 5 * platform.memPoolAlignment<dfloat>());

  deviceMemory<dfloat> o_rhsq0 = platform.reserve<dfloat>(N);
  deviceMemory<dfloat> o_rhsq = platform.reserve<dfloat>((Nstages-1)*N);

  deviceMemory<dfloat> o_rhspmlq0 = platform.reserve<dfloat>(Npml);
  deviceMemory<dfloat> o_rhspmlq = platform.reserve<dfloat>((Nstages-1)*Npml);

  deviceMemory<dfloat> o_fQM = platform.reserve<dfloat>(Ntraces);

  dfloat time = start;

  //set timesteps and shifting index
  for (int lev=0;lev<Nlevels;lev++) {
    mrdt[lev] = dt*(1 << lev);
    h_shiftIndex[lev] = 0;
  }
  o_mrdt.copyFrom(mrdt);
  h_shiftIndex.copyTo(o_shiftIndex);

  solver.Report(time,0);

  dfloat outputInterval=0.0;
  solver.settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat outputTime = time + outputInterval;

  // Populate Trace Buffer
  traceUpdateKernel(mesh.mrNelements[Nlevels-1],
                    mesh.o_mrElements[Nlevels-1],
                    mesh.o_mrLevel,
                    mesh.o_vmapM,
                    N,
                    o_shiftIndex,
                    o_mrdt,
                    o_zeros,
                    o_rhsq0,
                    o_rhsq,
                    o_q,
                    o_fQM);

  dfloat DT = dt*(1 << (Nlevels-1));

  int tstep=0;
  int order=0;
  while (time < end) {
    Step(solver,
         o_q, o_pmlq, o_rhsq0, o_rhsq,
         o_rhspmlq0, o_rhspmlq, o_fQM,
         time, dt, order);
    time += DT;
    tstep++;
    if (order<Nstages-1) order++;

    if (time>outputTime) {
      //report state
      solver.Report(outputTime,tstep);
      outputTime += outputInterval;
    }
  }
}

void mrab3::Step(solver_t& solver,
                 deviceMemory<dfloat> o_q,
                 std::optional<deviceMemory<dfloat>> o_pmlq,
                 deviceMemory<dfloat> o_rhsq0,
                 deviceMemory<dfloat> o_rhsq,
                 deviceMemory<dfloat> o_rhspmlq0,
                 deviceMemory<dfloat> o_rhspmlq,
                 deviceMemory<dfloat> o_fQM,
                 dfloat time, dfloat _dt, int order) {

  deviceMemory<dfloat> o_A = o_ab_a+order*Nstages;
  deviceMemory<dfloat> o_B = o_ab_b+order*Nstages;

  for (int Ntick=0; Ntick < (1 << (Nlevels-1));Ntick++) {

    // intermediate stage time
    dfloat currentTime = time + dt*Ntick;

    int lev=0;
    for (;lev<Nlevels-1;lev++)
      if (Ntick % (1<<(lev+1)) != 0) break; //find the max lev to compute rhs

    //evaluate ODE rhs = f(q,t)
    if (o_pmlq.has_value()) {
      solver.rhsf_MR_pml(o_q, o_pmlq.value(),
                         o_rhsq0, o_rhspmlq0,
                         o_fQM, currentTime, lev);
    } else {
      solver.rhsf_MR(o_q, o_rhsq0, o_fQM, currentTime, lev);
    }

    for (lev=0;lev<Nlevels-1;lev++)
      if ((Ntick+1) % (1<<(lev+1)) !=0) break; //find the max lev to update

    // update all elements of level <= lev
    if (mesh.mrNelements[lev])
      updateKernel(mesh.mrNelements[lev],
                   mesh.o_mrElements[lev],
                   mesh.o_mrLevel,
                   mesh.o_vmapM,
                   N,
                   o_shiftIndex,
                   o_mrdt,
                   o_A,
                   o_rhsq0,
                   o_rhsq,
                   o_fQM,
                   o_q);

    if (o_pmlq.has_value()) {
      if (mesh.mrNpmlElements[lev])
        pmlUpdateKernel(mesh.mrNpmlElements[lev],
                       mesh.o_mrPmlElements[lev],
                       mesh.o_mrPmlIds[lev],
                       mesh.o_mrLevel,
                       Npml,
                       Npmlfields,
                       o_shiftIndex,
                       o_mrdt,
                       o_A,
                       o_rhspmlq0,
                       o_rhspmlq,
                       o_pmlq.value());
    }

    //rotate index
    if (Nstages>2)
      for (int l=0; l<=lev; l++)
        h_shiftIndex[l] = (h_shiftIndex[l]+Nstages-2)%(Nstages-1);

    //compute intermediate trace values on lev+1 / lev interface
    if (lev+1<Nlevels && mesh.mrInterfaceNelements[lev+1])
      traceUpdateKernel(mesh.mrInterfaceNelements[lev+1],
                        mesh.o_mrInterfaceElements[lev+1],
                        mesh.o_mrLevel,
                        mesh.o_vmapM,
                        N,
                        o_shiftIndex,
                        o_mrdt,
                        o_B,
                        o_rhsq0,
                        o_rhsq,
                        o_q,
                        o_fQM);

    h_shiftIndex.copyTo(o_shiftIndex, properties_t("async", true));
  }
}

} //namespace TimeStepper

} //namespace libp
