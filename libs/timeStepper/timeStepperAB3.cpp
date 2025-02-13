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

/* Adams Bashforth, order 3 */
ab3::ab3(dlong Nelements, dlong NhaloElements,
         int Np, int Nfields,
         platform_t& _platform, comm_t _comm):
  ab3(Nelements, 0, NhaloElements,
      Np, Nfields, 0,
      _platform, _comm) {}

/* Adams Bashforth, order 3 */
ab3::ab3(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
         int Np, int Nfields, int Npmlfields,
         platform_t& _platform, comm_t _comm):
  timeStepperBase_t(Nelements, NpmlElements, NhaloElements,
                    Np, Nfields, Npmlfields,
                    _platform, _comm) {

  //Nstages = 3;
  shiftIndex = 0;

  properties_t kernelInfo = platform.props(); //copy base occa properties from solver

  const int blocksize=256;

  kernelInfo["defines/" "p_blockSize"] = blocksize;
  kernelInfo["defines/" "p_Nstages"] = Nstages;

  updateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperAB.okl",
                                    "abUpdate",
                                    kernelInfo);

  // initialize AB time stepping coefficients
  dfloat _ab_a[Nstages*Nstages] = {
                           1.0,      0.0,    0.0,
                         3./2.,   -1./2.,    0.0,
                       23./12., -16./12., 5./12.};

  ab_a.malloc(Nstages*Nstages);
  ab_a.copyFrom(_ab_a);

  o_ab_a = platform.malloc<dfloat>(ab_a);
}

void ab3::Run(solver_t& solver,
              deviceMemory<dfloat> o_q,
              std::optional<deviceMemory<dfloat>> o_pmlq,
              dfloat start, dfloat end) {

  /*Pre-reserve memory pool space to avoid some unnecessary re-sizing*/
  platform.reserve<dfloat>(Nstages*N + Nstages*Npml
                           + 2 * platform.memPoolAlignment<dfloat>());

  deviceMemory<dfloat> o_rhsq = platform.reserve<dfloat>(Nstages*N);
  deviceMemory<dfloat> o_rhspmlq = platform.reserve<dfloat>(Nstages*Npml);

  dfloat time = start;

  solver.Report(time,0);

  dfloat outputInterval=0.0;
  solver.settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat outputTime = time + outputInterval;

  int tstep=0;
  int order=0;
  while (time < end) {
    Step(solver, o_q, o_rhsq, o_pmlq, o_rhspmlq, time, dt, order);
    time += dt;
    tstep++;
    if (order<Nstages-1) order++;

    if (time>outputTime) {
      //report state
      solver.Report(time,tstep);
      outputTime += outputInterval;
    }
  }
}

void ab3::Step(solver_t& solver,
               deviceMemory<dfloat> o_q,
               deviceMemory<dfloat> o_rhsq,
               std::optional<deviceMemory<dfloat>> o_pmlq,
               deviceMemory<dfloat> o_rhspmlq,
               dfloat time, dfloat _dt, int order) {

  //rhs at current index
  deviceMemory<dfloat> o_rhsq0 = o_rhsq + shiftIndex*N;
  deviceMemory<dfloat> o_rhspmlq0 = o_rhspmlq + shiftIndex*Npml;

  //A coefficients at current order
  deviceMemory<dfloat> o_A = o_ab_a + order*Nstages;

  //evaluate ODE rhs = f(q,t)
  if (o_pmlq.has_value()) {
    solver.rhsf_pml(o_q, o_pmlq.value(), o_rhsq0, o_rhspmlq0, time);
  } else {
    solver.rhsf(o_q, o_rhsq0, time);
  }

  //update q
  updateKernel(N,
               dt,
               shiftIndex,
               o_A,
               o_rhsq,
               o_q);

  //update pmlq
  if (o_pmlq.has_value()) {
    updateKernel(Npml,
                 dt,
                 shiftIndex,
                 o_A,
                 o_rhspmlq,
                 o_pmlq.value());
  }

  //rotate index
  shiftIndex = (shiftIndex+Nstages-1)%Nstages;
}

} //namespace TimeStepper

} //namespace libp
