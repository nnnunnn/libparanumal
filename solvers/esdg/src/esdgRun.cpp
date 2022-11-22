/*

  The MIT License (MIT)

  Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "esdg.hpp"

void esdg_t::Run(){

  int useCheckpointLoad= settings.compareSetting("CHECKPOINT LOAD", "TRUE");
  
  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  if(useCheckpointLoad){
    loadCheckpoint(q, startTime);
    o_q.copyFrom(q);
  }
  
  std::cout << "Testing stability of (E=" << mesh.Nelements << ", N=" << mesh.N
	    << ", Nhat=" << fluxN << ") integrating from T=" << startTime << " to "
	    << finalTime << std::endl;
  
  platform.device.finish();

  if(!useCheckpointLoad){
      
      // do a projection onto nodes
      initialConditionKernel(mesh.Nelements,
			     gamma,
			     startTime,
			     mesh.o_x,
			     mesh.o_y,
			     mesh.o_z,
			     o_cx,
			     o_cy,
			     o_cz,
			     o_q);
  }
  
  timeStepper.Run(*this, o_q, startTime, finalTime);
}