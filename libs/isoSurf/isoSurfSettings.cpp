/*

The MIT License (MIT)

Copyright (c) 2017-2023 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "isoSurf.hpp"

using namespace libp;

isoSurfSettings_t::isoSurfSettings_t(const comm_t& _comm) :
  settings_t(_comm) {

  newSetting("ISOSURFACE MODE", 
            "SURFACE", 
            "Extract 3D surface or planar slices");

  newSetting("ISOSURFACE FILE NAME", 
            "iso", 
            "Output name for isosurface plots");

  newSetting("ISOSURFACE FIELD ID", 
             "3", 
             "Field for isosurfing");

  newSetting("ISOSURFACE COLOR ID", 
             "4", 
             "Field for coloring isosurface");
  
  newSetting("ISOSURFACE NUM LEVELS", 
             "1", 
             "Number of isosurface levels");
  
  newSetting("ISOSURFACE MAX NUMTRIS", 
             "5e6", 
             "Max num tris from isosurface kernel");

  newSetting("ISOSURFACE CONTOUR MIN", 
             "1.0", 
             "Isosurface min value");
  
  newSetting("ISOSURFACE CONTOUR MAX", 
             "1.0", 
             "Isosurface max value");
  
  newSetting("ISOSURFACE VOXEL SIZE", 
             "0.001", 
             "grid voxel size for merging small triangles");
  
  newSetting("ISOSURFACE VERTEX TOL", 
             "0.001", 
             "tolerance for merging close vertices");
}



void isoSurfSettings_t::report() {
  if (comm.rank() == 0) {
    std::cerr << "isoSurf Settings:\n\n";

    reportSetting("ISOSURFACE MODE");
    reportSetting("ISOSURFACE FILE NAME");

    reportSetting("ISOSURFACE FIELD ID");
    reportSetting("ISOSURFACE COLOR ID");
    reportSetting("ISOSURFACE NUM LEVELS");
    reportSetting("ISOSURFACE MAX NUMTRIS");

    reportSetting("ISOSURFACE CONTOUR MIN");
    reportSetting("ISOSURFACE CONTOUR MAX");
    reportSetting("ISOSURFACE VOXEL SIZE");
    reportSetting("ISOSURFACE VERTEX TOL");
  }
}
