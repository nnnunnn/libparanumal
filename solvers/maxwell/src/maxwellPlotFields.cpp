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

#include "maxwell.hpp"

// interpolate data to plot nodes and save to file (one per process
void maxwell_t::PlotFields(memory<dfloat> Q, const std::string fileName){

  FILE *fp;

  fp = fopen(fileName.c_str(), "w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
          mesh.Nelements*mesh.plotNp,
          mesh.Nelements*mesh.plotNelements);

  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");

  //scratch space for interpolation
  size_t Nscratch = std::max(mesh.Np, mesh.plotNp);
  memory<dfloat> scratch(2*Nscratch);

  memory<dfloat> Ix(mesh.plotNp);
  memory<dfloat> Iy(mesh.plotNp);
  memory<dfloat> Iz(mesh.plotNp);

  // compute plot node coordinates on the fly
  for(dlong e=0;e<mesh.Nelements;++e){
    mesh.PlotInterp(mesh.x + e*mesh.Np, Ix, scratch);
    mesh.PlotInterp(mesh.y + e*mesh.Np, Iy, scratch);
    if(mesh.dim==3)
      mesh.PlotInterp(mesh.z + e*mesh.Np, Iz, scratch);

    if (mesh.dim==2) {
      for(int n=0;n<mesh.plotNp;++n){
        fprintf(fp, "       ");
        fprintf(fp, "%g %g %g\n", Ix[n],Iy[n],0.0);
      }
    } else {
      for(int n=0;n<mesh.plotNp;++n){
        fprintf(fp, "       ");
        fprintf(fp, "%g %g %g\n", Ix[n],Iy[n],Iz[n]);
      }
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");

  memory<dfloat> IHx(mesh.plotNp);
  memory<dfloat> IHy(mesh.plotNp);
  memory<dfloat> IHz(mesh.plotNp);

  // write out magnetic field
  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"MagneticField\" NumberOfComponents=\"%d\" Format=\"ascii\">\n", mesh.dim);
  for(dlong e=0;e<mesh.Nelements;++e){
    mesh.PlotInterp(Q + 0*mesh.Np + e*mesh.Np*Nfields, IHx, scratch);
    mesh.PlotInterp(Q + 1*mesh.Np + e*mesh.Np*Nfields, IHy, scratch);
    if(mesh.dim==3)
      mesh.PlotInterp(Q + 2*mesh.Np + e*mesh.Np*Nfields, IHz, scratch);

    for(int n=0;n<mesh.plotNp;++n){
      fprintf(fp, "       ");
      fprintf(fp, "       ");
      if (mesh.dim==2)
        fprintf(fp, "%f %f\n", IHx[n], IHy[n]);
      else
        fprintf(fp, "%f %f %f\n", IHx[n], IHy[n], IHz[n]);
    }
  }
  fprintf(fp, "       </DataArray>\n");

  memory<dfloat> IEx(mesh.plotNp);
  memory<dfloat> IEy(mesh.plotNp);
  memory<dfloat> IEz(mesh.plotNp);
  
  // write out electric field
  if(mesh.dim==2){
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"ElectricField\" NumberOfComponents=\"%d\" Format=\"ascii\">\n", 1);
  }else{
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"ElectricField\" NumberOfComponents=\"%d\" Format=\"ascii\">\n", mesh.dim);
  }

  for(dlong e=0;e<mesh.Nelements;++e){
    if(mesh.dim==2){
      mesh.PlotInterp(Q + 2*mesh.Np + e*mesh.Np*Nfields, IEz, scratch);
    }
    else{
      mesh.PlotInterp(Q + 3*mesh.Np + e*mesh.Np*Nfields, IEx, scratch);
      mesh.PlotInterp(Q + 4*mesh.Np + e*mesh.Np*Nfields, IEy, scratch);
      mesh.PlotInterp(Q + 5*mesh.Np + e*mesh.Np*Nfields, IEz, scratch);
    }

    for(int n=0;n<mesh.plotNp;++n){
      fprintf(fp, "       ");
      fprintf(fp, "       ");
      if (mesh.dim==2)
        fprintf(fp, "%f\n", IEz[n]);
      else
        fprintf(fp, "%f %f %f\n", IEx[n], IEy[n], IEz[n]);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  fprintf(fp, "     </PointData>\n");

  
  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");

  for(dlong e=0;e<mesh.Nelements;++e){
    for(int n=0;n<mesh.plotNelements;++n){
      fprintf(fp, "       ");
      for(int m=0;m<mesh.plotNverts;++m){
        fprintf(fp, "%d ", e*mesh.plotNp + mesh.plotEToV[n*mesh.plotNverts+m]);
      }
      fprintf(fp, "\n");
    }
  }
  fprintf(fp, "        </DataArray>\n");

  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  dlong cnt = 0;
  for(dlong e=0;e<mesh.Nelements;++e){
    for(int n=0;n<mesh.plotNelements;++n){
      cnt += mesh.plotNverts;
      fprintf(fp, "       ");
      fprintf(fp, "%d\n", cnt);
    }
  }
  fprintf(fp, "       </DataArray>\n");

  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh.Nelements;++e){
    for(int n=0;n<mesh.plotNelements;++n){
      if(mesh.dim==2)
        fprintf(fp, "5\n");
      else
        fprintf(fp, "10\n");
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}
