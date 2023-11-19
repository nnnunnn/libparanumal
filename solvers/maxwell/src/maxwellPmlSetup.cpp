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

// report ranges for scale factors
#define NBN_SHOW_SCALE_RANGES 1

//#########################################################
// NBN: toggle use of {beta_a} vs {beta_ab}
#define NBN_USE_BETA_AB_FOR_3D 1
// from page 8 of pdf notes:
// beta_ab = alpha_a + sigma_a/(kappa_a*kappa_b)

// NBN: switch pmlWeights and pmlInvWeights
#define NBN_USE_INVERTED_WEIGHTS 1
//#########################################################



void maxwell_t::PmlSetup(){

  //build pml element lists
  mesh.PmlSetup();

  pmlcubature=0;

  int sigmaOrder = 1;
  int kappaOrder = 1;
  int alphaOrder = 1;

  dfloat sigmaXmax = 0, sigmaYmax = 0, sigmaZmax = 0;
  dfloat kappaXmax = 0, kappaYmax = 0, kappaZmax = 0;
  dfloat alphaXmax = 0, alphaYmax = 0, alphaZmax = 0;

  //set up the pml
  if (mesh.NpmlElements) {

    //get settings from solver
    settings.getSetting("PML SIGMA PROFILE ORDER", sigmaOrder);
    settings.getSetting("PML SIGMAX MAX", sigmaXmax);
    settings.getSetting("PML SIGMAY MAX", sigmaYmax);
    settings.getSetting("PML SIGMAZ MAX", sigmaZmax);

    settings.getSetting("PML KAPPA PROFILE ORDER", kappaOrder);
    settings.getSetting("PML KAPPAX MAX", kappaXmax);
    settings.getSetting("PML KAPPAY MAX", kappaYmax);
    settings.getSetting("PML KAPPAZ MAX", kappaZmax);

    settings.getSetting("PML ALPHA PROFILE ORDER", alphaOrder);
    settings.getSetting("PML ALPHAX MAX", alphaXmax);
    settings.getSetting("PML ALPHAY MAX", alphaYmax);
    settings.getSetting("PML ALPHAZ MAX", alphaZmax);

    pmlcubature = (settings.compareSetting("PML INTEGRATION", "CUBATURE")) ? 1:0;

    //find the bounding box of the whole domain and interior domain
    dfloat xmin = 1e9, xmax =-1e9;
    dfloat ymin = 1e9, ymax =-1e9;
    dfloat zmin = 1e9, zmax =-1e9;
    dfloat pmlxmin = 1e9, pmlxmax =-1e9;
    dfloat pmlymin = 1e9, pmlymax =-1e9;
    dfloat pmlzmin = 1e9, pmlzmax =-1e9;

    //for all elements
    for (dlong e=0;e<mesh.Nelements;++e) {
      for (int n=0;n<mesh.Nverts;++n) {
        dfloat x =0, y=0, z = 0;
        x  = mesh.EX[e*mesh.Nverts+n];
        y = mesh.EY[e * mesh.Nverts + n];

        pmlxmin = std::min(pmlxmin, x);
        pmlxmax = std::max(pmlxmax, x);
        pmlymin = std::min(pmlymin, y);
        pmlymax = std::max(pmlymax, y);

        if (mesh.dim == 3) {
          z = mesh.EZ[e * mesh.Nverts + n];
          pmlzmin = std::min(pmlzmin, z);
          pmlzmax = std::max(pmlzmax, z);
        }
      }
    }

    //for interior elements
    for (dlong m=0;m<mesh.NnonPmlElements;++m) {
      dlong e = mesh.nonPmlElements[m];
      for (int n=0;n<mesh.Nverts;++n) {
        dfloat x=0., y=0., z=0.;
        x = mesh.EX[e*mesh.Nverts+n];
        y = mesh.EY[e*mesh.Nverts+n];

        xmin = std::min(xmin, x);
        xmax = std::max(xmax, x);
        ymin = std::min(ymin, y);
        ymax = std::max(ymax, y);

        if (mesh.dim == 3) {
          z = mesh.EZ[e*mesh.Nverts+n];
          zmin = std::min(zmin, z);
          zmax = std::max(zmax, z);
        }
      }
    }

    dfloat sigmaXmaxScale = pow(fabs(pmlxmax-xmax),sigmaOrder);
    dfloat sigmaXminScale = pow(fabs(pmlxmin-xmin),sigmaOrder);
    dfloat sigmaYmaxScale = pow(fabs(pmlymax-ymax),sigmaOrder);
    dfloat sigmaYminScale = pow(fabs(pmlymin-ymin),sigmaOrder);
    dfloat sigmaZmaxScale = 0.f, sigmaZminScale = 0.f;
    if(mesh.dim==3){
      sigmaZmaxScale = pow(fabs(pmlzmax-zmax),sigmaOrder);
      sigmaZminScale = pow(fabs(pmlzmin-zmin),sigmaOrder);
    }

    dfloat kappaXmaxScale = pow(fabs(pmlxmax-xmax),kappaOrder);
    dfloat kappaXminScale = pow(fabs(pmlxmin-xmin),kappaOrder);
    dfloat kappaYmaxScale = pow(fabs(pmlymax-ymax),kappaOrder);
    dfloat kappaYminScale = pow(fabs(pmlymin-ymin),kappaOrder);
    dfloat kappaZmaxScale = 0.f, kappaZminScale = 0.f;
    if(mesh.dim==3){
      kappaZmaxScale = pow(fabs(pmlzmax-zmax),kappaOrder);
      kappaZminScale = pow(fabs(pmlzmin-zmin),kappaOrder);
    }

    dfloat alphaXmaxScale = pow(fabs(pmlxmax-xmax),alphaOrder);
    dfloat alphaXminScale = pow(fabs(pmlxmin-xmin),alphaOrder);
    dfloat alphaYmaxScale = pow(fabs(pmlymax-ymax),alphaOrder);
    dfloat alphaYminScale = pow(fabs(pmlymin-ymin),alphaOrder);
    dfloat alphaZmaxScale = 0.f, alphaZminScale = 0.f;
    if(mesh.dim==3){
      alphaZmaxScale = pow(fabs(pmlzmax-zmax),alphaOrder);
      alphaZminScale = pow(fabs(pmlzmin-zmin),alphaOrder);
    }

    // Set the size of pml nodes
    int pmlNp = (pmlcubature) ? mesh.cubNp : mesh.Np;
    int pmlNq = (pmlcubature) ? mesh.cubNq : mesh.Nq;

    //###################################
    // Note: for Tet3D, pmlNq is not set
    //###################################

    memory<dfloat> pmlr, pmls, pmlt;
    if(pmlcubature){
      pmlr = mesh.cubr;
      pmls = mesh.cubs;
      pmlt = mesh.cubt;
    } else {
      pmlr = mesh.r;
      pmls = mesh.s;
      pmlt = mesh.t;
    }

    // printf("Setting PML Coefficient \n");
    //set up damping parameter
    pmlSigma.malloc(mesh.dim*mesh.NpmlElements*pmlNp, 0.0);
    pmlKappa.malloc(mesh.dim*mesh.NpmlElements*pmlNp, 1.0);
    pmlAlpha.malloc(mesh.dim*mesh.NpmlElements*pmlNp, 1.0);

    if (2 == mesh.dim) {
      // NBN: beta_a = alpha_a + sigma_a/kappa_a
      pmlBeta.malloc(mesh.dim*mesh.NpmlElements*pmlNp, 0.0);
    } 
    else {
#if (NBN_USE_BETA_AB_FOR_3D)
      // NBN: from page 8 of pdf notes:
      // NBN: beta_ab = alpha_a + sigma_a/(kappa_a*kappa_b)
      pmlBeta.malloc(2*mesh.dim * mesh.NpmlElements*pmlNp, 0.0);
#else
      pmlBeta.malloc(  mesh.dim * mesh.NpmlElements*pmlNp, 0.0);
#endif
    }

    pmlInvWeights.malloc(3*mesh.NpmlElements*pmlNp, 0.0);
    pmlWeights.malloc(3*mesh.NpmlElements*pmlNp, 0.0);

    for (dlong m=0;m<mesh.NpmlElements;++m){
      dlong e     = mesh.pmlElements[m];
      hlong type  = mesh.elementInfo[e];

      //element vertices
      memory<dfloat> xe = mesh.EX + e*mesh.Nverts;
      memory<dfloat> ye = mesh.EY + e*mesh.Nverts;
      memory<dfloat> ze = mesh.EZ + e*mesh.Nverts;

      for(int n=0;n<pmlNp;++n){  // for each pml node
        dfloat x=0.f, y=0.f, z=0.f;
        dfloat rn=0.f, sn=0.f, tn=0.f;
        if(mesh.elementType==Mesh::TRIANGLES){
          rn = pmlr[n];
          sn = pmls[n];

          x = -0.5*(rn+sn)*xe[0] + 0.5*(1+rn)*xe[1] + 0.5*(1+sn)*xe[2];
          y = -0.5*(rn+sn)*ye[0] + 0.5*(1+rn)*ye[1] + 0.5*(1+sn)*ye[2];
        } else if(mesh.elementType==Mesh::QUADRILATERALS){
          const int i = n%pmlNq;
          const int j = n/pmlNq;
          rn = pmlr[i];
          sn = pmlr[j];

          x =  0.25*( (1.0-rn)*(1-sn)*xe[0]+(1.0-rn)*(1+sn)*xe[1]+(1.0+rn)*(1+sn)*xe[2]+(1.0+rn)*(1-sn)*xe[3]);
          y =  0.25*( (1.0-rn)*(1-sn)*ye[0]+(1.0-rn)*(1+sn)*ye[1]+(1.0+rn)*(1+sn)*ye[2]+(1.0+rn)*(1-sn)*ye[3]);
        } else if(mesh.elementType==Mesh::TETRAHEDRA){
          rn = pmlr[n];
          sn = pmls[n];
          tn = pmlt[n];

          x = -0.5*(rn+sn+tn+1)*xe[0] + 0.5*(1+rn)*xe[1] + 0.5*(1+sn)*xe[2] + 0.5*(tn+1)*xe[3];
          y = -0.5*(rn+sn+tn+1)*ye[0] + 0.5*(1+rn)*ye[1] + 0.5*(1+sn)*ye[2] + 0.5*(tn+1)*ye[3];
          z = -0.5*(rn+sn+tn+1)*ze[0] + 0.5*(1+rn)*ze[1] + 0.5*(1+sn)*ze[2] + 0.5*(tn+1)*ze[3];
        } else if(mesh.elementType==Mesh::HEXAHEDRA){
          const int i = n%pmlNq;
          const int j = (n/pmlNq)%pmlNq;
          const int k = (n/pmlNq)/pmlNq;
          rn = pmlr[i];
          sn = pmlr[j];
          tn = pmlr[k];

          x =
            +0.125*(1-rn)*(1-sn)*(1-tn)*xe[0]
            +0.125*(1+rn)*(1-sn)*(1-tn)*xe[1]
            +0.125*(1+rn)*(1+sn)*(1-tn)*xe[2]
            +0.125*(1-rn)*(1+sn)*(1-tn)*xe[3]
            +0.125*(1-rn)*(1-sn)*(1+tn)*xe[4]
            +0.125*(1+rn)*(1-sn)*(1+tn)*xe[5]
            +0.125*(1+rn)*(1+sn)*(1+tn)*xe[6]
            +0.125*(1-rn)*(1+sn)*(1+tn)*xe[7];

          y =
            +0.125*(1-rn)*(1-sn)*(1-tn)*ye[0]
            +0.125*(1+rn)*(1-sn)*(1-tn)*ye[1]
            +0.125*(1+rn)*(1+sn)*(1-tn)*ye[2]
            +0.125*(1-rn)*(1+sn)*(1-tn)*ye[3]
            +0.125*(1-rn)*(1-sn)*(1+tn)*ye[4]
            +0.125*(1+rn)*(1-sn)*(1+tn)*ye[5]
            +0.125*(1+rn)*(1+sn)*(1+tn)*ye[6]
            +0.125*(1-rn)*(1+sn)*(1+tn)*ye[7];

          z =
            +0.125*(1-rn)*(1-sn)*(1-tn)*ze[0]
            +0.125*(1+rn)*(1-sn)*(1-tn)*ze[1]
            +0.125*(1+rn)*(1+sn)*(1-tn)*ze[2]
            +0.125*(1-rn)*(1+sn)*(1-tn)*ze[3]
            +0.125*(1-rn)*(1-sn)*(1+tn)*ze[4]
            +0.125*(1+rn)*(1-sn)*(1+tn)*ze[5]
            +0.125*(1+rn)*(1+sn)*(1+tn)*ze[6]
            +0.125*(1-rn)*(1+sn)*(1+tn)*ze[7];
        }

        dlong baseid = mesh.dim*pmlNp*m + n;

        if (type == 100 || type == 300 || type == 500 || type == 700) {
          if (x >= xmax) {
            pmlSigma[baseid + 0*pmlNp] =     sigmaXmax * pow(fabs(x-xmax), sigmaOrder) / sigmaXmaxScale;
            pmlKappa[baseid + 0*pmlNp] = 1 + kappaXmax * pow(fabs(x-xmax), kappaOrder) / kappaXmaxScale;
            pmlAlpha[baseid + 0*pmlNp] = 1 + alphaXmax * pow(fabs(x-xmax), alphaOrder) / alphaXmaxScale;
          } else if (x <= xmin) {
            pmlSigma[baseid + 0*pmlNp] =     sigmaXmax * pow(fabs(xmin-x), sigmaOrder) / sigmaXminScale;
            pmlKappa[baseid + 0*pmlNp] = 1 + kappaXmax * pow(fabs(xmin-x), kappaOrder) / kappaXminScale;
            pmlAlpha[baseid + 0*pmlNp] = 1 + alphaXmax * pow(fabs(xmin-x), alphaOrder) / alphaXminScale;
          }
        }
        else if (type == 200 || type == 300 || type == 600 || type == 700) {
          if (y >= ymax) {
            pmlSigma[baseid + 1*pmlNp] =     sigmaYmax * pow(fabs(y-ymax), sigmaOrder) / sigmaYmaxScale;
            pmlKappa[baseid + 1*pmlNp] = 1 + kappaYmax * pow(fabs(y-ymax), kappaOrder) / kappaYmaxScale;
            pmlAlpha[baseid + 1*pmlNp] = 1 + alphaYmax * pow(fabs(y-ymax), alphaOrder) / alphaYmaxScale;
          } else if (y <= ymin) {
            pmlSigma[baseid + 1*pmlNp] =     sigmaYmax * pow(fabs(ymin-y), sigmaOrder) / sigmaYminScale;
            pmlKappa[baseid + 1*pmlNp] = 1 + kappaYmax * pow(fabs(ymin-y), kappaOrder) / kappaYminScale;
            pmlAlpha[baseid + 1*pmlNp] = 1 + alphaYmax * pow(fabs(ymin-y), alphaOrder) / alphaYminScale;
          }
        }

        if (mesh.dim == 3) {
          if (type == 400 || type == 500 || type == 600 || type == 700) {
            if (z >= zmax) {
              pmlSigma[baseid + 2*pmlNp] =     sigmaZmax * pow(fabs(z-zmax), sigmaOrder) / sigmaZmaxScale;
              pmlKappa[baseid + 2*pmlNp] = 1 + kappaZmax * pow(fabs(z-zmax), kappaOrder) / kappaZmaxScale;
              pmlAlpha[baseid + 2*pmlNp] = 1 + alphaZmax * pow(fabs(z-zmax), alphaOrder) / alphaZmaxScale;
            } else if (z <= zmin) {
              pmlSigma[baseid + 2*pmlNp] =     sigmaZmax * pow(fabs(zmin-z), sigmaOrder) / sigmaZminScale;
              pmlKappa[baseid + 2*pmlNp] = 1 + kappaZmax * pow(fabs(zmin-z), kappaOrder) / kappaZminScale;
              pmlAlpha[baseid + 2*pmlNp] = 1 + alphaZmax * pow(fabs(zmin-z), alphaOrder) / alphaZminScale;
            }
          }
        }
      }
    }

#if (NBN_SHOW_SCALE_RANGES)
    dfloat minSigmaX = 1e9, maxSigmaX = -1e9;
    dfloat minSigmaY = 1e9, maxSigmaY = -1e9;
    dfloat minSigmaZ = 1e9, maxSigmaZ = -1e9;
    for (dlong m=0;m<mesh.NpmlElements;++m) {
      for (int n=0;n<pmlNp; ++n) {  // for each node
        minSigmaX = std::min(minSigmaX, pmlSigma[mesh.dim*pmlNp*m + 0*pmlNp + n]);
        maxSigmaX = std::max(maxSigmaX, pmlSigma[mesh.dim*pmlNp*m + 0*pmlNp + n]);
        minSigmaY = std::min(minSigmaY, pmlSigma[mesh.dim*pmlNp*m + 1*pmlNp + n]);
        maxSigmaY = std::max(maxSigmaY, pmlSigma[mesh.dim*pmlNp*m + 1*pmlNp + n]);
        if (3 == mesh.dim) {
          minSigmaZ = std::min(minSigmaZ, pmlSigma[mesh.dim*pmlNp*m + 2*pmlNp + n]);
          maxSigmaZ = std::max(maxSigmaZ, pmlSigma[mesh.dim*pmlNp*m + 2*pmlNp + n]);
        }
      }
    }
    nnLOG(1, "SigmaX in range [%0.2f, %0.2f]\n", minSigmaX, maxSigmaX);
    nnLOG(1, "SigmaY in range [%0.2f, %0.2f]\n", minSigmaY, maxSigmaY);
    if (3 == mesh.dim) { nnLOG(1, "SigmaZ in range [%0.2f, %0.2f]\n", minSigmaZ, maxSigmaZ); }

    dfloat minKappaX = 1e9, maxKappaX = -1e9;
    dfloat minKappaY = 1e9, maxKappaY = -1e9;
    dfloat minKappaZ = 1e9, maxKappaZ = -1e9;
    for (dlong m=0;m<mesh.NpmlElements;++m) {
      for (int n=0;n<pmlNp;++n) {   // for each node
        minKappaX = std::min(minKappaX, pmlKappa[mesh.dim*pmlNp*m + 0*pmlNp + n]);
        maxKappaX = std::max(maxKappaX, pmlKappa[mesh.dim*pmlNp*m + 0*pmlNp + n]);
        minKappaY = std::min(minKappaY, pmlKappa[mesh.dim*pmlNp*m + 1*pmlNp + n]);
        maxKappaY = std::max(maxKappaY, pmlKappa[mesh.dim*pmlNp*m + 1*pmlNp + n]);
        if (3 == mesh.dim) {
          minKappaZ = std::min(minKappaZ, pmlKappa[mesh.dim*pmlNp*m + 2*pmlNp + n]);
          maxKappaZ = std::max(maxKappaZ, pmlKappa[mesh.dim*pmlNp*m + 2*pmlNp + n]);
        }
      }
    }
    nnLOG(1, "KappaX in range [%0.2f, %0.2f]\n", minKappaX, maxKappaX);
    nnLOG(1, "KappaY in range [%0.2f, %0.2f]\n", minKappaY, maxKappaY);
    if (3 == mesh.dim) { nnLOG(1, "KappaZ in range [%0.2f, %0.2f]\n", minKappaZ, maxKappaZ); }

    dfloat minAlphaX = 1e9, maxAlphaX = -1e9;
    dfloat minAlphaY = 1e9, maxAlphaY = -1e9;
    dfloat minAlphaZ = 1e9, maxAlphaZ = -1e9;
    for (dlong m=0;m<mesh.NpmlElements;++m) {
      for (int n=0;n<pmlNp;++n) {   // for each node
        minAlphaX = std::min(minAlphaX, pmlAlpha[mesh.dim*pmlNp*m + 0*pmlNp + n]);
        maxAlphaX = std::max(maxAlphaX, pmlAlpha[mesh.dim*pmlNp*m + 0*pmlNp + n]);
        minAlphaY = std::min(minAlphaY, pmlAlpha[mesh.dim*pmlNp*m + 1*pmlNp + n]);
        maxAlphaY = std::max(maxAlphaY, pmlAlpha[mesh.dim*pmlNp*m + 1*pmlNp + n]);
        if (3 == mesh.dim) {
          minAlphaZ = std::min(minAlphaZ, pmlAlpha[mesh.dim*pmlNp*m + 2*pmlNp + n]);
          maxAlphaZ = std::max(maxAlphaZ, pmlAlpha[mesh.dim*pmlNp*m + 2*pmlNp + n]);
        }
      }
    }
    nnLOG(1, "AlphaX in range [%0.2f, %0.2f]\n", minAlphaX, maxAlphaX);
    nnLOG(1, "AlphaY in range [%0.2f, %0.2f]\n", minAlphaY, maxAlphaY);
    if (3 == mesh.dim) { nnLOG(1, "AlphaZ in range [%0.2f, %0.2f]\n", minAlphaZ, maxAlphaZ); }
#endif // (NBN_SHOW_SCALE_RANGES)
    
    for(dlong m=0;m<mesh.NpmlElements;++m){

      if(m < 1000) nnTRC(1, "\n *** pmlElement %6d\n", (int)m);

      for(int n=0;n<pmlNp;++n){   // for each pml node
        dlong base = mesh.dim*pmlNp*m + n;
        dfloat kappax = pmlKappa[base + 0*pmlNp];
        dfloat kappay = pmlKappa[base + 1*pmlNp];
        dfloat kappaz = (mesh.dim == 3) ? pmlKappa[base + 2*pmlNp] : 1.;

        dfloat alphax = pmlAlpha[base + 0*pmlNp];
        dfloat alphay = pmlAlpha[base + 1*pmlNp];
        dfloat alphaz = (mesh.dim == 3) ? pmlAlpha[base + 2*pmlNp] : 1.;

        dfloat sigmax = pmlSigma[base + 0*pmlNp];
        dfloat sigmay = pmlSigma[base + 1*pmlNp];
        dfloat sigmaz = (mesh.dim == 3) ? pmlSigma[base + 2*pmlNp] : 1.;

        dlong wbase = 3*pmlNp*m + n;

#if (NBN_USE_INVERTED_WEIGHTS)
        pmlWeights[wbase+0*pmlNp] = (kappax*kappax) / (kappax*kappay*kappaz);
        pmlWeights[wbase+1*pmlNp] = (kappay*kappay) / (kappax*kappay*kappaz);
        pmlWeights[wbase+2*pmlNp] = (kappaz*kappaz) / (kappax*kappay*kappaz);

        pmlInvWeights[wbase+0*pmlNp] = 1.0 / pmlWeights[wbase+0*pmlNp];
        pmlInvWeights[wbase+1*pmlNp] = 1.0 / pmlWeights[wbase+1*pmlNp];
        pmlInvWeights[wbase+2*pmlNp] = 1.0 / pmlWeights[wbase+2*pmlNp];
#else
        pmlInvWeights[wbase + 0*pmlNp] = (kappax*kappax) / (kappax*kappay*kappaz);
        pmlInvWeights[wbase + 1*pmlNp] = (kappay*kappay) / (kappax*kappay*kappaz);
        pmlInvWeights[wbase + 2*pmlNp] = (kappaz*kappaz) / (kappax*kappay*kappaz);

        pmlWeights[wbase + 0*pmlNp] = 1.0 / pmlInvWeights[wbase + 0*pmlNp];
        pmlWeights[wbase + 1*pmlNp] = 1.0 / pmlInvWeights[wbase + 1*pmlNp];
        pmlWeights[wbase + 2*pmlNp] = 1.0 / pmlInvWeights[wbase + 2*pmlNp];
#endif

        if (2==mesh.dim){
          dlong bbase = mesh.dim*pmlNp*m + n;
          pmlBeta[bbase + 0 * pmlNp] = alphax + sigmax / kappax;
          pmlBeta[bbase + 1 * pmlNp] = alphay + sigmay / kappay;
        }
        else {
#if (NBN_USE_BETA_AB_FOR_3D)
          //################################################################
          // NBN: [from pdf p.8]
          // beta_ab = alpha_a + sigma_a/(kappa_a*kappa_b)
          dlong bbase = 6*pmlNp*m + n;                            // beta_ab
          pmlBeta[bbase + 0*pmlNp] = alphax + sigmax/(kappax*kappay); //  xy
          pmlBeta[bbase + 1*pmlNp] = alphax + sigmax/(kappax*kappaz); //  xz

          pmlBeta[bbase + 2*pmlNp] = alphay + sigmay/(kappay*kappax); //  yx
          pmlBeta[bbase + 3*pmlNp] = alphay + sigmay/(kappay*kappaz); //  yz

          pmlBeta[bbase + 4*pmlNp] = alphaz + sigmaz/(kappaz*kappax); //  zx
          pmlBeta[bbase + 5*pmlNp] = alphaz + sigmaz/(kappaz*kappay); //  zy
          //################################################################
#else
          dlong bbase = 3*pmlNp*m + n;
          pmlBeta[bbase + 0*pmlNp] = alphax + sigmax/kappax;
          pmlBeta[bbase + 1*pmlNp] = alphay + sigmay/kappay;
          pmlBeta[bbase + 2*pmlNp] = alphaz + sigmaz/kappaz;
#endif
        }

#if (1)
        if (m<1000) {
          nnTRC(1, "kappa=(%5.1f, %5.1f, %5.1f), alpha=(%5.1f, %5.1f, %5.1f), sigma=(%5.1f, %5.1f, %5.1f) ... ",
            kappax, kappay, kappaz, alphax, alphay, alphaz, sigmax, sigmay, sigmaz);

          nnTRC(1, "Wts=(%6.2f, %6.2f, %6.2f) ... ", 
            pmlWeights[wbase+0*pmlNp],pmlWeights[wbase+1*pmlNp],pmlWeights[wbase+2*pmlNp]);

          if (2 == mesh.dim) {
            dlong bbase = 2*pmlNp*m + n;
            nnTRC(1, "beta=(%5.1f, %5.1f)\n", pmlBeta[bbase+0*pmlNp], pmlBeta[bbase+1*pmlNp]);
          }else{
            dlong bbase = 6*pmlNp*m + n;
            nnTRC(1, "beta=(%5.1f, %5.1f) (%5.1f, %5.1f) (%5.1f, %5.1f)\n",
                      pmlBeta[bbase+0*pmlNp], pmlBeta[bbase+1*pmlNp], pmlBeta[bbase+2*pmlNp],
                      pmlBeta[bbase+3*pmlNp], pmlBeta[bbase+4*pmlNp], pmlBeta[bbase+5*pmlNp]);
          }
        }
#endif
      }
    }

    if (mesh.NpmlElements){      
      o_pmlSigma = platform.malloc<dfloat>(pmlSigma);
      o_pmlBeta  = platform.malloc<dfloat>(pmlBeta);
      o_pmlInvWeights = platform.malloc<dfloat>(pmlInvWeights);
      o_pmlWeights = platform.malloc<dfloat>(pmlWeights);
    }
  }
}
