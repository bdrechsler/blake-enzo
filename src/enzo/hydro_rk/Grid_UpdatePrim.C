/***********************************************************************
/
/  GRID CLASS (UPDATE HYDRO VARIABLES)
/
/  written by: Peng Wang
/  date:       April, 2007
/  modified1:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"
#include "EOS.h"

// extern "C" void FORTRAN_NAME(krome_renormref)(float* specabund, float* ref);
extern "C" void FORTRAN_NAME(krome_initref)(float* d,
   float* De, float* CHI, float* OI, float* HNCI, 
   float* HCNI, float* H2I, float* CI, float* HI, 
   float* H2OI, float* OHI, float* O2I, float* CH2I, 
   float* H2COI, float* HCOI, float* MGI, float* NH3I, 
   float* NOI, float* CNI, float* COI, float* N2I, 
   float* NH2I, float* CH3I, float* CH4I, float* NI, 
   float* NHI, float* HeI, float* HNOI, float* CH3OHI, 
   float* CO2I, float* H2CNI, float* HNCOI, float* NO2I, 
   float* O2HI, float* OCNI, float* CH3OH_DUSTI, float* HNCO_DUSTI, 
   float* H2CO_DUSTI, float* CH4_DUSTI, float* CO_DUSTI, float* H2O_DUSTI, 
   float* NO_DUSTI, float* CO2_DUSTI, float* N2_DUSTI, float* HCN_DUSTI, 
   float* NH3_DUSTI, float* O2_DUSTI, float* NO2_DUSTI, float* HNO_DUSTI, 
   float* O2H_DUSTI, float* H2CN_DUSTI, float* MG_DUSTI, float* HNC_DUSTI, 
   float* E_DUSTI, float* HCOII, float* HII, float* HOCII, 
   float* CII, float* CH2II, float* CHII, float* H2COII, 
   float* MGII, float* NH3II, float* NOII, float* CNII, 
   float* COII, float* N2II, float* O2II, float* H2OII, 
   float* NH2II, float* OII, float* OHII, float* CH3II, 
   float* CH4II, float* NII, float* HCNII, float* NHII, 
   float* H2II, float* HeII, float* HNOII, float* H2NOII, 
   float* H3II, float* H3COII, float* H3OII, float* HCNHII, 
   float* HCO2II, float* HeHII, float* N2HII, float* O2HII, 
   int* in, int* jn, int* kn, int* idim, 
   int* is, int* js, int* ks, int* ie, int* je, int* ke, 
   float* ATOM_C, float* ATOM_E, float* ATOM_H, float* ATOM_MG,
   float* ATOM_O, float* ATOM_N, float* ATOM_HE);

extern "C" void FORTRAN_NAME(krome_renormref)(float* d,
   float* De, float* CHI, float* OI, float* HNCI, 
   float* HCNI, float* H2I, float* CI, float* HI, 
   float* H2OI, float* OHI, float* O2I, float* CH2I, 
   float* H2COI, float* HCOI, float* MGI, float* NH3I, 
   float* NOI, float* CNI, float* COI, float* N2I, 
   float* NH2I, float* CH3I, float* CH4I, float* NI, 
   float* NHI, float* HeI, float* HNOI, float* CH3OHI, 
   float* CO2I, float* H2CNI, float* HNCOI, float* NO2I, 
   float* O2HI, float* OCNI, float* CH3OH_DUSTI, float* HNCO_DUSTI, 
   float* H2CO_DUSTI, float* CH4_DUSTI, float* CO_DUSTI, float* H2O_DUSTI, 
   float* NO_DUSTI, float* CO2_DUSTI, float* N2_DUSTI, float* HCN_DUSTI, 
   float* NH3_DUSTI, float* O2_DUSTI, float* NO2_DUSTI, float* HNO_DUSTI, 
   float* O2H_DUSTI, float* H2CN_DUSTI, float* MG_DUSTI, float* HNC_DUSTI, 
   float* E_DUSTI, float* HCOII, float* HII, float* HOCII, 
   float* CII, float* CH2II, float* CHII, float* H2COII, 
   float* MGII, float* NH3II, float* NOII, float* CNII, 
   float* COII, float* N2II, float* O2II, float* H2OII, 
   float* NH2II, float* OII, float* OHII, float* CH3II, 
   float* CH4II, float* NII, float* HCNII, float* NHII, 
   float* H2II, float* HeII, float* HNOII, float* H2NOII, 
   float* H3II, float* H3COII, float* H3OII, float* HCNHII, 
   float* HCO2II, float* HeHII, float* N2HII, float* O2HII, 
   int* in, int* jn, int* kn, int* idim, 
   int* is, int* js, int* ks, int* ie, int* je, int* ke, 
   float* ATOM_C, float* ATOM_E, float* ATOM_H, float* ATOM_MG,
   float* ATOM_O, float* ATOM_N, float* ATOM_HE);

int grid::UpdatePrim(float **dU, float c1, float c2)
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int i, j, k, n, n_dU, dim, igrid, field, size, activesize;
  for (dim = 0, size = 1; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  for (dim = 0, activesize = 1; dim < GridRank; dim++) {
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);
  }

  float *D, *sum;
  float SmallX = 1e-20;

  /*
  if ( (NSpecies+NColor) > 0) {
    D = new float[size];
    sum = new float[size];
    for (i = 0; i < size; i++) {
      D[i] = 0.0;
      sum[i] = 0.0;
    }
  }
  */

  // ORIGINAL
  if ( (NSpecies+NColor) > 0) {
    D = new float[activesize];
    sum = new float[activesize];
    for (i = 0; i < activesize; i++) {
      D[i] = 0.0;
      sum[i] = 0.0;
    }
  }

  float *Prim[NEQ_HYDRO+NSpecies+NColor];
  float *OldPrim[NEQ_HYDRO+NSpecies+NColor];
  this->ReturnHydroRKPointers(Prim, false);         //##### species in Prim[] are already fractions fractionalized in Grid_RK2_[12]Step -> ReturnHydroRKPointer, thus no need for ReturnMassFraction
  this->ReturnOldHydroRKPointers(OldPrim, false);   //##### whereas species in OldPrim[] are not

  //##### Want to mix species and colors for renormalization?  Normally you don't
  int NSpecies_renorm;

  if (MixSpeciesAndColors) 
    NSpecies_renorm = NSpecies+NColor;
  else  if (NoMultiSpeciesButColors) {
    NSpecies_renorm = NSpecies;
  }
  else
    switch (MultiSpecies) {  //update pure species! not colours!
    case 0:  NSpecies_renorm = 0;  break;
    case 1:  NSpecies_renorm = 5;  break;
    case 2:  NSpecies_renorm = 8;  break;
    case 3:  NSpecies_renorm = 11; break;
#ifdef USE_KROME
    case KROMESPECIES:  NSpecies_renorm = NKROMESPECIES; break;
#endif
    default: NSpecies_renorm = 0;  break;
    }

#ifdef USE_KROME
  float *atomabund[NKROMEATOMS];
  for (i=0; i<NKROMEATOMS; i++) {
    atomabund[i] = new float[size];
  }

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum;
  int CHINum, OINum, HNCINum, HCNINum, CINum, H2OINum,
    OHINum, O2INum, CH2INum, H2COINum, HCOINum,
    MGINum, NH3INum, NOINum, CNINum, COINum,
    N2INum, NH2INum, CH3INum, CH4INum, NINum,
    NHINum, HNOINum, CH3OHINum, CO2INum, H2CNINum,
    HNCOINum, NO2INum, O2HINum, OCNINum, CH3OH_DUSTINum,
    HNCO_DUSTINum, H2CO_DUSTINum, CH4_DUSTINum,
    CO_DUSTINum, H2O_DUSTINum, NO_DUSTINum, CO2_DUSTINum,
    N2_DUSTINum, HCN_DUSTINum, NH3_DUSTINum,
    O2_DUSTINum, NO2_DUSTINum, HNO_DUSTINum,
    O2H_DUSTINum, H2CN_DUSTINum, MG_DUSTINum,
    HNC_DUSTINum, E_DUSTINum, HCOIINum, HOCIINum,
    CIINum, CH2IINum, CHIINum, H2COIINum, MGIINum,
    NH3IINum, NOIINum, CNIINum, COIINum, N2IINum,
    O2IINum, H2OIINum, NH2IINum, OIINum, OHIINum,
    CH3IINum, CH4IINum, NIINum, HCNIINum, NHIINum,
    HNOIINum, H2NOIINum, H3IINum, H3COIINum,
    H3OIINum, HCNHIINum, HCO2IINum, HeHIINum,
    N2HIINum, O2HIINum;

  if (use_kromeconserve) {
    if (IdentifySpeciesFieldsKrome(
      DeNum, CHINum, OINum, HNCINum, HCNINum, H2INum,
      CINum, HINum, H2OINum, OHINum, O2INum, CH2INum,
      H2COINum, HCOINum, MGINum, NH3INum, NOINum,
      CNINum, COINum, N2INum, NH2INum, CH3INum,
      CH4INum, NINum, NHINum, HeINum, HNOINum,
      CH3OHINum, CO2INum, H2CNINum, HNCOINum, NO2INum,
      O2HINum, OCNINum, CH3OH_DUSTINum, HNCO_DUSTINum,
      H2CO_DUSTINum, CH4_DUSTINum, CO_DUSTINum,
      H2O_DUSTINum, NO_DUSTINum, CO2_DUSTINum,
      N2_DUSTINum, HCN_DUSTINum, NH3_DUSTINum,
      O2_DUSTINum, NO2_DUSTINum, HNO_DUSTINum,
      O2H_DUSTINum, H2CN_DUSTINum, MG_DUSTINum,
      HNC_DUSTINum, E_DUSTINum, HCOIINum, HIINum,
      HOCIINum, CIINum, CH2IINum, CHIINum, H2COIINum,
      MGIINum, NH3IINum, NOIINum, CNIINum, COIINum,
      N2IINum, O2IINum, H2OIINum, NH2IINum, OIINum,
      OHIINum, CH3IINum, CH4IINum, NIINum, HCNIINum,
      NHIINum, H2IINum, HeIINum, HNOIINum, H2NOIINum,
      H3IINum, H3COIINum, H3OIINum, HCNHIINum,
      HCO2IINum, HeHIINum, N2HIINum, O2HIINum
      ) == FAIL) {
            ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }

    FORTRAN_NAME(krome_initref)(
      BaryonField[iden], 
      BaryonField[DeNum], BaryonField[CHINum], BaryonField[OINum],
      BaryonField[HNCINum], BaryonField[HCNINum],
      BaryonField[H2INum], BaryonField[CINum],
      BaryonField[HINum], BaryonField[H2OINum],
      BaryonField[OHINum], BaryonField[O2INum],
      BaryonField[CH2INum], BaryonField[H2COINum],
      BaryonField[HCOINum], BaryonField[MGINum],
      BaryonField[NH3INum], BaryonField[NOINum],
      BaryonField[CNINum], BaryonField[COINum],
      BaryonField[N2INum], BaryonField[NH2INum],
      BaryonField[CH3INum], BaryonField[CH4INum],
      BaryonField[NINum], BaryonField[NHINum],
      BaryonField[HeINum], BaryonField[HNOINum],
      BaryonField[CH3OHINum], BaryonField[CO2INum],
      BaryonField[H2CNINum], BaryonField[HNCOINum],
      BaryonField[NO2INum], BaryonField[O2HINum],
      BaryonField[OCNINum], BaryonField[CH3OH_DUSTINum],
      BaryonField[HNCO_DUSTINum], BaryonField[H2CO_DUSTINum],
      BaryonField[CH4_DUSTINum], BaryonField[CO_DUSTINum],
      BaryonField[H2O_DUSTINum], BaryonField[NO_DUSTINum],
      BaryonField[CO2_DUSTINum], BaryonField[N2_DUSTINum],
      BaryonField[HCN_DUSTINum], BaryonField[NH3_DUSTINum],
      BaryonField[O2_DUSTINum], BaryonField[NO2_DUSTINum],
      BaryonField[HNO_DUSTINum], BaryonField[O2H_DUSTINum],
      BaryonField[H2CN_DUSTINum], BaryonField[MG_DUSTINum],
      BaryonField[HNC_DUSTINum], BaryonField[E_DUSTINum],
      BaryonField[HCOIINum], BaryonField[HIINum],
      BaryonField[HOCIINum], BaryonField[CIINum],
      BaryonField[CH2IINum], BaryonField[CHIINum],
      BaryonField[H2COIINum], BaryonField[MGIINum],
      BaryonField[NH3IINum], BaryonField[NOIINum],
      BaryonField[CNIINum], BaryonField[COIINum],
      BaryonField[N2IINum], BaryonField[O2IINum],
      BaryonField[H2OIINum], BaryonField[NH2IINum],
      BaryonField[OIINum], BaryonField[OHIINum],
      BaryonField[CH3IINum], BaryonField[CH4IINum],
      BaryonField[NIINum], BaryonField[HCNIINum],
      BaryonField[NHIINum], BaryonField[H2IINum],
      BaryonField[HeIINum], BaryonField[HNOIINum],
      BaryonField[H2NOIINum], BaryonField[H3IINum],
      BaryonField[H3COIINum], BaryonField[H3OIINum],
      BaryonField[HCNHIINum], BaryonField[HCO2IINum],
      BaryonField[HeHIINum], BaryonField[N2HIINum],
      BaryonField[O2HIINum], 
      GridDimension, GridDimension+1, GridDimension+2, 
      &GridRank, GridStartIndex, GridStartIndex+1, GridStartIndex+2, 
      GridEndIndex, GridEndIndex+1, GridEndIndex+2,
      atomabund[0], atomabund[1], atomabund[2], atomabund[3],
      atomabund[4], atomabund[5], atomabund[6]);
  }
#endif

  // update species
  /*
  for (field = NEQ_HYDRO; field < NEQ_HYDRO+NSpecies_renorm; field++) { 
    n = 0;
    n_dU = 0;
    for (k = 0; k < GridDimension[2]; k++) {
      for (j = 0; j < GridDimension[1]; j++) {
	igrid = (k * GridDimension[1] + j) * GridDimension[0];
        for (i = 0; i < GridDimension[0]; i++, n++, igrid++) {
	  // dU exists only for active region
          if (i >= GridStartIndex[0] && i <= GridEndIndex[0] &&
	      j >= GridStartIndex[1] && j <= GridEndIndex[1] &&
	      k >= GridStartIndex[2] && k <= GridEndIndex[2]) 
	    Prim[field][igrid] = c1*OldPrim[field][igrid] +
	      (1-c1)*Prim[field][igrid]*Prim[iden][igrid] + c2*dU[field][n_dU++];
          D[n] += Prim[field][igrid];
        }
      }
    }
  }

  // renormalize species

  for (field = NEQ_HYDRO; field < NEQ_HYDRO+NSpecies_renorm; field++) {
    n = 0;
    for (k = 0; k < GridDimension[2]; k++) {
      for (j = 0; j < GridDimension[1]; j++) {
	igrid = (k * GridDimension[1] + j) * GridDimension[0];
        for (i = 0; i < GridDimension[0]; i++, n++, igrid++) {
          Prim[field][igrid] = min(1.0, max((Prim[field][igrid]/D[n]), SmallX));
	  Prim[field][igrid] = Prim[field][igrid]/D[n];
          sum[n] += Prim[field][igrid];
        }
      }
    }
  }

  for (field = NEQ_HYDRO; field < NEQ_HYDRO+NSpecies_renorm; field++) {
    n = 0;
    for (k = 0; k < GridDimension[2]; k++) {
      for (j = 0; j < GridDimension[1]; j++) {
	igrid = (k * GridDimension[1] + j) * GridDimension[0];
        for (i = 0; i < GridDimension[0]; i++, n++, igrid++)
          Prim[field][igrid] /= sum[n];
      }
    }
  }
  */

  // ORIGINAL   //#####
  
  // update species

  for (field = NEQ_HYDRO; field < NEQ_HYDRO+NSpecies_renorm; field++) { 
    n = 0;
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
        for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++, igrid++) {
          Prim[field][igrid] = c1*OldPrim[field][igrid] +
            (1-c1)*Prim[field][igrid]*Prim[iden][igrid] + c2*dU[field][n];
          if (NoMultiSpeciesButColors != 1)
	    D[n] += Prim[field][igrid];
        }
      }
    }
  }

#ifdef USE_KROME
  if (use_kromeconserve) {
    FORTRAN_NAME(krome_renormref)(
      BaryonField[iden], 
      BaryonField[DeNum], BaryonField[CHINum], BaryonField[OINum],
      BaryonField[HNCINum], BaryonField[HCNINum],
      BaryonField[H2INum], BaryonField[CINum],
      BaryonField[HINum], BaryonField[H2OINum],
      BaryonField[OHINum], BaryonField[O2INum],
      BaryonField[CH2INum], BaryonField[H2COINum],
      BaryonField[HCOINum], BaryonField[MGINum],
      BaryonField[NH3INum], BaryonField[NOINum],
      BaryonField[CNINum], BaryonField[COINum],
      BaryonField[N2INum], BaryonField[NH2INum],
      BaryonField[CH3INum], BaryonField[CH4INum],
      BaryonField[NINum], BaryonField[NHINum],
      BaryonField[HeINum], BaryonField[HNOINum],
      BaryonField[CH3OHINum], BaryonField[CO2INum],
      BaryonField[H2CNINum], BaryonField[HNCOINum],
      BaryonField[NO2INum], BaryonField[O2HINum],
      BaryonField[OCNINum], BaryonField[CH3OH_DUSTINum],
      BaryonField[HNCO_DUSTINum], BaryonField[H2CO_DUSTINum],
      BaryonField[CH4_DUSTINum], BaryonField[CO_DUSTINum],
      BaryonField[H2O_DUSTINum], BaryonField[NO_DUSTINum],
      BaryonField[CO2_DUSTINum], BaryonField[N2_DUSTINum],
      BaryonField[HCN_DUSTINum], BaryonField[NH3_DUSTINum],
      BaryonField[O2_DUSTINum], BaryonField[NO2_DUSTINum],
      BaryonField[HNO_DUSTINum], BaryonField[O2H_DUSTINum],
      BaryonField[H2CN_DUSTINum], BaryonField[MG_DUSTINum],
      BaryonField[HNC_DUSTINum], BaryonField[E_DUSTINum],
      BaryonField[HCOIINum], BaryonField[HIINum],
      BaryonField[HOCIINum], BaryonField[CIINum],
      BaryonField[CH2IINum], BaryonField[CHIINum],
      BaryonField[H2COIINum], BaryonField[MGIINum],
      BaryonField[NH3IINum], BaryonField[NOIINum],
      BaryonField[CNIINum], BaryonField[COIINum],
      BaryonField[N2IINum], BaryonField[O2IINum],
      BaryonField[H2OIINum], BaryonField[NH2IINum],
      BaryonField[OIINum], BaryonField[OHIINum],
      BaryonField[CH3IINum], BaryonField[CH4IINum],
      BaryonField[NIINum], BaryonField[HCNIINum],
      BaryonField[NHIINum], BaryonField[H2IINum],
      BaryonField[HeIINum], BaryonField[HNOIINum],
      BaryonField[H2NOIINum], BaryonField[H3IINum],
      BaryonField[H3COIINum], BaryonField[H3OIINum],
      BaryonField[HCNHIINum], BaryonField[HCO2IINum],
      BaryonField[HeHIINum], BaryonField[N2HIINum],
      BaryonField[O2HIINum], 
      GridDimension, GridDimension+1, GridDimension+2, 
      &GridRank, GridStartIndex, GridStartIndex+1, GridStartIndex+2, 
      GridEndIndex, GridEndIndex+1, GridEndIndex+2,
      atomabund[0], atomabund[1], atomabund[2], atomabund[3],
      atomabund[4], atomabund[5], atomabund[6]);
  }
#endif

  // renormalize species
  if (NoMultiSpeciesButColors != 1) {
    for (field = NEQ_HYDRO; field < NEQ_HYDRO+NSpecies_renorm; field++) {
      n = 0;
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++, igrid++) {
	    Prim[field][igrid] = min(1.0, max((Prim[field][igrid]/D[n]), SmallX));
	    sum[n] += Prim[field][igrid];
	  }
	}
      }
    }
  

    for (field = NEQ_HYDRO; field < NEQ_HYDRO+NSpecies_renorm; field++) {
      n = 0;
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++, igrid++)
	    Prim[field][igrid] /= sum[n];
	}
      }
    }
  }


  // update conserved variables
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;

  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum);

  float rho_old, vx_old, vy_old, vz_old, e_old, etot_old, Tau_old, eint_old,
    rho, vx, vy, vz, e, etot, Tau, eint, p, v2,
    D_new, S1_new, S2_new, S3_new, Tau_new, h, cs, dpdrho, dpde, Eint_new;

  n = 0;
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	// first convert to conserved variables to do the update
	igrid = (k * GridDimension[1] + j) * GridDimension[0] + i;
	rho_old  = OldBaryonField[DensNum][igrid];
	vx_old   = OldBaryonField[Vel1Num][igrid];
	vy_old   = OldBaryonField[Vel2Num][igrid];
	vz_old   = OldBaryonField[Vel3Num][igrid];
	etot_old = OldBaryonField[TENum][igrid];
	if (DualEnergyFormalism) {
	  eint_old = OldBaryonField[GENum][igrid];
	}
	Tau_old = rho_old*etot_old;

	rho  = BaryonField[DensNum][igrid];
	vx   = BaryonField[Vel1Num][igrid];
	vy   = BaryonField[Vel2Num][igrid];
	vz   = BaryonField[Vel3Num][igrid];
	etot = BaryonField[TENum][igrid];
	if (DualEnergyFormalism) {
	  eint = BaryonField[GENum][igrid];
	}
	Tau = rho*etot;

	D_new   = c1*rho_old + (1.0-c1)*rho + c2*dU[iD][n];
	S1_new  = c1*rho_old*vx_old + (1.0-c1)*rho*vx + c2*dU[iS1][n];
	S2_new  = c1*rho_old*vy_old + (1.0-c1)*rho*vy + c2*dU[iS2][n];
	S3_new  = c1*rho_old*vz_old + (1.0-c1)*rho*vz + c2*dU[iS3][n];
	Tau_new = c1*Tau_old + (1.0-c1)*Tau + c2*dU[iEtot][n];


	
	if (DualEnergyFormalism) {
	  Eint_new = c1*rho_old*eint_old + (1.0-c1)*rho*eint + c2*dU[iEint][n];
	  /*if (Eint_new < 0) {
	    printf("UpdatePrim: eint < 0 in dual energy update\n");
	    printf("eint_old=%"GSYM",eint=%"GSYM",eint_new=%"GSYM", dU=%"GSYM"\n",
		   eint_old, eint, Eint_new, dU[iEint][n]);
	    return FAIL;
	    }*/
	}

	if (D_new < 0 || isnan(D_new)) {
	  PrintToScreenBoundaries(BaryonField[0], "Density", 1, j);

	  printf("UpdatePrim: rho <0 at %"ISYM" %"ISYM" %"ISYM": rho_old=%"GSYM", rho=%"GSYM", rho_new=%"GSYM", dU[iD]=%"GSYM"\n", 
		 i, j, k, rho_old, rho, D_new, dU[iD][n]);
	  //D_new = max(D_new, SmallRho);
	  D_new = rho;
	  //D_new = rho;
	  return FAIL;
	}

	//D_new = max(D_new, SmallRho);

	// convert back to primitives
	vx = S1_new/D_new;
	vy = S2_new/D_new;
	vz = S3_new/D_new;
	etot = Tau_new/D_new;



	v2 = vx*vx + vy*vy + vz*vz;
	// If using polytropic EOS, calcuate etot using density
	if (EOSType > 0) { 
	  EOS(p, D_new, eint, h, cs, dpdrho, dpde, EOSType, 0);
	  etot = eint + 0.5*v2;
	}


	if (etot < 0 && EOSType == 0) {
	  float v2_old = vx_old*vx_old + vy_old*vy_old + vz_old*vz_old;
	  printf("UpdatePrim: tau < 0. etot_old=%"GSYM", etot=%"GSYM", etot_new=%"GSYM", v2=%"GSYM", v2old=%"GSYM", dU[iTau] = %"GSYM"\n", 
		 Tau_old/rho_old, Tau/rho, Tau_new/D_new, v2, v2_old, dU[iEtot][n]*CellWidth[0][0]/dtFixed);
	  printf("rho_new=%"GSYM", rho=%"GSYM", rho_old=%"GSYM"\n", D_new, rho, rho_old);
	  printf("i=%"ISYM", j=%"ISYM", k=%"ISYM"\n", i,j,k);

	  //return FAIL;
	}
	
	BaryonField[DensNum][igrid] = D_new;
	BaryonField[Vel1Num][igrid] = vx;
	BaryonField[Vel2Num][igrid] = vy;
	BaryonField[Vel3Num][igrid] = vz;
	BaryonField[TENum][igrid] = etot;




	if (DualEnergyFormalism) {
	  v2 = vx*vx + vy*vy + vz*vz;
	  eint = Eint_new/D_new;
	  float eint_du = eint;
	  float eint1 = etot - 0.5*v2;
          float emin = SmallT/(Mu*(Gamma-1.0));
	  //float emin = 4.0*0.48999*D_new*pow(CellWidth[0][0],2)/(Gamma*(Gamma-1.0));

	  if (eint1 > 0) {
	    EOS(p, D_new, eint1, h, cs, dpdrho, dpde, EOSType, 2);
	  }
	  else {
	    cs = 0.0;
	  }
	  if (cs*cs > DualEnergyFormalismEta1*v2 && eint1 > 0.5*eint) {
	    eint = eint1;
	  }
	  eint = max(eint, emin);
	  BaryonField[GENum][igrid] = eint;
	  BaryonField[TENum][igrid] = eint + 0.5*v2;
	  
	  if (eint < 0.0) {
	    printf("UpdatePrim:eint < 0, cs2=%"GSYM", eta*v2=%"GSYM", eint=%"GSYM", etot=%"GSYM", v2=%"GSYM", p=%"GSYM", rho=%"GSYM",eint1=%"GSYM"\n", 
		   cs*cs, DualEnergyFormalismEta1*v2, eint_du, etot, v2, p, D_new, eint1);
	    printf("old rho=%"GSYM", old v:%"GSYM" %"GSYM" %"GSYM", old etot=%"GSYM" oldeint=%"GSYM" \n", 
		   rho_old, vx_old, vy_old, vz_old,
		   etot_old, eint_old);
	    printf("dU=%"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM", Tau_old=%"GSYM"\n", 
		   dU[iD][n], dU[iS1][n], dU[iS2][n], dU[iS3][n], dU[iEtot][n], Tau_old);
	    //return FAIL;
	    BaryonField[GENum][igrid] = OldBaryonField[GENum][igrid];
	    BaryonField[TENum][igrid] = OldBaryonField[TENum][igrid];
	  }
	}
      }
    }
  }

  // convert species from mass fraction to density (this reverts what Grid_ReturnHydroRKPointers did in Grid_RungeKutta_[12]Step)
  if (NoMultiSpeciesButColors != 1)
    for (field = NEQ_HYDRO; field < NEQ_HYDRO+NSpecies+NColor; field++)   
      for (n = 0; n < size; n++) 
	Prim[field][n] *= Prim[iden][n];
  
#ifdef USE_KROME
  for (i=0; i<NKROMEATOMS; i++) {
    delete [] atomabund[i];
  }
#else 
  this->UpdateElectronDensity();
#endif
 
  if ( (NSpecies+NColor) > 0) {
    delete [] D;
    delete [] sum;
  }
  
  return SUCCESS;
}


