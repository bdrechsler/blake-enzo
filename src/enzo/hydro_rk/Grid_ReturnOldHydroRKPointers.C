/***********************************************************************
/
/  GRID CLASS (RETURNS AN ARRAY OF POINTERS THAT ARE COMPATIBLE WITH 
/              THE HYDRO_RK SOLVERS)
/
/  written by: John Wise
/  date:       July, 2009
/  modified1:
/
/
************************************************************************/

#include <stdio.h>
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

int FindField(int field, int farray[], int numfields);

int grid::ReturnOldHydroRKPointers(float **Prim, bool ReturnMassFractions)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  int i, n, dim, size, nfield, n0;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
#ifdef USE_KROME
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
#endif

  /* Add the physical quantities */

  if (HydroMethod == HD_RK) {
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				     Vel3Num, TENum);
    nfield = n0 = NEQ_HYDRO;
  }

  else if (HydroMethod == MHD_RK) {
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				     Vel3Num, TENum, B1Num, B2Num, B3Num, 
				     PhiNum);
    nfield = n0 = NEQ_MHD;
  }
  
  Prim[iden] = OldBaryonField[DensNum];
  Prim[ivx] = OldBaryonField[Vel1Num];
  if (GridRank > 1)
    Prim[ivy] = OldBaryonField[Vel2Num];
  if (GridRank > 2)
    Prim[ivz] = OldBaryonField[Vel3Num];
  Prim[ietot] = OldBaryonField[TENum];
  if (DualEnergyFormalism)
    Prim[ieint] = OldBaryonField[GENum];

  if (HydroMethod == MHD_RK) {
    Prim[iBx] = OldBaryonField[B1Num];
    Prim[iBy] = OldBaryonField[B2Num];
    Prim[iBz] = OldBaryonField[B3Num];
    Prim[iPhi] = OldBaryonField[PhiNum];
  }

  /* Add the species */

  if (MultiSpecies) {

    this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
				HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);

#ifdef USE_KROME
    if (MultiSpecies == KROMESPECIES)
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
#endif


    //Prim[nfield++] = OldBaryonField[DeNum];
    Prim[nfield++] = OldBaryonField[HINum];
    Prim[nfield++] = OldBaryonField[HIINum];
    Prim[nfield++] = OldBaryonField[HeINum];
    Prim[nfield++] = OldBaryonField[HeIINum];
    Prim[nfield++] = OldBaryonField[HeIIINum];

    if (MultiSpecies > 1) {
      Prim[nfield++] = OldBaryonField[HMNum];
      Prim[nfield++] = OldBaryonField[H2INum];
      Prim[nfield++] = OldBaryonField[H2IINum];
    }

    if (MultiSpecies > 2) {
      Prim[nfield++] = OldBaryonField[DINum];
      Prim[nfield++] = OldBaryonField[DIINum];
      Prim[nfield++] = OldBaryonField[HDINum];
    }

#ifdef USE_KROME
    if (MultiSpecies == KROMESPECIES) {
      Prim[nfield++] =  OldBaryonField[CHINum];
      Prim[nfield++] =  OldBaryonField[OINum];
      Prim[nfield++] =  OldBaryonField[HNCINum];
      Prim[nfield++] =  OldBaryonField[HCNINum];
      Prim[nfield++] =  OldBaryonField[CINum];
      Prim[nfield++] =  OldBaryonField[H2OINum];
      Prim[nfield++] =  OldBaryonField[OHINum];
      Prim[nfield++] =  OldBaryonField[O2INum];
      Prim[nfield++] =  OldBaryonField[CH2INum];
      Prim[nfield++] =  OldBaryonField[H2COINum];
      Prim[nfield++] =  OldBaryonField[HCOINum];
      Prim[nfield++] =  OldBaryonField[MGINum];
      Prim[nfield++] =  OldBaryonField[NH3INum];
      Prim[nfield++] =  OldBaryonField[NOINum];
      Prim[nfield++] =  OldBaryonField[CNINum];
      Prim[nfield++] =  OldBaryonField[COINum];
      Prim[nfield++] =  OldBaryonField[N2INum];
      Prim[nfield++] =  OldBaryonField[NH2INum];
      Prim[nfield++] =  OldBaryonField[CH3INum];
      Prim[nfield++] =  OldBaryonField[CH4INum];
      Prim[nfield++] =  OldBaryonField[NINum];
      Prim[nfield++] =  OldBaryonField[NHINum];
      Prim[nfield++] =  OldBaryonField[HNOINum];
      Prim[nfield++] =  OldBaryonField[CH3OHINum];
      Prim[nfield++] =  OldBaryonField[CO2INum];
      Prim[nfield++] =  OldBaryonField[H2CNINum];
      Prim[nfield++] =  OldBaryonField[HNCOINum];
      Prim[nfield++] =  OldBaryonField[NO2INum];
      Prim[nfield++] =  OldBaryonField[O2HINum];
      Prim[nfield++] =  OldBaryonField[OCNINum];
      Prim[nfield++] =  OldBaryonField[CH3OH_DUSTINum];
      Prim[nfield++] =  OldBaryonField[HNCO_DUSTINum];
      Prim[nfield++] =  OldBaryonField[H2CO_DUSTINum];
      Prim[nfield++] =  OldBaryonField[CH4_DUSTINum];
      Prim[nfield++] =  OldBaryonField[CO_DUSTINum];
      Prim[nfield++] =  OldBaryonField[H2O_DUSTINum];
      Prim[nfield++] =  OldBaryonField[NO_DUSTINum];
      Prim[nfield++] =  OldBaryonField[CO2_DUSTINum];
      Prim[nfield++] =  OldBaryonField[N2_DUSTINum];
      Prim[nfield++] =  OldBaryonField[HCN_DUSTINum];
      Prim[nfield++] =  OldBaryonField[NH3_DUSTINum];
      Prim[nfield++] =  OldBaryonField[O2_DUSTINum];
      Prim[nfield++] =  OldBaryonField[NO2_DUSTINum];
      Prim[nfield++] =  OldBaryonField[HNO_DUSTINum];
      Prim[nfield++] =  OldBaryonField[O2H_DUSTINum];
      Prim[nfield++] =  OldBaryonField[H2CN_DUSTINum];
      Prim[nfield++] =  OldBaryonField[MG_DUSTINum];
      Prim[nfield++] =  OldBaryonField[HNC_DUSTINum];
      Prim[nfield++] =  OldBaryonField[E_DUSTINum];
      Prim[nfield++] =  OldBaryonField[HCOIINum];
      Prim[nfield++] =  OldBaryonField[HOCIINum];
      Prim[nfield++] =  OldBaryonField[CIINum];
      Prim[nfield++] =  OldBaryonField[CH2IINum];
      Prim[nfield++] =  OldBaryonField[CHIINum];
      Prim[nfield++] =  OldBaryonField[H2COIINum];
      Prim[nfield++] =  OldBaryonField[MGIINum];
      Prim[nfield++] =  OldBaryonField[NH3IINum];
      Prim[nfield++] =  OldBaryonField[NOIINum];
      Prim[nfield++] =  OldBaryonField[CNIINum];
      Prim[nfield++] =  OldBaryonField[COIINum];
      Prim[nfield++] =  OldBaryonField[N2IINum];
      Prim[nfield++] =  OldBaryonField[O2IINum];
      Prim[nfield++] =  OldBaryonField[H2OIINum];
      Prim[nfield++] =  OldBaryonField[NH2IINum];
      Prim[nfield++] =  OldBaryonField[OIINum];
      Prim[nfield++] =  OldBaryonField[OHIINum];
      Prim[nfield++] =  OldBaryonField[CH3IINum];
      Prim[nfield++] =  OldBaryonField[CH4IINum];
      Prim[nfield++] =  OldBaryonField[NIINum];
      Prim[nfield++] =  OldBaryonField[HCNIINum];
      Prim[nfield++] =  OldBaryonField[NHIINum];
      Prim[nfield++] =  OldBaryonField[HNOIINum];
      Prim[nfield++] =  OldBaryonField[H2NOIINum];
      Prim[nfield++] =  OldBaryonField[H3IINum];
      Prim[nfield++] =  OldBaryonField[H3COIINum];
      Prim[nfield++] =  OldBaryonField[H3OIINum];
      Prim[nfield++] =  OldBaryonField[HCNHIINum];
      Prim[nfield++] =  OldBaryonField[HCO2IINum];
      Prim[nfield++] =  OldBaryonField[HeHIINum];
      Prim[nfield++] =  OldBaryonField[N2HIINum];
      Prim[nfield++] =  OldBaryonField[O2HIINum];

    }
#endif

  } // ENDIF MultiSpecies

  /* Add the colours (treat them as species) */

  int SNColourNum, MetalNum, MetalIaNum, MetalIINum, MBHColourNum, Galaxy1ColourNum, 
    Galaxy2ColourNum; 

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum, MBHColourNum, 
				 Galaxy1ColourNum, Galaxy2ColourNum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyColourFields.\n");
    return FAIL;
  }
  
  if (MetalNum != -1) {
    Prim[nfield++] = OldBaryonField[MetalNum];
    if (StarMakerTypeIaSNe)
      Prim[nfield++] = OldBaryonField[MetalIaNum];
    if (StarMakerTypeIISNeMetalField)
      Prim[nfield++] = OldBaryonField[MetalIINum];
    if (MultiMetals || TestProblemData.MultiMetals) {
      Prim[nfield++] = OldBaryonField[MetalNum+1];
      Prim[nfield++] = OldBaryonField[MetalNum+2];
    }
  }

  if (SNColourNum      != -1) Prim[nfield++] = OldBaryonField[SNColourNum];  
  /*   //##### These fields are currently not being used and only causing interpolation problems
  if (MBHColourNum     != -1) Prim[nfield++] = OldBaryonField[MBHColourNum];
  if (Galaxy1ColourNum != -1) Prim[nfield++] = OldBaryonField[Galaxy1ColourNum];
  if (Galaxy2ColourNum != -1) Prim[nfield++] = OldBaryonField[Galaxy2ColourNum];
  */

  /* Convert the species and color fields into mass fractions */

  for (dim = 0, size = 1; dim < GridRank; dim++)
    size *= GridDimension[dim];

  if (ReturnMassFractions)  
    for (n = n0; n < nfield; n++)
      for (i = 0; i < size; i++)
	Prim[n][i] /= Prim[iden][i];  

  return SUCCESS;

}
