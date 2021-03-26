/***********************************************************************
/
/  GRID CLASS (UPDATE MHD VARIABLES)
/
/  written by: Peng Wang
/  date:       June, 2007
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

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

extern "C" void FORTRAN_NAME(krome_initref)(float *d, float *urho,
                                            float *De, float *CHI, float *OI, float *HNCI,
                                            float *HCNI, float *H2I, float *CI, float *HI,
                                            float *H2OI, float *OHI, float *O2I, float *CH2I,
                                            float *H2COI, float *HCOI, float *MGI, float *NH3I,
                                            float *NOI, float *SII, float *SIC2I, float *SIC3I,
                                            float *SICI, float *SIH2I, float *SIH3I,
                                            float *CNI, float *COI, float *N2I, float *NH2I,
                                            float *CH3I, float *CH4I, float *NI, float *NHI,
                                            float *SIH4I, float *SIHI, float *SIOI, float *HeI,
                                            float *HNOI, float *CH3OHI, float *CO2I,
                                            float *H2CNI, float *H2SIOI, float *HNCOI,
                                            float *NO2I, float *O2HI, float *OCNI, float *CH3OH_DUSTI,
                                            float *HNCO_DUSTI, float *H2CO_DUSTI, float *SIH4_DUSTI,
                                            float *H2SIO_DUSTI, float *SIC_DUSTI, float *SIC2_DUSTI,
                                            float *SIC3_DUSTI, float *CH4_DUSTI, float *CO_DUSTI,
                                            float *H2O_DUSTI, float *NO_DUSTI, float *CO2_DUSTI,
                                            float *N2_DUSTI, float *HCN_DUSTI, float *NH3_DUSTI,
                                            float *O2_DUSTI, float *NO2_DUSTI, float *HNO_DUSTI,
                                            float *O2H_DUSTI, float *H2CN_DUSTI, float *MG_DUSTI,
                                            float *HNC_DUSTI, float *E_DUSTI, float *SIO_DUSTI,
                                            float *HCOII, float *HII, float *HOCII, float *CII,
                                            float *CH2II, float *CHII, float *H2COII,
                                            float *MGII, float *NH3II, float *NOII, float *SIII,
                                            float *SIC2II, float *SIC3II, float *SICII,
                                            float *SIH2II, float *SIH3II, float *CNII,
                                            float *COII, float *N2II, float *O2II, float *H2OII,
                                            float *NH2II, float *OII, float *OHII, float *CH3II,
                                            float *CH4II, float *NII, float *HCNII, float *NHII,
                                            float *SIH4II, float *SIHII, float *SIOII,
                                            float *H2II, float *HeII, float *HNOII, float *H2NOII,
                                            float *H3II, float *H3COII, float *H3OII,
                                            float *HCNHII, float *HCO2II, float *HeHII,
                                            float *N2HII, float *O2HII, float *SIH5II,
                                            float *SIOHII, int *in, int *jn, int *kn, int *idim,
                                            int *is, int *js, int *ks, int *ie, int *je, int *ke,
                                            float *ATOM_C, float *ATOM_H, float *ATOM_MG,
                                            float *ATOM_O, float *ATOM_N, float *ATOM_SI, float *ATOM_HE);

extern "C" void FORTRAN_NAME(krome_renormref)(float *d, float *urho,
                                              float *De, float *CHI, float *OI, float *HNCI,
                                              float *HCNI, float *H2I, float *CI, float *HI,
                                              float *H2OI, float *OHI, float *O2I, float *CH2I,
                                              float *H2COI, float *HCOI, float *MGI, float *NH3I,
                                              float *NOI, float *SII, float *SIC2I, float *SIC3I,
                                              float *SICI, float *SIH2I, float *SIH3I,
                                              float *CNI, float *COI, float *N2I, float *NH2I,
                                              float *CH3I, float *CH4I, float *NI, float *NHI,
                                              float *SIH4I, float *SIHI, float *SIOI, float *HeI,
                                              float *HNOI, float *CH3OHI, float *CO2I,
                                              float *H2CNI, float *H2SIOI, float *HNCOI,
                                              float *NO2I, float *O2HI, float *OCNI, float *CH3OH_DUSTI,
                                              float *HNCO_DUSTI, float *H2CO_DUSTI, float *SIH4_DUSTI,
                                              float *H2SIO_DUSTI, float *SIC_DUSTI, float *SIC2_DUSTI,
                                              float *SIC3_DUSTI, float *CH4_DUSTI, float *CO_DUSTI,
                                              float *H2O_DUSTI, float *NO_DUSTI, float *CO2_DUSTI,
                                              float *N2_DUSTI, float *HCN_DUSTI, float *NH3_DUSTI,
                                              float *O2_DUSTI, float *NO2_DUSTI, float *HNO_DUSTI,
                                              float *O2H_DUSTI, float *H2CN_DUSTI, float *MG_DUSTI,
                                              float *HNC_DUSTI, float *E_DUSTI, float *SIO_DUSTI,
                                              float *HCOII, float *HII, float *HOCII, float *CII,
                                              float *CH2II, float *CHII, float *H2COII,
                                              float *MGII, float *NH3II, float *NOII, float *SIII,
                                              float *SIC2II, float *SIC3II, float *SICII,
                                              float *SIH2II, float *SIH3II, float *CNII,
                                              float *COII, float *N2II, float *O2II, float *H2OII,
                                              float *NH2II, float *OII, float *OHII, float *CH3II,
                                              float *CH4II, float *NII, float *HCNII, float *NHII,
                                              float *SIH4II, float *SIHII, float *SIOII,
                                              float *H2II, float *HeII, float *HNOII, float *H2NOII,
                                              float *H3II, float *H3COII, float *H3OII,
                                              float *HCNHII, float *HCO2II, float *HeHII,
                                              float *N2HII, float *O2HII, float *SIH5II,
                                              float *SIOHII, int *in, int *jn, int *kn, int *idim,
                                              int *is, int *js, int *ks, int *ie, int *je, int *ke,
                                              float *ATOM_C, float *ATOM_H, float *ATOM_MG,
                                              float *ATOM_O, float *ATOM_N, float *ATOM_SI, float *ATOM_HE);

int grid::UpdateMHDPrim(float **dU, float c1, float c2)
{

  if (ProcessorNumber != MyProcessorNumber)
  {
    return SUCCESS;
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num,
      B1Num, B2Num, B3Num, PhiNum;

  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                   Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum);

  int i, j, k, n, n_dU, dim, igrid, field, size, activesize;
  for (dim = 0, size = 1; dim < GridRank; dim++)
  {
    size *= GridDimension[dim];
  }

  for (dim = 0, activesize = 1; dim < GridRank; dim++)
  {
    activesize *= (GridDimension[dim] - 2 * NumberOfGhostZones);
  }

  float *D, *sum;
  float SmallX = 1e-20;

  // ORIGINAL
  if ((NSpecies + NColor) > 0)
  {
    D = new float[activesize];
    sum = new float[activesize];
    for (i = 0; i < activesize; i++)
    {
      D[i] = 0.0;
      sum[i] = 0.0;
    }
  }

  float *Prim[NEQ_MHD + NSpecies + NColor];
  float *OldPrim[NEQ_MHD + NSpecies + NColor];
  this->ReturnHydroRKPointers(Prim, false);
  this->ReturnOldHydroRKPointers(OldPrim, false);

  //##### Want to mix species and colors for renormalization?  Normally you don't
  int NSpecies_renorm;
  if (MixSpeciesAndColors)
    NSpecies_renorm = NSpecies + NColor;
  else if (NoMultiSpeciesButColors)
  {
    NSpecies_renorm = NSpecies;
  }
  else
    switch (MultiSpecies)
    { //update pure species! not colours!
    case 0:
      NSpecies_renorm = 0;
      break;
    case 1:
      NSpecies_renorm = 5;
      break;
    case 2:
      NSpecies_renorm = 8;
      break;
    case 3:
      NSpecies_renorm = 11;
      break;
#ifdef USE_KROME
    case KROMESPECIES:
      NSpecies_renorm = NKROMESPECIES;
      break;
#endif
    default:
      NSpecies_renorm = 0;
      break;
    }

#ifdef USE_KROME
  float *atomabund[NKROMEATOMS];
  for (i = 0; i < NKROMEATOMS; i++)
  {
    atomabund[i] = new float[size];
  }

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  int CHINum, OINum, HNCINum, HCNINum, CINum, H2OINum,
      OHINum, O2INum, CH2INum, H2COINum, HCOINum,
      MGINum, NH3INum, NOINum, SIINum, SIC2INum,
      SIC3INum, SICINum, SIH2INum, SIH3INum, CNINum,
      COINum, N2INum, NH2INum, CH3INum, CH4INum,
      NINum, NHINum, SIH4INum, SIHINum, SIOINum,
      HNOINum, CH3OHINum, CO2INum, H2CNINum, H2SIOINum,
      HNCOINum, NO2INum, O2HINum, OCNINum, CH3OH_DUSTINum,
      HNCO_DUSTINum, H2CO_DUSTINum, SIH4_DUSTINum,
      H2SIO_DUSTINum, SIC_DUSTINum, SIC2_DUSTINum,
      SIC3_DUSTINum, CH4_DUSTINum, CO_DUSTINum,
      H2O_DUSTINum, NO_DUSTINum, CO2_DUSTINum,
      N2_DUSTINum, HCN_DUSTINum, NH3_DUSTINum,
      O2_DUSTINum, NO2_DUSTINum, HNO_DUSTINum,
      O2H_DUSTINum, H2CN_DUSTINum, MG_DUSTINum,
      HNC_DUSTINum, E_DUSTINum, SIO_DUSTINum, HCOIINum,
      HOCIINum, CIINum, CH2IINum, CHIINum, H2COIINum,
      MGIINum, NH3IINum, NOIINum, SIIINum, SIC2IINum,
      SIC3IINum, SICIINum, SIH2IINum, SIH3IINum,
      CNIINum, COIINum, N2IINum, O2IINum, H2OIINum,
      NH2IINum, OIINum, OHIINum, CH3IINum, CH4IINum,
      NIINum, HCNIINum, NHIINum, SIH4IINum, SIHIINum,
      SIOIINum, HNOIINum, H2NOIINum, H3IINum, H3COIINum,
      H3OIINum, HCNHIINum, HCO2IINum, HeHIINum,
      N2HIINum, O2HIINum, SIH5IINum, SIOHIINum;

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
        VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, Time) == FAIL)
  {
    ENZO_FAIL("Error in GetUnits.");
  }
  if (use_kromeconserve)
  {
    if (IdentifySpeciesFieldsKrome(
            DeNum, CHINum, OINum, HNCINum, HCNINum, H2INum,
            CINum, HINum, H2OINum, OHINum, O2INum, CH2INum,
            H2COINum, HCOINum, MGINum, NH3INum, NOINum,
            SIINum, SIC2INum, SIC3INum, SICINum, SIH2INum,
            SIH3INum, CNINum, COINum, N2INum, NH2INum,
            CH3INum, CH4INum, NINum, NHINum, SIH4INum,
            SIHINum, SIOINum, HeINum, HNOINum, CH3OHINum,
            CO2INum, H2CNINum, H2SIOINum, HNCOINum, NO2INum,
            O2HINum, OCNINum, CH3OH_DUSTINum, HNCO_DUSTINum,
            H2CO_DUSTINum, SIH4_DUSTINum, H2SIO_DUSTINum,
            SIC_DUSTINum, SIC2_DUSTINum, SIC3_DUSTINum,
            CH4_DUSTINum, CO_DUSTINum, H2O_DUSTINum,
            NO_DUSTINum, CO2_DUSTINum, N2_DUSTINum, HCN_DUSTINum,
            NH3_DUSTINum, O2_DUSTINum, NO2_DUSTINum,
            HNO_DUSTINum, O2H_DUSTINum, H2CN_DUSTINum,
            MG_DUSTINum, HNC_DUSTINum, E_DUSTINum, SIO_DUSTINum,
            HCOIINum, HIINum, HOCIINum, CIINum, CH2IINum,
            CHIINum, H2COIINum, MGIINum, NH3IINum, NOIINum,
            SIIINum, SIC2IINum, SIC3IINum, SICIINum,
            SIH2IINum, SIH3IINum, CNIINum, COIINum, N2IINum,
            O2IINum, H2OIINum, NH2IINum, OIINum, OHIINum,
            CH3IINum, CH4IINum, NIINum, HCNIINum, NHIINum,
            SIH4IINum, SIHIINum, SIOIINum, H2IINum, HeIINum,
            HNOIINum, H2NOIINum, H3IINum, H3COIINum,
            H3OIINum, HCNHIINum, HCO2IINum, HeHIINum,
            N2HIINum, O2HIINum, SIH5IINum, SIOHIINum) == FAIL)
    {
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }

    FORTRAN_NAME(krome_initref)
    (
        BaryonField[iden], &DensityUnits,
        BaryonField[DeNum], BaryonField[CHINum], BaryonField[OINum],
        BaryonField[HNCINum], BaryonField[HCNINum],
        BaryonField[H2INum], BaryonField[CINum],
        BaryonField[HINum], BaryonField[H2OINum],
        BaryonField[OHINum], BaryonField[O2INum],
        BaryonField[CH2INum], BaryonField[H2COINum],
        BaryonField[HCOINum], BaryonField[MGINum],
        BaryonField[NH3INum], BaryonField[NOINum],
        BaryonField[SIINum], BaryonField[SIC2INum],
        BaryonField[SIC3INum], BaryonField[SICINum],
        BaryonField[SIH2INum], BaryonField[SIH3INum],
        BaryonField[CNINum], BaryonField[COINum],
        BaryonField[N2INum], BaryonField[NH2INum],
        BaryonField[CH3INum], BaryonField[CH4INum],
        BaryonField[NINum], BaryonField[NHINum],
        BaryonField[SIH4INum], BaryonField[SIHINum],
        BaryonField[SIOINum], BaryonField[HeINum],
        BaryonField[HNOINum], BaryonField[CH3OHINum],
        BaryonField[CO2INum], BaryonField[H2CNINum],
        BaryonField[H2SIOINum], BaryonField[HNCOINum],
        BaryonField[NO2INum], BaryonField[O2HINum],
        BaryonField[OCNINum], BaryonField[CH3OH_DUSTINum],
        BaryonField[HNCO_DUSTINum], BaryonField[H2CO_DUSTINum],
        BaryonField[SIH4_DUSTINum], BaryonField[H2SIO_DUSTINum],
        BaryonField[SIC_DUSTINum], BaryonField[SIC2_DUSTINum],
        BaryonField[SIC3_DUSTINum], BaryonField[CH4_DUSTINum],
        BaryonField[CO_DUSTINum], BaryonField[H2O_DUSTINum],
        BaryonField[NO_DUSTINum], BaryonField[CO2_DUSTINum],
        BaryonField[N2_DUSTINum], BaryonField[HCN_DUSTINum],
        BaryonField[NH3_DUSTINum], BaryonField[O2_DUSTINum],
        BaryonField[NO2_DUSTINum], BaryonField[HNO_DUSTINum],
        BaryonField[O2H_DUSTINum], BaryonField[H2CN_DUSTINum],
        BaryonField[MG_DUSTINum], BaryonField[HNC_DUSTINum],
        BaryonField[E_DUSTINum], BaryonField[SIO_DUSTINum],
        BaryonField[HCOIINum], BaryonField[HIINum],
        BaryonField[HOCIINum], BaryonField[CIINum],
        BaryonField[CH2IINum], BaryonField[CHIINum],
        BaryonField[H2COIINum], BaryonField[MGIINum],
        BaryonField[NH3IINum], BaryonField[NOIINum],
        BaryonField[SIIINum], BaryonField[SIC2IINum],
        BaryonField[SIC3IINum], BaryonField[SICIINum],
        BaryonField[SIH2IINum], BaryonField[SIH3IINum],
        BaryonField[CNIINum], BaryonField[COIINum],
        BaryonField[N2IINum], BaryonField[O2IINum],
        BaryonField[H2OIINum], BaryonField[NH2IINum],
        BaryonField[OIINum], BaryonField[OHIINum],
        BaryonField[CH3IINum], BaryonField[CH4IINum],
        BaryonField[NIINum], BaryonField[HCNIINum],
        BaryonField[NHIINum], BaryonField[SIH4IINum],
        BaryonField[SIHIINum], BaryonField[SIOIINum],
        BaryonField[H2IINum], BaryonField[HeIINum],
        BaryonField[HNOIINum], BaryonField[H2NOIINum],
        BaryonField[H3IINum], BaryonField[H3COIINum],
        BaryonField[H3OIINum], BaryonField[HCNHIINum],
        BaryonField[HCO2IINum], BaryonField[HeHIINum],
        BaryonField[N2HIINum], BaryonField[O2HIINum],
        BaryonField[SIH5IINum], BaryonField[SIOHIINum],
        GridDimension, GridDimension + 1, GridDimension + 2,
        &GridRank, GridStartIndex, GridStartIndex + 1, GridStartIndex + 2,
        GridEndIndex, GridEndIndex + 1, GridEndIndex + 2,
        atomabund[0], atomabund[1], atomabund[2], atomabund[3],
        atomabund[4], atomabund[5], atomabund[6]);
  }
#endif

  // update species
  for (field = NEQ_MHD; field < NEQ_MHD + NSpecies_renorm; field++)
  {
    n = 0;
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    {
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
      {
        igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
        for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++, igrid++)
        {
          Prim[field][igrid] = c1 * OldPrim[field][igrid] +
                               (1 - c1) * Prim[field][igrid] * Prim[iden][igrid] + c2 * dU[field][n];
          if (NoMultiSpeciesButColors != 1)
            D[n] += Prim[field][igrid];
        }
      }
    }
  }

#ifdef USE_KROME
  if (use_kromeconserve)
  {
    FORTRAN_NAME(krome_renormref)
    (
        BaryonField[iden], &DensityUnits,
        BaryonField[DeNum], BaryonField[CHINum], BaryonField[OINum],
        BaryonField[HNCINum], BaryonField[HCNINum],
        BaryonField[H2INum], BaryonField[CINum],
        BaryonField[HINum], BaryonField[H2OINum],
        BaryonField[OHINum], BaryonField[O2INum],
        BaryonField[CH2INum], BaryonField[H2COINum],
        BaryonField[HCOINum], BaryonField[MGINum],
        BaryonField[NH3INum], BaryonField[NOINum],
        BaryonField[SIINum], BaryonField[SIC2INum],
        BaryonField[SIC3INum], BaryonField[SICINum],
        BaryonField[SIH2INum], BaryonField[SIH3INum],
        BaryonField[CNINum], BaryonField[COINum],
        BaryonField[N2INum], BaryonField[NH2INum],
        BaryonField[CH3INum], BaryonField[CH4INum],
        BaryonField[NINum], BaryonField[NHINum],
        BaryonField[SIH4INum], BaryonField[SIHINum],
        BaryonField[SIOINum], BaryonField[HeINum],
        BaryonField[HNOINum], BaryonField[CH3OHINum],
        BaryonField[CO2INum], BaryonField[H2CNINum],
        BaryonField[H2SIOINum], BaryonField[HNCOINum],
        BaryonField[NO2INum], BaryonField[O2HINum],
        BaryonField[OCNINum], BaryonField[CH3OH_DUSTINum],
        BaryonField[HNCO_DUSTINum], BaryonField[H2CO_DUSTINum],
        BaryonField[SIH4_DUSTINum], BaryonField[H2SIO_DUSTINum],
        BaryonField[SIC_DUSTINum], BaryonField[SIC2_DUSTINum],
        BaryonField[SIC3_DUSTINum], BaryonField[CH4_DUSTINum],
        BaryonField[CO_DUSTINum], BaryonField[H2O_DUSTINum],
        BaryonField[NO_DUSTINum], BaryonField[CO2_DUSTINum],
        BaryonField[N2_DUSTINum], BaryonField[HCN_DUSTINum],
        BaryonField[NH3_DUSTINum], BaryonField[O2_DUSTINum],
        BaryonField[NO2_DUSTINum], BaryonField[HNO_DUSTINum],
        BaryonField[O2H_DUSTINum], BaryonField[H2CN_DUSTINum],
        BaryonField[MG_DUSTINum], BaryonField[HNC_DUSTINum],
        BaryonField[E_DUSTINum], BaryonField[SIO_DUSTINum],
        BaryonField[HCOIINum], BaryonField[HIINum],
        BaryonField[HOCIINum], BaryonField[CIINum],
        BaryonField[CH2IINum], BaryonField[CHIINum],
        BaryonField[H2COIINum], BaryonField[MGIINum],
        BaryonField[NH3IINum], BaryonField[NOIINum],
        BaryonField[SIIINum], BaryonField[SIC2IINum],
        BaryonField[SIC3IINum], BaryonField[SICIINum],
        BaryonField[SIH2IINum], BaryonField[SIH3IINum],
        BaryonField[CNIINum], BaryonField[COIINum],
        BaryonField[N2IINum], BaryonField[O2IINum],
        BaryonField[H2OIINum], BaryonField[NH2IINum],
        BaryonField[OIINum], BaryonField[OHIINum],
        BaryonField[CH3IINum], BaryonField[CH4IINum],
        BaryonField[NIINum], BaryonField[HCNIINum],
        BaryonField[NHIINum], BaryonField[SIH4IINum],
        BaryonField[SIHIINum], BaryonField[SIOIINum],
        BaryonField[H2IINum], BaryonField[HeIINum],
        BaryonField[HNOIINum], BaryonField[H2NOIINum],
        BaryonField[H3IINum], BaryonField[H3COIINum],
        BaryonField[H3OIINum], BaryonField[HCNHIINum],
        BaryonField[HCO2IINum], BaryonField[HeHIINum],
        BaryonField[N2HIINum], BaryonField[O2HIINum],
        BaryonField[SIH5IINum], BaryonField[SIOHIINum],
        GridDimension, GridDimension + 1, GridDimension + 2,
        &GridRank, GridStartIndex, GridStartIndex + 1, GridStartIndex + 2,
        GridEndIndex, GridEndIndex + 1, GridEndIndex + 2,
        atomabund[0], atomabund[1], atomabund[2], atomabund[3],
        atomabund[4], atomabund[5], atomabund[6]);
  }
#endif

  // renormalize species
  if (NoMultiSpeciesButColors != 1)
  {
    for (field = NEQ_MHD; field < NEQ_MHD + NSpecies_renorm; field++)
    {
      n = 0;
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      {
        for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
        {
          igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
          for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++, igrid++)
          {
            Prim[field][igrid] = min(1.0, max((Prim[field][igrid] / D[n]), SmallX));
            sum[n] += Prim[field][igrid];
          }
        }
      }
    }

    for (field = NEQ_MHD; field < NEQ_MHD + NSpecies_renorm; field++)
    {
      n = 0;
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      {
        for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
        {
          igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
          for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++, igrid++)
            Prim[field][igrid] /= sum[n];
        }
      }
    }
  } //close if (NoMultiSpeciesButColors == 1)
  /* Update conserved variables */

  float rho_old, vx_old, vy_old, vz_old, e_old, etot_old, Tau_old, eint_old,
      rho, vx, vy, vz, e, etot, Tau, eint, p, v2,
      D_new, S1_new, S2_new, S3_new, Tau_new, h, cs, dpdrho, dpde, Eint_new,
      Bx_old, By_old, Bz_old, Bx, By, Bz, Bx_new, By_new, Bz_new,
      Phi_old, Phi, Phi_new, B2;

  float rhou, lenu, tempu, tu, velu;
  GetUnits(&rhou, &lenu, &tempu, &tu, &velu, Time);

  n = 0;
  FLOAT x, y, z, r;

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
  {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
    {
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++)
      {
        // first convert to conserved variables to do the update
        igrid = (k * GridDimension[1] + j) * GridDimension[0] + i;
        r = sqrt(x * x + y * y + z * z);

        rho_old = OldBaryonField[DensNum][igrid];
        vx_old = OldBaryonField[Vel1Num][igrid];
        vy_old = OldBaryonField[Vel2Num][igrid];
        vz_old = OldBaryonField[Vel3Num][igrid];
        etot_old = OldBaryonField[TENum][igrid];
        Tau_old = rho_old * etot_old;
        if (DualEnergyFormalism)
        {
          eint_old = OldBaryonField[GENum][igrid];
        }
        Bx_old = OldBaryonField[B1Num][igrid];
        By_old = OldBaryonField[B2Num][igrid];
        Bz_old = OldBaryonField[B3Num][igrid];
        Phi_old = OldBaryonField[PhiNum][igrid];

        rho = BaryonField[DensNum][igrid];
        vx = BaryonField[Vel1Num][igrid];
        vy = BaryonField[Vel2Num][igrid];
        vz = BaryonField[Vel3Num][igrid];
        etot = BaryonField[TENum][igrid];
        Tau = rho * etot;
        if (DualEnergyFormalism)
        {
          eint = BaryonField[GENum][igrid];
        }
        Bx = BaryonField[B1Num][igrid];
        By = BaryonField[B2Num][igrid];
        Bz = BaryonField[B3Num][igrid];
        Phi = BaryonField[PhiNum][igrid];

        D_new = c1 * rho_old + (1.0 - c1) * rho + c2 * dU[iD][n];
        S1_new = c1 * rho_old * vx_old + (1.0 - c1) * rho * vx + c2 * dU[iS1][n];
        S2_new = c1 * rho_old * vy_old + (1.0 - c1) * rho * vy + c2 * dU[iS2][n];
        S3_new = c1 * rho_old * vz_old + (1.0 - c1) * rho * vz + c2 * dU[iS3][n];
        Tau_new = c1 * Tau_old + (1.0 - c1) * Tau + c2 * dU[iEtot][n];
        Bx_new = c1 * Bx_old + (1.0 - c1) * Bx + c2 * dU[iBx][n];
        By_new = c1 * By_old + (1.0 - c1) * By + c2 * dU[iBy][n];
        Bz_new = c1 * Bz_old + (1.0 - c1) * Bz + c2 * dU[iBz][n];
        Phi_new = c1 * Phi_old + (1.0 - c1) * Phi + c2 * dU[iPhi][n];
        if (DualEnergyFormalism)
        {
          Eint_new = c1 * rho_old * eint_old + (1.0 - c1) * rho * eint + c2 * dU[iEint][n];
        }

        if (D_new < 0 || isnan(D_new))
        {

          printf("UpdateMHDPrim: rho <0 at %" ISYM " %" ISYM " %" ISYM ": rho_old=%" FSYM ", rho=%" FSYM ", rho_new=%" FSYM ", dU[iD]=%" FSYM "\n",
                 i, j, k, rho_old, rho, D_new, dU[iD][n]);
          printf("%f %f %f %f %f %f\n", Bx_old, By_old, Bz_old, Bx_new, By_new, Bz_new);
          D_new = max(rho, SmallRho);
          printf("UpdateMHDPrim: use rho: %" FSYM "\n", D_new);
          //	  D_new = rho;

          return FAIL;
        }

        // convert back to primitives
        vx = S1_new / D_new;
        vy = S2_new / D_new;
        vz = S3_new / D_new;
        etot = Tau_new / D_new;

        if (etot < 0 && EOSType == 0)
        {
          float v2_old = vx_old * vx_old + vy_old * vy_old + vz_old * vz_old;
          float B2_old = Bx_old * vx_old + By_old * By_old + Bz_old * Bz_old;
          printf("UpdateMHDPrim: tau < 0. etot_old=%" GSYM ", etot=%" GSYM ", etot_new=%" GSYM ", v2=%" GSYM ", v2old=%" GSYM ", dU[iTau] = %" GSYM ", dtFixed = %" GSYM "\n",
                 Tau_old / rho_old, Tau / rho, Tau_new / D_new, v2, v2_old, dU[iEtot][n] * CellWidth[0][0] / dtFixed, dtFixed);
          printf("rho_new=%" GSYM ", rho=%" GSYM ", rho_old=%" GSYM ", B2_old/rho_old=%" GSYM "\n", D_new, rho, rho_old, B2_old / rho_old);
          printf("location: %" GSYM ", %" GSYM ", %" GSYM "\n", CellLeftEdge[0][i], CellLeftEdge[1][j], CellLeftEdge[2][k]);
          //return FAIL;
        }

        // if using polytropic EOS, calculate etot directly from density
        if (EOSType > 0)
        {
          EOS(p, D_new, eint, h, cs, dpdrho, dpde, EOSType, 0);
          v2 = vx * vx + vy * vy + vz * vz;
          B2 = Bx_new * Bx_new + By_new * By_new + Bz_new * Bz_new;
          etot = eint + 0.5 * v2 + 0.5 * B2 / D_new;
        }

        BaryonField[DensNum][igrid] = D_new;

        BaryonField[Vel1Num][igrid] = vx;
        BaryonField[Vel2Num][igrid] = vy;
        BaryonField[Vel3Num][igrid] = vz;
        BaryonField[TENum][igrid] = etot;
        BaryonField[B1Num][igrid] = Bx_new;
        BaryonField[B2Num][igrid] = By_new;
        BaryonField[B3Num][igrid] = Bz_new;
        BaryonField[PhiNum][igrid] = Phi_new * exp(-c1 * dtFixed * pow(C_h / C_p, 2));

        if (DualEnergyFormalism)
        {
          v2 = vx * vx + vy * vy + vz * vz;
          B2 = Bx_new * Bx_new + By_new * By_new + Bz_new * Bz_new;
          eint = Eint_new / D_new;
          float emin = SmallT / (Mu * (Gamma - 1.0));

          float eint1 = etot - 0.5 * v2 - 0.5 * B2 / D_new;
          if (eint1 > 0)
          {
            EOS(p, D_new, eint1, h, cs, dpdrho, dpde, EOSType, 2);
          }
          else
          {
            cs = 0.0;
          }
          if (cs * cs > DualEnergyFormalismEta1 * v2 &&
              cs * cs > DualEnergyFormalismEta1 * B2 / D_new &&
              eint1 > 0.5 * eint)
          {
            eint = eint1;
          }
          eint = max(eint, emin);
          BaryonField[GENum][igrid] = eint;
          BaryonField[TENum][igrid] = eint + 0.5 * v2 + 0.5 * B2 / D_new;
          if (BaryonField[GENum][igrid] < 0.0)
          {
            printf("UpdateMHDPrim: eint < 0, cs2=%" GSYM ", eta*v2=%" GSYM ", eint=%" GSYM ", etot=%" GSYM ", 0.5*v2=%" GSYM ", p=%" GSYM ", rho=%" GSYM ",0.5*B2/rho=%" GSYM "\n",
                   cs * cs, DualEnergyFormalismEta1 * v2, eint, etot, 0.5 * v2, p, D_new, 0.5 * B2 / rho);
            printf("dU[%" ISYM "]=%" GSYM ", dU[ieint]=%" GSYM ", eint_old=%" GSYM ",eint1=%" GSYM "\n", iEtot, dU[iEtot][n], dU[iEint][n], eint_old, eint1);
            return FAIL;
          }
        }
      }
    }
  }

  // Convert species from mass fraction to density
  if (NoMultiSpeciesButColors != 1)
    for (field = NEQ_MHD; field < NEQ_MHD + NSpecies + NColor; field++)
      for (n = 0; n < size; n++)
        Prim[field][n] *= Prim[iden][n];

#ifdef USE_KROME
  for (i = 0; i < NKROMEATOMS; i++)
  {
    delete[] atomabund[i];
  }
#else
  this->UpdateElectronDensity();
#endif

  if ((NSpecies + NColor) > 0)
  {
    delete[] D;
    delete[] sum;
  }

  return SUCCESS;
}
