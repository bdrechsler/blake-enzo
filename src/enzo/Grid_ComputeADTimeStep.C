/***********************************************************************
/
/  GRID CLASS (Compute thermal conduction time-scale)
/
/  written by:  Duncan Christie (based on Grid_ComputeConductionTimeStep.C)
/  date:        June 2016
/
/  PURPOSE:  Calculates the shortest time scale for ambipolar diffusion on a
/  grid
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "phys_constants.h"
#include "CosmologyParameters.h"

// Function prototypes
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, double *MassUnits, FLOAT Time);
float GetResistivity(float rho, float B, float T, float Time,
                     float rho_HI, float rho_H2I, float rho_HeI,
                     float rho_E, float rho_HII, float rho_HeII, float rho_H3II,
                     float rho_MII, float rho_AII);

// Member functions
int grid::ComputeADTimeStep(float &dt)
{

  dt = huge_number;

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
  if (NumberOfBaryonFields == 0)
    return SUCCESS;
  this->DebugCheck("ComputeADTimeStep");

  // Some locals
  float rho, Bx, By, Bz, B2, eta, Bmag;
  int i, j, k;

  if (!UseAmbipolarDiffusion)
  {
    return SUCCESS;
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num,
      B1Num, B2Num, B3Num, PhiNum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num,
                                   TENum, B1Num, B2Num, B3Num, PhiNum);

#ifdef USE_KROME
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, CHINum, OINum, HNCINum, HCNINum, CINum, H2OINum,
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
            HCO2IINum, HeHIINum, N2HIINum, O2HIINum) == FAIL)
    {
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }
#endif

  /* Compute the field size. */
  int size = 1;
  int dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  float *Temp = new float[size];
  float T;
  /* The ambipolar diffusion module requires the temperature*/
  if (this->ComputeTemperatureField(Temp) == FAIL)
  {
    ENZO_FAIL("ComputeADTimeStep(): Error when calling ComputeTemperatureField()\n");
  }

  /*Grid sizes*/
  FLOAT dxinv = 1.0 / CellWidth[0][0];
  FLOAT dyinv = (GridRank > 1) ? 1.0 / CellWidth[1][0] : 0.0;
  FLOAT dzinv = (GridRank > 2) ? 1.0 / CellWidth[2][0] : 0.0;

  dxinv = max(dxinv, max(dyinv, dzinv));

  float dtinv_ad = 0.;
  int n = 0;
  float rho_HI = 0.0, rho_H2I = 0.0, rho_HeI = 0.0,
        rho_E = 0.0, rho_HII = 0.0, rho_HeII = 0.0,
        rho_H3II = 0.0, rho_MII = 0.0, rho_AII = 0.0;
  for (k = 0; k < GridDimension[2]; k++)
  {
    for (j = 0; j < GridDimension[1]; j++)
    {
      for (i = 0; i < GridDimension[0]; i++, n++)
      {
        rho = BaryonField[DensNum][n];
        Bx = BaryonField[B1Num][n];
        By = BaryonField[B2Num][n];
        Bz = BaryonField[B3Num][n];

        B2 = Bx * Bx + By * By + Bz * Bz;

        Bmag = sqrt(B2);
        T = Temp[n];

#ifdef USE_KROME
        rho_HI = BaryonField[HINum][n];
        rho_H2I = BaryonField[H2INum][n];
        rho_HeI = BaryonField[HeINum][n];
        rho_E = BaryonField[DeNum][n];
        rho_HII = BaryonField[HIINum][n];
        rho_HeII = BaryonField[HeIINum][n];
        rho_H3II = BaryonField[H3IINum][n];
        rho_MII = BaryonField[HCOIINum][n] + BaryonField[HOCIINum][n] + BaryonField[CH2IINum][n] + BaryonField[CHIINum][n] + BaryonField[H2COIINum][n] + BaryonField[NH3IINum][n] + BaryonField[NOIINum][n] + BaryonField[CNIINum][n] + BaryonField[COIINum][n] + BaryonField[N2IINum][n] + BaryonField[O2IINum][n] + BaryonField[H2OIINum][n] + BaryonField[NH2IINum][n] + BaryonField[OHIINum][n] + BaryonField[CH3IINum][n] + BaryonField[CH3IINum][n] + BaryonField[HCNIINum][n] + BaryonField[NHIINum][n] + BaryonField[HNOIINum][n] + BaryonField[H2NOIINum][n] + BaryonField[H3COIINum][n] + BaryonField[H3OIINum][n] + BaryonField[HCNHIINum][n] + BaryonField[HCO2IINum][n] + BaryonField[HeHIINum][n] + BaryonField[N2HIINum][n] + BaryonField[O2HIINum][n];
        rho_AII = BaryonField[CIINum][n] + BaryonField[MGIINum][n] + BaryonField[OIINum][n] + BaryonField[NIINum][n];
#endif

        eta = GetResistivity(rho, Bmag, T, Time, rho_HI, rho_H2I, rho_HeI,
                             rho_E, rho_HII, rho_HeII, rho_H3II,
                             rho_MII, rho_AII);

        dtinv_ad = max(dtinv_ad, 10.e0 * eta * dxinv * dxinv);
      }
    }
  }

  delete[] Temp;

  if (dtinv_ad != 0)
    dt = 1.0 / dtinv_ad;

  if (dt == 0.)
    fprintf(stderr, "AD Error: %f\n", dtinv_ad);

  return SUCCESS;
}
