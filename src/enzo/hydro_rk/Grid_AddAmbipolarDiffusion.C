/***********************************************************************
/
/  GRID CLASS (ADD AMBIPOLAR DIFFUSION TERM)
/
/  written by: Peng Wang
/  date:       July, 2008
/  modified1:  Duncan Christie, Spring 2016 
/              Rewrote much of routine to fix timestep errors within the code.
/
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "EOS.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

float GetResistivity(float rho, float B, float T, float Time,
                     float rho_HI, float rho_H2I, float rho_HeI,
                     float rho_E, float rho_HII, float rho_HeII, float rho_H3II,
                     float rho_MII, float rho_AII);

int grid::AddAmbipolarDiffusion()
{

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum;
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

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1,
        TimeUnits = 1.0, VelocityUnits = 1.0;
  float A;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
           &TimeUnits, &VelocityUnits, Time);

  float T;

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
  {
    size *= GridDimension[dim];
  }

  double *D[3];
  for (int dim = 0; dim < 3; dim++)
  {
    D[dim] = new double[size];
  }

  FLOAT dt_ad = dtFixed, dt1;
  FLOAT dx = CellWidth[0][0];
  FLOAT dy = (GridRank > 1) ? CellWidth[1][0] : 1.;
  FLOAT dz = (GridRank > 2) ? CellWidth[2][0] : 1.;

  double rho, Bx, By, Bz, B2, Bmag;
  int nx = (GridRank > 1) ? GridDimension[0] : 0, ny = (GridRank > 2) ? GridDimension[1] : 0;
  int igrid;
  int kmp = (GridRank > 2) ? 1 : 0,
      jmp = (GridRank > 1) ? 1 : 0;

  float *Temp = new float[size]; // Gas temperature

  // Get the temperature
  if (this->ComputeTemperatureField(Temp) == FAIL)
  {
    ENZO_FAIL("Error in grid->ComputeTemperatureField.");
  }

  /* First calculate D */

  double dBxdy, dBxdz, dBydx, dBydz, dBzdx, dBzdy;
  double j2;
  double clight_code = clight / VelocityUnits;
  float rho_HI = 0.0, rho_H2I = 0.0, rho_HeI = 0.0,
        rho_E = 0.0, rho_HII = 0.0, rho_HeII = 0.0,
        rho_H3II = 0.0, rho_MII = 0.0, rho_AII = 0.0;

  for (int k = GridStartIndex[2] - kmp; k <= GridEndIndex[2] + kmp; k++)
  {
    for (int j = GridStartIndex[1] - jmp; j <= GridEndIndex[1] + jmp; j++)
    {
      for (int i = GridStartIndex[0] - 1; i <= GridEndIndex[0] + 1; i++)
      {

        igrid = i + (j + k * ny) * nx;

        rho = BaryonField[DensNum][igrid];
        Bx = BaryonField[B1Num][igrid];
        By = BaryonField[B2Num][igrid];
        Bz = BaryonField[B3Num][igrid];

        B2 = Bx * Bx + By * By + Bz * Bz;
        Bmag = sqrt(B2);

        T = Temp[igrid];

        dBxdy = (GridRank > 1) ? ((BaryonField[B1Num][igrid + nx] - BaryonField[B1Num][igrid - nx])) / (2.0 * dy) : 0;
        dBxdz = (GridRank > 2) ? ((BaryonField[B1Num][igrid + nx * ny] - BaryonField[B1Num][igrid - nx * ny])) / (2.0 * dz) : 0;
        dBydx = ((BaryonField[B2Num][igrid + 1] - BaryonField[B2Num][igrid - 1])) / (2.0 * dx);
        dBydz = (GridRank > 2) ? ((BaryonField[B2Num][igrid + nx * ny] - BaryonField[B2Num][igrid - nx * ny])) / (2.0 * dz) : 0;
        dBzdx = ((BaryonField[B3Num][igrid + 1] - BaryonField[B3Num][igrid - 1])) / (2.0 * dx);
        dBzdy = (GridRank > 1) ? ((BaryonField[B3Num][igrid + nx] - BaryonField[B3Num][igrid - nx])) / (2.0 * dy) : 0;

        D[0][igrid] = (Bx * Bz * (dBydx - dBxdy) - (By * By + Bz * Bz) * (dBzdy - dBydz) + Bx * By * (dBxdz - dBzdx)) / B2;
        D[1][igrid] = (Bx * By * (dBzdy - dBydz) - (Bx * Bx + Bz * Bz) * (dBxdz - dBzdx) + By * Bz * (dBydx - dBxdy)) / B2;
        D[2][igrid] = (By * Bz * (dBxdz - dBzdx) - (Bx * Bx + By * By) * (dBydx - dBxdy) + Bx * Bz * (dBzdy - dBydz)) / B2;

#ifdef USE_KROME
        rho_HI = BaryonField[HINum][igrid];
        rho_H2I = BaryonField[H2INum][igrid];
        rho_HeI = BaryonField[HeINum][igrid];
        rho_E = BaryonField[DeNum][igrid];
        rho_HII = BaryonField[HIINum][igrid];
        rho_HeII = BaryonField[HeIINum][igrid];
        rho_H3II = BaryonField[H3IINum][igrid];
        rho_MII = BaryonField[HCOIINum][igrid] + BaryonField[HOCIINum][igrid] + BaryonField[CH2IINum][igrid] + BaryonField[CHIINum][igrid] + BaryonField[H2COIINum][igrid] + BaryonField[NH3IINum][igrid] + BaryonField[NOIINum][igrid] + BaryonField[CNIINum][igrid] + BaryonField[COIINum][igrid] + BaryonField[N2IINum][igrid] + BaryonField[O2IINum][igrid] + BaryonField[H2OIINum][igrid] + BaryonField[NH2IINum][igrid] + BaryonField[OHIINum][igrid] + BaryonField[CH3IINum][igrid] + BaryonField[CH3IINum][igrid] + BaryonField[HCNIINum][igrid] + BaryonField[NHIINum][igrid] + BaryonField[HNOIINum][igrid] + BaryonField[H2NOIINum][igrid] + BaryonField[H3COIINum][igrid] + BaryonField[H3OIINum][igrid] + BaryonField[HCNHIINum][igrid] + BaryonField[HCO2IINum][igrid] + BaryonField[HeHIINum][igrid] + BaryonField[N2HIINum][igrid] + BaryonField[O2HIINum][igrid];
        rho_AII = BaryonField[CIINum][igrid] + BaryonField[MGIINum][igrid] + BaryonField[OIINum][igrid] + BaryonField[NIINum][igrid];
#endif

        A = GetResistivity(rho, Bmag, T, Time, rho_HI, rho_H2I, rho_HeI,
                           rho_E, rho_HII, rho_HeII, rho_H3II,
                           rho_MII, rho_AII);

        /*Add in the Joule heating*/

        if (ADJouleHeating)
        {
          j2 = D[0][igrid] * D[0][igrid] + D[1][igrid] * D[1][igrid] + D[2][igrid] * D[2][igrid];
          j2 = dtFixed * A * j2 / rho;
          BaryonField[TENum][igrid] += j2;
          if (DualEnergyFormalism)
          {
            BaryonField[GENum][igrid] += j2;
          }
        }

        D[0][igrid] *= A;
        D[1][igrid] *= A;
        D[2][igrid] *= A;
      }
    }
  }

  /* Then calculte the AD term and update B field */

  double AD[3], Bx_old, By_old, Bz_old, Bx_new, By_new, Bz_new, B2_old, B2_new, B2_err, x, y, z;
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
  {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
    {
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
      {

        igrid = i + (j + k * ny) * nx;

        x = CellLeftEdge[0][i] + 0.5 * CellWidth[0][i];
        y = CellLeftEdge[1][j] + 0.5 * CellWidth[1][j];
        z = CellLeftEdge[2][k] + 0.5 * CellWidth[2][k];

        AD[0] = ((GridRank > 1) ? (D[2][igrid + nx] - D[2][igrid - nx]) / (2.0 * dy) : 0) - ((GridRank > 2) ? (D[1][igrid + nx * ny] - D[1][igrid - nx * ny]) / (2.0 * dz) : 0);
        AD[1] = ((GridRank > 2) ? (D[0][igrid + nx * ny] - D[0][igrid - nx * ny]) / (2.0 * dz) : 0) - (D[2][igrid + 1] - D[2][igrid - 1]) / (2.0 * dx);
        AD[2] = (D[1][igrid + 1] - D[1][igrid - 1]) / (2.0 * dx) - ((GridRank > 1) ? (D[0][igrid + nx] - D[0][igrid - nx]) / (2.0 * dy) : 0);

        /* Subtract out the magnetic energy*/
        B2 = pow(BaryonField[B1Num][igrid], 2) + pow(BaryonField[B2Num][igrid], 2) + pow(BaryonField[B3Num][igrid], 2);
        rho = BaryonField[DensNum][igrid];
        BaryonField[TENum][igrid] -= 0.5 * B2 / rho;

        /* Update the field*/

        BaryonField[B1Num][igrid] += dtFixed * AD[0];
        BaryonField[B2Num][igrid] += dtFixed * AD[1];
        BaryonField[B3Num][igrid] += dtFixed * AD[2];

        /* Add back in the magnetic energy*/
        B2 = pow(BaryonField[B1Num][igrid], 2) + pow(BaryonField[B2Num][igrid], 2) + pow(BaryonField[B3Num][igrid], 2);
        BaryonField[TENum][igrid] += 0.5 * B2 / rho;
      }
    }
  }

  for (int dim = 0; dim < 3; dim++)
  {
    delete D[dim];
  }

  delete[] Temp;

  return SUCCESS;
}
