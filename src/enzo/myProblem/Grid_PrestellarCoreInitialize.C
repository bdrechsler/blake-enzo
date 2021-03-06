/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR PRESTELLAR CORE) 
/
/  written by: Chia-Jung Hsu
/  date:       March, 2019
/  modified1:  
/
/  PURPOSE: Sets the density and turbulence in a prestellar core.
/
/  RETURNS: FAIL or SUCCESS
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
#include "Grid.h"

void TurbulenceGenerator_FFTW(const int field_type, const int gridsize, const int randomSeed,
        const float power_turb, float* vfield);

void Turbulence_Generator(float **vel, int dim0, int dim1, int dim2, int ind, 
        float kmin, float kmax, float dk,
        FLOAT **LeftEdge, FLOAT **CellWidth, int seed);

int grid::PrestellarCoreInitializeGrid(
            float PrestellarCoreRadius,
            float PrestellarCoreDensity,
            float PrestellarCoreSurfaceDensity,
            float PrestellarCoreDensityJump,
            float PrestellarCoreInternalEnergy,
            float PrestellarCoreAngularVelocity,
            float PrestellarCoreBzField,
            float PrestellarCoreAmbientBzField,
            float PrestellarCoreVelocityDispersion,
            float PrestellarCoreTurbulenceKStart,
            float PrestellarCoreTurbulenceKEnd,
            float PrestellarCoreOPR,
            float PrestellarCoreNDeplete,
            float PrestellarCoreCODeplete,
            int PrestellarCoreRandomSeed,
            int level,
            int* baseDims,
            float* PrestellarCoreInitAbundance,
            float* PrestellarCoreMoleMass,
            float* Turbulence,
            bool SetBaryonField)
{

  /* declarations */

  int size = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* set fields in the rotating core region: x^2+y^2+z^2 < dr^2. */

  // static int count = 0, count2 = 0;
  int index, jndex, kndex, i, j, k, n;
  float zonex, zoney, zonez, radius;
  float r_f = 0.15 * PrestellarCoreRadius,
        r_s = 0.05 * PrestellarCoreRadius,
        r_c = PrestellarCoreRadius;
  float r2 = PrestellarCoreRadius*PrestellarCoreRadius;
  float xcenter = (DomainRightEdge[0] + DomainLeftEdge[0])/2.0;
  float ycenter = (DomainRightEdge[1] + DomainLeftEdge[1])/2.0;
  float zcenter = (DomainRightEdge[2] + DomainLeftEdge[2])/2.0;
  float dSize   = (DomainRightEdge[0] - DomainLeftEdge[0]);
  float PrestellarCoreOuterDensity = PrestellarCoreDensityJump * PrestellarCoreSurfaceDensity;
  float ambientInternalEnergy = PrestellarCoreInternalEnergy / PrestellarCoreDensityJump;

  int DensNum, TENum, GENum, V1Num, V2Num, V3Num;
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

  int *speciesMap = new int [NSpecies+1];
  for (int i=0; i<NSpecies+1; i++) speciesMap[i] = -1;

  NumberOfBaryonFields = 0;
  FieldType[DensNum = NumberOfBaryonFields++] = Density;
  FieldType[TENum   = NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[GENum   = NumberOfBaryonFields++] = InternalEnergy;
  int vel = NumberOfBaryonFields;
  FieldType[V1Num   = NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1 || HydroMethod > 2)
    FieldType[V2Num   = NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2 || HydroMethod > 2)
    FieldType[V3Num   = NumberOfBaryonFields++] = Velocity3;
  if ( UseMHD ) {
    FieldType[B1Num = NumberOfBaryonFields++] = Bfield1;
    FieldType[B2Num = NumberOfBaryonFields++] = Bfield2;
    FieldType[B3Num = NumberOfBaryonFields++] = Bfield3;
    if( HydroMethod == MHD_RK ){
        FieldType[PhiNum = NumberOfBaryonFields++] = PhiField;
    }
    // if (UseDivergenceCleaning) {
    //   FieldType[NumberOfBaryonFields++] = Phi_pField;
    // }
  }

  if (WritePotential)
    FieldType[NumberOfBaryonFields++] = GravPotential;

  int colorfields = NumberOfBaryonFields;

  if (MultiSpecies) {
    FieldType[DeNum     = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum     = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum    = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum    = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum   = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum  = NumberOfBaryonFields++] = HeIIIDensity;
    if (MultiSpecies > 1) {
      FieldType[HMNum   = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum  = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }

#ifdef USE_KROME
    if (MultiSpecies == KROMESPECIES) {
      FieldType[CHINum          = NumberOfBaryonFields++] = CHIDensity;
      FieldType[OINum           = NumberOfBaryonFields++] = OIDensity;
      FieldType[HNCINum         = NumberOfBaryonFields++] = HNCIDensity;
      FieldType[HCNINum         = NumberOfBaryonFields++] = HCNIDensity;
      FieldType[CINum           = NumberOfBaryonFields++] = CIDensity;
      FieldType[H2OINum         = NumberOfBaryonFields++] = H2OIDensity;
      FieldType[OHINum          = NumberOfBaryonFields++] = OHIDensity;
      FieldType[O2INum          = NumberOfBaryonFields++] = O2IDensity;
      FieldType[CH2INum         = NumberOfBaryonFields++] = CH2IDensity;
      FieldType[H2COINum        = NumberOfBaryonFields++] = H2COIDensity;
      FieldType[HCOINum         = NumberOfBaryonFields++] = HCOIDensity;
      FieldType[MGINum          = NumberOfBaryonFields++] = MGIDensity;
      FieldType[NH3INum         = NumberOfBaryonFields++] = NH3IDensity;
      FieldType[NOINum          = NumberOfBaryonFields++] = NOIDensity;
      FieldType[CNINum          = NumberOfBaryonFields++] = CNIDensity;
      FieldType[COINum          = NumberOfBaryonFields++] = COIDensity;
      FieldType[N2INum          = NumberOfBaryonFields++] = N2IDensity;
      FieldType[NH2INum         = NumberOfBaryonFields++] = NH2IDensity;
      FieldType[CH3INum         = NumberOfBaryonFields++] = CH3IDensity;
      FieldType[CH4INum         = NumberOfBaryonFields++] = CH4IDensity;
      FieldType[NINum           = NumberOfBaryonFields++] = NIDensity;
      FieldType[NHINum          = NumberOfBaryonFields++] = NHIDensity;
      FieldType[HNOINum         = NumberOfBaryonFields++] = HNOIDensity;
      FieldType[CH3OHINum       = NumberOfBaryonFields++] = CH3OHIDensity;
      FieldType[CO2INum         = NumberOfBaryonFields++] = CO2IDensity;
      FieldType[H2CNINum        = NumberOfBaryonFields++] = H2CNIDensity;
      FieldType[HNCOINum        = NumberOfBaryonFields++] = HNCOIDensity;
      FieldType[NO2INum         = NumberOfBaryonFields++] = NO2IDensity;
      FieldType[O2HINum         = NumberOfBaryonFields++] = O2HIDensity;
      FieldType[OCNINum         = NumberOfBaryonFields++] = OCNIDensity;
      FieldType[CH3OH_DUSTINum  = NumberOfBaryonFields++] = CH3OH_DUSTIDensity;
      FieldType[HNCO_DUSTINum   = NumberOfBaryonFields++] = HNCO_DUSTIDensity;
      FieldType[H2CO_DUSTINum   = NumberOfBaryonFields++] = H2CO_DUSTIDensity;
      FieldType[CH4_DUSTINum    = NumberOfBaryonFields++] = CH4_DUSTIDensity;
      FieldType[CO_DUSTINum     = NumberOfBaryonFields++] = CO_DUSTIDensity;
      FieldType[H2O_DUSTINum    = NumberOfBaryonFields++] = H2O_DUSTIDensity;
      FieldType[NO_DUSTINum     = NumberOfBaryonFields++] = NO_DUSTIDensity;
      FieldType[CO2_DUSTINum    = NumberOfBaryonFields++] = CO2_DUSTIDensity;
      FieldType[N2_DUSTINum     = NumberOfBaryonFields++] = N2_DUSTIDensity;
      FieldType[HCN_DUSTINum    = NumberOfBaryonFields++] = HCN_DUSTIDensity;
      FieldType[NH3_DUSTINum    = NumberOfBaryonFields++] = NH3_DUSTIDensity;
      FieldType[O2_DUSTINum     = NumberOfBaryonFields++] = O2_DUSTIDensity;
      FieldType[NO2_DUSTINum    = NumberOfBaryonFields++] = NO2_DUSTIDensity;
      FieldType[HNO_DUSTINum    = NumberOfBaryonFields++] = HNO_DUSTIDensity;
      FieldType[O2H_DUSTINum    = NumberOfBaryonFields++] = O2H_DUSTIDensity;
      FieldType[H2CN_DUSTINum   = NumberOfBaryonFields++] = H2CN_DUSTIDensity;
      FieldType[MG_DUSTINum     = NumberOfBaryonFields++] = MG_DUSTIDensity;
      FieldType[HNC_DUSTINum    = NumberOfBaryonFields++] = HNC_DUSTIDensity;
      FieldType[E_DUSTINum      = NumberOfBaryonFields++] = E_DUSTIDensity;
      FieldType[HCOIINum        = NumberOfBaryonFields++] = HCOIIDensity;
      FieldType[HOCIINum        = NumberOfBaryonFields++] = HOCIIDensity;
      FieldType[CIINum          = NumberOfBaryonFields++] = CIIDensity;
      FieldType[CH2IINum        = NumberOfBaryonFields++] = CH2IIDensity;
      FieldType[CHIINum         = NumberOfBaryonFields++] = CHIIDensity;
      FieldType[H2COIINum       = NumberOfBaryonFields++] = H2COIIDensity;
      FieldType[MGIINum         = NumberOfBaryonFields++] = MGIIDensity;
      FieldType[NH3IINum        = NumberOfBaryonFields++] = NH3IIDensity;
      FieldType[NOIINum         = NumberOfBaryonFields++] = NOIIDensity;
      FieldType[CNIINum         = NumberOfBaryonFields++] = CNIIDensity;
      FieldType[COIINum         = NumberOfBaryonFields++] = COIIDensity;
      FieldType[N2IINum         = NumberOfBaryonFields++] = N2IIDensity;
      FieldType[O2IINum         = NumberOfBaryonFields++] = O2IIDensity;
      FieldType[H2OIINum        = NumberOfBaryonFields++] = H2OIIDensity;
      FieldType[NH2IINum        = NumberOfBaryonFields++] = NH2IIDensity;
      FieldType[OIINum          = NumberOfBaryonFields++] = OIIDensity;
      FieldType[OHIINum         = NumberOfBaryonFields++] = OHIIDensity;
      FieldType[CH3IINum        = NumberOfBaryonFields++] = CH3IIDensity;
      FieldType[CH4IINum        = NumberOfBaryonFields++] = CH4IIDensity;
      FieldType[NIINum          = NumberOfBaryonFields++] = NIIDensity;
      FieldType[HCNIINum        = NumberOfBaryonFields++] = HCNIIDensity;
      FieldType[NHIINum         = NumberOfBaryonFields++] = NHIIDensity;
      FieldType[HNOIINum        = NumberOfBaryonFields++] = HNOIIDensity;
      FieldType[H2NOIINum       = NumberOfBaryonFields++] = H2NOIIDensity;
      FieldType[H3IINum         = NumberOfBaryonFields++] = H3IIDensity;
      FieldType[H3COIINum       = NumberOfBaryonFields++] = H3COIIDensity;
      FieldType[H3OIINum        = NumberOfBaryonFields++] = H3OIIDensity;
      FieldType[HCNHIINum       = NumberOfBaryonFields++] = HCNHIIDensity;
      FieldType[HCO2IINum       = NumberOfBaryonFields++] = HCO2IIDensity;
      FieldType[HeHIINum        = NumberOfBaryonFields++] = HeHIIDensity;
      FieldType[N2HIINum        = NumberOfBaryonFields++] = N2HIIDensity;
      FieldType[O2HIINum        = NumberOfBaryonFields++] = O2HIIDensity;


      speciesMap[0] = DeNum;
      speciesMap[1] = CHINum;
      speciesMap[2] = OINum;
      speciesMap[3] = HNCINum;
      speciesMap[4] = HCNINum;
      speciesMap[5] = H2INum;
      speciesMap[6] = CINum;
      speciesMap[7] = HINum;
      speciesMap[8] = H2OINum;
      speciesMap[9] = OHINum;
      speciesMap[10] = O2INum;
      speciesMap[11] = CH2INum;
      speciesMap[12] = H2COINum;
      speciesMap[13] = HCOINum;
      speciesMap[14] = MGINum;
      speciesMap[15] = NH3INum;
      speciesMap[16] = NOINum;
      speciesMap[17] = CNINum;
      speciesMap[18] = COINum;
      speciesMap[19] = N2INum;
      speciesMap[20] = NH2INum;
      speciesMap[21] = CH3INum;
      speciesMap[22] = CH4INum;
      speciesMap[23] = NINum;
      speciesMap[24] = NHINum;
      speciesMap[25] = HeINum;
      speciesMap[26] = HNOINum;
      speciesMap[27] = CH3OHINum;
      speciesMap[28] = CO2INum;
      speciesMap[29] = H2CNINum;
      speciesMap[30] = HNCOINum;
      speciesMap[31] = NO2INum;
      speciesMap[32] = O2HINum;
      speciesMap[33] = OCNINum;
      speciesMap[34] = CH3OH_DUSTINum;
      speciesMap[35] = HNCO_DUSTINum;
      speciesMap[36] = H2CO_DUSTINum;
      speciesMap[37] = CH4_DUSTINum;
      speciesMap[38] = CO_DUSTINum;
      speciesMap[39] = H2O_DUSTINum;
      speciesMap[40] = NO_DUSTINum;
      speciesMap[41] = CO2_DUSTINum;
      speciesMap[42] = N2_DUSTINum;
      speciesMap[43] = HCN_DUSTINum;
      speciesMap[44] = NH3_DUSTINum;
      speciesMap[45] = O2_DUSTINum;
      speciesMap[46] = NO2_DUSTINum;
      speciesMap[47] = HNO_DUSTINum;
      speciesMap[48] = O2H_DUSTINum;
      speciesMap[49] = H2CN_DUSTINum;
      speciesMap[50] = MG_DUSTINum;
      speciesMap[51] = HNC_DUSTINum;
      speciesMap[52] = E_DUSTINum;
      speciesMap[53] = HCOIINum;
      speciesMap[54] = HIINum;
      speciesMap[55] = HOCIINum;
      speciesMap[56] = CIINum;
      speciesMap[57] = CH2IINum;
      speciesMap[58] = CHIINum;
      speciesMap[59] = H2COIINum;
      speciesMap[60] = MGIINum;
      speciesMap[61] = NH3IINum;
      speciesMap[62] = NOIINum;
      speciesMap[63] = CNIINum;
      speciesMap[64] = COIINum;
      speciesMap[65] = N2IINum;
      speciesMap[66] = O2IINum;
      speciesMap[67] = H2OIINum;
      speciesMap[68] = NH2IINum;
      speciesMap[69] = OIINum;
      speciesMap[70] = OHIINum;
      speciesMap[71] = CH3IINum;
      speciesMap[72] = CH4IINum;
      speciesMap[73] = NIINum;
      speciesMap[74] = HCNIINum;
      speciesMap[75] = NHIINum;
      speciesMap[76] = H2IINum;
      speciesMap[77] = HeIINum;
      speciesMap[78] = HNOIINum;
      speciesMap[79] = H2NOIINum;
      speciesMap[80] = H3IINum;
      speciesMap[81] = H3COIINum;
      speciesMap[82] = H3OIINum;
      speciesMap[83] = HCNHIINum;
      speciesMap[84] = HCO2IINum;
      speciesMap[85] = HeHIINum;
      speciesMap[86] = N2HIINum;
      speciesMap[87] = O2HIINum;

    }
#endif
  }

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (!SetBaryonField)
    return SUCCESS;

  this->AllocateGrids(); 

  float coreMass1 = 0.0, coreMass2 = 0.0, coreMass3;
  float xlen, ylen, zlen, Vol;
  for (i = 0; i < size; i++) {

    index = i % GridDimension[0];
    jndex = ((i-index) % (GridDimension[0]*GridDimension[1]))/GridDimension[0];
    kndex = (i-index-jndex*GridDimension[0])/GridDimension[0]/GridDimension[1];

    zonex = *(CellLeftEdge[0] + index) + 0.5*(*(CellWidth[0] + index)) - xcenter;
    zoney = *(CellLeftEdge[1] + jndex) + 0.5*(*(CellWidth[1] + jndex)) - ycenter;
    zonez = *(CellLeftEdge[2] + kndex) + 0.5*(*(CellWidth[2] + kndex)) - zcenter;
    radius = sqrt(zonex*zonex + zoney*zoney + zonez*zonez);

    // if (count<1) 
    //     printf("coreDens: %13.7e, outerDens: %13.7e, radius: %13.7e, coreRad: %13.7e\n", 
    //         PrestellarCoreDensity, PrestellarCoreOuterDensity, 
    //         radius, r_c);

    BaryonField[DensNum][i] = PrestellarCoreOuterDensity
        + (PrestellarCoreDensity - PrestellarCoreOuterDensity)
        / (1.0 + pow( (radius/r_f), 1.5))
        * (0.5 - 0.5*tanh( (radius-r_c)/ r_s ));

    xlen = *(CellWidth[0] + index);
    ylen = *(CellWidth[1] + jndex);
    zlen = *(CellWidth[2] + kndex);
    Vol = xlen * ylen * zlen;
    if (!ParallelRootGridIO){
      if (radius < r_c) coreMass1 += BaryonField[DensNum][i] * Vol;
      if (BaryonField[DensNum][i] > PrestellarCoreOuterDensity*10.0) coreMass2 += BaryonField[DensNum][i] * Vol;
      if (BaryonField[DensNum][i] > PrestellarCoreOuterDensity*1.1) coreMass3 += BaryonField[DensNum][i] * Vol;
    }

    // ambient gas temperature is 10 times larger than core
    float transfer = 0.5*tanh( (radius-r_c)/ r_s );
    BaryonField[TENum][i] = ambientInternalEnergy * (0.5 + transfer)
                          + PrestellarCoreInternalEnergy * (0.5 - transfer);
    // BaryonField[TENum][i] = PrestellarCoreInternalEnergy * 10.0;

    if (DualEnergyFormalism) {
      BaryonField[GENum][i] = ambientInternalEnergy * (0.5 + transfer)
                            + PrestellarCoreInternalEnergy * (0.5 - transfer);
    }
    
    if (zonex*zonex + zoney*zoney + zonez*zonez < r2) {
      // BaryonField[TENum][i] /= 10.0;
      BaryonField[V1Num][i] =   PrestellarCoreAngularVelocity*zoney;
      BaryonField[V2Num][i] = - PrestellarCoreAngularVelocity*zonex;
      BaryonField[V3Num][i] = 0.0; // the rotation axis is || to z.

      // if there is no turbulence, calculate kinetic energy now
      if (PrestellarCoreTurbulenceKStart > PrestellarCoreTurbulenceKEnd)
        BaryonField[TENum][i] += (BaryonField[V1Num][i]*BaryonField[V1Num][i] 
                               +  BaryonField[V2Num][i]*BaryonField[V2Num][i]
                               +  BaryonField[V3Num][i]*BaryonField[V3Num][i])/2.0;
    }

    if (UseMHD) {
      radius = sqrt(zonex*zonex + zoney*zoney);
      BaryonField[B1Num][i] = 0.0;
      BaryonField[B2Num][i] = 0.0;
      BaryonField[B3Num][i] = PrestellarCoreAmbientBzField + 
        (PrestellarCoreBzField - PrestellarCoreAmbientBzField)
        / (1.0 + pow( (radius/r_f), 0.5))
        * (0.5 - 0.5*tanh( (radius-r_c)/ r_s ));
      BaryonField[TENum][i] += 0.5*BaryonField[B3Num][i]*BaryonField[B3Num][i]/BaryonField[DensNum][i];
      // printf("zonex: %13.7e, zoney: %13.7e, zonez: %13.7e, Bz: %13.7e \n", zonex, zoney, zonez, BaryonField[B3Num][i]);
    }

#ifdef USE_KROME
    if (MultiSpecies == KROMESPECIES) {
      // for (int abNum=0; abNum<NSpecies+1; abNum++) {
      //   int speciesNum = speciesMap[abNum];
      //   if (speciesNum != -1) {
      //     BaryonField[speciesNum][i] = 1.0e-20 * PrestellarCoreMoleMass[abNum] 
      //                                * BaryonField[0][i] / 1.40045;
      //   }
      // }
      for (int speciesNum = DeNum; speciesNum <= O2HIINum; speciesNum ++) {
        BaryonField[speciesNum][i] = 1e-20*BaryonField[0][i];
      }
      if (PrestellarCoreOPR < 999.0){
        /* set your preferable initial abundances */
        BaryonField[H2INum][i]          = 1.00e+0*BaryonField[0][i] / 1.4;
        BaryonField[HINum][i]           = 5.00e-5*BaryonField[0][i] / 1.4;
        BaryonField[HeINum][i]          = 3.90e-1*BaryonField[0][i] / 1.4;
        BaryonField[NINum][i]           = 14.0*7.5e-5*BaryonField[0][i] / 1.4;
        BaryonField[OINum][i]           = 16.0*1.8e-4*BaryonField[0][i] / 1.4;
        BaryonField[COINum][i]          = 28.0*1.4e-4*BaryonField[0][i] / 1.4;
        BaryonField[MGINum][i]          = 24.0*7.0e-9*BaryonField[0][i] / 1.4;
      }
      else {
        // read in table
        for (int abNum=0; abNum<NSpecies+1; abNum++) {
          int speciesNum = speciesMap[abNum];
          if (speciesNum != -1) {
            BaryonField[speciesNum][i] = PrestellarCoreInitAbundance[abNum] 
                                       * PrestellarCoreMoleMass[abNum] 
                                       * BaryonField[0][i] / 1.40045;
          }
        }
      }
    }
#else
    if (MultiSpecies) {
      for (int speciesNum = DeNum; speciesNum <= HeIIINum; speciesNum ++) {
        BaryonField[speciesNum][i] = 1e-20*BaryonField[0][i];
      }
      BaryonField[H2INum][i]          = 5.00e-1*BaryonField[0][i] / 1.412;
      BaryonField[HINum][i]           = 5.00e-1*BaryonField[0][i] / 1.412;
      BaryonField[HeINum][i]          = 4.00e-1*BaryonField[0][i] / 1.412;
    }
#endif


    // if (count < 1){
    //     printf("zonex: %13.7e, zoney: %13.7e\n", zonex, zoney);
    //     printf("Eint: %13.7e, xvel: %13.7e\n", BaryonField[1][i], BaryonField[2][i]);
    // }
    // count++;
  }

  if (!ParallelRootGridIO){
    printf("Total core mass inside Rc: %13.7e\n", coreMass1);
    printf("Total core mass where rho > 10.0 rho0: %13.7e\n", coreMass2);
    printf("Total core mass where rho > 1.1 rho0: %13.7e\n", coreMass3);
  }

  // float k1 = 65.0/dSize, 
  //       k2 = 73.0/dSize,
  //       //dk = 1.0/dSize;
  //       //k2 = float(int(GridDimension[0]/10.0)), 
  //       dk = max(1.0/dSize,(k2-k1)/13.0);
  float k1 = PrestellarCoreTurbulenceKStart, k2 = PrestellarCoreTurbulenceKEnd, 
        dk = 1.0;
        //dk = max(1.0, (k2-k1)/10.0);


  if (k1 < k2){

    // increase grid dimensional with level to initialize turbulence in higher resolution
    int activesize = 1, turbsize=1;
    for (dim = 0; dim < GridRank; dim++) {
      activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);
      turbsize *= baseDims[dim] * pow(RefineBy, level);
    }
    int gridsize = baseDims[0] * pow(RefineBy, level);
    float *TurbulenceVelocity[3];
    for (dim = 0; dim < GridRank; dim++) {
      TurbulenceVelocity[dim] = new float[activesize];
    }

    // look for the sub-region from the grids of turbulence initialization
    int TurbGridStartIdx[3] = {0},
        TurbGridEndIdx[3] = {0};

    for (dim = 0; dim < GridRank; dim++) {
      TurbGridStartIdx[dim] = (GridLeftEdge[dim] - DomainLeftEdge[dim])/CellWidth[dim][0];
      TurbGridEndIdx[dim]   = TurbGridStartIdx[dim] + (GridDimension[dim]-2*NumberOfGhostZones) - 1;
      printf("StartIndex=%d, EndIndex=%d\n", TurbGridStartIdx[dim], TurbGridEndIdx[dim]);
    }
 
    for (dim = 0; dim < GridRank; dim++)
    for (i = TurbGridStartIdx[0], n=0; i <= TurbGridEndIdx[0]; i++)
    for (j = TurbGridStartIdx[1]; j <= TurbGridEndIdx[1]; j++)
    for (k = TurbGridStartIdx[2]; k <= TurbGridEndIdx[2]; k++, n++){
      int ref = dim + 3*i + 3*gridsize*j + 3*gridsize*gridsize*k;
      TurbulenceVelocity[dim][n] = Turbulence[ref];
    }

    // Old turbulence generator
    // Turbulence_Generator(TurbulenceVelocity, 
    //         GridDimension[0]-2*NumberOfGhostZones,
    //         GridDimension[1]-2*NumberOfGhostZones,
    //         GridDimension[2]-2*NumberOfGhostZones,
    //         4.0, k1, k2, dk,
    //         CellLeftEdge, CellWidth, 65536);    

    n = 0;
    for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++, n++) {
      int igrid = i + GridDimension[0]*(j+k*GridDimension[1]);
      // zonex = *(CellLeftEdge[0] + i) + 0.5*(*(CellWidth[0] + i)) - xcenter;
      // zoney = *(CellLeftEdge[1] + j) + 0.5*(*(CellWidth[1] + j)) - ycenter;
      // zonez = *(CellLeftEdge[2] + k) + 0.5*(*(CellWidth[2] + k)) - zcenter;
      // radius = sqrt(zonex*zonex + zoney*zoney + zonez*zonez);

      BaryonField[V1Num][igrid] += TurbulenceVelocity[0][n];
                                 //* (0.5 - 0.5*tanh( (radius-r_c)/ r_s ));
      BaryonField[V2Num][igrid] += TurbulenceVelocity[1][n];
                                 //* (0.5 - 0.5*tanh( (radius-r_c)/ r_s ));
      BaryonField[V3Num][igrid] += TurbulenceVelocity[2][n];
                                 //* (0.5 - 0.5*tanh( (radius-r_c)/ r_s ));
      // total energy
      BaryonField[TENum][igrid] += (BaryonField[V1Num][igrid]*BaryonField[V1Num][igrid] 
                                 +  BaryonField[V2Num][igrid]*BaryonField[V2Num][igrid]
                                 +  BaryonField[V3Num][igrid]*BaryonField[V3Num][igrid])/2.0;

      // float vel = sqrt(BaryonField[2][igrid]*BaryonField[2][igrid] 
      //           +      BaryonField[3][igrid]*BaryonField[3][igrid]
      //           +      BaryonField[4][igrid]*BaryonField[4][igrid]);
      // printf("xvel: %13.7e, yvel: %13.7e, zvel: %13.7e, vel: %13.7e \n",
      //         BaryonField[2][igrid], BaryonField[3][igrid], BaryonField[4][igrid], vel);
    }}}
  
    int idx = GridStartIndex[0] + GridDimension[0]*(GridStartIndex[1]+GridStartIndex[2]*GridDimension[1]);
    // if(count2<1)
    //   printf("Eint: %13.7e, xvel: %13.7e\n", BaryonField[1][idx], BaryonField[2][idx]);
    // count2 ++;
  
    for (dim = 0; dim < GridRank; dim++)
      delete [] TurbulenceVelocity[dim];

  } // if (k1 < k2), set turbulence
 
  delete [] speciesMap;

  return SUCCESS;
}

