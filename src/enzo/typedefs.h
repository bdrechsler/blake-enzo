#ifndef __typedefs_h_
#define __typedefs_h_
/***********************************************************************
/
/  MISCELANEOUS TYPEDEFS AND ENUMERATIONS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include "CloudyCoolingData.h"
#include "CoolData.h"
#include "RateData.h"
#include "RadiationFieldData.h"
#include "TestProblemData.h"
#include "ResistivityTable.h"

/* These are the different types of baryon fields. */

#ifdef SMALL_INTS
typedef int field_type;
typedef int boundary_type;
typedef int gravity_boundary_type;
typedef int interpolation_type;
typedef int hydro_method;
typedef int star_type;
typedef int enum_type;
typedef int staggering;
typedef int fieldtype;
typedef int mhd_ct_method;
typedef int forcing_type;
#endif

#ifdef LARGE_INTS
typedef long_int field_type;
typedef long_int boundary_type;
typedef long_int gravity_boundary_type;
typedef long_int interpolation_type;
typedef long_int hydro_method;
typedef long_int star_type;
typedef long_int enum_type;
typedef long_int staggering;
typedef long_int fieldtype;
typedef int mhd_ct_method;
typedef long_int forcing_type;
#endif

const field_type
    Density = 0,
    TotalEnergy = 1,
    InternalEnergy = 2,
    Pressure = 3,
    Velocity1 = 4,
    Velocity2 = 5,
    Velocity3 = 6,
    ElectronDensity = 7,
    HIDensity = 8,
    HIIDensity = 9,
    HeIDensity = 10,
    HeIIDensity = 11,
    HeIIIDensity = 12,
    HMDensity = 13,
    H2IDensity = 14,
    H2IIDensity = 15,
    DIDensity = 16,
    DIIDensity = 17,
    HDIDensity = 18,
    SNColour = 19,
    Metallicity = 20,
    ExtraType0 = 21,
    ExtraType1 = 22,
    kphHI = 23,
    PhotoGamma = 24,
    kphHeI = 25,
    gammaHeI = 26,
    kphHeII = 27,
    gammaHeII = 28,
    kdissH2I = 29,
    GravPotential = 30,
    Acceleration0 = 31,
    Acceleration1 = 32,
    Acceleration2 = 33,
    RadPressure0 = 34,
    RadPressure1 = 35,
    RadPressure2 = 36,
    Emissivity0 = 37,

    /* these pseudo-fields are used to access grid data 
   the "g" prefix is to avoid namespace conflict */

    gParticlePosition = 37,
    gParticleVelocity = 38,
    gParticleMass = 39,
    gParticleAcceleration = 40,
    gParticleNumber = 41,
    gParticleType = 42,
    gParticleAttribute = 43,
    gPotentialField = 44,
    gAccelerationField = 45,
    gGravitatingMassField = 46,
    gFlaggingField = 47,
    gVelocity = 48,

    Bfield1 = 49,
    Bfield2 = 50,
    Bfield3 = 51,
    PhiField = 52,
    Phi_pField = 53,
    DebugField = 54,

    DrivingField1 = 55,
    DrivingField2 = 56,
    DrivingField3 = 57,

    AccelerationField1 = 58,
    AccelerationField2 = 59,
    AccelerationField3 = 60,

    Galaxy1Colour = 61,
    Galaxy2Colour = 62,
    /* these are required for Sam Skillman's Shock/Cosmic ray models. */
    Mach = 63,
    PreShockTemperature = 64,
    PreShockDensity = 65,

    /* these are required for Simon Glover's chemistry (which also needs some of the
   other fields, which are used for MultiSpecies) */
    CIDensity = 66,
    CIIDensity = 67,
    OIDensity = 68,
    OIIDensity = 69,
    SiIDensity = 70,
    SiIIDensity = 71,
    SiIIIDensity = 72,
    CHIDensity = 73,
    CH2IDensity = 74,
    CH3IIDensity = 75,
    C2IDensity = 76,
    COIDensity = 77,
    HCOIIDensity = 78,
    OHIDensity = 79,
    H2OIDensity = 80,
    O2IDensity = 81,

    MBHColour = 82,
    ForbiddenRefinement = 83,

    /* FLD radiation module stuff (D. Reynolds) */
    RadiationFreq0 = 84,
    RadiationFreq1 = 85,
    RadiationFreq2 = 86,
    RadiationFreq3 = 87,
    RadiationFreq4 = 88,
    RadiationFreq5 = 89,
    RadiationFreq6 = 90,
    RadiationFreq7 = 91,
    RadiationFreq8 = 92,
    RadiationFreq9 = 93,

    /* Number of ray segments for ray tracing load balancing */
    RaySegments = 94,

    /* Metals from Type Ia SNe */
    MetalSNIaDensity = 95,
    MetalSNIIDensity = 96,

    /* Cosmic Ray Energy Density */
    CRDensity = 97,

    /* IR photodetachment fields */
    kdissH2II = 98,
    kphHM = 99,

    /* Real and Imag of Wave Function */
    RePsi = 101,
    ImPsi = 102,
    FDMDensity = 103,

    HNCIDensity = 104,
    HCNIDensity = 105,
    H2COIDensity = 106,
    HCOIDensity = 107,
    MGIDensity = 108,
    NH3IDensity = 109,
    NOIDensity = 110,
    SIIDensity = 111,
    SIC2IDensity = 112,
    SIC3IDensity = 113,
    SICIDensity = 114,
    SIH2IDensity = 115,
    SIH3IDensity = 116,
    CNIDensity = 117,
    N2IDensity = 118,
    NH2IDensity = 119,
    CH3IDensity = 120,
    CH4IDensity = 121,
    NIDensity = 122,
    NHIDensity = 123,
    SIH4IDensity = 124,
    SIHIDensity = 125,
    SIOIDensity = 126,
    HNOIDensity = 127,
    CH3OHIDensity = 128,
    CO2IDensity = 129,
    H2CNIDensity = 130,
    H2SIOIDensity = 131,
    HNCOIDensity = 132,
    NO2IDensity = 133,
    O2HIDensity = 134,
    OCNIDensity = 135,
    CH3OH_DUSTIDensity = 136,
    HNCO_DUSTIDensity = 137,
    H2CO_DUSTIDensity = 138,
    SIH4_DUSTIDensity = 139,
    H2SIO_DUSTIDensity = 140,
    SIC_DUSTIDensity = 141,
    SIC2_DUSTIDensity = 142,
    SIC3_DUSTIDensity = 143,
    CH4_DUSTIDensity = 144,
    CO_DUSTIDensity = 145,
    H2O_DUSTIDensity = 146,
    NO_DUSTIDensity = 147,
    CO2_DUSTIDensity = 148,
    N2_DUSTIDensity = 149,
    HCN_DUSTIDensity = 150,
    NH3_DUSTIDensity = 151,
    O2_DUSTIDensity = 152,
    NO2_DUSTIDensity = 153,
    HNO_DUSTIDensity = 154,
    O2H_DUSTIDensity = 155,
    H2CN_DUSTIDensity = 156,
    MG_DUSTIDensity = 157,
    HNC_DUSTIDensity = 158,
    E_DUSTIDensity = 159,
    SIO_DUSTIDensity = 160,
    HOCIIDensity = 161,
    CH2IIDensity = 162,
    CHIIDensity = 163,
    H2COIIDensity = 164,
    MGIIDensity = 165,
    NH3IIDensity = 166,
    NOIIDensity = 167,
    SIIIDensity = 168,
    SIC2IIDensity = 169,
    SIC3IIDensity = 170,
    SICIIDensity = 171,
    SIH2IIDensity = 172,
    SIH3IIDensity = 173,
    CNIIDensity = 174,
    COIIDensity = 175,
    N2IIDensity = 176,
    O2IIDensity = 177,
    H2OIIDensity = 178,
    NH2IIDensity = 179,
    OHIIDensity = 180,
    CH4IIDensity = 181,
    NIIDensity = 182,
    HCNIIDensity = 183,
    NHIIDensity = 184,
    SIH4IIDensity = 185,
    SIHIIDensity = 186,
    SIOIIDensity = 187,
    HNOIIDensity = 188,
    H2NOIIDensity = 189,
    H3IIDensity = 190,
    H3COIIDensity = 191,
    H3OIIDensity = 192,
    HCNHIIDensity = 193,
    HCO2IIDensity = 194,
    HeHIIDensity = 195,
    N2HIIDensity = 196,
    O2HIIDensity = 197,
    SIH5IIDensity = 198,
    SIOHIIDensity = 199,

    FieldUndefined = 200;

/*
enum field_type {Density, TotalEnergy, InternalEnergy, Pressure,
		 Velocity1, Velocity2, Velocity3, 
		 ElectronDensity, HIDensity, HIIDensity,  HeIDensity, 
		 HeIIDensity, HeIIIDensity, HMDensity, H2IDensity, 
		 H2IIDensity, DIDensity, DIIDensity, HDIDensity,
                 Metallicity, ExtraType0, ExtraType1, GravPotential,
		 Acceleration0, Acceleration1,Acceleration2,
                 FieldUndefined};
*/

#define FieldTypeIsDensity(A) ((((A) >= TotalEnergy && (A) <= Velocity3) || ((A) >= kphHI && (A) <= kdissH2I) || ((A) >= RadiationFreq0 && (A) <= RaySegments) || ((A) >= Bfield1 && (A) <= AccelerationField3)) ? FALSE : TRUE)
#define FieldTypeIsRadiation(A) ((((A) >= kphHI && (A) <= kdissH2I) || ((A) >= RadiationFreq0 && (A) <= RadiationFreq9)) ? TRUE : FALSE)
#define FieldTypeNoInterpolate(A) (((((A) >= Mach) && ((A) <= PreShockDensity)) || ((A) == GravPotential) || ((A) == RaySegments)) ? TRUE : FALSE)

/* Different stochastic forcing types */
const forcing_type
    None = 0,
    Peak = 1,
    Parabolic = 2,
    Band = 3;

const enum_type
    /* indices used for vectors/Jacobians in SGS model */
    SGSX = 0,
    SGSY = 1,
    SGSZ = 2,
    /* indices used for symmetric tensors */
    SGSXX = 0,
    SGSYY = 1,
    SGSZZ = 2,
    SGSXY = 3,
    SGSYZ = 4,
    SGSXZ = 5,
    SGSYX = 3,
    SGSZY = 4,
    SGSZX = 5;

/* These are the different types of fluid boundary conditions. */

const boundary_type
    reflecting = 0,
    outflow = 1,
    inflow = 2,
    periodic = 3,
    shearing = 4,
    BoundaryUndefined = 5;

// enum boundary_type {reflecting, outflow, inflow, periodic, BoundaryUndefined};

/* These are the different types of gravity boundary conditions. */

const gravity_boundary_type
    TopGridPeriodic = 0,
    TopGridIsolated = 1,
    SubGridIsolated = 2,
    GravityUndefined = 3;

// enum gravity_boundary_type {TopGridPeriodic, TopGridIsolated,
// 				    SubGridIsolated, GravityUndefined};

/* Interpolation types. */

const interpolation_type
    ThirdOrderA = 0,
    SecondOrderA = 1,
    SecondOrderB = 2,
    SecondOrderC = 3,
    FirstOrderA = 4,
    InterpolationUndefined = 5;

// enum interpolation_type {ThirdOrderA, SecondOrderA, SecondOrderB, SecondOrderC,
// 			 FirstOrderA, InterpolationUndefined};

/* Hydrodynamics methods. */

const hydro_method
    PPM_DirectEuler = 0,
    PPM_LagrangeRemap = 1,
    Zeus_Hydro = 2,
    HD_RK = 3,
    MHD_RK = 4,
    NoHydro = 5,
    MHD_Li = 6,
    HydroMethodUndefined = 7;

// enum hydro_method {PPM_DirectEuler, PPM_LagrangeRemap, Zeus_Hydro};

const enum_type iHI = 0, iHeI = 1, iHeII = 2, LW = 3, IR = 4, XRAYS = 5,
                TRACINGSPECTRUM = 6, H2II = 7;
const enum_type Cartesian = 0, Spherical = 1, Cylindrical = 2;
const enum_type PLM = 0, PPM = 1, CENO = 2, WENO3 = 3, WENO5 = 4, ZERO = 5;
const enum_type FluxReconstruction = 0, HLL = 1, Marquina = 2,
                LLF = 3, HLLC = 4, TwoShock = 5, HLLD = 6;
const enum_type Neumann = 0, Dirichlet = 1;
const enum_type Isotropic = 1, Beamed = -2, Episodic = -3;

/* Stanford RK MUSCL solvers support */
//enum {Cartesian, Spherical, Cylindrical};
//enum {PLM, PPM, CENO, WENO3, WENO5};
//enum {FluxReconstruction, HLL, Marquina, LLF, HLLC};

/* These are the different types of poisson cleaining boundary conditions. */
//enum{Neumann, Dirichlet};

const mhd_ct_method
    CT_None = 0,
    CT_BalsaraSpicer = 1,
    CT_Athena_LF = 2,
    CT_Athena_Switch = 3,
    CT_Biermann = 4;

/* Definitions for streaming format */

const staggering VERTEX_CENTERED = 0, CELL_CENTERED = 1;
const fieldtype SCALAR = 1, VECTOR = 3;

/* Star particle types */

const star_type
    PopIII = PARTICLE_TYPE_SINGLE_STAR,
    PopII = PARTICLE_TYPE_CLUSTER,
    NormalStar = PARTICLE_TYPE_STAR,
    SimpleSource = PARTICLE_TYPE_SIMPLE_SOURCE,
    BlackHole = PARTICLE_TYPE_BLACK_HOLE,
    PopIII_CF = PARTICLE_TYPE_COLOR_STAR, // Non-radiating PopIII
    MBH = PARTICLE_TYPE_MBH,
    RadSource = PARTICLE_TYPE_RAD,
    Kravtsov = PARTICLE_TYPE_STAR,
    CenOstriker = PARTICLE_TYPE_STAR,
    AccretingParticle = PARTICLE_TYPE_MUST_REFINE;

/* Define a float/int union. */

union float_int
{
  long_int ival;
  PINT IVAL;
  float fval;
  FLOAT FVAL;
};

//struct ParticleMoveList {
//  int FromGrid;
//  int ToGrid[6];
//  int NumberToMove[6];
//  float_int *Pointer[6];
//};

struct particle_data
{
  FLOAT pos[MAX_DIMENSION];
  float vel[MAX_DIMENSION];
  float mass;
  float attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
  PINT id;
  int type;
  int grid;
  int proc;
};

#include "StarBuffer.h"
struct star_data
{
  StarBuffer data;
  int grid;
  int proc;
};

struct hilbert_data
{
  double hkey;
  int grid_num;
};

// Used in DepositParticleMassFlaggingField.C
struct two_int
{
  int grid;
  int proc;
};

#endif
