/***********************************************************************
/
/  GetResistivity()
/
/  written by: Duncan Christie
/  date:       April, 2016
/  PURPOSE: Returns the ambipolar diffusion resistivity coefficient.
/
/  NOTE:  The value returned includes all the factors of clight and 4*Pi
/  absorbed into the coefficient.  Some papers (see any papers from
/  the Mouschovias group, for example) keep those constants out of the
/  resistivity.  This makes more sense aesthetically, but it is preferential
/  to not have all those factors of clight floating around the code.
/ 
/            
/  INPUT:
/    rho - density
/    T - temperature
/    B - magnetic field
/    Time - time
/  RETURNS:
/    
/    The resistivity (in code units)
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
#include "phys_constants.h"

#define SMALL_LOG_VALUE -99.0

/**************************** Functions Prototypes ******************************/

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, float Time);

float ResistivityTableLookup(float D, float T, float B);

int GetConductivity(float B, float T, float Time, float rho_HI, float rho_H2I, float rho_HeI,
                    float rho_E, float rho_HII, float rho_HeII, float rho_H3II,
                    float rho_MII, float rho_AII, float cond[3])
{
  // momentum exchange timescale / conductivity of e-, H+, He+, H3+, M+, and A+
  // reference: Christie+2017 Table 1.
  float tau_s[6], sigma_s[6];
  // temporary variables store the exchange coefficient/timescale of x species and H, He, H2
  float coeff_sH, coeff_sHe, coeff_sH2;
  float tau_sH, tau_sHe, tau_sH2;

  // polarizabilities of H, He,  H2
  float p_He = 0.207, p_H = 0.667, p_H2 = 0.804;

  float sqrtT = sqrt(T);
  float theta = log(T);
  // printf("Tmperature = %13.7e\n", T);

  // collisions with e-
  coeff_sH = 1e-9 * sqrtT * (2.841 + 0.093 * theta + 0.245 * pow(theta, 2) - 0.089 * pow(theta, 3));
  coeff_sHe = 1e-9 * sqrtT * 0.428;
  coeff_sH2 = 1e-9 * sqrtT * (0.535 + 0.203 * theta - 0.163 * pow(theta, 2) + 0.050 * pow(theta, 3));
  // momentum exchange
  tau_sH = (me + mh) / (rho_HI * coeff_sH);
  tau_sHe = (me + 4.0 * mh) / (rho_HeI * coeff_sH);
  tau_sH2 = (me + 2.0 * mh) / (rho_H2I * coeff_sH);
  tau_s[0] = 1.0 / (1.0 / tau_sH + 1.0 / tau_sHe + 1.0 / tau_sH2);

  // collisions with H+
  coeff_sH = 1e-9 * 0.649 * pow(T, 0.375);
  coeff_sHe = 1e-9 * (1.424 + 7.438e-6 * T - 6.734e-9 * pow(T, 2));
  coeff_sH2 = 1e-9 * (1.003 + 0.05 * theta + 0.136 * pow(theta, 2) - 0.014 * pow(theta, 3));
  tau_sH = (mh + mh) / (rho_HI * coeff_sH);
  tau_sHe = (mh + 4.0 * mh) / (rho_HeI * coeff_sH);
  tau_sH2 = (mh + 2.0 * mh) / (rho_H2I * coeff_sH);
  tau_s[1] = 1.0 / (1.0 / tau_sH + 1.0 / tau_sHe + 1.0 / tau_sH2);

  // collisions with He+
  coeff_sH = 1e-9 * 4.71e-1;
  coeff_sHe = 1e-9 * 2.0 * 8.73e-2 * pow(1.0 - 0.093 * theta, 2);
  // Langevin approximation H -> H2
  coeff_sH2 = sqrt((4.0 + 2.0) * 1.0 * p_H2 / (4.0 + 1.0) * 2.0 * p_H) * tau_sH;
  tau_sH = (4.0 * mh + mh) / (rho_HI * coeff_sH);
  tau_sHe = (4.0 * mh + 4.0 * mh) / (rho_HeI * coeff_sH);
  tau_sH2 = (4.0 * mh + 2.0 * mh) / (rho_H2I * coeff_sH);
  tau_s[2] = 1.0 / (1.0 / tau_sH + 1.0 / tau_sHe + 1.0 / tau_sH2);

  // collisions with H3+
  coeff_sH2 = 1e-9 * (2.693 - 1.238 * theta + 0.663 * pow(theta, 2) - 0.089 * pow(theta, 3));
  // Langevin approximation H2 -> H
  coeff_sH = sqrt((3.0 + 1.0) * 2.0 * p_H / (3.0 + 2.0) * 1.0 * p_H2) * tau_sH2;
  // Langevin approximation H2 -> He
  coeff_sHe = sqrt((3.0 + 4.0) * 2.0 * p_He / (3.0 + 2.0) * 4.0 * p_H2) * tau_sH2;
  tau_sH = (3.0 * mh + mh) / (rho_HI * coeff_sH);
  tau_sHe = (3.0 * mh + 4.0 * mh) / (rho_HeI * coeff_sH);
  tau_sH2 = (3.0 * mh + 2.0 * mh) / (rho_H2I * coeff_sH);
  tau_s[3] = 1.0 / (1.0 / tau_sH + 1.0 / tau_sHe + 1.0 / tau_sH2);

  // collisions with molecular ions (approximation from HCO+)
  coeff_sH2 = 1e-9 * sqrtT * (1.476 - 1.409 * theta + 0.555 * pow(theta, 2) - 0.0775 * pow(theta, 3));
  // Langevin approximation H2 -> H
  coeff_sH = sqrt((29.0 + 1.0) * 2.0 * p_H / (29.0 + 2.0) * 1.0 * p_H2) * tau_sH2;
  // Langevin approximation H2 -> He
  coeff_sHe = sqrt((29.0 + 4.0) * 2.0 * p_He / (29.0 + 2.0) * 4.0 * p_H2) * tau_sH2;
  tau_sH = (29.0 * mh + mh) / (rho_HI * coeff_sH);
  tau_sHe = (29,0 * mh + 4.0 * mh) / (rho_HeI * coeff_sH);
  tau_sH2 = (29.0 * mh + 2.0 * mh) / (rho_H2I * coeff_sH);
  tau_s[4] = 1.0 / (1.0 / tau_sH + 1.0 / tau_sHe + 1.0 / tau_sH2);

  // collisions with atomic ions (approximation from C+)
  coeff_sH = 1e-9 * (1.983 + 0.425 * theta - 0.431 * pow(theta, 2) + 0.114 * pow(theta, 3));
  // Langevin approximation H -> H2
  coeff_sH2 = sqrt((12.0 + 2.0) * 1.0 * p_H2 / (12.0 + 1.0) * 2.0 * p_H) * tau_sH;
  // Langevin approximation H -> He
  coeff_sHe = sqrt((12.0 + 4.0) * 1.0 * p_He / (12.0 + 1.0) * 4.0 * p_H) * tau_sH;
  tau_sH = (12.0 * mh + mh) / (rho_HI * coeff_sH);
  tau_sHe = (12.0 * mh + 4.0 * mh) / (rho_HeI * coeff_sH);
  tau_sH2 = (12.0 * mh + 2.0 * mh) / (rho_H2I * coeff_sH);
  tau_s[5] = 1.0 / (1.0 / tau_sH + 1.0 / tau_sHe + 1.0 / tau_sH2);

  // printf("tau e: %13.7e, H+: %13.7e, He+: %13.7e, H3+: %13.7e, M+: %13.7e, A+: %13.7e\n", 
  //         tau_s[0], tau_s[1], tau_s[2], tau_s[3], tau_s[4], tau_s[5]);

  const float ucharge = 4.80320425e-10; // statcoulombs = cm^(3/2) g^(1/2) s^-1

  sigma_s[0] = (rho_E / mh) * tau_s[0] * ucharge * ucharge / me;                     // conductivity of e-, n_e = rho_e/mh
  sigma_s[1] = (rho_HII / mh) * tau_s[1] * ucharge * ucharge / mh;                   // H+
  sigma_s[2] = (rho_HeII / (4.0 * mh)) * tau_s[2] * ucharge * ucharge / (4.0 * mh);  // He+
  sigma_s[3] = (rho_H3II / (3.0 * mh)) * tau_s[3] * ucharge * ucharge / (3.0 * mh);  // H3+
  sigma_s[4] = (rho_MII / (29.0 * mh)) * tau_s[4] * ucharge * ucharge / (29.0 * mh); // M+ := HCO+
  sigma_s[5] = (rho_AII / (12.0 * mh)) * tau_s[5] * ucharge * ucharge / (12.0 * mh); // A+ := C+

  // printf("sigma e: %13.7e, H+: %13.7e, He+: %13.7e, H3+: %13.7e, M+: %13.7e, A+: %13.7e\n", 
  //         sigma_s[0], sigma_s[1], sigma_s[2], sigma_s[3], sigma_s[4], sigma_s[5]);

  float omega_s[6];
  omega_s[0] = ucharge * B / (me * clight);
  omega_s[1] = ucharge * B / (mh * clight);
  omega_s[2] = ucharge * B / (4.0 * mh * clight);
  omega_s[3] = ucharge * B / (3.0 * mh * clight);
  omega_s[4] = ucharge * B / (29.0 * mh * clight);
  omega_s[5] = ucharge * B / (12.0 * mh * clight);

  // printf("omega e: %13.7e, H+: %13.7e, He+: %13.7e, H3+: %13.7e, M+: %13.7e, A+: %13.7e\n", 
  //         omega_s[0], omega_s[1], omega_s[2], omega_s[3], omega_s[4], omega_s[5]);

  cond[0] = cond[1] = cond[2] = 0.0;
  for (int i = 0; i < 6; i++)
  {
    cond[0] += sigma_s[i];
    cond[1] += sigma_s[i] / (1.0 + pow(omega_s[i] * tau_s[i], 2));
    cond[2] -= sigma_s[i] * omega_s[i] * tau_s[i] / (1.0 + pow(omega_s[i] * tau_s[i], 2));
  }

  return SUCCESS;
}

// float GetResistivity(float rho, float B, float T, float Time)
float GetResistivity(float rho, float B, float T, float Time,
                     // neutral species
                     float rho_HI, float rho_H2I, float rho_HeI,
                     // ion species
                     float rho_E, float rho_HII, float rho_HeII, float rho_H3II,
                     float rho_MII, float rho_AII)
{
  float eta_perp = 0.;
  float rhoi = 0.;
  float D = 0.;
  float alpha;
  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1,
        TimeUnits = 1.0, VelocityUnits = 1.0;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
           &TimeUnits, &VelocityUnits, Time);
  float MagneticUnits = sqrt(DensityUnits * 4.0 * M_PI) * VelocityUnits;

  float clight_code = clight / VelocityUnits;

  rho *= DensityUnits;
  B *= MagneticUnits;

  rho_HI *= DensityUnits;
  rho_H2I *= DensityUnits;
  rho_HeI *= DensityUnits;
  rho_E *= DensityUnits;
  rho_HII *= DensityUnits;
  rho_HeII *= DensityUnits;
  rho_H3II *= DensityUnits;
  rho_MII *= DensityUnits;
  rho_AII *= DensityUnits;

  //T *= TemperatureUnits;

  if (debug)
  {
    if (T != T)
    {
      fprintf(stderr, "GetResistivity(): T is NaN!\n");
      fprintf(stderr, "GetResistivity(): D = %e, T = %e, B= %e\n", D, T, B);
    }
  }

  switch (ADResistivityType)
  {
  case 1:
    return ADResistivityScale * TimeUnits / (LengthUnits * LengthUnits);
  case 2:
    eta_perp = 0.;
    eta_perp = ADResistivityScale * B * B / (4. * pi * pow(rho, 1.5));
    eta_perp *= TimeUnits / (LengthUnits * LengthUnits);
    return eta_perp;
  case 3:
    eta_perp = clight_code * clight_code * ResistivityTableLookup(rho, T, B) / (4. * pi);
    eta_perp /= TimeUnits;
    // printf("ADResistivityType = 3, eta = %13.7e\n", eta_perp);
    return eta_perp;
  case 4:
    float ni, nn;
    nn = rho / (2.33 * 1.6733e-24);
    ni = 1e-3 * pow(nn / 1e5, 0.5) + 1e-3 * pow(nn / 1e3, -2);
    rhoi = 29.0 * 1.6733e-24 * ni;
    alpha = 3.7e13;
    eta_perp = B * B / (4.0 * pi * alpha * rhoi * rho);
    eta_perp *= TimeUnits / (LengthUnits * LengthUnits);
    // printf("ADResistivityType = 4, eta = %13.7e\n", eta_perp);
    return eta_perp;
  case 5:
    float cond[3], eta[3], eta_AD;
    GetConductivity(B, T, Time, rho_HI, rho_H2I, rho_HeI,
                    rho_E, rho_HII, rho_HeII, rho_H3II,
                    rho_MII, rho_AII, &cond[0]);

    eta[0] = 1.0 / cond[0];                                     // eta_parallel
    eta[1] = cond[1] / (cond[1] * cond[1] + cond[2] * cond[2]); // eta_perpendicular
    // eta[2] = cond[2] / (cond[1] * cond[1] + cond[2] * cond[2]); // eta_H
    eta_AD = eta[1] - eta[0];
    eta_AD *= clight * clight / (4.0 * pi);
    eta_AD *= TimeUnits / (LengthUnits * LengthUnits);
    // printf("ADResistivityType = 5, cond_para = %13.7e, cond_perp = %13.7e, cond_H = %13.7e\n", cond[0], cond[1], cond[2]);
    // printf("ADResistivityType = 5, eta_para = %13.7e, eta_perp = %13.7e, eta = %13.7e\n", eta[0], eta[1], eta_AD);
    return eta_AD;
  default:
    return 0.;
  }
}
