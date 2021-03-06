/***********************************************************************
/
/  GRID CLASS (HANDLE CALLING AND SOLVING COOLING/CHEMISTRY)
/
/  written by: Matthew Turk
/  date:       June, 2009
/  modified1:
/
/  PURPOSE: Move logic for chemistry/cooling module selection here
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include "preincludes.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::MultiSpeciesHandler()
{
  if ((!MultiSpecies) && (!RadiativeCooling)) return SUCCESS; 
  if (GadgetEquilibriumCooling != 0) return SUCCESS;

  LCAPERF_START("grid_MultiSpeciesHandler");

#ifdef USE_GRACKLE
  if (grackle_data->use_grackle == TRUE) {
    grackle_data->radiative_transfer_intermediate_step = FALSE;
    if (this->GrackleWrapper() == FAIL) {
      ENZO_FAIL("Error in GrackleWrapper.\n");
    }
#ifdef USE_KROME
    // In case of using krome and grackle at the same time.
    if (!use_krome) return SUCCESS;
#else
    return SUCCESS;
#endif
  }
#endif

//#ifdef USE_KROME
//  if (use_krome) {
//    KromeSolver();
//  }
//  else if (MultiSpecies && RadiativeCooling ) {
//#else
//  if (MultiSpecies && RadiativeCooling ) {
//#endif
#ifdef USE_KROME
  if ((MultiSpecies && RadiativeCooling) || use_krome ) {
    if (use_krome && use_kromestep == 2) return SUCCESS;
#else
  if (MultiSpecies && RadiativeCooling) {
#endif
    int RTCoupledSolverIntermediateStep = FALSE;
    this->SolveRateAndCoolEquations(RTCoupledSolverIntermediateStep);
  } 
  else {
    if (MultiSpecies)
      this->SolveRateEquations();
    if (RadiativeCooling)
      this->SolveRadiativeCooling();
  }

  if (ProblemType == 62)
    this->CoolingTestResetEnergies();

  LCAPERF_STOP("grid_MultiSpeciesHandler");
  return SUCCESS;
}
