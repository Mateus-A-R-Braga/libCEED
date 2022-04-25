#ifndef setup_fe_h
#define setup_fe_h

#include <petsc.h>
#include <petscdmplex.h>
#include <petscsys.h>
#include <ceed.h>
#include "structs.h"

// ---------------------------------------------------------------------------
// Set-up FE
// ---------------------------------------------------------------------------
PetscErrorCode SetupFE(ProblemData problem_data,
                       AppCtx app_ctx, DM dm);

#endif // setup_fe_h
