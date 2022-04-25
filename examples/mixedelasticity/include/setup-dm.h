#ifndef setup_dm_h
#define setup_dm_h

#include <petsc.h>
#include <petscdmplex.h>
#include <petscsys.h>
#include <ceed.h>
#include "structs.h"

// ---------------------------------------------------------------------------
// Set-up DM
// ---------------------------------------------------------------------------
PetscErrorCode CreateDM(MPI_Comm comm, VecType vec_type, DM *dm);

#endif // setup_dm_h
