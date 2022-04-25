#ifndef register_problem_h
#define register_problem_h

#include "structs.h"

// Register problems to be available on the command line
PetscErrorCode RegisterProblems_MixedElasticity(AppCtx app_ctx);
// -----------------------------------------------------------------------------
// Set up problems function prototype
// -----------------------------------------------------------------------------
// 1) linear2d
PetscErrorCode MixedElasticity_LINEAR2D(Ceed ceed, ProblemData problem_data,
                                        void *ctx);

// 2) linear3d
PetscErrorCode MixedElasticity_LINEAR3D(Ceed ceed, ProblemData problem_data,
                                        void *ctx);

// 3) ...

// 4) ...

#endif // register_problem_h
