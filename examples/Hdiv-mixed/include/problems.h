#ifndef problems_h
#define problems_h

#include "../include/structs.h"

// -----------------------------------------------------------------------------
// Set up problems function prototype
// -----------------------------------------------------------------------------
// 1) darcy2d
PetscErrorCode Hdiv_DARCY2D(ProblemData *problem_data, void *ctx);

// 2) darcy3d
PetscErrorCode Hdiv_DARCY3D(ProblemData *problem_data, void *ctx);

// 3) darcy3dprism

// 4) richard

#endif // problems_h
