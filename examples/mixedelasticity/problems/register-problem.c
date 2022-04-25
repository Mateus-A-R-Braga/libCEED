// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
// reserved. See files LICENSE and NOTICE for details.
//
// This file is part of CEED, a collection of benchmarks, miniapps, software
// libraries and APIs for efficient high-order finite element and spectral
// element discretizations for exascale applications. For more information and
// source code availability see http://github.com/ceed.
//
// The CEED research is supported by the Exascale Computing Project 17-SC-20-SC,
// a collaborative effort of two U.S. Department of Energy organizations (Office
// of Science and the National Nuclear Security Administration) responsible for
// the planning and preparation of a capable exascale ecosystem, including
// software, applications, hardware, advanced system engineering and early
// testbed platforms, in support of the nation's exascale computing imperative.

/// @file
/// Command line option processing for mixed-elasticity example using PETSc

#include "../include/register-problem.h"

// Register problems to be available on the command line
PetscErrorCode RegisterProblems_MixedElasticity(AppCtx app_ctx) {
  app_ctx->problems = NULL;
  PetscFunctionBeginUser;
  // 1) linear2d (MixedElasticity_LINEAR2D is created in linear2d.c)
  PetscCall( PetscFunctionListAdd(&app_ctx->problems, "linear2d",
                                  MixedElasticity_LINEAR2D) );
  // 2) linear3d (MixedElasticity_LINEAR3D is created in linear3d.c)
  PetscCall( PetscFunctionListAdd(&app_ctx->problems, "linear3d",
                                  MixedElasticity_LINEAR3D) );
  // 3) ...

  // 4) ...

  PetscFunctionReturn(0);
}