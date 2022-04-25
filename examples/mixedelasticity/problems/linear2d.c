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
/// Utility functions for setting up mixed-elasticity problem in 2D

#include "../include/register-problem.h"
#include "../qfunctions/volumetric-geometry2d.h"
#include "../qfunctions/linear-force2d.h"
#include "../qfunctions/linear-system2d.h"
#include "../qfunctions/linear-error2d.h"
#include "../qfunctions/pressure-boundary2d.h"

PetscErrorCode MixedElasticity_LINEAR2D(Ceed ceed, ProblemData problem_data,
                                        void *ctx) {
  AppCtx               app_ctx = *(AppCtx *)ctx;
  PetscFunctionBeginUser;

  // ------------------------------------------------------
  //               SET UP POISSON_QUAD2D
  // ------------------------------------------------------
  problem_data->dim                     = 2;
  problem_data->q_data_size             = 5;
  problem_data->q_data_size_face        = 3;
  problem_data->quadrature_mode         = CEED_GAUSS;
  problem_data->setup_geo               = SetupVolumeGeometry2D;
  problem_data->setup_geo_loc           = SetupVolumeGeometry2D_loc;
  problem_data->force                   = LinearForce2D;
  problem_data->force_loc               = LinearForce2D_loc;
  problem_data->residual                = LinearSystem2D;
  problem_data->residual_loc            = LinearSystem2D_loc;
  problem_data->jacobian                = JacobianLinearSystem2D;
  problem_data->jacobian_loc            = JacobianLinearSystem2D_loc;
  problem_data->error                   = LinearError2D;
  problem_data->error_loc               = LinearError2D_loc;
  problem_data->bc_pressure             = BCPressure2D;
  problem_data->bc_pressure_loc         = BCPressure2D_loc;

  // ------------------------------------------------------
  //              Command line Options
  // ------------------------------------------------------
  PetscOptionsBegin(app_ctx->comm, NULL, "Options for mixed-elasticity problem",
                    NULL);

  PetscOptionsEnd();
  PetscFunctionReturn(0);
}
