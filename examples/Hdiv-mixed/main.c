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

//                        libCEED + PETSc Example: Mixed-Poisson in H(div) space
//
// This example demonstrates a simple usage of libCEED with PETSc to solve
//   elasticity problems.
//
// The code uses higher level communication protocols in DMPlex.
//
// Build with: make
// Run with:
//   ./main -pc_type svd -problem darcy2d -dm_plex_dim 2 -dm_plex_box_faces 4,4
//   ./main -pc_type svd -problem darcy3d -dm_plex_dim 3 -dm_plex_box_faces 4,4,4
//   ./main -pc_type svd -problem darcy3d -dm_plex_filename /path to the mesh file
//   ./main -pc_type svd -problem darcy2d -dm_plex_dim 2 -dm_plex_box_faces 4,4 -bc_pressure 1
//   ./main -pc_type svd -problem darcy2d -dm_plex_dim 2 -dm_plex_box_faces 4,4 -bc_pressure 1,2,3,4
//   ./main -pc_type svd -problem richard2d -dm_plex_dim 2 -dm_plex_box_faces 4,4

#include "ceed/ceed.h"
#include <stdio.h>
const char help[] = "Solve H(div)-mixed problem using PETSc and libCEED\n";

#include "main.h"

int main(int argc, char **argv) {
  // ---------------------------------------------------------------------------
  // Initialize PETSc
  // ---------------------------------------------------------------------------
  PetscCall( PetscInitialize(&argc, &argv, NULL, help) );

  // ---------------------------------------------------------------------------
  // Initialize libCEED
  // ---------------------------------------------------------------------------
  // -- Initialize backend
  Ceed ceed;
  CeedInit("/cpu/self/ref/serial", &ceed);
  CeedMemType mem_type_backend;
  CeedGetPreferredMemType(ceed, &mem_type_backend);

  VecType        vec_type = NULL;
  switch (mem_type_backend) {
  case CEED_MEM_HOST: vec_type = VECSTANDARD; break;
  case CEED_MEM_DEVICE: {
    const char *resolved;
    CeedGetResource(ceed, &resolved);
    if (strstr(resolved, "/gpu/cuda")) vec_type = VECCUDA;
    else if (strstr(resolved, "/gpu/hip/occa"))
      vec_type = VECSTANDARD; // https://github.com/CEED/libCEED/issues/678
    else if (strstr(resolved, "/gpu/hip")) vec_type = VECHIP;
    else vec_type = VECSTANDARD;
  }
  }

  // ---------------------------------------------------------------------------
  // Create structs
  // ---------------------------------------------------------------------------
  AppCtx app_ctx;
  PetscCall( PetscCalloc1(1, &app_ctx) );

  ProblemData problem_data = NULL;
  PetscCall( PetscCalloc1(1, &problem_data) );

  CeedData ceed_data;
  PetscCall( PetscCalloc1(1, &ceed_data) );

  Physics phys_ctx;
  PetscCall( PetscCalloc1(1, &phys_ctx) );

  OperatorApplyContext ctx_residual_ut, ctx_initial;
  PetscCall( PetscCalloc1(1, &ctx_residual_ut) );
  PetscCall( PetscCalloc1(1, &ctx_initial) );

  // ---------------------------------------------------------------------------
  // Process command line options
  // ---------------------------------------------------------------------------
  // -- Register problems to be available on the command line
  PetscCall( RegisterProblems_Hdiv(app_ctx) );

  // -- Process general command line options
  MPI_Comm comm = PETSC_COMM_WORLD;
  app_ctx->comm = comm;
  PetscCall( ProcessCommandLineOptions(app_ctx) );

  // ---------------------------------------------------------------------------
  // Choose the problem from the list of registered problems
  // ---------------------------------------------------------------------------
  {
    PetscErrorCode (*p)(Ceed, ProblemData, void *);
    PetscCall( PetscFunctionListFind(app_ctx->problems, app_ctx->problem_name,
                                     &p) );
    if (!p) SETERRQ(PETSC_COMM_SELF, 1, "Problem '%s' not found",
                      app_ctx->problem_name);
    PetscCall( (*p)(ceed, problem_data, &app_ctx) );
  }

  // ---------------------------------------------------------------------------
  // Create DM
  // ---------------------------------------------------------------------------
  DM             dm;
  PetscCall( CreateDM(comm, vec_type, &dm) );
  // TODO: add mesh option
  // perturb to have smooth random mesh
  // PetscCall( PerturbVerticesSmooth(dm) );

  // ---------------------------------------------------------------------------
  // Setup FE
  // ---------------------------------------------------------------------------
  SetupFE(comm, dm);

  // ---------------------------------------------------------------------------
  // Setup libCEED
  // ---------------------------------------------------------------------------
  // -- Set up libCEED objects
  PetscCall( SetupLibceed(dm, ceed, app_ctx, ctx_residual_ut,
                          problem_data, ceed_data) );
  //CeedVectorView(force_ceed, "%12.8f", stdout);
  //PetscCall( DMAddBoundariesPressure(ceed, ceed_data, app_ctx, problem_data, dm,
  //                                   bc_pressure) );


  // ---------------------------------------------------------------------------
  // Create Global Solution
  // ---------------------------------------------------------------------------
  Vec U; // U = [p,u]
  PetscCall( DMCreateGlobalVector(dm, &U) );

  // ---------------------------------------------------------------------------
  // Setup TSSolve for Richard problem
  // ---------------------------------------------------------------------------
  TS ts;
  SNES snes;
  KSP ksp;
  if (problem_data->has_ts) {
    // ---------------------------------------------------------------------------
    // Create global initial conditions
    // ---------------------------------------------------------------------------

    SetupResidualOperatorCtx_U0(dm, ceed, ceed_data, ctx_initial);
    CreateInitialConditions(ceed_data, U, ctx_initial);
    VecView(U, PETSC_VIEWER_STDOUT_WORLD);
    SetupResidualOperatorCtx_Ut(dm, ceed, ceed_data, ctx_residual_ut);
    PetscCall( VecZeroEntries(ctx_residual_ut->X_t_loc) );
    PetscCall( TSSolveRichard(dm, ceed, ceed_data, app_ctx, ctx_residual_ut,
                              &U, &ts) );
  }

  if (!problem_data->has_ts) {
    // ---------------------------------------------------------------------------
    // Solve PDE
    // ---------------------------------------------------------------------------
    // Create SNES
    PetscCall( SNESCreate(comm, &snes) );
    PetscCall( SNESGetKSP(snes, &ksp) );
    PetscCall( PDESolver(comm, dm, ceed, ceed_data, vec_type, snes, ksp, &U) );
    //VecView(U, PETSC_VIEWER_STDOUT_WORLD);
  }

  // ---------------------------------------------------------------------------
  // Compute L2 error of mms problem
  // ---------------------------------------------------------------------------
  CeedScalar l2_error_u, l2_error_p;
  PetscCall( ComputeL2Error(dm, ceed,ceed_data, U, &l2_error_u,
                            &l2_error_p) );
  // ---------------------------------------------------------------------------
  // Print output results
  // ---------------------------------------------------------------------------
  PetscCall( PrintOutput(ceed, mem_type_backend, ts,
                         snes, ksp, U, l2_error_u, l2_error_p, app_ctx, problem_data->has_ts) );
  // ---------------------------------------------------------------------------
  // Save solution (paraview)
  // ---------------------------------------------------------------------------
  PetscViewer viewer;

  PetscCall( PetscViewerVTKOpen(comm,"solution.vtu",FILE_MODE_WRITE,&viewer) );
  PetscCall( VecView(U, viewer) );
  PetscCall( PetscViewerDestroy(&viewer) );

  // ---------------------------------------------------------------------------
  // Free objects
  // ---------------------------------------------------------------------------

  // Free PETSc objects
  PetscCall( DMDestroy(&dm) );
  PetscCall( VecDestroy(&U) );
  if (problem_data->has_ts) {
    PetscCall( TSDestroy(&ts) );
  } else {
    PetscCall( SNESDestroy(&snes) );
  }

  // -- Function list
  PetscCall( PetscFunctionListDestroy(&app_ctx->problems) );

  // Free libCEED objects
  //CeedVectorDestroy(&bc_pressure);
  PetscCall( CeedDataDestroy(ceed_data, problem_data) );
  // -- Structs
  PetscCall( PetscFree(app_ctx) );
  PetscCall( PetscFree(problem_data) );
  PetscCall( PetscFree(phys_ctx) );
  PetscCall( PetscFree(ctx_residual_ut) );
  PetscCall( PetscFree(ctx_initial) );
  CeedDestroy(&ceed);

  return PetscFinalize();
}
