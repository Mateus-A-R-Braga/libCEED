#include "../include/post-processing.h"
#include "../include/setup-solvers.h"
// -----------------------------------------------------------------------------
// This function print the output
// -----------------------------------------------------------------------------
PetscErrorCode PrintOutput(Ceed ceed, AppCtx app_ctx, PetscBool has_ts,
                           CeedMemType mem_type_backend,
                           TS ts, SNES snes, KSP ksp,
                           Vec U, CeedScalar l2_error_u, CeedScalar l2_error_p) {

  PetscFunctionBeginUser;

  const char *used_resource;
  CeedGetResource(ceed, &used_resource);
  char hostname[PETSC_MAX_PATH_LEN];
  PetscCall( PetscGetHostName(hostname, sizeof hostname) );
  PetscInt comm_size;
  PetscCall( MPI_Comm_size(app_ctx->comm, &comm_size) );
  PetscCall( PetscPrintf(app_ctx->comm,
                         "\n-- Mixed H(div) Example - libCEED + PETSc --\n"
                         "  MPI:\n"
                         "    Hostname                           : %s\n"
                         "    Total ranks                        : %d\n"
                         "  libCEED:\n"
                         "    libCEED Backend                    : %s\n"
                         "    libCEED Backend MemType            : %s\n",
                         hostname, comm_size, used_resource, CeedMemTypes[mem_type_backend]) );

  VecType vecType;
  PetscCall( VecGetType(U, &vecType) );
  PetscCall( PetscPrintf(app_ctx->comm,
                         "  PETSc:\n"
                         "    PETSc Vec Type                     : %s\n",
                         vecType) );

  PetscInt       U_l_size, U_g_size;
  PetscCall( VecGetSize(U, &U_g_size) );
  PetscCall( VecGetLocalSize(U, &U_l_size) );
  PetscCall( PetscPrintf(app_ctx->comm,
                         "  Problem:\n"
                         "    Problem Name                       : %s\n"
                         "    Global nodes (u + p)               : %" PetscInt_FMT "\n"
                         "    Owned nodes (u + p)                : %" PetscInt_FMT "\n",
                         app_ctx->problem_name, U_g_size, U_l_size
                        ) );
  // --TS
  if (has_ts) {
    PetscInt ts_steps;
    TSType ts_type;
    TSConvergedReason ts_reason;
    PetscCall( TSGetStepNumber(ts, &ts_steps) );
    PetscCall( TSGetType(ts, &ts_type) );
    PetscCall( TSGetConvergedReason(ts, &ts_reason) );
    PetscCall( PetscPrintf(app_ctx->comm,
                           "  TS:\n"
                           "    TS Type                            : %s\n"
                           "    TS Convergence                     : %s\n"
                           "    Number of TS steps                 : %" PetscInt_FMT "\n"
                           "    Final time                         : %g\n",
                           ts_type, TSConvergedReasons[ts_reason],
                           ts_steps, (double)app_ctx->t_final) );

    PetscCall( TSGetSNES(ts, &snes) );
  }
  // -- SNES
  PetscInt its, snes_its = 0;
  PetscCall( SNESGetIterationNumber(snes, &its) );
  snes_its += its;
  SNESType snes_type;
  SNESConvergedReason snes_reason;
  PetscReal snes_rnorm;
  PetscCall( SNESGetType(snes, &snes_type) );
  PetscCall( SNESGetConvergedReason(snes, &snes_reason) );
  PetscCall( SNESGetFunctionNorm(snes, &snes_rnorm) );
  PetscCall( PetscPrintf(app_ctx->comm,
                         "  SNES:\n"
                         "    SNES Type                          : %s\n"
                         "    SNES Convergence                   : %s\n"
                         "    Total SNES Iterations              : %" PetscInt_FMT "\n"
                         "    Final rnorm                        : %e\n",
                         snes_type, SNESConvergedReasons[snes_reason],
                         snes_its, (double)snes_rnorm) );
  if (!has_ts) {
    PetscInt ksp_its = 0;
    PetscCall( SNESGetLinearSolveIterations(snes, &its) );
    ksp_its += its;
    KSPType ksp_type;
    KSPConvergedReason ksp_reason;
    PetscReal ksp_rnorm;
    PC pc;
    PCType pc_type;
    PetscCall( KSPGetPC(ksp, &pc) );
    PetscCall( PCGetType(pc, &pc_type) );
    PetscCall( KSPGetType(ksp, &ksp_type) );
    PetscCall( KSPGetConvergedReason(ksp, &ksp_reason) );
    PetscCall( KSPGetIterationNumber(ksp, &ksp_its) );
    PetscCall( KSPGetResidualNorm(ksp, &ksp_rnorm) );
    PetscCall( PetscPrintf(app_ctx->comm,
                           "  KSP:\n"
                           "    KSP Type                           : %s\n"
                           "    PC Type                            : %s\n"
                           "    KSP Convergence                    : %s\n"
                           "    Total KSP Iterations               : %" PetscInt_FMT "\n"
                           "    Final rnorm                        : %e\n",
                           ksp_type, pc_type, KSPConvergedReasons[ksp_reason], ksp_its,
                           (double)ksp_rnorm ) );
  }

  PetscCall( PetscPrintf(app_ctx->comm,
                         "  L2 Error (MMS):\n"
                         "    L2 error of u and p                : %e, %e\n",
                         (double)l2_error_u,
                         (double)l2_error_p) );
  PetscFunctionReturn(0);
};

// -----------------------------------------------------------------------------
// Setup operator context data for initial condition, u field
// -----------------------------------------------------------------------------
PetscErrorCode SetupProjectVelocityCtx_Hdiv(MPI_Comm comm, DM dm, Ceed ceed,
    CeedData ceed_data,
    OperatorApplyContext ctx_post_Hdiv) {
  PetscFunctionBeginUser;

  ctx_post_Hdiv->comm = comm;
  ctx_post_Hdiv->dm = dm;
  PetscCall( DMCreateLocalVector(dm, &ctx_post_Hdiv->X_loc) );
  ctx_post_Hdiv->x_ceed = ceed_data->u_ceed;
  //ctx_project_velocity->y_ceed = ceed_data->v0_ceed;
  ctx_post_Hdiv->ceed = ceed;
  //ctx_project_velocity->op_apply = ceed_data->op_ics_u;

  PetscFunctionReturn(0);
}

PetscErrorCode SetupProjectVelocityCtx_H1(MPI_Comm comm, DM dm_H1, Ceed ceed,
    CeedData ceed_data,
    OperatorApplyContext ctx_post_H1) {
  PetscFunctionBeginUser;

  ctx_post_H1->comm = comm;
  ctx_post_H1->dm = dm_H1;
  PetscCall( DMCreateLocalVector(dm_H1, &ctx_post_H1->X_loc) );
  PetscCall( VecDuplicate(ctx_post_H1->X_loc, &ctx_post_H1->Y_loc) );
  ctx_post_H1->x_ceed = ceed_data->up_ceed;
  ctx_post_H1->y_ceed = ceed_data->vp_ceed;
  ctx_post_H1->ceed = ceed;
  ctx_post_H1->op_apply = ceed_data->op_post_mass;

  PetscFunctionReturn(0);
}
// -----------------------------------------------------------------------------
// This function print the output
// -----------------------------------------------------------------------------
PetscErrorCode ProjectVelocity(CeedData ceed_data,
                               Vec U, VecType vec_type, Vec *U_H1,
                               OperatorApplyContext ctx_post_Hdiv,
                               OperatorApplyContext ctx_post_H1) {

  PetscFunctionBeginUser;
  const PetscScalar *x;
  PetscMemType  x_mem_type;

  // ----------------------------------------------
  // Create local rhs for u field
  // ----------------------------------------------
  Vec post_rhs_loc;
  PetscScalar *ru;
  PetscMemType ru_mem_type;
  PetscCall( DMCreateLocalVector(ctx_post_H1->dm, &post_rhs_loc) );
  PetscCall( VecZeroEntries(post_rhs_loc) );
  PetscCall( VecGetArrayAndMemType(post_rhs_loc, &ru, &ru_mem_type) );
  CeedElemRestrictionCreateVector(ceed_data->elem_restr_u_post,
                                  &ceed_data->post_rhs_ceed,
                                  NULL);
  CeedVectorSetArray(ceed_data->post_rhs_ceed, MemTypeP2C(ru_mem_type),
                     CEED_USE_POINTER, ru);


  // Global-to-local: map final U in Hdiv space to local
  PetscCall( DMGlobalToLocal(ctx_post_Hdiv->dm,
                             U, INSERT_VALUES, ctx_post_Hdiv->X_loc) );
  // Place Hdiv PETSc vectors in CEED vectors
  PetscCall( VecGetArrayReadAndMemType(ctx_post_Hdiv->X_loc,
                                       &x, &x_mem_type) );
  CeedVectorSetArray(ctx_post_Hdiv->x_ceed, MemTypeP2C(x_mem_type),
                     CEED_USE_POINTER, (PetscScalar *)x);

  // Apply operator to create RHS for u field
  CeedOperatorApply(ceed_data->op_post_rhs, ceed_data->x_coord,
                    ceed_data->post_rhs_ceed, CEED_REQUEST_IMMEDIATE);

  // Restore vectors Hdiv vector
  CeedVectorTakeArray(ctx_post_Hdiv->x_ceed,
                      MemTypeP2C(x_mem_type), NULL);
  PetscCall( VecRestoreArrayReadAndMemType(ctx_post_Hdiv->X_loc, &x) );

  // ----------------------------------------------
  // Create global rhs for u field
  // ----------------------------------------------
  Vec post_rhs;
  CeedVectorTakeArray(ceed_data->post_rhs_ceed, MemTypeP2C(ru_mem_type), NULL);
  PetscCall( VecRestoreArrayAndMemType(post_rhs_loc, &ru) );
  PetscCall( DMCreateGlobalVector(ctx_post_H1->dm, &post_rhs) );
  PetscCall( VecZeroEntries(post_rhs) );
  PetscCall( DMLocalToGlobal(ctx_post_H1->dm, post_rhs_loc, ADD_VALUES,
                             post_rhs) );

  // ----------------------------------------------
  // Solve for U_H1, M*U_H1 = post_rhs
  // ----------------------------------------------
  PetscInt UH1_g_size, UH1_l_size;
  PetscCall( VecGetSize(*U_H1, &UH1_g_size) );
  // Local size for matShell
  PetscCall( VecGetLocalSize(*U_H1, &UH1_l_size) );

  // Operator
  Mat mat_ksp_projection;
  // -- Form Action of residual on u
  PetscCall( MatCreateShell(ctx_post_H1->comm, UH1_l_size, UH1_l_size, UH1_g_size,
                            UH1_g_size, ceed_data->ctx_post_H1, &mat_ksp_projection) );
  PetscCall( MatShellSetOperation(mat_ksp_projection, MATOP_MULT,
                                  (void (*)(void))ApplyMatOp) );
  PetscCall( MatShellSetVecType(mat_ksp_projection, vec_type) );

  KSP ksp_projection;
  PetscCall( KSPCreate(ctx_post_H1->comm, &ksp_projection) );
  PetscCall( KSPSetOperators(ksp_projection, mat_ksp_projection,
                             mat_ksp_projection) );
  PetscCall( KSPSetFromOptions(ksp_projection) );
  PetscCall( KSPSetUp(ksp_projection) );
  PetscCall( KSPSolve(ksp_projection, post_rhs, *U_H1) );

  // Clean up
  PetscCall( VecDestroy(&post_rhs_loc) );
  PetscCall( VecDestroy(&post_rhs) );
  PetscCall( VecDestroy(&ctx_post_H1->X_loc) );
  PetscCall( VecDestroy(&ctx_post_H1->Y_loc) );
  PetscCall( VecDestroy(&ctx_post_Hdiv->X_loc) );
  PetscCall( MatDestroy(&mat_ksp_projection) );
  PetscCall( KSPDestroy(&ksp_projection) );

  PetscFunctionReturn(0);
};
// -----------------------------------------------------------------------------
