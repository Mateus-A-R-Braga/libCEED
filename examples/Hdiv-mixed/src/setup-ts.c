#include "../include/setup-ts.h"
#include "../include/setup-matops.h"
#include "../include/setup-libceed.h"
#include "ceed/ceed.h"
#include "petscerror.h"
#include <stdio.h>


// -----------------------------------------------------------------------------
// Setup operator context data for Richard problem
// -----------------------------------------------------------------------------
PetscErrorCode SetupResidualOperatorCtx_Ut(DM dm, Ceed ceed, CeedData ceed_data,
    OperatorApplyContext ctx_residual_ut) {
  PetscFunctionBeginUser;

  ctx_residual_ut->dm = dm;
  PetscCall( DMCreateLocalVector(dm, &ctx_residual_ut->X_loc) );
  PetscCall( DMCreateLocalVector(dm, &ctx_residual_ut->X_t_loc) );
  PetscCall( DMCreateLocalVector(dm, &ctx_residual_ut->Y_loc) );
  ctx_residual_ut->x_ceed = ceed_data->x_ceed;
  ctx_residual_ut->x_t_ceed = ceed_data->x_t_ceed;
  ctx_residual_ut->y_ceed = ceed_data->y_ceed;
  ctx_residual_ut->ceed = ceed;
  ctx_residual_ut->op_apply = ceed_data->op_residual;
  ctx_residual_ut->op_true = ceed_data->op_true;

  PetscFunctionReturn(0);
}

PetscErrorCode SetupResidualOperatorCtx_U0(DM dm, Ceed ceed, CeedData ceed_data,
    OperatorApplyContext ctx_initial) {
  PetscFunctionBeginUser;

  ctx_initial->dm = dm;
  PetscCall( DMCreateLocalVector(dm, &ctx_initial->X_loc) );
  PetscCall( DMCreateLocalVector(dm, &ctx_initial->X_t_loc) );
  PetscCall( DMCreateLocalVector(dm, &ctx_initial->Y_loc) );
  ctx_initial->x_ceed = ceed_data->x_ceed;
  ctx_initial->x_t_ceed = ceed_data->x_t_ceed;
  ctx_initial->y_ceed = ceed_data->y_ceed;
  ctx_initial->ceed = ceed;
  ctx_initial->op_apply = ceed_data->op_ics;

  PetscFunctionReturn(0);
}
// -----------------------------------------------------------------------------
// Create global initial conditions vector
// -----------------------------------------------------------------------------
PetscErrorCode CreateInitialConditions(CeedData ceed_data, Vec U0,
                                       OperatorApplyContext ctx_initial) {
  PetscScalar *u0;
  PetscMemType u0_mem_type;

  PetscFunctionBeginUser;

  PetscCall( VecGetArrayAndMemType(ctx_initial->Y_loc, &u0, &u0_mem_type) );
  CeedVectorSetArray(ctx_initial->y_ceed, MemTypeP2C(u0_mem_type),
                     CEED_USE_POINTER, u0);

  // Apply libCEED operator
  CeedOperatorApply(ctx_initial->op_apply, ceed_data->x_coord,
                    ctx_initial->y_ceed, CEED_REQUEST_IMMEDIATE);
  // Restore PETSc vectors
  CeedVectorTakeArray(ctx_initial->y_ceed, MemTypeP2C(u0_mem_type), NULL);
  PetscCall( VecRestoreArrayAndMemType(ctx_initial->Y_loc, &u0) );

  // Create global initial conditions
  PetscCall( VecZeroEntries(U0) );
  // Local-to-global
  PetscCall( DMLocalToGlobal(ctx_initial->dm, ctx_initial->Y_loc, ADD_VALUES,
                             U0) );
  PetscFunctionReturn(0);

}

// -----------------------------------------------------------------------------
// Create IResidual for Richard's problem
// -----------------------------------------------------------------------------
PetscErrorCode TSFormIResidual(TS ts, PetscReal time, Vec X, Vec X_t, Vec Y,
                               void *ctx_residual_ut) {
  OperatorApplyContext ctx   = (OperatorApplyContext)ctx_residual_ut;
  PetscScalar  *x, *x_t, *y;
  PetscMemType x_mem_type, x_t_mem_type, y_mem_type;
  PetscFunctionBeginUser;

  if(ctx->t != time) {
    CeedOperatorContextSetDouble(ctx->op_apply,
                                 ctx->solution_time_label, &time);
    ctx->t = time;
  }
  // Global-to-local
  PetscCall( DMGlobalToLocal(ctx->dm, X, INSERT_VALUES, ctx->X_loc) );
  PetscCall( DMGlobalToLocal(ctx->dm, X_t, INSERT_VALUES, ctx->X_t_loc) );

  // Place PETSc vectors in CEED vectors
  PetscCall( VecGetArrayReadAndMemType(ctx->X_loc,
                                       (const PetscScalar **)&x,
                                       &x_mem_type) );
  PetscCall( VecGetArrayReadAndMemType(ctx->X_t_loc,
                                       (const PetscScalar **)&x_t,
                                       &x_t_mem_type) );
  PetscCall( VecGetArrayAndMemType(ctx->Y_loc, &y, &y_mem_type) );
  CeedVectorSetArray(ctx->x_ceed, MemTypeP2C(x_mem_type),
                     CEED_USE_POINTER, x);
  CeedVectorSetArray(ctx->x_t_ceed, MemTypeP2C(x_t_mem_type),
                     CEED_USE_POINTER, x_t);
  CeedVectorSetArray(ctx->y_ceed, MemTypeP2C(y_mem_type),
                     CEED_USE_POINTER, y);

  // Apply CEED operator
  CeedOperatorApply(ctx->op_apply, ctx->x_ceed, ctx->y_ceed,
                    CEED_REQUEST_IMMEDIATE);

  // Restore vectors
  CeedVectorTakeArray(ctx->x_ceed, MemTypeP2C(x_mem_type), NULL);
  CeedVectorTakeArray(ctx->x_t_ceed, MemTypeP2C(x_t_mem_type), NULL);
  CeedVectorTakeArray(ctx->y_ceed, MemTypeP2C(y_mem_type), NULL);
  PetscCall( VecRestoreArrayReadAndMemType(ctx->X_loc,
             (const PetscScalar **)&x) );
  PetscCall( VecRestoreArrayReadAndMemType(ctx->X_t_loc,
             (const PetscScalar **)&x_t) );
  PetscCall( VecRestoreArrayAndMemType(ctx->Y_loc, &y) );

  // Local-to-Global
  PetscCall( VecZeroEntries(Y) );
  PetscCall( DMLocalToGlobal(ctx->dm, ctx->Y_loc, ADD_VALUES, Y) );

  PetscFunctionReturn(0);
}

// -----------------------------------------------------------------------------
// TS: Create, setup, and solve
// -----------------------------------------------------------------------------
PetscErrorCode TSSolveRichard(DM dm, Ceed ceed, CeedData ceed_data,
                              AppCtx app_ctx, OperatorApplyContext ctx_residual_ut,
                              Vec *U, TS *ts) {
  MPI_Comm       comm = app_ctx->comm;
  TSAdapt        adapt;
  PetscFunctionBeginUser;

  PetscCall( TSCreate(comm, ts) );
  PetscCall( TSSetDM(*ts, dm) );
  PetscCall( TSSetType(*ts, TSBDF) );
  PetscCall( TSSetIFunction(*ts, NULL, TSFormIResidual, ctx_residual_ut) );

  PetscCall( TSSetMaxTime(*ts, app_ctx->t_final) );
  PetscCall( TSSetExactFinalTime(*ts, TS_EXACTFINALTIME_STEPOVER) );
  PetscCall( TSSetTimeStep(*ts, 2.e-1) );
  PetscCall( TSGetAdapt(*ts, &adapt) );
  PetscCall( TSAdaptSetStepLimits(adapt, 1.e-12, 1.e2) );
  PetscCall( TSSetFromOptions(*ts) );
  ctx_residual_ut->t = -1.0;

  // Solve
  PetscScalar start_time;
  PetscCall( TSGetTime(*ts, &start_time) );

  PetscCall(TSSetTime(*ts, start_time));
  PetscCall(TSSetStepNumber(*ts, 0));

  PetscCall( TSSolve(*ts, *U) );

  PetscScalar    final_time;
  PetscCall( TSGetSolveTime(*ts, &final_time) );
  app_ctx->t_final = final_time;

  PetscFunctionReturn(0);
}


// -----------------------------------------------------------------------------
