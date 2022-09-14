#ifndef setup_ts_h
#define setup_ts_h

#include <ceed.h>
#include <petsc.h>

#include "structs.h"

PetscErrorCode CreateInitialConditions(CeedData ceed_data, AppCtx app_ctx,
                                       VecType vec_type, Vec U);
PetscErrorCode SetupResidualOperatorCtx_Ut(MPI_Comm comm, DM dm, Ceed ceed,
    CeedData ceed_data, OperatorApplyContext ctx_residual_ut);
PetscErrorCode SetupResidualOperatorCtx_U0(MPI_Comm comm, DM dm, Ceed ceed,
    CeedData ceed_data, OperatorApplyContext ctx_initial_u0);
PetscErrorCode SetupResidualOperatorCtx_P0(MPI_Comm comm, DM dm, Ceed ceed,
    CeedData ceed_data, OperatorApplyContext ctx_initial_p0);
PetscErrorCode TSFormIResidual(TS ts, PetscReal time, Vec X, Vec X_t, Vec Y,
                               void *ctx_residual_ut);
PetscErrorCode TSSolveRichard(CeedData ceed_data, AppCtx app_ctx,
                              TS ts, Vec *U);

#endif // setup_ts_h