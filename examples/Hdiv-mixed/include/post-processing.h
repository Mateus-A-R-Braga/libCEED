#ifndef post_processing_h
#define post_processing_h

#include <ceed.h>
#include <petsc.h>

#include "structs.h"
#include "../include/setup-libceed.h"
PetscErrorCode PrintOutput(Ceed ceed, AppCtx app_ctx, PetscBool has_ts,
                           CeedMemType mem_type_backend,
                           TS ts, SNES snes, KSP ksp,
                           Vec U, CeedScalar l2_error_u,
                           CeedScalar l2_error_p);
PetscErrorCode SetupProjectVelocityCtx_Hdiv(MPI_Comm comm, DM dm, Ceed ceed,
    CeedData ceed_data, OperatorApplyContext ctx_post_Hdiv);
PetscErrorCode SetupProjectVelocityCtx_H1(MPI_Comm comm, DM dm_H1, Ceed ceed,
    CeedData ceed_data, OperatorApplyContext ctx_post_H1);
PetscErrorCode ProjectVelocity(CeedData ceed_data,
                               Vec U, VecType vec_type, Vec *U_H1,
                               OperatorApplyContext ctx_post_Hdiv,
                               OperatorApplyContext ctx_post_H1);
#endif // post_processing_h
