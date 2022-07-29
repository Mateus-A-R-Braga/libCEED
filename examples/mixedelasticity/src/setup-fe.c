#include "../include/setup-fe.h"
#include "../include/setup-boundary.h"

// ---------------------------------------------------------------------------
// Set-up FE
// ---------------------------------------------------------------------------
PetscErrorCode SetupFE(ProblemData problem_data,
                       AppCtx app_ctx, DM dm) {
  // Two FE space for displacement and pressure
  PetscFE   fe[2];
  // number of quadrature points
  PetscInt  q_degree = app_ctx->u_degree + app_ctx->q_extra;
  PetscBool is_simplex = PETSC_TRUE;
  MPI_Comm  comm;
  PetscFunctionBeginUser;

  // Check if simplex or tensor-product element
  PetscCall( DMPlexIsSimplex(dm, &is_simplex) );
  // Create FE space
  PetscCall( PetscObjectGetComm((PetscObject)dm, &comm) );
  PetscCall( PetscFECreateLagrange(comm, problem_data->dim,
                                   problem_data->dim, is_simplex,
                                   app_ctx->u_degree, q_degree, &fe[0]) );
  PetscCall( PetscFECreateLagrange(comm, problem_data->dim, 1, is_simplex,
                                   app_ctx->p_degree, q_degree, &fe[1]) );
  PetscCall( PetscFECopyQuadrature(fe[0], fe[1]) );
  PetscCall( PetscObjectSetName((PetscObject)fe[0], "Displacement") );
  PetscCall( PetscObjectSetName((PetscObject)fe[1], "Pressure") );
  PetscCall( DMAddField(dm, NULL, (PetscObject)fe[0]) );
  PetscCall( DMAddField(dm, NULL, (PetscObject)fe[1]) );
  PetscCall( DMCreateDS(dm) );

  {
    // create FE field for coordinates
    PetscFE fe_coords;
    PetscInt num_comp_coord;
    PetscCall( DMGetCoordinateDim(dm, &num_comp_coord) );
    PetscCall( PetscFECreateLagrange(comm, problem_data->dim, num_comp_coord,
                                     is_simplex, 1, q_degree,
                                     &fe_coords) );
    PetscCall( DMProjectCoordinates(dm, fe_coords) );
    PetscCall( PetscFEDestroy(&fe_coords) );
  }

  if (app_ctx->setup_dirichlet) {
    // Dirichlet boundaries
    PetscCall(DMAddBoundariesDirichlet(dm));
  }

  if (!is_simplex) {
    DM dm_coord;
    PetscCall(DMGetCoordinateDM(dm, &dm_coord));
    PetscCall(DMPlexSetClosurePermutationTensor(dm, PETSC_DETERMINE, NULL));
    PetscCall(DMPlexSetClosurePermutationTensor(dm_coord, PETSC_DETERMINE, NULL));
  }

  // Cleanup
  PetscCall( PetscFEDestroy(&fe[0]) );
  PetscCall( PetscFEDestroy(&fe[1]) );

  PetscFunctionReturn(0);
};