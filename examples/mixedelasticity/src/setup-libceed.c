#include "../include/setup-libceed.h"
#include "../include/setup-boundary.h"
#include "../include/petsc-macros.h"

// -----------------------------------------------------------------------------
// Convert PETSc MemType to libCEED MemType
// -----------------------------------------------------------------------------
CeedMemType MemTypeP2C(PetscMemType mem_type) {
  return PetscMemTypeDevice(mem_type) ? CEED_MEM_DEVICE : CEED_MEM_HOST;
}
// -----------------------------------------------------------------------------
// Destroy libCEED objects
// -----------------------------------------------------------------------------
PetscErrorCode CeedDataDestroy(CeedData ceed_data) {

  PetscFunctionBegin;

  // Vectors
  CeedVectorDestroy(&ceed_data->x_ceed);
  CeedVectorDestroy(&ceed_data->y_ceed);
  CeedVectorDestroy(&ceed_data->x_coord);
  CeedVectorDestroy(&ceed_data->q_data);
  // Restrictions
  CeedElemRestrictionDestroy(&ceed_data->elem_restr_x);
  CeedElemRestrictionDestroy(&ceed_data->elem_restr_u);
  CeedElemRestrictionDestroy(&ceed_data->elem_restr_U_i);
  CeedElemRestrictionDestroy(&ceed_data->elem_restr_f_i);
  CeedElemRestrictionDestroy(&ceed_data->elem_restr_p);
  // Bases
  CeedBasisDestroy(&ceed_data->basis_x);
  CeedBasisDestroy(&ceed_data->basis_u);
  CeedBasisDestroy(&ceed_data->basis_p);
  CeedBasisDestroy(&ceed_data->basis_u_face);
  // QFunctions
  CeedQFunctionDestroy(&ceed_data->qf_true);
  CeedQFunctionDestroy(&ceed_data->qf_residual);
  CeedQFunctionDestroy(&ceed_data->qf_jacobian);
  CeedQFunctionDestroy(&ceed_data->qf_error);
  // Operators
  CeedOperatorDestroy(&ceed_data->op_true);
  CeedOperatorDestroy(&ceed_data->op_residual);
  CeedOperatorDestroy(&ceed_data->op_jacobian);
  CeedOperatorDestroy(&ceed_data->op_error);
  PetscCall( PetscFree(ceed_data) );

  PetscFunctionReturn(0);
};

// -----------------------------------------------------------------------------
// Get CEED restriction data from DMPlex
// -----------------------------------------------------------------------------
PetscErrorCode CreateRestrictionForDomain(Ceed ceed, DM dm, CeedInt height,
    DMLabel domain_label, CeedInt value, CeedElemRestriction *elem_restr_x,
    CeedElemRestriction *elem_restr_u, CeedElemRestriction *elem_restr_p) {
  PetscInt num_elem, elem_size, num_dof, num_comp, *elem_restr_offsets;
  DM       dm_coord;

  PetscFunctionBeginUser;
  // Create coordinate restriction
  PetscCall( DMGetCoordinateDM(dm, &dm_coord) );
  PetscCall( DMPlexGetLocalOffsets(dm_coord, domain_label, value, height, 0,
                                   &num_elem,
                                   &elem_size, &num_comp, &num_dof, &elem_restr_offsets) );
  CeedElemRestrictionCreate(ceed, num_elem, elem_size, num_comp,
                            1, num_dof, CEED_MEM_HOST, CEED_COPY_VALUES,
                            elem_restr_offsets, elem_restr_x);
  // Create solution restriction- displacement field
  PetscCall( DMPlexGetLocalOffsets(dm, domain_label, value, height, 0, &num_elem,
                                   &elem_size, &num_comp, &num_dof, &elem_restr_offsets) );
  CeedElemRestrictionCreate(ceed, num_elem, elem_size, num_comp,
                            1, num_dof, CEED_MEM_HOST, CEED_COPY_VALUES,
                            elem_restr_offsets, elem_restr_u);
  // Create solution restriction- displacement field
  PetscCall( DMPlexGetLocalOffsets(dm, domain_label, value, height, 1, &num_elem,
                                   &elem_size, &num_comp, &num_dof, &elem_restr_offsets) );
  CeedElemRestrictionCreate(ceed, num_elem, elem_size, num_comp,
                            1, num_dof, CEED_MEM_HOST, CEED_COPY_VALUES,
                            elem_restr_offsets, elem_restr_p);
  PetscCall( PetscFree(elem_restr_offsets) );

  PetscFunctionReturn(0);
};

// -----------------------------------------------------------------------------
// Set up libCEED on the fine grid for a given degree
// -----------------------------------------------------------------------------
PetscErrorCode SetupLibceed(DM dm, Ceed ceed, AppCtx app_ctx,
                            ProblemData problem_data,
                            CeedData ceed_data) {
  CeedInt       P_u = app_ctx->u_degree + 1, P_p = app_ctx->p_degree + 1;
  // Number of quadratures in 1D, q_extra is set in cl-options.c
  CeedInt       Q = P_u + 1 + app_ctx->q_extra;
  CeedInt       dim, num_comp_x, num_comp_u, num_comp_p, num_qpts;
  CeedInt       q_data_size = problem_data->q_data_size;
  Vec           coords;
  PetscInt      c_start, c_end, num_elem;
  const PetscScalar *coordArray;
  CeedQFunction qf_setup_geo;
  CeedOperator  op_setup_geo;

  PetscFunctionBeginUser;
  // ---------------------------------------------------------------------------
  // libCEED bases
  // ---------------------------------------------------------------------------
  dim = problem_data->dim;
  num_comp_x = dim;
  num_comp_u = dim;
  num_comp_p = 1;   // one scalar dof
  // Number of quadratures per element

  CeedBasisCreateTensorH1Lagrange(ceed, dim, num_comp_u, P_u, Q,
                                  problem_data->quadrature_mode, &ceed_data->basis_u);
  CeedBasisCreateTensorH1Lagrange(ceed, dim, num_comp_p, P_p, Q,
                                  problem_data->quadrature_mode, &ceed_data->basis_p);
  CeedBasisCreateTensorH1Lagrange(ceed, dim, num_comp_x, 2, Q,
                                  problem_data->quadrature_mode, &ceed_data->basis_x);

  // ---------------------------------------------------------------------------
  // libCEED restrictions
  // ---------------------------------------------------------------------------
  CeedInt height = 0; // 0 means no boundary conditions
  DMLabel domain_label = 0;
  PetscInt value = 0;
  CreateRestrictionForDomain(ceed, dm, height, domain_label, value,
                             &ceed_data->elem_restr_x, &ceed_data->elem_restr_u,
                             &ceed_data->elem_restr_p);
  // -- Geometric ceed_data restriction
  CeedBasisGetNumQuadraturePoints(ceed_data->basis_u, &num_qpts);
  PetscCall( DMPlexGetHeightStratum(dm, 0, &c_start, &c_end) );
  num_elem = c_end - c_start;
  CeedElemRestrictionCreateStrided(ceed, num_elem, num_qpts, q_data_size,
                                   num_elem * num_qpts * q_data_size,
                                   CEED_STRIDES_BACKEND, &ceed_data->elem_restr_q_data);
  // -- restriction for true solution U=[p,u]
  CeedElemRestrictionCreateStrided(ceed, num_elem, num_qpts, (dim+1),
                                   (dim+1)*num_elem*num_qpts,
                                   CEED_STRIDES_BACKEND, &ceed_data->elem_restr_U_i);
  CeedElemRestrictionCreateStrided(ceed, num_elem, num_qpts, dim,
                                   dim*num_elem*num_qpts,
                                   CEED_STRIDES_BACKEND, &ceed_data->elem_restr_f_i);
  // ---------------------------------------------------------------------------
  // Get physical coordinates
  // ---------------------------------------------------------------------------
  PetscCall( DMGetCoordinatesLocal(dm, &coords) );
  PetscCall( VecGetArrayRead(coords, &coordArray) );
  CeedElemRestrictionCreateVector(ceed_data->elem_restr_x, &ceed_data->x_coord,
                                  NULL);
  CeedVectorSetArray(ceed_data->x_coord, CEED_MEM_HOST, CEED_COPY_VALUES,
                     (PetscScalar *)coordArray);
  PetscCall( VecRestoreArrayRead(coords, &coordArray) );

  // ---------------------------------------------------------------------------
  // Persistent libCEED vectors
  // ---------------------------------------------------------------------------
  // -- Operator action variables: we use them in setup-solvers.c
  CeedElemRestrictionCreateVector(ceed_data->elem_restr_u, &ceed_data->x_ceed,
                                  NULL);
  CeedElemRestrictionCreateVector(ceed_data->elem_restr_u, &ceed_data->y_ceed,
                                  NULL);
  // ---------------------------------------------------------------------------
  // Setup volumetric geometric data
  // ---------------------------------------------------------------------------
  // -- Geometric data vector
  CeedVectorCreate(ceed, num_elem*num_qpts*q_data_size,
                   &ceed_data->q_data);

  // -- QFunction
  CeedQFunctionCreateInterior(ceed, 1, problem_data->setup_geo,
                              problem_data->setup_geo_loc, &qf_setup_geo);
  CeedQFunctionAddInput(qf_setup_geo, "dx", num_comp_x*dim, CEED_EVAL_GRAD);
  CeedQFunctionAddInput(qf_setup_geo, "weight", 1, CEED_EVAL_WEIGHT);
  CeedQFunctionAddOutput(qf_setup_geo, "q_data", q_data_size, CEED_EVAL_NONE);
  // -- Operator
  CeedOperatorCreate(ceed, qf_setup_geo, CEED_QFUNCTION_NONE,
                     CEED_QFUNCTION_NONE, &op_setup_geo);
  CeedOperatorSetField(op_setup_geo, "dx", ceed_data->elem_restr_x,
                       ceed_data->basis_x, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(op_setup_geo, "weight", CEED_ELEMRESTRICTION_NONE,
                       ceed_data->basis_x, CEED_VECTOR_NONE);
  CeedOperatorSetField(op_setup_geo, "q_data",
                       ceed_data->elem_restr_q_data,
                       CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
  // -- Compute the quadrature data
  CeedOperatorApply(op_setup_geo, ceed_data->x_coord, ceed_data->q_data,
                    CEED_REQUEST_IMMEDIATE);
  // -- Cleanup
  CeedQFunctionDestroy(&qf_setup_geo);
  CeedOperatorDestroy(&op_setup_geo);

  // ---------------------------------------------------------------------------
  // Setup true solution and force
  // ---------------------------------------------------------------------------
  CeedVector true_vec, true_force;
  CeedVectorCreate(ceed, num_elem*num_qpts*(dim+1), &true_vec);
  CeedVectorCreate(ceed, num_elem*num_qpts*dim, &true_force);
  // Create the q-function that sets up the RHS and true solution
  CeedQFunctionCreateInterior(ceed, 1, problem_data->true_solution,
                              problem_data->true_solution_loc, &ceed_data->qf_true);
  CeedQFunctionSetContext(ceed_data->qf_true, problem_data->qfunction_context);
  CeedQFunctionAddInput(ceed_data->qf_true, "x", num_comp_x, CEED_EVAL_INTERP);
  CeedQFunctionAddOutput(ceed_data->qf_true, "true force", dim, CEED_EVAL_NONE);
  CeedQFunctionAddOutput(ceed_data->qf_true, "true solution", dim+1,
                         CEED_EVAL_NONE);
  // Create the operator that builds the RHS and true solution
  CeedOperatorCreate(ceed, ceed_data->qf_true, CEED_QFUNCTION_NONE,
                     CEED_QFUNCTION_NONE,
                     &ceed_data->op_true);
  CeedOperatorSetField(ceed_data->op_true, "x", ceed_data->elem_restr_x,
                       ceed_data->basis_x, ceed_data->x_coord);
  CeedOperatorSetField(ceed_data->op_true, "true force",
                       ceed_data->elem_restr_f_i,
                       CEED_BASIS_COLLOCATED, true_force);
  CeedOperatorSetField(ceed_data->op_true, "true solution",
                       ceed_data->elem_restr_U_i,
                       CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
  // Setup RHS and true solution
  CeedOperatorApply(ceed_data->op_true, ceed_data->x_coord, true_vec,
                    CEED_REQUEST_IMMEDIATE);

  // Local residual evaluator
  // ---------------------------------------------------------------------------
  // Create the QFunction and Operator that computes the residual of the PDE.
  // ---------------------------------------------------------------------------
  // -- QFunction
  CeedQFunctionCreateInterior(ceed, 1, problem_data->residual,
                              problem_data->residual_loc, &ceed_data->qf_residual);
  CeedQFunctionSetContext(ceed_data->qf_residual,
                          problem_data->qfunction_context);
  CeedQFunctionAddInput(ceed_data->qf_residual, "dudX", num_comp_u*dim,
                        CEED_EVAL_GRAD);
  CeedQFunctionAddInput(ceed_data->qf_residual, "p", 1, CEED_EVAL_INTERP);
  CeedQFunctionAddInput(ceed_data->qf_residual, "q_data", q_data_size,
                        CEED_EVAL_NONE);
  CeedQFunctionAddInput(ceed_data->qf_residual, "true force", dim,
                        CEED_EVAL_NONE);
  CeedQFunctionAddOutput(ceed_data->qf_residual, "v", dim, CEED_EVAL_INTERP);
  CeedQFunctionAddOutput(ceed_data->qf_residual, "dvdX", num_comp_u*dim,
                         CEED_EVAL_GRAD);
  CeedQFunctionAddOutput(ceed_data->qf_residual, "q", 1, CEED_EVAL_INTERP);

  // -- Operator
  CeedOperatorCreate(ceed, ceed_data->qf_residual, CEED_QFUNCTION_NONE,
                     CEED_QFUNCTION_NONE,
                     &ceed_data->op_residual);
  CeedOperatorSetField(ceed_data->op_residual, "dudX", ceed_data->elem_restr_u,
                       ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_residual, "p", ceed_data->elem_restr_p,
                       ceed_data->basis_p, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_residual, "q_data",
                       ceed_data->elem_restr_q_data,
                       CEED_BASIS_COLLOCATED, ceed_data->q_data);
  CeedOperatorSetField(ceed_data->op_residual, "true force",
                       ceed_data->elem_restr_f_i,
                       CEED_BASIS_COLLOCATED, true_force);
  CeedOperatorSetField(ceed_data->op_residual, "v", ceed_data->elem_restr_u,
                       ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_residual, "dvdX", ceed_data->elem_restr_u,
                       ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_residual, "q", ceed_data->elem_restr_p,
                       ceed_data->basis_p, CEED_VECTOR_ACTIVE);
  // -- Save libCEED data to apply operator in matops.c

  // ---------------------------------------------------------------------------
  // Add Pressure boundary condition. See setup-boundary.c
  // ---------------------------------------------------------------------------
  //DMAddBoundariesPressure(ceed, ceed_data, app_ctx, problem_data, dm);

  // Local jacobian evaluator
  // ---------------------------------------------------------------------------
  // Create the QFunction and Operator that computes the jacobian of the PDE.
  // ---------------------------------------------------------------------------
  // -- QFunction
  CeedQFunctionCreateInterior(ceed, 1, problem_data->jacobian,
                              problem_data->jacobian_loc, &ceed_data->qf_jacobian);
  CeedQFunctionSetContext(ceed_data->qf_jacobian,
                          problem_data->qfunction_context);
  CeedQFunctionAddInput(ceed_data->qf_jacobian, "ddudX", num_comp_u*dim,
                        CEED_EVAL_GRAD);
  CeedQFunctionAddInput(ceed_data->qf_jacobian, "dp", 1, CEED_EVAL_INTERP);
  CeedQFunctionAddInput(ceed_data->qf_jacobian, "q_data", q_data_size,
                        CEED_EVAL_NONE);
  CeedQFunctionAddOutput(ceed_data->qf_jacobian, "dv", dim, CEED_EVAL_INTERP);
  CeedQFunctionAddOutput(ceed_data->qf_jacobian, "ddvdX", num_comp_u*dim,
                         CEED_EVAL_GRAD);
  CeedQFunctionAddOutput(ceed_data->qf_jacobian, "dq", 1, CEED_EVAL_INTERP);
  // -- Operator
  CeedOperatorCreate(ceed, ceed_data->qf_jacobian, CEED_QFUNCTION_NONE,
                     CEED_QFUNCTION_NONE,
                     &ceed_data->op_jacobian);
  CeedOperatorSetField(ceed_data->op_jacobian, "ddudX", ceed_data->elem_restr_u,
                       ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_jacobian, "dp", ceed_data->elem_restr_p,
                       ceed_data->basis_p, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_jacobian, "q_data",
                       ceed_data->elem_restr_q_data,
                       CEED_BASIS_COLLOCATED, ceed_data->q_data);
  CeedOperatorSetField(ceed_data->op_jacobian, "dv", ceed_data->elem_restr_u,
                       ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_jacobian, "ddvdX", ceed_data->elem_restr_u,
                       ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_jacobian, "dq", ceed_data->elem_restr_p,
                       ceed_data->basis_p, CEED_VECTOR_ACTIVE);
  // -- Save libCEED data to apply operator in matops.c

  // ---------------------------------------------------------------------------
  // Setup Error Qfunction
  // ---------------------------------------------------------------------------
  // Create the q-function that sets up the error
  CeedQFunctionCreateInterior(ceed, 1, problem_data->error,
                              problem_data->error_loc, &ceed_data->qf_error);
  CeedQFunctionSetContext(ceed_data->qf_error, problem_data->qfunction_context);
  CeedQFunctionAddInput(ceed_data->qf_error, "q_data", q_data_size,
                        CEED_EVAL_NONE);
  CeedQFunctionAddInput(ceed_data->qf_error, "u", dim, CEED_EVAL_INTERP);
  CeedQFunctionAddInput(ceed_data->qf_error, "p", 1, CEED_EVAL_INTERP);
  CeedQFunctionAddInput(ceed_data->qf_error, "true solution", dim+1,
                        CEED_EVAL_NONE);
  CeedQFunctionAddOutput(ceed_data->qf_error, "error", dim+1, CEED_EVAL_NONE);
  // Create the operator that builds the error
  CeedOperatorCreate(ceed, ceed_data->qf_error, CEED_QFUNCTION_NONE,
                     CEED_QFUNCTION_NONE,
                     &ceed_data->op_error);
  CeedOperatorSetField(ceed_data->op_error, "q_data",
                       ceed_data->elem_restr_q_data,
                       CEED_BASIS_COLLOCATED, ceed_data->q_data);
  CeedOperatorSetField(ceed_data->op_error, "u", ceed_data->elem_restr_u,
                       ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_error, "p", ceed_data->elem_restr_p,
                       ceed_data->basis_p, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_error, "true solution",
                       ceed_data->elem_restr_U_i,
                       CEED_BASIS_COLLOCATED, true_vec);
  CeedOperatorSetField(ceed_data->op_error, "error", ceed_data->elem_restr_U_i,
                       CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
  // -- Save libCEED data to apply operator in matops.c

  PetscFunctionReturn(0);
};
// -----------------------------------------------------------------------------