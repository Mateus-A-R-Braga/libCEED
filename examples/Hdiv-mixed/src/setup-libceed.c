#include "../include/setup-libceed.h"
#include "../include/setup-boundary.h"
#include "../include/petsc-macros.h"
#include "../basis/Hdiv-quad.h"
#include "../basis/Hdiv-hex.h"
#include "../basis/L2-P0.h"
#include "ceed/ceed.h"
#include <stdio.h>

// -----------------------------------------------------------------------------
// Convert PETSc MemType to libCEED MemType
// -----------------------------------------------------------------------------
CeedMemType MemTypeP2C(PetscMemType mem_type) {
  return PetscMemTypeDevice(mem_type) ? CEED_MEM_DEVICE : CEED_MEM_HOST;
}
// -----------------------------------------------------------------------------
// Destroy libCEED objects
// -----------------------------------------------------------------------------
PetscErrorCode CeedDataDestroy(CeedData ceed_data, ProblemData problem_data) {

  PetscFunctionBegin;

  // Vectors
  CeedVectorDestroy(&ceed_data->x_ceed);
  CeedVectorDestroy(&ceed_data->y_ceed);
  CeedVectorDestroy(&ceed_data->x_coord);
  CeedVectorDestroy(&ceed_data->x_t_ceed);
  // Restrictions
  CeedElemRestrictionDestroy(&ceed_data->elem_restr_x);
  CeedElemRestrictionDestroy(&ceed_data->elem_restr_u);
  CeedElemRestrictionDestroy(&ceed_data->elem_restr_U_i); // U = [p,u]
  CeedElemRestrictionDestroy(&ceed_data->elem_restr_p);
  CeedElemRestrictionDestroy(&ceed_data->elem_restr_p_i);
  // Bases
  CeedBasisDestroy(&ceed_data->basis_x);
  CeedBasisDestroy(&ceed_data->basis_u);
  CeedBasisDestroy(&ceed_data->basis_p);
  CeedBasisDestroy(&ceed_data->basis_u_face);
  // QFunctions
  CeedQFunctionDestroy(&ceed_data->qf_residual);
  CeedQFunctionDestroy(&ceed_data->qf_jacobian);
  CeedQFunctionDestroy(&ceed_data->qf_error);
  CeedQFunctionDestroy(&ceed_data->qf_true);
  // Operators
  CeedOperatorDestroy(&ceed_data->op_residual);
  CeedOperatorDestroy(&ceed_data->op_jacobian);
  CeedOperatorDestroy(&ceed_data->op_error);
  CeedOperatorDestroy(&ceed_data->op_true);
  if (problem_data->has_ts) {
    // -- Cleanup
    CeedQFunctionDestroy(&ceed_data->qf_ics);
    CeedOperatorDestroy(&ceed_data->op_ics);
  }
  PetscCall( PetscFree(ceed_data) );

  PetscFunctionReturn(0);
};

// -----------------------------------------------------------------------------
// Utility function - essential BC dofs are encoded in closure indices as -(i+1)
// -----------------------------------------------------------------------------
PetscInt Involute(PetscInt i) {
  return i >= 0 ? i : -(i + 1);
};

// -----------------------------------------------------------------------------
// Get CEED restriction data from DMPlex
// -----------------------------------------------------------------------------
PetscErrorCode CreateRestrictionFromPlex(Ceed ceed, DM dm, CeedInt height,
    DMLabel domain_label, CeedInt value, CeedElemRestriction *elem_restr) {
  PetscInt num_elem, elem_size, num_dof, num_comp, *elem_restr_offsets;

  PetscFunctionBeginUser;

  PetscCall( DMPlexGetLocalOffsets(dm, domain_label, value, height, 0, &num_elem,
                                   &elem_size, &num_comp, &num_dof, &elem_restr_offsets) );

  CeedElemRestrictionCreate(ceed, num_elem, elem_size, num_comp,
                            1, num_dof, CEED_MEM_HOST, CEED_COPY_VALUES,
                            elem_restr_offsets, elem_restr);
  PetscCall( PetscFree(elem_restr_offsets) );

  PetscFunctionReturn(0);
};

// -----------------------------------------------------------------------------
// Get Oriented CEED restriction data from DMPlex
// -----------------------------------------------------------------------------
PetscErrorCode CreateRestrictionFromPlexOriented(Ceed ceed, DM dm,
    CeedInt P, CeedElemRestriction *elem_restr_oriented,
    CeedElemRestriction *elem_restr) {
  PetscSection section;
  PetscInt p, num_elem, num_dof, *restr_indices_u, *restr_indices_p,
           elem_offset, num_fields, dim, c_start, c_end;
  Vec U_loc;
  const PetscInt *ornt; // this is for orientation of dof
  PetscFunctionBeginUser;
  PetscCall( DMGetDimension(dm, &dim) );
  PetscCall( DMGetLocalSection(dm, &section) );
  PetscCall( PetscSectionGetNumFields(section, &num_fields) );
  PetscInt num_comp[num_fields], field_offsets[num_fields+1];
  field_offsets[0] = 0;
  for (PetscInt f = 0; f < num_fields; f++) {
    PetscCall( PetscSectionGetFieldComponents(section, f, &num_comp[f]) );
    field_offsets[f+1] = field_offsets[f] + num_comp[f];
  }
  PetscCall( DMPlexGetHeightStratum(dm, 0, &c_start, &c_end) );
  num_elem = c_end - c_start;
  PetscCall( PetscMalloc1(num_elem*dim*PetscPowInt(P, dim),
                          &restr_indices_u) );
  PetscCall( PetscMalloc1(num_elem,&restr_indices_p) );
  bool *orient_indices_u; // to flip the dof
  PetscCall( PetscMalloc1(num_elem*dim*PetscPowInt(P, dim), &orient_indices_u) );
  for (p = 0, elem_offset = 0; p < num_elem; p++) {
    PetscInt num_indices, *indices, faces_per_elem, dofs_per_face;
    PetscCall( DMPlexGetClosureIndices(dm, section, section, p, PETSC_TRUE,
                                       &num_indices, &indices, NULL, NULL) );

    restr_indices_p[p] = indices[num_indices - 1];
    PetscCall( DMPlexGetConeOrientation(dm, p, &ornt) );
    // Get number of faces per element
    PetscCall( DMPlexGetConeSize(dm, p, &faces_per_elem) );
    dofs_per_face = faces_per_elem - 2;
    for (PetscInt f = 0; f < faces_per_elem; f++) {
      for (PetscInt i = 0; i < dofs_per_face; i++) {
        PetscInt ii = dofs_per_face*f + i;
        // Essential boundary conditions are encoded as -(loc+1), but we don't care so we decode.
        PetscInt loc = Involute(indices[ii*num_comp[0]]);
        restr_indices_u[elem_offset] = loc;
        // Set orientation
        orient_indices_u[elem_offset] = ornt[f] < 0;
        elem_offset++;
      }
    }
    PetscCall( DMPlexRestoreClosureIndices(dm, section, section, p, PETSC_TRUE,
                                           &num_indices, &indices, NULL, NULL) );
  }
  //if (elem_offset != num_elem*dim*PetscPowInt(P, dim))
  //  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_LIB,
  //          "ElemRestriction of size (%" PetscInt_FMT ", %" PetscInt_FMT" )
  //          initialized %" PetscInt_FMT " nodes", num_elem,
  //          dim*PetscPowInt(P, dim),elem_offset);

  PetscCall( DMGetLocalVector(dm, &U_loc) );
  PetscCall( VecGetLocalSize(U_loc, &num_dof) );
  PetscCall( DMRestoreLocalVector(dm, &U_loc) );
  // dof per element in Hdiv is dim*P^dim, for linear element P=2
  CeedElemRestrictionCreateOriented(ceed, num_elem, dim*PetscPowInt(P, dim),
                                    1, 1, num_dof, CEED_MEM_HOST, CEED_COPY_VALUES,
                                    restr_indices_u, orient_indices_u,
                                    elem_restr_oriented);
  CeedElemRestrictionCreate(ceed, num_elem, 1,
                            1, 1, num_dof, CEED_MEM_HOST, CEED_COPY_VALUES,
                            restr_indices_p, elem_restr);
  PetscCall( PetscFree(restr_indices_p) );
  PetscCall( PetscFree(restr_indices_u) );
  PetscCall( PetscFree(orient_indices_u) );
  PetscFunctionReturn(0);
};

// -----------------------------------------------------------------------------
// Set up libCEED on the fine grid for a given degree
// -----------------------------------------------------------------------------
PetscErrorCode SetupLibceed(DM dm, Ceed ceed, AppCtx app_ctx,
                            OperatorApplyContext ctx_residual_ut,
                            ProblemData problem_data, CeedData ceed_data) {
  CeedInt       P = app_ctx->degree + 1;
  // Number of quadratures in 1D, q_extra is set in cl-options.c
  CeedInt       Q = P + 1 + app_ctx->q_extra;
  CeedInt       dim, num_comp_x, num_comp_u, num_comp_p;
  DM            dm_coord;
  Vec           coords;
  PetscInt      c_start, c_end, num_elem;
  const PetscScalar *coordArray;

  PetscFunctionBeginUser;
  // ---------------------------------------------------------------------------
  // libCEED bases:Hdiv basis_u and Lagrange basis_x
  // ---------------------------------------------------------------------------
  dim = problem_data->dim;
  num_comp_x = dim;
  num_comp_u = 1;   // one vector dof
  num_comp_p = 1;   // one scalar dof
  // Number of quadratures per element
  CeedInt       num_qpts = PetscPowInt(Q, dim);
  // Pressure and velocity dof per element
  CeedInt       P_p = 1, P_u = dim*PetscPowInt(P, dim);
  CeedScalar    q_ref[dim*num_qpts], q_weights[num_qpts];
  CeedScalar    div[P_u*num_qpts], interp_u[dim*P_u*num_qpts],
                interp_p[P_p*num_qpts], *grad=NULL;
  if (dim == 2) {
    HdivBasisQuad(Q, q_ref, q_weights, interp_u, div,
                  problem_data->quadrature_mode);
    CeedBasisCreateHdiv(ceed, CEED_TOPOLOGY_QUAD, num_comp_u, P_u, num_qpts,
                        interp_u, div, q_ref, q_weights, &ceed_data->basis_u);
    L2BasisP0(dim, Q, q_ref, q_weights, interp_p, problem_data->quadrature_mode);
    CeedBasisCreateH1(ceed, CEED_TOPOLOGY_QUAD, num_comp_p, 1, num_qpts, interp_p,
                      grad, q_ref,q_weights, &ceed_data->basis_p);
    HdivBasisQuad(Q, q_ref, q_weights, interp_u, div,
                  CEED_GAUSS_LOBATTO);
    CeedBasisCreateHdiv(ceed, CEED_TOPOLOGY_QUAD, num_comp_u, P_u, num_qpts,
                        interp_u, div, q_ref, q_weights, &ceed_data->basis_u_face);
  } else {
    HdivBasisHex(Q, q_ref, q_weights, interp_u, div, problem_data->quadrature_mode);
    CeedBasisCreateHdiv(ceed, CEED_TOPOLOGY_HEX, num_comp_u, P_u, num_qpts,
                        interp_u, div, q_ref, q_weights, &ceed_data->basis_u);
    L2BasisP0(dim, Q, q_ref, q_weights, interp_p, problem_data->quadrature_mode);
    CeedBasisCreateH1(ceed, CEED_TOPOLOGY_HEX, num_comp_p, 1, num_qpts, interp_p,
                      grad, q_ref,q_weights, &ceed_data->basis_p);
    HdivBasisHex(Q, q_ref, q_weights, interp_u, div, CEED_GAUSS_LOBATTO);
    CeedBasisCreateHdiv(ceed, CEED_TOPOLOGY_HEX, num_comp_u, P_u, num_qpts,
                        interp_u, div, q_ref, q_weights, &ceed_data->basis_u_face);
  }

  CeedBasisCreateTensorH1Lagrange(ceed, dim, num_comp_x, 2, Q,
                                  problem_data->quadrature_mode, &ceed_data->basis_x);

  // ---------------------------------------------------------------------------
  // libCEED restrictions
  // ---------------------------------------------------------------------------
  PetscCall( DMGetCoordinateDM(dm, &dm_coord) );
  PetscCall( DMPlexSetClosurePermutationTensor(dm_coord, PETSC_DETERMINE, NULL) );
  CeedInt height = 0; // 0 means no boundary conditions
  DMLabel domain_label = 0;
  PetscInt value = 0;
  // -- Coordinate restriction
  PetscCall( CreateRestrictionFromPlex(ceed, dm_coord, height, domain_label,
                                       value, &ceed_data->elem_restr_x) );
  // -- Solution restriction
  PetscCall( CreateRestrictionFromPlexOriented(ceed, dm, P,
             &ceed_data->elem_restr_u, &ceed_data->elem_restr_p) );
  // -- Geometric ceed_data restriction
  PetscCall( DMPlexGetHeightStratum(dm, 0, &c_start, &c_end) );
  num_elem = c_end - c_start;

  CeedElemRestrictionCreateStrided(ceed, num_elem, num_qpts, (dim+1),
                                   (dim+1)*num_elem*num_qpts,
                                   CEED_STRIDES_BACKEND, &ceed_data->elem_restr_U_i);
  CeedElemRestrictionCreateStrided(ceed, num_elem, num_qpts, 1,
                                   1*num_elem*num_qpts,
                                   CEED_STRIDES_BACKEND, &ceed_data->elem_restr_p_i);
  // ---------------------------------------------------------------------------
  //  Get physical coordinates
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
  // -- Operator action variables: we use them in setup-ts.c
  CeedElemRestrictionCreateVector(ceed_data->elem_restr_u, &ceed_data->x_t_ceed,
                                  NULL);
  if (problem_data->has_ts) {
    CeedBasis basis_xc;
    CeedBasisCreateTensorH1Lagrange(ceed, dim, num_comp_x, 2, 2,
                                    CEED_GAUSS_LOBATTO, &basis_xc);
    // ---------------------------------------------------------------------------
    // Setup qfunction for initial conditions
    // ---------------------------------------------------------------------------
    CeedQFunctionCreateInterior(ceed, 1, problem_data->ics,
                                problem_data->ics_loc, &ceed_data->qf_ics);
    CeedQFunctionSetContext(ceed_data->qf_ics, problem_data->qfunction_context);
    CeedQFunctionAddInput(ceed_data->qf_ics, "x", num_comp_x, CEED_EVAL_INTERP);
    CeedQFunctionAddInput(ceed_data->qf_ics, "dx", dim*dim, CEED_EVAL_GRAD);
    CeedQFunctionAddOutput(ceed_data->qf_ics, "u_0", dim, CEED_EVAL_NONE);
    CeedQFunctionAddOutput(ceed_data->qf_ics, "p_0", 1, CEED_EVAL_NONE);
    // Create the operator that builds the initial conditions
    CeedOperatorCreate(ceed, ceed_data->qf_ics, CEED_QFUNCTION_NONE,
                       CEED_QFUNCTION_NONE,
                       &ceed_data->op_ics);
    CeedOperatorSetField(ceed_data->op_ics, "x", ceed_data->elem_restr_x,
                         basis_xc, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(ceed_data->op_ics, "dx", ceed_data->elem_restr_x,
                         basis_xc, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(ceed_data->op_ics, "u_0", ceed_data->elem_restr_u,
                         CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(ceed_data->op_ics, "p_0", ceed_data->elem_restr_p,
                         CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
    // -- Save libCEED data to apply operator in setup-ts.c
    //CeedVector U0;
    //CeedElemRestrictionCreateVector(ceed_data->elem_restr_u, &U0,
    //                              NULL);
    //CeedOperatorApply(ceed_data->op_ics, ceed_data->x_coord,
    //                U0,
    //                CEED_REQUEST_IMMEDIATE);
    //CeedVectorView(U0, "%12.8f", stdout);
    CeedBasisDestroy(&basis_xc);
  }

  // ---------------------------------------------------------------------------
  // Setup true solution and force
  // ---------------------------------------------------------------------------
  CeedVector true_vec, true_force;
  CeedVectorCreate(ceed, num_elem*num_qpts*(dim+1), &true_vec);
  CeedVectorCreate(ceed, num_elem*num_qpts*1, &true_force);
  // Create the q-function that sets up the RHS and true solution
  CeedQFunctionCreateInterior(ceed, 1, problem_data->true_solution,
                              problem_data->true_solution_loc, &ceed_data->qf_true);
  CeedQFunctionSetContext(ceed_data->qf_true, problem_data->qfunction_context);
  CeedQFunctionAddInput(ceed_data->qf_true, "x", num_comp_x, CEED_EVAL_INTERP);
  CeedQFunctionAddOutput(ceed_data->qf_true, "true force", 1, CEED_EVAL_NONE);
  CeedQFunctionAddOutput(ceed_data->qf_true, "true solution", dim+1,
                         CEED_EVAL_NONE);
  // Create the operator that builds the RHS and true solution
  CeedOperatorCreate(ceed, ceed_data->qf_true, CEED_QFUNCTION_NONE,
                     CEED_QFUNCTION_NONE,
                     &ceed_data->op_true);
  CeedOperatorSetField(ceed_data->op_true, "x", ceed_data->elem_restr_x,
                       ceed_data->basis_x, ceed_data->x_coord);
  CeedOperatorSetField(ceed_data->op_true, "true force",
                       ceed_data->elem_restr_p_i,
                       CEED_BASIS_COLLOCATED, true_force);
  CeedOperatorSetField(ceed_data->op_true, "true solution",
                       ceed_data->elem_restr_U_i,
                       CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
  if (problem_data->has_ts) {
    double final_time = app_ctx->t_final;
    CeedOperatorContextGetFieldLabel(ceed_data->op_true, "final_time",
                                     &ctx_residual_ut->final_time_label);
    CeedOperatorContextSetDouble(ceed_data->op_true,
                                 ctx_residual_ut->final_time_label, &final_time);
  }
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
  //CeedQFunctionContextDestroy(&problem_data->qfunction_context);
  CeedQFunctionAddInput(ceed_data->qf_residual, "weight", 1, CEED_EVAL_WEIGHT);
  CeedQFunctionAddInput(ceed_data->qf_residual, "dx", dim*dim, CEED_EVAL_GRAD);
  CeedQFunctionAddInput(ceed_data->qf_residual, "u", dim, CEED_EVAL_INTERP);
  CeedQFunctionAddInput(ceed_data->qf_residual, "div_u", 1, CEED_EVAL_DIV);
  CeedQFunctionAddInput(ceed_data->qf_residual, "p", 1, CEED_EVAL_INTERP);
  CeedQFunctionAddInput(ceed_data->qf_residual, "true force", 1, CEED_EVAL_NONE);
  CeedQFunctionAddInput(ceed_data->qf_residual, "x", num_comp_x,
                        CEED_EVAL_INTERP);
  if (problem_data->has_ts) {
    CeedQFunctionAddInput(ceed_data->qf_residual, "p_t", 1, CEED_EVAL_INTERP);
  }
  CeedQFunctionAddOutput(ceed_data->qf_residual, "v", dim, CEED_EVAL_INTERP);
  CeedQFunctionAddOutput(ceed_data->qf_residual, "div_v", 1, CEED_EVAL_DIV);
  CeedQFunctionAddOutput(ceed_data->qf_residual, "q", 1, CEED_EVAL_INTERP);

  // -- Operator
  CeedOperatorCreate(ceed, ceed_data->qf_residual, CEED_QFUNCTION_NONE,
                     CEED_QFUNCTION_NONE,
                     &ceed_data->op_residual);
  CeedOperatorSetField(ceed_data->op_residual, "weight",
                       CEED_ELEMRESTRICTION_NONE,
                       ceed_data->basis_x, CEED_VECTOR_NONE);
  CeedOperatorSetField(ceed_data->op_residual, "dx", ceed_data->elem_restr_x,
                       ceed_data->basis_x, ceed_data->x_coord);
  CeedOperatorSetField(ceed_data->op_residual, "u", ceed_data->elem_restr_u,
                       ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_residual, "div_u", ceed_data->elem_restr_u,
                       ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_residual, "p", ceed_data->elem_restr_p,
                       ceed_data->basis_p, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_residual, "true force",
                       ceed_data->elem_restr_p_i,
                       CEED_BASIS_COLLOCATED, true_force);
  CeedOperatorSetField(ceed_data->op_residual, "x", ceed_data->elem_restr_x,
                       ceed_data->basis_x, ceed_data->x_coord);
  if (problem_data->has_ts) {
    CeedOperatorSetField(ceed_data->op_residual, "p_t", ceed_data->elem_restr_p,
                         ceed_data->basis_p, ceed_data->x_t_ceed);
  }
  CeedOperatorSetField(ceed_data->op_residual, "v", ceed_data->elem_restr_u,
                       ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_residual, "div_v", ceed_data->elem_restr_u,
                       ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_residual, "q", ceed_data->elem_restr_p,
                       ceed_data->basis_p, CEED_VECTOR_ACTIVE);
  // -- Save libCEED data to apply operator in matops.c
  CeedOperatorContextGetFieldLabel(ceed_data->op_residual, "time",
                                   &ctx_residual_ut->solution_time_label);
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
  //CeedQFunctionContextDestroy(&problem_data->qfunction_context);
  CeedQFunctionAddInput(ceed_data->qf_jacobian, "weight", 1, CEED_EVAL_WEIGHT);
  CeedQFunctionAddInput(ceed_data->qf_jacobian, "dx", dim*dim, CEED_EVAL_GRAD);
  CeedQFunctionAddInput(ceed_data->qf_jacobian, "du", dim, CEED_EVAL_INTERP);
  CeedQFunctionAddInput(ceed_data->qf_jacobian, "div_du", 1, CEED_EVAL_DIV);
  CeedQFunctionAddInput(ceed_data->qf_jacobian, "dp", 1, CEED_EVAL_INTERP);
  CeedQFunctionAddInput(ceed_data->qf_jacobian, "x", num_comp_x,
                        CEED_EVAL_INTERP);
  //CeedQFunctionAddInput(qf_jacobian, "u", dim, CEED_EVAL_INTERP);
  //CeedQFunctionAddInput(qf_jacobian, "p", 1, CEED_EVAL_INTERP);
  CeedQFunctionAddOutput(ceed_data->qf_jacobian, "dv", dim, CEED_EVAL_INTERP);
  CeedQFunctionAddOutput(ceed_data->qf_jacobian, "div_dv", 1, CEED_EVAL_DIV);
  CeedQFunctionAddOutput(ceed_data->qf_jacobian, "dq", 1, CEED_EVAL_INTERP);
  // -- Operator
  CeedOperatorCreate(ceed, ceed_data->qf_jacobian, CEED_QFUNCTION_NONE,
                     CEED_QFUNCTION_NONE,
                     &ceed_data->op_jacobian);
  CeedOperatorSetField(ceed_data->op_jacobian, "weight",
                       CEED_ELEMRESTRICTION_NONE,
                       ceed_data->basis_x, CEED_VECTOR_NONE);
  CeedOperatorSetField(ceed_data->op_jacobian, "dx", ceed_data->elem_restr_x,
                       ceed_data->basis_x, ceed_data->x_coord);
  CeedOperatorSetField(ceed_data->op_jacobian, "du", ceed_data->elem_restr_u,
                       ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_jacobian, "div_du", ceed_data->elem_restr_u,
                       ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_jacobian, "dp", ceed_data->elem_restr_p,
                       ceed_data->basis_p, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_jacobian, "x", ceed_data->elem_restr_x,
                       ceed_data->basis_x, ceed_data->x_coord);
  //CeedOperatorSetField(op_jacobian, "u", ceed_data->elem_restr_u,
  //                     ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  //CeedOperatorSetField(op_jacobian, "p", ceed_data->elem_restr_p,
  //                     ceed_data->basis_p, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_jacobian, "dv", ceed_data->elem_restr_u,
                       ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_jacobian, "div_dv", ceed_data->elem_restr_u,
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
  CeedQFunctionAddInput(ceed_data->qf_error, "weight", 1, CEED_EVAL_WEIGHT);
  CeedQFunctionAddInput(ceed_data->qf_error, "dx", dim*dim, CEED_EVAL_GRAD);
  CeedQFunctionAddInput(ceed_data->qf_error, "u", dim, CEED_EVAL_INTERP);
  CeedQFunctionAddInput(ceed_data->qf_error, "p", 1, CEED_EVAL_INTERP);
  CeedQFunctionAddInput(ceed_data->qf_error, "true solution", dim+1,
                        CEED_EVAL_NONE);
  CeedQFunctionAddOutput(ceed_data->qf_error, "error", dim+1, CEED_EVAL_NONE);
  // Create the operator that builds the error
  CeedOperatorCreate(ceed, ceed_data->qf_error, CEED_QFUNCTION_NONE,
                     CEED_QFUNCTION_NONE,
                     &ceed_data->op_error);
  CeedOperatorSetField(ceed_data->op_error, "weight", CEED_ELEMRESTRICTION_NONE,
                       ceed_data->basis_x, CEED_VECTOR_NONE);
  CeedOperatorSetField(ceed_data->op_error, "dx", ceed_data->elem_restr_x,
                       ceed_data->basis_x, ceed_data->x_coord);
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

  // -- Cleanup
  CeedVectorDestroy(&true_vec);
  CeedVectorDestroy(&true_force);
  PetscFunctionReturn(0);
};
// -----------------------------------------------------------------------------