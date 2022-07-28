#ifndef structs_h
#define structs_h

#include <ceed.h>
#include <petsc.h>

// Application context from user command line options
typedef struct AppCtx_ *AppCtx;
struct AppCtx_ {
  MPI_Comm          comm;
  // Order of fe space, extra quadrature pts
  PetscInt          u_degree, p_degree;
  PetscInt          q_extra;
  // For applying traction BCs
  PetscInt          bc_pressure_count;
  PetscInt          bc_pressure_faces[16];
  PetscScalar       bc_pressure_vector[16][3];
  // Problem type arguments
  PetscFunctionList problems;
  char              problem_name[PETSC_MAX_PATH_LEN];

};

// PETSc operator contexts
typedef struct OperatorApplyContext_ *OperatorApplyContext;
struct OperatorApplyContext_ {
  MPI_Comm        comm;
  Vec             X_loc, Y_loc;
  CeedVector      x_ceed, y_ceed;
  CeedOperator    op_apply;
  DM              dm;
  Ceed            ceed;
};

// libCEED data struct
typedef struct CeedData_ *CeedData;
struct CeedData_ {
  CeedBasis            basis_x, basis_u, basis_p, basis_u_face;
  CeedElemRestriction  elem_restr_x, elem_restr_u, elem_restr_U_i,
                       elem_restr_p, elem_restr_q_data, elem_restr_f_i;
  CeedQFunction        qf_true, qf_residual, qf_jacobian, qf_error;
  CeedOperator         op_true, op_residual, op_jacobian, op_error;
  CeedVector           x_ceed, y_ceed, x_coord, q_data;
};

// Problem specific data
typedef struct ProblemData_ *ProblemData;
struct ProblemData_ {
  CeedQFunctionUser setup_geo, true_solution, residual, jacobian, error,
                    bc_pressure;
  const char        *setup_geo_loc, *true_solution_loc, *residual_loc,
        *jacobian_loc,
        *error_loc, *bc_pressure_loc;
  CeedQuadMode      quadrature_mode;
  CeedInt           dim, q_data_size, q_data_size_face;
  CeedQFunctionContext qfunction_context;
};

#endif // structs_h