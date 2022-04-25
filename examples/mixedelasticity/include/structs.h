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

// libCEED data struct
typedef struct CeedData_ *CeedData;
struct CeedData_ {
  CeedBasis            basis_x, basis_u, basis_p, basis_u_face;
  CeedElemRestriction  elem_restr_x, elem_restr_u, elem_restr_true_soln,
                       elem_restr_p, elem_restr_q_data;
  CeedQFunction        qf_residual, qf_jacobian, qf_error;
  CeedOperator         op_residual, op_jacobian, op_error;
  CeedVector           x_ceed, y_ceed, x_coord, q_data;
  CeedQFunctionContext pq2d_context;
};

// 1) linear2d
#ifndef PHYSICS_LINEAR2D_STRUCT
#define PHYSICS_LINEAR2D_STRUCT
typedef struct LINEAR2DContext_ *LINEAR2DContext;
struct LINEAR2DContext_ {
  CeedScalar kappa;
};
#endif

// 2) linear3d
#ifndef PHYSICS_LINEAR3D_STRUCT
#define PHYSICS_LINEAR3D_STRUCT
typedef struct LINEAR3DContext_ *LINEAR3DContext;
struct LINEAR3DContext_ {
  CeedScalar kappa;
};
#endif

// 4) richard

// Struct that contains all enums and structs used for the physics of all problems
typedef struct Physics_ *Physics;
struct Physics_ {
  LINEAR2DContext            linear2d_ctx;
  LINEAR3DContext            linear3d_ctx;
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

// Problem specific data
typedef struct ProblemData_ *ProblemData;
struct ProblemData_ {
  CeedQFunctionUser force, residual, jacobian, error,
                    setup_true, bc_pressure, setup_geo;
  const char        *force_loc, *residual_loc, *jacobian_loc,
        *error_loc, *setup_true_loc, *bc_pressure_loc, *setup_geo_loc;
  CeedQuadMode      quadrature_mode;
  CeedInt           dim, q_data_size, q_data_size_face;

};

#endif // structs_h