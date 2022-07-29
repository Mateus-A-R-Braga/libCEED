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

/// @file
/// Force of linear mixed-elasticity problem 2D (quad element) using PETSc

#ifndef LINEAR_TRUE2D_H
#define LINEAR_TRUE2D_H

#include <math.h>
#include "utils.h"

// -----------------------------------------------------------------------------
// Strong form:
//  \div(\sigma) + f    = 0      in \Omega
//  \div(u )     - p/k  = 0      in \Omega
//  u                   = u_b    on \Gamma_D
//  \sigma.n            = t_b    on \Gamma_N
//
// Constitutive Equations:
//  \sigma = p*I + 2*mu*e_d
// where
// e_d = e - tr(e)/3 * I,
// p = k*tr(e) = k*\div(u)
// e = 0.5*( \grad(u) + \grad(u)^T )
// k = E/3*(1-2*nu) is bulk modulus
// mu = G = E/2*(1+nu)
//
// This QFunction setup the true solution and forcing f of the above equation
// Inputs:
//   coords: physical coordinate
//
// Output:
//   true_force     : = -\div(\sigma)
//   true_solution  : = [p, u] where p, u are the exact solution solution
// -----------------------------------------------------------------------------
#ifndef PHYSICS_CTX
#define PHYSICS_CTX
typedef struct PhysicsCtx_ *PhysicsCtx;
struct PhysicsCtx_ {
  CeedScalar   nu;      // Poisson's ratio
  CeedScalar   E;       // Young's Modulus
};
#endif
CEED_QFUNCTION(LinearTrue2D)(void *ctx, const CeedInt Q,
                             const CeedScalar *const *in,
                             CeedScalar *const *out) {
  // *INDENT-OFF*
  // Inputs
  const CeedScalar (*coords) = in[0]; 
  // Outputs
  CeedScalar (*true_force) = out[0], (*true_solution) = out[1];
  // Context
  PhysicsCtx  context = (PhysicsCtx)ctx;
  const CeedScalar E        = context->E;
  const CeedScalar nu       = context->nu;
  const CeedScalar mu       = E / (2*(1 + nu));
  const CeedScalar k        = E / (3*(1 - 2*nu)); // Bulk Modulus
  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    CeedScalar x = coords[i+0*Q], y = coords[i+1*Q];  
    CeedScalar u1     = exp(2*x) * sin(3*y);
    CeedScalar u1_1   = 2*u1;
    CeedScalar u1_12  = 6*exp(2*x) * cos(3*y);
    CeedScalar u1_11  = 4*u1;
    //CeedScalar u1_2   = 3*exp(2*x) * cos(3*y);
    CeedScalar u1_21  = 6*exp(2*x) * cos(3*y);
    CeedScalar u1_22  = -9*u1;

    CeedScalar u2     = exp(3*y) * sin(4*x);
    //CeedScalar u2_1   = 4*exp(3*y) * cos(4*x);
    CeedScalar u2_12  = 12*exp(3*y) * cos(4*x);
    CeedScalar u2_11  = -16*u2;
    CeedScalar u2_2   = 3*u2;
    CeedScalar u2_21  = 12*exp(3*y) * cos(4*x);
    CeedScalar u2_22  = 9*u2;

    CeedScalar p  = k*(u1_1 + u2_2); 
    CeedScalar f1 = -(2*mu/3)*(2*u1_11 - u2_21) - k*(u1_11 + u2_21) - mu*(u1_22 + u2_12);
    CeedScalar f2 = -(2*mu/3)*(2*u2_22 - u1_12) - k*(u1_12 + u2_22) - mu*(u1_21 + u2_11);
    // True Force: f = \div(u)
    true_force[i+0*Q] = f1;
    true_force[i+1*Q] = f2;
    // True Solution
    true_solution[i+0*Q] = p;
    true_solution[i+1*Q] = u1;
    true_solution[i+2*Q] = u2;
  } // End of Quadrature Point Loop
  return 0;
}
// -----------------------------------------------------------------------------

#endif //End of LINEAR_TRUE2D_H
