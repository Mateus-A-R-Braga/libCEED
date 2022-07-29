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
/// Linear mixed-elasticity problem 2D (quad element) using PETSc

#ifndef LINEAR_SYSTEM2D_H
#define LINEAR_SYSTEM2D_H

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
// e = 0.5*( \grad(u) + \grad(u)^T )
// k = E/3*(1-2*nu) is bulk modulus
// mu = G = E/2*(1+nu)
//
// Weak form: Find (u,p) \in VxQ (V=H^1, Q=L^2) on \Omega
//  (\grad(v), \sigma)   - (v, f)     = <v, t_b>_{\Gamma_N}
//  (q, \div(u))         - (q, p/k)   = 0
// This QFunction setup the mixed form of the above equation
// Inputs:
// dudX    : du/dX, u is basis and X is reference coordinate [-1,1]^dim
//   p     : basis_p at quadrature points
// q_data  : w*det_J and dX/dx (inverse of the Jacobian)
//   f     : force vector created in true qfunction
//
// Output:
//   v     :-(v,f) = -\int( v * f * w*det_J )dX
// dvdX    : (\grad(v), \sigma) = \int (\grad(v):sigma*w*det_J)dX
//   q     : (q, \div(u))-(q, p/k) = \int (q*( \div(u)-p/k )*w*det_J)dX
// -----------------------------------------------------------------------------
#ifndef PHYSICS_CTX
#define PHYSICS_CTX
typedef struct PhysicsCtx_ *PhysicsCtx;
struct PhysicsCtx_ {
  CeedScalar   nu;      // Poisson's ratio
  CeedScalar   E;       // Young's Modulus
};
#endif
// -----------------------------------------------------------------------------
// Residual evaluation for linear mixed-elasticity problem
// -----------------------------------------------------------------------------
CEED_QFUNCTION(LinearSystem2D)(void *ctx, CeedInt Q,
                               const CeedScalar *const *in,
                               CeedScalar *const *out) {
  // *INDENT-OFF*
  // Inputs
  const CeedScalar (*du_dX)[2][CEED_Q_VLA] = (const CeedScalar(*)[2][CEED_Q_VLA])in[0],
                   (*p) = (const CeedScalar(*))in[1],
                   (*q_data)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[2],
                   (*f) = in[3];

  // Outputs
  CeedScalar (*v)[CEED_Q_VLA] = (CeedScalar(*)[CEED_Q_VLA])out[0],
             (*dvdX)[2][CEED_Q_VLA] = (CeedScalar(*)[2][CEED_Q_VLA])out[1],
             (*q) = (CeedScalar(*))out[2];
  // Context
  PhysicsCtx  context = (PhysicsCtx)ctx;
  const CeedScalar E        = context->E;
  const CeedScalar nu       = context->nu;
  const CeedScalar mu       = E / (2*(1 + nu));
  const CeedScalar k        = E / (3*(1 - 2*nu)); // Bulk Modulus

  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    // *INDENT-OFF*
    const CeedScalar dudX[2][2]= {{du_dX[0][0][i], du_dX[1][0][i]},
                                  {du_dX[0][1][i], du_dX[1][1][i]}
                                 };
    // -- Qdata
    const CeedScalar wdetJ      =   q_data[0][i];
    const CeedScalar dXdx[2][2] = {{q_data[1][i],
                                    q_data[2][i]},
                                   {q_data[3][i],
                                    q_data[4][i]}
                                  };

    CeedScalar grad_u[2][2];
    AlphaMatMatMult2x2(1., dudX, dXdx, grad_u);
    // e = 0.5*( \grad(u) + \grad(u)^T )
    const CeedScalar e[2][2] = {{0.5*(grad_u[0][0] + grad_u[0][0]),
                                 0.5*(grad_u[0][1] + grad_u[1][0])},
                                {0.5*(grad_u[1][0] + grad_u[0][1]),
                                 0.5*(grad_u[1][1] + grad_u[1][1])}
                               };
    //AlphaMatPlusMatTranspose2x2(0.5, grad_u, grad_u, e);
    CeedScalar tr_e = e[0][0] + e[1][1];
    // e_d = e - tr(e)/3 * I,
    const CeedScalar e_d[2][2] = {{e[0][0]-tr_e/3., e[0][1]},
                                  {e[1][0], e[1][1]-tr_e/3.}
                                 }; 
    // \sigma = p*I + 2*mu*e_d
    const CeedScalar sigma[2][2] = {{p[i] + 2*mu*e_d[0][0], 2*mu*e_d[0][1]},
                                    {2*mu*e_d[1][0], p[i] + 2*mu*e_d[0][0]}};

    // *INDENT-ON*
    // Apply dXdx^T and weight to sigma
    for (CeedInt j = 0; j < 2; j++)     // Component
      for (CeedInt k = 0; k < 2; k++) { // Derivative
        dvdX[k][j][i] = 0;
        for (CeedInt m = 0; m < 2; m++)
          dvdX[k][j][i] += dXdx[k][m] * sigma[j][m] * wdetJ;
      }

    for (CeedInt k = 0; k < 2; k++) {
      v[k][i] = -f[i+k*Q]*wdetJ;
    }

    CeedScalar div_u = grad_u[0][0] + grad_u[1][1];
    q[i] = (div_u -p[i]/k)* wdetJ;
  } // End of Quadrature Point Loop

  return 0;
}

// -----------------------------------------------------------------------------
// Jacobian evaluation for linear mixed-elasticity problem
// -----------------------------------------------------------------------------
CEED_QFUNCTION(JacobianLinearSystem2D)(void *ctx, CeedInt Q,
                                       const CeedScalar *const *in,
                                       CeedScalar *const *out) {
  // *INDENT-OFF*
  // Inputs
  const CeedScalar (*ddu_dX)[2][CEED_Q_VLA] = (const CeedScalar(*)[2][CEED_Q_VLA])in[0],
                   (*dp) = (const CeedScalar(*))in[1],
                   (*q_data)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[2];
  // Outputs
  CeedScalar (*dv)[CEED_Q_VLA] = (CeedScalar(*)[CEED_Q_VLA])out[0],
             (*ddvdX)[2][CEED_Q_VLA] = (CeedScalar(*)[2][CEED_Q_VLA])out[1],
             (*dq) = (CeedScalar(*))out[2];
  // Context
  PhysicsCtx  context = (PhysicsCtx)ctx;
  const CeedScalar E        = context->E;
  const CeedScalar nu       = context->nu;
  const CeedScalar mu       = E / (2*(1 + nu));
  const CeedScalar k        = E / (3*(1 - 2*nu)); // Bulk Modulus

  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    // *INDENT-OFF*
    const CeedScalar ddudX[2][2]= {{ddu_dX[0][0][i], ddu_dX[1][0][i]},
                                   {ddu_dX[0][1][i], ddu_dX[1][1][i]}
                                  };
    // -- Qdata
    const CeedScalar wdetJ      =   q_data[0][i];
    const CeedScalar dXdx[2][2] = {{q_data[1][i],
                                    q_data[2][i]},
                                   {q_data[3][i],
                                    q_data[4][i]}
                                  };

    CeedScalar grad_du[2][2];
    AlphaMatMatMult2x2(1., ddudX, dXdx, grad_du);
    // e = 0.5*( \grad(u) + \grad(u)^T )
    const CeedScalar de[2][2] = {{0.5*(grad_du[0][0] + grad_du[0][0]),
                                  0.5*(grad_du[0][1] + grad_du[1][0])},
                                 {0.5*(grad_du[1][0] + grad_du[0][1]),
                                  0.5*(grad_du[1][1] + grad_du[1][1])}
                                };
    //AlphaMatPlusMatTranspose2x2(0.5, grad_u, grad_u, e);
    CeedScalar tr_de = de[0][0] + de[1][1];
    // e_d = e - tr(e)/3 * I,
    const CeedScalar de_d[2][2] = {{de[0][0]-tr_de/3., de[0][1]},
                                   {de[1][0], de[1][1]-tr_de/3.}
                                  }; 
    // \sigma = p*I + 2*mu*e_d
    const CeedScalar dsigma[2][2] = {{dp[i] + 2*mu*de_d[0][0], 2*mu*de_d[0][1]},
                                     {2*mu*de_d[1][0], dp[i] + 2*mu*de_d[0][0]}};

    // *INDENT-ON*
    // Apply dXdx^T and weight to sigma
    for (CeedInt j = 0; j < 2; j++)     // Component
      for (CeedInt k = 0; k < 2; k++) { // Derivative
        ddvdX[k][j][i] = 0;
        for (CeedInt m = 0; m < 2; m++)
          ddvdX[k][j][i] += dXdx[k][m] * dsigma[j][m] * wdetJ;
      }

    for (CeedInt k = 0; k < 2; k++) {
      dv[k][i] = 0.;
    }

    CeedScalar div_du = grad_du[0][0] + grad_du[1][1];
    dq[i] = (div_du -dp[i]/k)* wdetJ;
  } // End of Quadrature Point Loop

  return 0;
}

// -----------------------------------------------------------------------------

#endif //End of LINEAR_SYSTEM2D_H
