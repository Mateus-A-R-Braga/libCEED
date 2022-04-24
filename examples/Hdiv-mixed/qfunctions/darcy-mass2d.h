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
/// Darcy problem 2D (quad element) using PETSc

#ifndef DARCY_MASS2D_H
#define DARCY_MASS2D_H

#include <math.h>

// -----------------------------------------------------------------------------
// Strong form:
//  u       = -\grad(p)
//  \div(u) = f
// Weak form: Find (u,p) \in VxQ (V=H(div), Q=L^2) on \Omega
//  (u, v) - (p, \div(v)) = -<p, v\cdot n>
// -(q, \div(u))          = -(q, f)
// This QFunction setup the mixed form of the above equation
// Inputs:
//   w     : weight of quadrature
//   J     : dx/dX. x physical coordinate, X reference coordinate [-1,1]^dim
//   u     : basis_u at quadrature points
// div(u)  : divergence of basis_u at quadrature points
//   p     : basis_p at quadrature points
//
// Output:
//   v     : (v,u) = \int (v^T * u detJ*w) ==> \int (v^T J^T*J*u*w/detJ)
// div(v)  : -(\div(v), p) = -\int (div(v)^T * p *w)
//   q     : -(q, \div(u)) = -\int (q^T * div(u) *w)
// which create the following coupled system
//                            D = [ M  B^T ]
//                                [ B   0  ]
// M = (v,u), B = -(q, \div(u))
// Note we need to apply Piola map on the basis_u, which is J*u/detJ
// So (v,u) = \int (v^T * u detJ*w) ==> \int (v^T J^T*J*u*w/detJ)
// -----------------------------------------------------------------------------
CEED_QFUNCTION(SetupDarcyMass2D)(void *ctx, CeedInt Q,
                                 const CeedScalar *const *in,
                                 CeedScalar *const *out) {
  // *INDENT-OFF*
  // Inputs
  const CeedScalar (*w) = in[0],
                   (*dxdX)[2][CEED_Q_VLA] = (const CeedScalar(*)[2][CEED_Q_VLA])in[1],
                   (*u)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[2],
                   (*div_u) = (const CeedScalar(*))in[3],
                   (*p) = (const CeedScalar(*))in[4];

  // Outputs
  CeedScalar (*v)[CEED_Q_VLA] = (CeedScalar(*)[CEED_Q_VLA])out[0],
             (*div_v) = (CeedScalar(*))out[1],
             (*q) = (CeedScalar(*))out[2];
  // *INDENT-ON*

  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    // *INDENT-OFF*
    // Setup, J = dx/dX
    const CeedScalar J[2][2] = {{dxdX[0][0][i], dxdX[1][0][i]},
                                {dxdX[0][1][i], dxdX[1][1][i]}};
    const CeedScalar detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];

    // *INDENT-ON*
    // Piola map: J^T*J*u*w/detJ
    // 1) Compute J^T * J
    CeedScalar JTJ[2][2];
    for (CeedInt j = 0; j < 2; j++) {
      for (CeedInt k = 0; k < 2; k++) {
        JTJ[j][k] = 0;
        for (CeedInt m = 0; m < 2; m++)
          JTJ[j][k] += J[m][j] * J[m][k];
      }
    }
    // 2) Compute J^T*J*u * w /detJ
    for (CeedInt k = 0; k < 2; k++) {
      v[k][i] = 0;
      for (CeedInt m = 0; m < 2; m++)
        v[k][i] += JTJ[k][m] * u[m][i] * w[i]/detJ;
    }

    div_v[i] = -p[i] * w[i];
    q[i] = -div_u[i] * w[i];
  } // End of Quadrature Point Loop

  return 0;
}

// -----------------------------------------------------------------------------
// Jacobian evaluation for Darcy problem
// -----------------------------------------------------------------------------
CEED_QFUNCTION(SetupJacobianDarcyMass2D)(void *ctx, CeedInt Q,
    const CeedScalar *const *in,
    CeedScalar *const *out) {
  // *INDENT-OFF*
  // Inputs
  const CeedScalar (*w) = in[0],
                   (*dxdX)[2][CEED_Q_VLA] = (const CeedScalar(*)[2][CEED_Q_VLA])in[1],
                   (*du)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[2],
                   (*div_du) = (const CeedScalar(*))in[3],
                   (*dp) = (const CeedScalar(*))in[4];

  // Outputs
  CeedScalar (*dv)[CEED_Q_VLA] = (CeedScalar(*)[CEED_Q_VLA])out[0],
             (*div_dv) = (CeedScalar(*))out[1],
             (*dq) = (CeedScalar(*))out[2];

  // *INDENT-ON*

  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    // *INDENT-OFF*
    // Setup, J = dx/dX
    const CeedScalar J[2][2] = {{dxdX[0][0][i], dxdX[1][0][i]},
                                {dxdX[0][1][i], dxdX[1][1][i]}};
    const CeedScalar detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];

    // *INDENT-ON*
    // Piola map: J^T*J*du*w/detJ
    // 1) Compute J^T * J
    CeedScalar JTJ[2][2];
    for (CeedInt j = 0; j < 2; j++) {
      for (CeedInt k = 0; k < 2; k++) {
        JTJ[j][k] = 0;
        for (CeedInt m = 0; m < 2; m++)
          JTJ[j][k] += J[m][j] * J[m][k];
      }
    }
    // 2) Compute J^T*J*du * w /detJ
    for (CeedInt k = 0; k < 2; k++) {
      dv[k][i] = 0;
      for (CeedInt m = 0; m < 2; m++)
        dv[k][i] += JTJ[k][m] * du[m][i] * w[i]/detJ;
    }

    div_dv[i] = -dp[i] * w[i];
    dq[i] = -div_du[i] * w[i];
  } // End of Quadrature Point Loop

  return 0;
}

// -----------------------------------------------------------------------------

#endif //End of DARCY_MASS2D_H