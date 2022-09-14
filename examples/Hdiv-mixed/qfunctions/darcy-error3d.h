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
/// Compute pointwise error of the H(div) example using PETSc

#ifndef DARCY_ERROR3D_H
#define DARCY_ERROR3D_H

#include <math.h>
#include <ceed.h>
#include "utils.h"


// -----------------------------------------------------------------------------
// Compuet error
// -----------------------------------------------------------------------------
#ifndef DARCY_CTX
#define DARCY_CTX
typedef struct DARCYContext_ *DARCYContext;
struct DARCYContext_ {
  CeedScalar kappa;
  CeedScalar g;
  CeedScalar rho_a0;
  CeedScalar alpha_a, b_a;
  CeedScalar lx, ly, lz;
};
#endif
CEED_QFUNCTION(DarcyError3D)(void *ctx, const CeedInt Q,
                             const CeedScalar *const *in,
                             CeedScalar *const *out) {
  // *INDENT-OFF*
  // Inputs
  const CeedScalar (*w) = in[0],
                   (*dxdX)[3][CEED_Q_VLA] = (const CeedScalar(*)[3][CEED_Q_VLA])in[1],
                   (*u)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[2],
                   (*p) = (const CeedScalar(*))in[3],
                   (*true_soln) = in[4];
  // Outputs
  CeedScalar (*error) = out[0];
  // Context
  DARCYContext  context = (DARCYContext)ctx;
  //const CeedScalar    kappa   = context->kappa;
  const CeedScalar rho_a0   = context->rho_a0;
  const CeedScalar g        = context->g;
  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    // Setup, J = dx/dX
    const CeedScalar J[3][3] = {{dxdX[0][0][i], dxdX[1][0][i], dxdX[2][0][i]},
                                {dxdX[0][1][i], dxdX[1][1][i], dxdX[2][1][i]},
                                {dxdX[0][2][i], dxdX[1][2][i], dxdX[2][2][i]}};
    const CeedScalar det_J = MatDet3x3(J);              
    // Compute Piola map:uh = J*u/detJ
    CeedScalar u1[3] = {u[0][i], u[1][i], u[2][i]}, uh[3];
    AlphaMatVecMult3x3(1/det_J, J, u1, uh);

    // Error
    CeedScalar psi = p[i] / (rho_a0 * g);
    error[i+0*Q] = (psi - true_soln[i+0*Q])*(psi - true_soln[i+0*Q])*w[i]*det_J;
    error[i+1*Q] = (uh[0] - true_soln[i+1*Q])*(uh[0] - true_soln[i+1*Q])*w[i]*det_J;
    error[i+2*Q] = (uh[1] - true_soln[i+2*Q])*(uh[1] - true_soln[i+2*Q])*w[i]*det_J;
    error[i+3*Q] = (uh[2] - true_soln[i+3*Q])*(uh[2] - true_soln[i+3*Q])*w[i]*det_J;
  } // End of Quadrature Point Loop

  return 0;
}
// -----------------------------------------------------------------------------

#endif // End DARCY_ERROR3D_H