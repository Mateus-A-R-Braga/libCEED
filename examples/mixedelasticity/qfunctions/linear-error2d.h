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
/// Compute pointwise error of the mixed-elasticity example using PETSc

#ifndef LINEAR_ERROR2D_H
#define LINEAR_ERROR2D_H

#include <math.h>

// -----------------------------------------------------------------------------
// Compuet error
// -----------------------------------------------------------------------------
CEED_QFUNCTION(LinearError2D)(void *ctx, const CeedInt Q,
                              const CeedScalar *const *in,
                              CeedScalar *const *out) {
  // *INDENT-OFF*
  // Inputs
  const CeedScalar (*q_data)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[0],
                   (*u)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[1],
                   (*p) = (const CeedScalar(*))in[2],
                   (*target) = in[3];
  // Outputs
  CeedScalar (*error) = out[0];
  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {          
    const CeedScalar wdetJ      =   q_data[0][i];
    // Error
    error[i+0*Q] = (p[i] - target[i+0*Q])*(p[i] - target[i+0*Q])*wdetJ;
    error[i+1*Q] = (u[0][i] - target[i+1*Q])*(u[0][i] - target[i+1*Q])*wdetJ;
    error[i+2*Q] = (u[1][i] - target[i+2*Q])*(u[1][i] - target[i+2*Q])*wdetJ;
  } // End of Quadrature Point Loop

  return 0;
}
// -----------------------------------------------------------------------------

#endif // End LINEAR_ERROR2D_H
