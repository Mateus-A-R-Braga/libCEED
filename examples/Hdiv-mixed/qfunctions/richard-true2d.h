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
/// Force of Richard problem 2D (quad element) using PETSc

#ifndef RICHARD_TRUE2D_H
#define RICHARD_TRUE2D_H

#include <math.h>
#include "utils.h"

// See Matthew Farthing, Christopher Kees, Cass Miller (2003)
// https://www.sciencedirect.com/science/article/pii/S0309170802001872
// -----------------------------------------------------------------------------
// Strong form:
//  u        = -rho*k_r*K *[grad(\psi) - rho*g_u]   in \Omega x [0,T]
//  -\div(u) = -f  + d (rho*theta)/dt              in \Omega x [0,T]
//  p        = p_b                                  on \Gamma_D x [0,T]
//  u.n      = u_b                                  on \Gamma_N x [0,T]
//  p        = p_0                                  in \Omega, t = 0
//
//  Where rho = rho_a/rho_a0, rho_a = rho_a0*exp(\beta * (p - p0)), p0 = 101325 Pa is atmospheric pressure
//  rho_a0 is the density at p_0, g_u = g/norm(g) where g is gravity.
//  k_r = b_a + alpha_a * (\psi - x2), where \psi = p / (rho_a0 * norm(g)) and x2 is vertical axis
//
// Weak form: Find (u, \psi) \in VxQ (V=H(div), Q=L^2) on \Omega
//  (v, K^{-1}/rho*k_r * u) -(v, rho*g_u) -(\div(v), \psi) = -<v, p_b*n>_{\Gamma_D}
// -(q, \div(u))  + (q, f) -(q, d (rho*\theta)/dt ) = 0
//
// We solve MMS for  K = kappa*I and beta=0 ==> rho=1 and \theta = alpha_a*\psi, so
// -(q, d (rho*\theta)/dt ) = -alpha_a*(q, d(\psi)/dt )
//
// This QFunction setup the true solution and forcing f of the above equation
// Inputs:
//   coords: physical coordinate
//
// Output:
//   true_force     : = div(u) + d (rho*theta)/dt
//   true_solution  : = [\psi, u] where \psi, u are the exact solution solution
// -----------------------------------------------------------------------------
// We have 3 experiment parameters as described in Table 1:P1, P2, P3
// Matthew Farthing, Christopher Kees, Cass Miller (2003)
// https://www.sciencedirect.com/science/article/pii/S0309170802001872
#ifndef RICHARD_CTX
#define RICHARD_CTX
typedef struct RICHARDContext_ *RICHARDContext;
struct RICHARDContext_ {
  CeedScalar kappa;
  CeedScalar g;
  CeedScalar rho_a0;
  CeedScalar alpha_a, b_a;
  CeedScalar beta, p0;
  CeedScalar t;
  CeedScalar gamma;
};
#endif
// -----------------------------------------------------------------------------
// Initial conditions for Richard problem
// -----------------------------------------------------------------------------
CEED_QFUNCTION(RichardICs2D)(void *ctx, const CeedInt Q,
                             const CeedScalar *const *in,
                             CeedScalar *const *out) {
  // *INDENT-OFF*
  // Inputs
  const CeedScalar (*coords) = in[0];
  // Outputs
  CeedScalar (*u_0) = out[0], (*psi_0) = out[1];
  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    // Setup, (x,y) and J = dx/dX
    CeedScalar x = coords[i+0*Q], y = coords[i+1*Q];

    // 1st eq: component 1
    u_0[i+0*Q] = 0.;
    // 1st eq: component 2
    u_0[i+1*Q] = 0.;
    // 2nd eq
    psi_0[i] = sin(PI_DOUBLE*x)*sin(PI_DOUBLE*y);
  } // End of Quadrature Point Loop
  return 0;
}

// -----------------------------------------------------------------------------
// True solution for Richard problem
// -----------------------------------------------------------------------------
CEED_QFUNCTION(RichardTrue2D)(void *ctx, const CeedInt Q,
                              const CeedScalar *const *in,
                              CeedScalar *const *out) {
  // *INDENT-OFF*
  // Inputs
  const CeedScalar (*coords) = in[0];
  // Outputs
  CeedScalar (*true_force) = out[0], (*true_solution) = out[1];
  // Context
  RICHARDContext  context = (RICHARDContext)ctx;
  const CeedScalar kappa    = context->kappa;
  const CeedScalar alpha_a  = context->alpha_a;
  const CeedScalar b_a      = context->b_a;
  const CeedScalar gamma   = context->gamma;
  CeedScalar t             = context->t;
  printf("time %f \n", t);
  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    CeedScalar x = coords[i+0*Q], y = coords[i+1*Q];
    CeedScalar psi    = exp(-gamma*t)*sin(PI_DOUBLE*x)*sin(PI_DOUBLE*y);
    CeedScalar psi_x  = PI_DOUBLE*exp(-gamma*t)*cos(PI_DOUBLE*x)*sin(PI_DOUBLE*y);
    CeedScalar psi_y  = PI_DOUBLE*exp(-gamma*t)*sin(PI_DOUBLE*x)*cos(PI_DOUBLE*y);

    // k_r = b_a + alpha_a * (1 - x*y)
    CeedScalar k_r = b_a + alpha_a*(1-x*y);
    // rho = rho_a/rho_a0
    CeedScalar rho = 1.;
    // u = -rho*k_r*K *[grad(\psi) - rho*g_u]
    CeedScalar u[2] = {-rho*k_r*kappa*psi_x, -rho*k_r*kappa*(psi_y-1)};
    CeedScalar div_u = -rho*kappa*(-alpha_a*y*psi_x - k_r*PI_DOUBLE*PI_DOUBLE*psi
                                   -alpha_a*x*psi_y - k_r*PI_DOUBLE*PI_DOUBLE*psi);

    // True Force: f = \div(u) + d (rho*theta)/dt
    true_force[i+0*Q] = div_u -alpha_a*gamma*psi;
    // True Solution
    true_solution[i+0*Q] = psi;
    true_solution[i+1*Q] = u[0];
    true_solution[i+2*Q] = u[1];
  } // End of Quadrature Point Loop
  return 0;
}
// -----------------------------------------------------------------------------

#endif //End of RICHARD_TRUE2D_H
