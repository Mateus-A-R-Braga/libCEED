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
  CeedScalar t, t_final;
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
  const CeedScalar (*coords) = in[0],
                   (*dxdX)[2][CEED_Q_VLA] = (const CeedScalar(*)[2][CEED_Q_VLA])in[1];
  // Outputs
  CeedScalar (*u_0) = out[0], (*psi_0) = out[1];
  // Quadrature Point Loop
  {
  // Setup, J = dx/dX
  const CeedScalar J0[2][2] = {{dxdX[0][0][0], dxdX[1][0][0]},
                               {dxdX[0][1][0], dxdX[1][1][0]}};
  const CeedScalar J1[2][2] = {{dxdX[0][0][1], dxdX[1][0][1]},
                               {dxdX[0][1][1], dxdX[1][1][1]}};
  const CeedScalar J2[2][2] = {{dxdX[0][0][2], dxdX[1][0][2]},
                               {dxdX[0][1][2], dxdX[1][1][2]}};
  const CeedScalar J3[2][2] = {{dxdX[0][0][3], dxdX[1][0][3]},
                               {dxdX[0][1][3], dxdX[1][1][3]}};
  CeedScalar x0 = coords[0+0*Q], y0 = coords[0+1*Q];
  CeedScalar x1 = coords[1+0*Q], y1 = coords[1+1*Q];
  CeedScalar x2 = coords[2+0*Q], y2 = coords[2+1*Q];
  CeedScalar x3 = coords[3+0*Q], y3 = coords[3+1*Q];
  CeedScalar xc = (x0+x1)/2., yc = (y0+y2)/2.;
  CeedScalar ue0[2] = {x0-y0, x0+y0};
  CeedScalar ue1[2] = {x1-y1, x1+y1};
  CeedScalar ue2[2] = {x2-y2, x2+y2};
  CeedScalar ue3[2] = {x3-y3, x3+y3};
  CeedScalar nl0[2] = {-J0[1][1],J0[0][1]};
  CeedScalar nb0[2] = {J0[1][0],-J0[0][0]};
  CeedScalar nr1[2] = {J1[1][1],-J1[0][1]};
  CeedScalar nb1[2] = {J1[1][0],-J1[0][0]};
  CeedScalar nl2[2] = {-J2[1][1],J2[0][1]};
  CeedScalar nt2[2] = {-J2[1][0],J2[0][0]};
  CeedScalar nr3[2] = {J3[1][1],-J3[0][1]};
  CeedScalar nt3[2] = {-J3[1][0],J3[0][0]};
  CeedScalar d0, d1, d2, d3, d4, d5, d6, d7;
  d0 = ue0[0]*nb0[0]+ue0[1]*nb0[1];
  d1 = ue1[0]*nb1[0]+ue1[1]*nb1[1];
  d2 = ue1[0]*nr1[0]+ue1[1]*nr1[1];
  d3 = ue3[0]*nr3[0]+ue3[1]*nr3[1];
  d4 = ue2[0]*nt2[0]+ue2[1]*nt2[1];
  d5 = ue3[0]*nt3[0]+ue3[1]*nt3[1];
  d6 = ue0[0]*nl0[0]+ue0[1]*nl0[1];
  d7 = ue2[0]*nl2[0]+ue2[1]*nl2[1];
  u_0[0] = d0;
  u_0[1] = d1;
  u_0[2] = d2;
  u_0[3] = d3;
  u_0[4] = d4;
  u_0[5] = d5;
  u_0[6] = d6;
  u_0[7] = d7;
  // 2nd eq
  psi_0[0] = sin(PI_DOUBLE*xc)*sin(PI_DOUBLE*yc);
  psi_0[1] = sin(PI_DOUBLE*xc)*sin(PI_DOUBLE*yc);
  psi_0[2] = sin(PI_DOUBLE*xc)*sin(PI_DOUBLE*yc);
  psi_0[3] = sin(PI_DOUBLE*xc)*sin(PI_DOUBLE*yc);
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
  const CeedScalar gamma    = context->gamma;
  CeedScalar t_final        = context->t_final;
  //printf("time final in True Qfunction %f\n", t_final);
  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    CeedScalar x = coords[i+0*Q], y = coords[i+1*Q];
    // psi = exp(-gamma*t)*sin(pi*x)*sin(pi*y)
    CeedScalar psi1 = sin(PI_DOUBLE*x)*sin(PI_DOUBLE*y);
    CeedScalar psi    = exp(-gamma*t_final)*psi1;
    CeedScalar psi1_x = PI_DOUBLE*cos(PI_DOUBLE*x)*sin(PI_DOUBLE*y);
    CeedScalar psi_x  = exp(-gamma*t_final)*psi1_x;
    CeedScalar psi1_y = PI_DOUBLE*sin(PI_DOUBLE*x)*cos(PI_DOUBLE*y);
    CeedScalar psi_y  = exp(-gamma*t_final)*psi1_y;

    // k_r = b_a + alpha_a * (1 - x*y)
    CeedScalar k_r = b_a + alpha_a*(1-x*y);
    // rho = rho_a/rho_a0
    CeedScalar rho = 1.;
    // u = -rho*k_r*K *[grad(\psi) - rho*g_u]
    CeedScalar u[2] = {-rho*k_r*kappa*psi_x, -rho*k_r*kappa*(psi_y-1)};
    // we factor out exp() term and applied in residual function richard-system2d.h
    CeedScalar div_u = -rho*kappa*(-alpha_a*y*psi1_x - k_r*PI_DOUBLE*PI_DOUBLE*psi1
                                   -alpha_a*x*psi1_y - k_r*PI_DOUBLE*PI_DOUBLE*psi1);

    // True Force: f = \div(u) + d (rho*theta)/dt
    true_force[i+0*Q] = div_u -alpha_a*gamma*psi1;
    // True Solution
    true_solution[i+0*Q] = psi;
    true_solution[i+1*Q] = u[0];
    true_solution[i+2*Q] = u[1];
  } // End of Quadrature Point Loop
  return 0;
}
// -----------------------------------------------------------------------------

#endif //End of RICHARD_TRUE2D_H
