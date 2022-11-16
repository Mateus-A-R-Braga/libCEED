#ifndef reynolds_stress_h
#define reynolds_stress_h

#include <ceed.h>
#include <math.h>
#include <stdlib.h>
#include "newtonian_state.h"
#include "newtonian_types.h"
#include "stabilization.h"
#include "utils.h"

// *****************************************************************************
// This QFunction calcualtes the Reynolds stress to be L2 projected back onto
// the FE space.
//
// Note, the mean of the product <U_iU_j> is what is projected and the Reynolds
// stress is computed from <u'_iu'j> = <U_iU_j> - <U_i><U_j>
// *****************************************************************************









CEED_QFUNCTION_HELPER int DivDiffusiveFlux(void *ctx, CeedInt Q,
    const CeedScalar *const *in, CeedScalar *const *out,
    StateFromQi_t StateFromQi, StateFromQi_fwd_t StateFromQi_fwd) {
  // *INDENT-OFF*
  // Inputs
  const CeedScalar (*q)[CEED_Q_VLA]         = (const CeedScalar(*)[CEED_Q_VLA])in[0],
                   (*Grad_q)[5][CEED_Q_VLA] = (const CeedScalar(*)[5][CEED_Q_VLA])in[1],
                   (*q_data)[CEED_Q_VLA]    = (const CeedScalar(*)[CEED_Q_VLA])in[2],
                   (*x)[CEED_Q_VLA]         = (const CeedScalar(*)[CEED_Q_VLA])in[3];
  // Outputs
  CeedScalar (*Grad_v)[5][CEED_Q_VLA] = (CeedScalar(*)[5][CEED_Q_VLA])out[0],
             (*jac_data)[CEED_Q_VLA]  = (CeedScalar(*)[CEED_Q_VLA])out[1];
  // *INDENT-ON*
  // Context
  const NewtonianIdealGasContext context        = (NewtonianIdealGasContext)ctx;
  const CeedScalar dt                           = context->dt;
  const StateConservative ZeroFlux              = {0., {0., 0., 0.}, 0.};
  const StateConservative ZeroInviscidFluxes[3] = {ZeroFlux, ZeroFlux, ZeroFlux};

  CeedPragmaSIMD
  // Quadrature Point Loop
  for(CeedInt i=0; i<Q; i++) {
    const CeedScalar x_i[3] = {x[0][i], x[1][i], x[2][i]};
    const CeedScalar qi[5]  = {q[0][i], q[1][i], q[2][i], q[3][i], q[4][i]};
    const State s = StateFromQi(context, qi, x_i);

    // -- Interp-to-Interp q_data
    const CeedScalar wdetJ      =   q_data[0][i];
    // -- Interp-to-Grad q_data
    // ---- Inverse of change of coordinate matrix: X_i,j
    // *INDENT-OFF*
    const CeedScalar dXdx[3][3] = {{q_data[1][i], q_data[2][i], q_data[3][i]},
                                   {q_data[4][i], q_data[5][i], q_data[6][i]},
                                   {q_data[7][i], q_data[8][i], q_data[9][i]}
                                  };
    // *INDENT-ON*
    State grad_s[3];
    for (CeedInt j=0; j<3; j++) {
      CeedScalar dx_i[3] = {0}, dqi[5];
      for (CeedInt k=0; k<5; k++)
        dqi[k] = Grad_q[0][k][i] * dXdx[0][j] +
                 Grad_q[1][k][i] * dXdx[1][j] +
                 Grad_q[2][k][i] * dXdx[2][j];
      dx_i[j] = 1.;
      grad_s[j] = StateFromQi_fwd(context, s, dqi, x_i, dx_i);
    }

    CeedScalar strain_rate[6], kmstress[6], stress[3][3], Fe[3];
    KMStrainRate(grad_s, strain_rate);
    NewtonianStress(context, strain_rate, kmstress);
    KMUnpack(kmstress, stress);
    ViscousEnergyFlux(context, s.Y, grad_s, stress, Fe);

    // Total flux
    CeedScalar DiffFlux[5][3];
    FluxTotal(ZeroInviscidFluxes, stress, Fe, DiffFlux);

    for (CeedInt j=0; j<3; j++) {
      for (CeedInt k=0; k<5; k++) {
        Grad_v[j][k][i] = -wdetJ * (dXdx[j][0] * DiffFlux[k][0] +
                                    dXdx[j][1] * DiffFlux[k][1] +
                                    dXdx[j][2] * DiffFlux[k][2]);
      }
    }

    CeedScalar Tau_d[5];
    Tau_diagPrim(context, s, dXdx, dt, Tau_d);

    for (CeedInt j=0; j<5; j++) jac_data[j][i]     = qi[j];
    for (CeedInt j=0; j<6; j++) jac_data[5+j][i]   = kmstress[j];
    for (CeedInt j=0; j<3; j++) jac_data[5+6+j][i] = Tau_d[j];
  } // End Quadrature Point Loop

  // Return
  return 0;
}

CEED_QFUNCTION(DivDiffusiveFlux_Conserv)(void *ctx, CeedInt Q,
    const CeedScalar *const *in, CeedScalar *const *out) {
  return DivDiffusiveFlux(ctx, Q, in, out, StateFromU, StateFromU_fwd);
}

CEED_QFUNCTION(DivDiffusiveFlux_Prim)(void *ctx, CeedInt Q,
                                      const CeedScalar *const *in, CeedScalar *const *out) {
  return DivDiffusiveFlux(ctx, Q, in, out, StateFromY, StateFromY_fwd);
}

#endif // reynolds_stress_h
