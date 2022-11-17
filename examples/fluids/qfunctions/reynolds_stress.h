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

CEED_QFUNCTION_HELPER int ReynoldsStress(void *ctx, CeedInt Q,
    const CeedScalar *const *in, CeedScalar *const *out,
    StateFromQi_t StateFromQi) {

  // *INDENT-OFF*
  // Inputs
  const CeedScalar (*q)[CEED_Q_VLA]         = (const CeedScalar(*)[CEED_Q_VLA])in[0],
                   (*q_data)[CEED_Q_VLA]    = (const CeedScalar(*)[CEED_Q_VLA])in[1],
                   (*x)[CEED_Q_VLA]         = (const CeedScalar(*)[CEED_Q_VLA])in[2];

  // Outputs
  CeedScalar (*U_prod)[CEED_Q_VLA] = (CeedScalar(*)[CEED_Q_VLA])out[0];

  // CEED_Q_VLA == number of quadrature points

  // *INDENT-ON*
  // Context
  const NewtonianIdealGasContext context        = (NewtonianIdealGasContext)ctx;
  const CeedScalar dt                           = context->dt;
//  const StateConservative ZeroFlux              = {0., {0., 0., 0.}, 0.};
//  const StateConservative ZeroInviscidFluxes[3] = {ZeroFlux, ZeroFlux, ZeroFlux};

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
//    const CeedScalar dXdx[3][3] = {{q_data[1][i], q_data[2][i], q_data[3][i]},
//                                   {q_data[4][i], q_data[5][i], q_data[6][i]},
//                                   {q_data[7][i], q_data[8][i], q_data[9][i]}
//                                  };

   // CeedScalar U_prod[6];
   // Using Voight notation (array ordering)
    U_prod[0][i] = s.Y.velocity[0] * qi[1] * wdetJ // U*U
    U_prod[1][i] = s.Y.velocity[1] * qi[2] * wdetJ // V*V
    U_prod[2][i] = qi[3] * qi[3] * wdetJ // W*W
    U_prod[3][i] = qi[2] * qi[3] * wdetJ // V*W
    U_prod[4][i] = qi[1] * qi[3] * wdetJ // U*W
    U_prod[5][i] = qi[1] * qi[2] * wdetJ // U*V
    
   // Where to multiply wdetJ ?



    // *INDENT-ON*
//    State grad_s[3];
//    for (CeedInt j=0; j<3; j++) {
//      CeedScalar dx_i[3] = {0}, dqi[5];
//      for (CeedInt k=0; k<5; k++)
//        dqi[k] = Grad_q[0][k][i] * dXdx[0][j] +
//                 Grad_q[1][k][i] * dXdx[1][j] +
//                 Grad_q[2][k][i] * dXdx[2][j];
//      dx_i[j] = 1.;
//      grad_s[j] = StateFromQi_fwd(context, s, dqi, x_i, dx_i);
//    }

//    CeedScalar strain_rate[6], kmstress[6], stress[3][3], Fe[3];
//    KMStrainRate(grad_s, strain_rate);
//    NewtonianStress(context, strain_rate, kmstress);
//    KMUnpack(kmstress, stress);
//    ViscousEnergyFlux(context, s.Y, grad_s, stress, Fe);

    // Total flux
//    CeedScalar DiffFlux[5][3];
//    FluxTotal(ZeroInviscidFluxes, stress, Fe, DiffFlux);

//    for (CeedInt j=0; j<3; j++) {
//      for (CeedInt k=0; k<5; k++) {
//        Grad_v[j][k][i] = -wdetJ * (dXdx[j][0] * DiffFlux[k][0] +
//                                    dXdx[j][1] * DiffFlux[k][1] +
//                                    dXdx[j][2] * DiffFlux[k][2]);
//      }
//    }

//    CeedScalar Tau_d[5];
//    Tau_diagPrim(context, s, dXdx, dt, Tau_d);






  } // End Quadrature Point Loop

  // Return
  return 0;
}

CEED_QFUNCTION(ReynoldsStress_Conserv)(void *ctx, CeedInt Q,
    const CeedScalar *const *in, CeedScalar *const *out) {
  return ReynoldsStress(ctx, Q, in, out, StateFromU);
}

CEED_QFUNCTION(ReynoldsStress_Prim)(void *ctx, CeedInt Q,
    const CeedScalar *const *in, CeedScalar *const *out) {
  return ReynoldsStress(ctx, Q, in, out, StateFromY);
}

#endif // reynolds_stress_h
