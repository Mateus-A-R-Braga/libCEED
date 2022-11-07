// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and other CEED contributors.
// All Rights Reserved. See the top-level LICENSE and NOTICE files for details.
//
// SPDX-License-Identifier: BSD-2-Clause
//
// This file is part of CEED:  http://github.com/ceed

/// @file
/// Utility functions for setting up vortex shedding problem behind a cylinder

#include "../navierstokes.h"
#include "../qfunctions/vortexshedding.h"
#include "../qfunctions/freestream_bc.h"
#include "ceed/ceed-f64.h"
#include "stg_shur14.h"

PetscErrorCode NS_VORTEXSHEDDING(ProblemData *problem, DM dm, void *ctx) {

  PetscInt ierr;
  User      user    = *(User *)ctx;
  MPI_Comm  comm    = PETSC_COMM_WORLD;
  VortexsheddingContext vortexshedding_ctx;
  FreestreamContext freestream_ctx;
  NewtonianIdealGasContext newtonian_ig_ctx;
  CeedQFunctionContext vortexshedding_context, freestream_context;

  PetscFunctionBeginUser;
  ierr = NS_NEWTONIAN_IG(problem, dm, ctx); CHKERRQ(ierr);
  ierr = PetscCalloc1(1, &vortexshedding_ctx); CHKERRQ(ierr);

  // ------------------------------------------------------
  //               SET UP Vortex Shedding
  // ------------------------------------------------------
  CeedQFunctionContextDestroy(&problem->ics.qfunction_context);

  switch (user->phys->state_var) {
  case STATEVAR_CONSERVATIVE:
    problem->ics.qfunction      = ICsVortexshedding_Conserv;
    problem->ics.qfunction_loc  = ICsVortexshedding_Conserv_loc;
    problem->apply_freestream.qfunction     = Freestream_Conserv;
    problem->apply_freestream.qfunction_loc = Freestream_Conserv_loc;
  case STATEVAR_PRIMITIVE:
    problem->ics.qfunction      = ICsVortexshedding_Prim;
    problem->ics.qfunction_loc  = ICsVortexshedding_Prim_loc;
    problem->apply_freestream.qfunction     = Freestream_Prim;
    problem->apply_freestream.qfunction_loc = Freestream_Prim_loc;
  }

  CeedScalar U_in[3]            = {1.0};       // m/s
  CeedScalar T_in               = 300.;        // K
  CeedScalar P0                 = 1.e5;        // Pa
  CeedScalar L;
  CeedScalar H;
  CeedScalar D;
  CeedScalar radius;

  PetscOptionsBegin(comm, NULL, "Options for VORTEX SHEDDING problem", NULL);
  PetscInt narray=3;
  ierr = PetscOptionsScalarArray("-velocity_inflow",
                            "Velocity at inflow",
                            NULL, U_in, &narray, NULL); CHKERRQ(ierr);
  PetscCheck(narray == 3, comm, PETSC_ERR_ARG_SIZ,
             "-velocity_inflow should recieve array of size 3, instead recieved size %"
             PetscInt_FMT".", narray);
  ierr = PetscOptionsScalar("-P0", "Pressure at outflow",
                            NULL, P0, &P0, NULL); CHKERRQ(ierr);
  ierr = PetscOptionsScalar("-temperature_inflow", "Temperature at inflow",
                            NULL, T_in, &T_in, NULL); CHKERRQ(ierr);
  ierr = PetscOptionsScalar("-L", "Length of the rectangular channel",
                            NULL, L, &L, NULL); CHKERRQ(ierr);
  ierr = PetscOptionsScalar("-H", "Heigth of the rectangular channel",
                            NULL, H, &H, NULL); CHKERRQ(ierr);
  ierr = PetscOptionsScalar("-D", "Cylinder diameter",
                            NULL, D, &D, NULL); CHKERRQ(ierr);
  ierr = PetscOptionsScalar("-radius", "Cylinder radius",
                            NULL, radius, &radius, NULL); CHKERRQ(ierr);
  PetscOptionsEnd();

  PetscScalar meter  = user->units->meter;
  PetscScalar second = user->units->second;
  PetscScalar Kelvin = user->units->Kelvin;
  PetscScalar Pascal = user->units->Pascal;

  T_in   *= Kelvin;
  P0     *= Pascal;
  for (int i=0; i<3; i++) {
    U_in[i]     *= meter / second;
  }
  L      *= meter;
  H      *= meter;
  D      *= meter;
  radius *= meter;

  CeedQFunctionContextGetData(problem->apply_vol_rhs.qfunction_context,
                              CEED_MEM_HOST, &newtonian_ig_ctx);
  State S_in;
  {
    CeedScalar Y[5] = {P0, U_in[0], U_in[1], U_in[2], T_in};
    CeedScalar x[3] = {0.};
    S_in = StateFromY(newtonian_ig_ctx, Y, x);
  }

  //-- Setup Problem information
  CeedScalar center;
  {
    PetscReal domain_min[3], domain_max[3], domain_size[3];
    ierr = DMGetBoundingBox(dm, domain_min, domain_max); CHKERRQ(ierr);
    // compute L and H
    for (PetscInt i=0; i<3; i++) domain_size[i] = domain_max[i] - domain_min[i];
    L = domain_size[0]*meter;
    H = domain_size[1]*meter;
    center = 0.5 * H; // this should be checked
  }
  // Some properties depend on parameters from NewtonianIdealGas
  CeedQFunctionContextGetData(problem->apply_vol_rhs.qfunction_context,
                              CEED_MEM_HOST, &newtonian_ig_ctx);

  vortexshedding_ctx->L        = L;
  vortexshedding_ctx->H        = H;
  vortexshedding_ctx->D        = D;
  vortexshedding_ctx->center   = center;
  vortexshedding_ctx->radius   = radius;
  vortexshedding_ctx->T_in     = T_in;
  vortexshedding_ctx->P0       = P0;
  vortexshedding_ctx->S_in     = S_in;
  vortexshedding_ctx->implicit = user->phys->implicit;
  vortexshedding_ctx->newtonian_ctx = *newtonian_ig_ctx;

  ierr = PetscCalloc1(1, &freestream_ctx); CHKERRQ(ierr);
  freestream_ctx->newtonian_ctx = *newtonian_ig_ctx;
  freestream_ctx->S_in       = S_in;

  CeedQFunctionContextRestoreData(problem->apply_vol_rhs.qfunction_context,
                                  &newtonian_ig_ctx);

  CeedQFunctionContextCreate(user->ceed, &vortexshedding_context);
  CeedQFunctionContextSetData(vortexshedding_context, CEED_MEM_HOST,
                              CEED_USE_POINTER,
                              sizeof(*vortexshedding_ctx), vortexshedding_ctx);
  CeedQFunctionContextSetDataDestroy(vortexshedding_context, CEED_MEM_HOST,
                                     FreeContextPetsc);
  CeedQFunctionContextDestroy(&problem->ics.qfunction_context);
  problem->ics.qfunction_context = vortexshedding_context;


  CeedQFunctionContextCreate(user->ceed, &freestream_context);
  CeedQFunctionContextSetData(freestream_context, CEED_MEM_HOST, CEED_USE_POINTER,
                              sizeof(*freestream_ctx), freestream_ctx);
  CeedQFunctionContextSetDataDestroy(freestream_context, CEED_MEM_HOST,
                                     FreeContextPetsc);
  problem->apply_freestream.qfunction_context = freestream_context;

  PetscFunctionReturn(0);
}