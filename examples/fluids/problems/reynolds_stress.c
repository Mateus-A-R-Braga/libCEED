// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and other CEED contributors.
// All Rights Reserved. See the top-level LICENSE and NOTICE files for details.
//
// SPDX-License-Identifier: BSD-2-Clause
//
// This file is part of CEED:  http://github.com/ceed

/// @file
/// Utility functions for setting up problems using the various statistics Qfunctions such as 
/// qfunctions/reynolds_stress.h

#include "../navierstokes.h"
#include "../qfunctions/setupgeo.h"
#include "../qfunctions/reynolds_stress.h"

// Need to make two functions in this file:
// 1) function that sets up the Qfunction structure that specifies what qfunciton needs to happen and context for it. This is called inside newtonian_ns (add a flag to turn stats collection on and off)
// 2) create the ceed operator itself


// -- Create QFunction for Reynolds stress
PetscErrorCode CreateStatsOperator(stats, ceed_data, ...........){
// stats could be somthing like the Reynolds stress and is an instance of the problem QFunction specifier

// setup restriction
// setup basis
// create operator


  int num_comp_q = 5;
  int q_data_size_vol = ; // What size should this be?
  int num_comp_x = 3;
  int num_comp_stats = 6;

  CeedQFunction qf_stats
  CeedQFunctionCreateInterior(ceed, 1, stats.qfunction, rs.qfunction_loc, &qf_stats);
  CeedQFunctionSetContext(&qf_stats, stats.qfunction_context);
  CeedQFunctionContextDestroy(&stats.qfunction_context);
  CeedQFunctionAddInput(qf_stats, "q", num_comp_q, CEED_EVAL_INTERP); // This sets the QFunction input 
  CeedQFunctionAddInput(qf_stats, "q_data", q_data_size_vol, CEED_EVAL_NONE); // This sets the QFunction input
  CeedQFunctionAddInput(qf_stats, "x", num_comp_x, CEED_EVAL_INTERP); // This sets the QFunction input
  CeedQFunctionAddOutput(qf_stats, "U_prod", num_comp_stats, CEED_EVAL_INTERP); // This sets the Qfunction output
  CeedBasisCreateTensorH1Lagrange(ceed, dim, num_comp_stats, P, Q, CEED_GAUSS, &basis_stats); //This creates the basis

// What should the 2nd and 3rd arguments of the CeedQFunctionAddInput or Output be?



// do the same proceedure for creating the ceed operator (lines 350 ish in src/setuplibceed.c)
// -- CEED operator for stats collection function
PetscErrorCode CreateCeedOperator(    ){

  CeedOperator op;
  CeedOperatorCreate(ceed, ceed_data->qf_stats, NULL, NULL, &op);
  CeedOperatorSetField(op, "q", ceed_data->elem_restr_q, ceed_data->basis_q, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(op, "q_data", ceed_data->elem_restr_qd_i, CEED_BASIS_COLLOCATED, ceed_data->q_data);
  CeedOperatorSetField(op, "x", ceed_data->elem_restr_x, ceed_data->basis_x, ceed_data->x_coord);
  CeedOperatorSetField(op, "U_prod", ceed_data->elem_restr_q, ceed_data->basis_q, CEED_VECTOR_ACTIVE);
  CeedQFunctionDestroy(&qf_ijacobian_vol);
}
