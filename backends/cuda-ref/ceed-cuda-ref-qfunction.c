// Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. LLNL-CODE-734707.
// All Rights reserved. See files LICENSE and NOTICE for details.
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

#include <ceed/ceed.h>
#include <ceed/backend.h>
#include <cuda.h>
#include <stdio.h>
#include <string.h>
#include "ceed-cuda-ref.h"
#include "ceed-cuda-ref-qfunction-load.h"
#include "../cuda/ceed-cuda-compile.h"

//------------------------------------------------------------------------------
// Apply QFunction
//------------------------------------------------------------------------------
static int CeedQFunctionApply_Cuda(CeedQFunction qf, CeedInt Q,
                                   CeedVector *U, CeedVector *V) {
  int ierr;
  Ceed ceed;
  ierr = CeedQFunctionGetCeed(qf, &ceed); CeedChkBackend(ierr);

  // Build and compile kernel, if not done
  ierr = CeedCudaBuildQFunction(qf); CeedChkBackend(ierr);

  CeedQFunction_Cuda *data;
  ierr = CeedQFunctionGetData(qf, &data); CeedChkBackend(ierr);
  Ceed_Cuda *ceed_Cuda;
  ierr = CeedGetData(ceed, &ceed_Cuda); CeedChkBackend(ierr);
  CeedInt num_input_fields, num_output_fields;
  ierr = CeedQFunctionGetNumArgs(qf, &num_input_fields, &num_output_fields);
  CeedChkBackend(ierr);

  // Read vectors
  for (CeedInt i = 0; i < num_input_fields; i++) {
    ierr = CeedVectorGetArrayRead(U[i], CEED_MEM_DEVICE, &data->fields.inputs[i]);
    CeedChkBackend(ierr);
  }
  for (CeedInt i = 0; i < num_output_fields; i++) {
    ierr = CeedVectorGetArrayWrite(V[i], CEED_MEM_DEVICE, &data->fields.outputs[i]);
    CeedChkBackend(ierr);
  }

  // Get context data
  CeedQFunctionContext ctx;
  ierr = CeedQFunctionGetInnerContext(qf, &ctx); CeedChkBackend(ierr);
  if (ctx) {
    ierr = CeedQFunctionContextGetData(ctx, CEED_MEM_DEVICE, &data->d_c);
    CeedChkBackend(ierr);
  }

  // Run kernel
  void *args[] = {&data->d_c, (void *) &Q, &data->fields};
  ierr = CeedRunKernelAutoblockCuda(ceed, data->QFunction, Q, args);
  CeedChkBackend(ierr);

  // Restore vectors
  for (CeedInt i = 0; i < num_input_fields; i++) {
    ierr = CeedVectorRestoreArrayRead(U[i], &data->fields.inputs[i]);
    CeedChkBackend(ierr);
  }
  for (CeedInt i = 0; i < num_output_fields; i++) {
    ierr = CeedVectorRestoreArray(V[i], &data->fields.outputs[i]);
    CeedChkBackend(ierr);
  }

  // Restore context
  if (ctx) {
    ierr = CeedQFunctionContextRestoreData(ctx, &data->d_c);
    CeedChkBackend(ierr);
  }
  return CEED_ERROR_SUCCESS;
}

//------------------------------------------------------------------------------
// Destroy QFunction
//------------------------------------------------------------------------------
static int CeedQFunctionDestroy_Cuda(CeedQFunction qf) {
  int ierr;
  CeedQFunction_Cuda *data;
  ierr = CeedQFunctionGetData(qf, &data); CeedChkBackend(ierr);
  Ceed ceed;
  ierr = CeedQFunctionGetCeed(qf, &ceed); CeedChkBackend(ierr);
  if (data->module)
    CeedChk_Cu(ceed, cuModuleUnload(data->module));
  ierr = CeedFree(&data); CeedChkBackend(ierr);

  return CEED_ERROR_SUCCESS;
}

//------------------------------------------------------------------------------
// Set User QFunction
//------------------------------------------------------------------------------
static int CeedQFunctionSetCUDAUserFunction_Cuda(CeedQFunction qf,
    CUfunction f) {
  int ierr;
  CeedQFunction_Cuda *data;
  ierr = CeedQFunctionGetData(qf, &data); CeedChkBackend(ierr);
  data->QFunction = f;

  return CEED_ERROR_SUCCESS;
}

//------------------------------------------------------------------------------
// Create QFunction
//------------------------------------------------------------------------------
int CeedQFunctionCreate_Cuda(CeedQFunction qf) {
  int ierr;
  Ceed ceed;
  CeedQFunctionGetCeed(qf, &ceed);
  CeedQFunction_Cuda *data;
  ierr = CeedCalloc(1, &data); CeedChkBackend(ierr);
  ierr = CeedQFunctionSetData(qf, data); CeedChkBackend(ierr);

  // Read QFunction source
  ierr = CeedQFunctionGetKernelName(qf, &data->qfunction_name);
  CeedChkBackend(ierr);
  CeedDebug256(ceed, 2, "----- Loading QFunction User Source -----\n");
  ierr = CeedQFunctionLoadSourceToBuffer(qf, &data->qfunction_source);
  CeedChkBackend(ierr);
  CeedDebug256(ceed, 2, "----- Loading QFunction User Source Complete! -----\n");

  // Register backend functions
  ierr = CeedSetBackendFunction(ceed, "QFunction", qf, "Apply",
                                CeedQFunctionApply_Cuda); CeedChkBackend(ierr);
  ierr = CeedSetBackendFunction(ceed, "QFunction", qf, "Destroy",
                                CeedQFunctionDestroy_Cuda); CeedChkBackend(ierr);
  ierr = CeedSetBackendFunction(ceed, "QFunction", qf, "SetCUDAUserFunction",
                                CeedQFunctionSetCUDAUserFunction_Cuda);
  CeedChkBackend(ierr);
  return CEED_ERROR_SUCCESS;
}
//------------------------------------------------------------------------------