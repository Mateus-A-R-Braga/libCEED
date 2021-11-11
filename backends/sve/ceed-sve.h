// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and other CEED contributors.
// All Rights Reserved. See the top-level LICENSE and NOTICE files for details.
//
// SPDX-License-Identifier: BSD-2-Clause
//
// This file is part of CEED:  http://github.com/ceed

#ifndef _ceed_sve_h
#define _ceed_sve_h

#include <ceed/ceed.h>
#include <ceed/backend.h>

CEED_INTERN int CeedTensorContractCreate_f32_Sve(CeedBasis basis,
    CeedTensorContract contract);
CEED_INTERN int CeedTensorContractCreate_f64_Sve(CeedBasis basis,
    CeedTensorContract contract);

#endif // _ceed_sve_h
