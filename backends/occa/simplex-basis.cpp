// Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

#include "simplex-basis.hpp"


namespace ceed {
  namespace occa {
    SimplexBasis::SimplexBasis() {
      OCCA_DEBUG_TRACE("simplex-basis: SimplexBasis");
    }

    SimplexBasis::~SimplexBasis() {
      OCCA_DEBUG_TRACE("simplex-basis: ~SimplexBasis");
    }

    int SimplexBasis::apply(const CeedInt elementCount,
                            CeedTransposeMode tmode,
                            CeedEvalMode emode,
                            Vector *u,
                            Vector *v) {
      OCCA_DEBUG_TRACE("simplex-basis: apply");

      return 0;
    }

    //---[ Ceed Callbacks ]-------------
    int SimplexBasis::ceedCreate(CeedElemTopology topo, CeedInt dim,
                                 CeedInt ndof, CeedInt nqpts,
                                 const CeedScalar *interp,
                                 const CeedScalar *grad,
                                 const CeedScalar *qref,
                                 const CeedScalar *qweight,
                                 CeedBasis basis) {
      OCCA_DEBUG_TRACE("simplex-basis: ceedCreate");

      int ierr;
      Ceed ceed;
      ierr = CeedBasisGetCeed(basis, &ceed); CeedChk(ierr);

      SimplexBasis *basis_ = new SimplexBasis();
      ierr = CeedBasisSetData(basis, (void**) &basis_); CeedChk(ierr);

      ierr = registerBasisFunction(ceed, basis, "Apply",
                                   (ceed::occa::ceedFunction) Basis::ceedApply);
      CeedChk(ierr);

      ierr = registerBasisFunction(ceed, basis, "Destroy",
                                   (ceed::occa::ceedFunction) Basis::ceedDestroy);
      CeedChk(ierr);

      return 0;
    }
  }
}
