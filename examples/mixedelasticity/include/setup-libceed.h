#ifndef setuplibceed_h
#define setuplibceed_h

#include "structs.h"

// Convert PETSc MemType to libCEED MemType
CeedMemType MemTypeP2C(PetscMemType mtype);
// Destroy libCEED objects
PetscErrorCode CeedDataDestroy(CeedData ceed_data);
// Utility function to create local CEED restriction from DMPlex
PetscErrorCode CreateRestrictionForDomain(Ceed ceed, DM dm, CeedInt height,
    DMLabel domain_label, CeedInt value, CeedElemRestriction *elem_restr_x,
    CeedElemRestriction *elem_restr_u, CeedElemRestriction *elem_restr_p);
// Set up libCEED for a given degree
PetscErrorCode SetupLibceed(DM dm, Ceed ceed, AppCtx app_ctx,
                            ProblemData problem_data,
                            CeedData ceed_data);
#endif // setuplibceed_h
