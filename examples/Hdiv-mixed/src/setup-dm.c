#include "../include/setup-dm.h"
#include "petscerror.h"

// ---------------------------------------------------------------------------
// Setup DM
// ---------------------------------------------------------------------------
PetscErrorCode CreateDM(MPI_Comm comm, VecType vec_type, DM *dm) {

  PetscFunctionBeginUser;

  // Create DMPLEX
  PetscCall( DMCreate(comm, dm) );
  PetscCall( DMSetType(*dm, DMPLEX) );
  PetscCall( DMSetVecType(*dm, vec_type) );
  // Set Tensor elements
  PetscCall( PetscOptionsSetValue(NULL, "-dm_plex_simplex", "0") );
  // Set CL options
  PetscCall( DMSetFromOptions(*dm) );
  PetscCall( DMViewFromOptions(*dm, NULL, "-dm_view") );

  PetscFunctionReturn(0);
};

PetscErrorCode PerturbVerticesSmooth(DM dm) {

  Vec          coordinates;
  PetscSection coordSection;
  PetscScalar *coords;
  PetscInt     v,vStart,vEnd,offset,dim;
  PetscReal    x,y,z;

  PetscFunctionBeginUser;

  PetscCall( DMGetDimension(dm, &dim) );
  PetscCall( DMGetCoordinateSection(dm, &coordSection) );
  PetscCall( DMGetCoordinatesLocal(dm, &coordinates) );
  PetscCall( DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd) );
  PetscCall( VecGetArray(coordinates,&coords) );
  for(v=vStart; v<vEnd; v++) {
    PetscCall( PetscSectionGetOffset(coordSection,v,&offset) );
    if(dim==2) {
      PetscReal domain_min[2], domain_max[2], domain_size[2];
      PetscCall( DMGetBoundingBox(dm, domain_min, domain_max) );
      for (PetscInt i=0; i<2; i++) domain_size[i] = domain_max[i] - domain_min[i];
      x = coords[offset]; y = coords[offset+1];
      coords[offset]   = x + (0.06*domain_size[0])*PetscSinReal(
                           2.0*PETSC_PI*x/domain_size[0])*PetscSinReal(
                           2.0*PETSC_PI*y/domain_size[1]);
      coords[offset+1] = y - (0.05*domain_size[1])*PetscSinReal(
                           2.0*PETSC_PI*x/domain_size[0])*PetscSinReal(
                           2.0*PETSC_PI*y/domain_size[1]);
    } else {
      PetscReal domain_min[3], domain_max[3], domain_size[3];
      PetscCall( DMGetBoundingBox(dm, domain_min, domain_max) );
      for (PetscInt i=0; i<3; i++) domain_size[i] = domain_max[i] - domain_min[i];
      x = coords[offset]; y = coords[offset+1]; z = coords[offset+2];
      coords[offset]   = x + (0.03*domain_size[0])*PetscSinReal(
                           3*PETSC_PI*x/domain_size[0])*PetscCosReal(
                           3*PETSC_PI*y/domain_size[1])*PetscCosReal(3*PETSC_PI*z/domain_size[2]);
      coords[offset+1] = y - (0.04*domain_size[1])*PetscCosReal(
                           3*PETSC_PI*x/domain_size[0])*PetscSinReal(
                           3*PETSC_PI*y/domain_size[1])*PetscCosReal(3*PETSC_PI*z/domain_size[2]);
      coords[offset+2] = z + (0.05*domain_size[2])*PetscCosReal(
                           3*PETSC_PI*x/domain_size[0])*PetscCosReal(
                           3*PETSC_PI*y/domain_size[1])*PetscSinReal(3*PETSC_PI*z/domain_size[2]);
    }
  }
  PetscCall( VecRestoreArray(coordinates,&coords) );
  PetscFunctionReturn(0);
}
