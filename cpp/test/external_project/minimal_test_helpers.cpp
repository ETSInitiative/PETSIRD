/*
  Copyright (C) 2026 University College London

  SPDX-License-Identifier: Apache-2.0
*/

// (un)comment if you want HDF5 or binary output
#define USE_HDF5

#ifdef USE_HDF5
#  include "petsird/hdf5/protocols.h"
using petsird::hdf5::PETSIRDReader;
#else
#  include "petsird/binary/protocols.h"
using petsird::binary::PETSIRDReader;
#endif
#include "petsird_helpers.h"
#include "petsird_helpers/geometry.h"
#include "petsird_helpers/create.h"

int
main(int, char const* argv[])
{
  PETSIRDReader reader(argv[1]);
  petsird::Header header;
  reader.ReadHeader(header);

  return 0;
}
