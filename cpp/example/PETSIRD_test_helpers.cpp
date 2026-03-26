#include <iostream>
#include "petsird/binary/protocols.h"
using petsird::binary::PETSIRDReader;
#include "petsird_helpers/create.h"

int
main(int argc, char const* argv[])
{
  auto prog_name = argv[0];

  std::string filename;
  // option processing
  while (argc > 1 && (strncmp(argv[1], "-", 1) == 0))
    {
      if (strcmp(argv[1], "--input") == 0 || strcmp(argv[1], "-i") == 0)
        {
          filename = argv[2];
          ++argv;
          --argc;
        }
      else
        {
          std::cerr << "Wrong options\n";
          return 1;
        }
      ++argv;
      --argc;
    }

  // Open the file
  PETSIRDReader reader(filename);
  petsird::Header header;
  reader.ReadHeader(header);

  return 0;
}
