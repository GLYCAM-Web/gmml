#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFUNCTIONS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFUNCTIONS_HPP_

#include <string>

namespace pdb
{
    int checkShiftFromSerialNumberOverrun(std::string line);
};     // namespace pdb
#endif /* INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFUNCTIONS_HPP_ */
