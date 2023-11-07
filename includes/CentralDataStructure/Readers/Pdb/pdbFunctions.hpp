#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFUNCTIONS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFUNCTIONS_HPP_

#include <string>
#include "includes/CentralDataStructure/coordinate.hpp"

namespace pdb
{
    int checkShiftFromSerialNumberOverrun(const std::string& line);
    int checkSecondShiftFromResidueNumberOverrun(const std::string& line, const int shift = 0);
    cds::Coordinate checkShiftsAndExtractCoordinate(const std::string& line);

};     // namespace pdb
#endif /* INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFUNCTIONS_HPP_ */
