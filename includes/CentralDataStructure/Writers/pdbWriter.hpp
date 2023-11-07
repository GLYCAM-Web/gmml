#ifndef INCLUDES_CENTRALDATASTRUCTURE_WRITERS_PDBWRITER_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_WRITERS_PDBWRITER_HPP_

#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"

#include <string>
#include <iostream>

// Not all the options are available now, like e.g. writing molecule chain numbers, but I will add them as I need them
namespace cds
{
    void writeEnsembleToPdb(std::ostream& stream, const std::vector<cds::Assembly*> molecules);
    void writeAssemblyToPdb(std::ostream& stream, const std::vector<cds::Molecule*> molecules);
    void writeMoleculeToPdb(std::ostream& stream, const std::vector<cds::Residue*> residues,
                            unsigned int coordinateSetNumber = 0);
    void writeResidueToPdb(std::ostream& stream, const cds::Residue* residue, const std::string recordName = "ATOM",
                           unsigned int coordinateSetNumber = 0);
    void writeAtomToPdb(std::ostream& stream, const cds::Atom* atom, const std::string recordName = "ATOM",
                        const std::string residueName = "", const int residueNumber = 1, const std::string chainId = "",
                        const std::string insertionCode = "", const double occupancy = 1.00,
                        const double temperatureFactor = 0.00);
    void writeConectCards(std::ostream& stream, const std::vector<cds::Residue*> residues);
    void writeTrajectoryToPdb(std::ostream& stream, const std::vector<cds::Molecule*> molecules);

} // namespace cds
#endif /* INCLUDES_CENTRALDATASTRUCTURE_WRITERS_PDBWRITER_HPP_ */
