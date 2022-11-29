#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_LIBRARYRESIDUE_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_LIBRARYRESIDUE_HPP_

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Readers/LibraryAtom.hpp"
#include "includes/CodeUtils/constants.hpp"
#include <sstream>

namespace lib
{
class LibraryResidue : public cds::Residue
{
public:
    LibraryResidue(std::stringstream& residueStream, const std::string& name);
private:
//  I don't think we ever use the box
//    double box_angle_;
//    double box_length_;
//    double box_width_;
//    double box_height_;
    int head_atom_index_ = constants::iNotSet;
    int tail_atom_index_ = constants::iNotSet;
};
}//namespace lib
#endif
