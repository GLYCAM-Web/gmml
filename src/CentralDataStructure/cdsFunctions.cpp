#include "includes/CentralDataStructure/cdsFunctions.hpp"

void cds::serializeAtomNumbers(std::vector<cds::Atom*> atoms)
{
    unsigned int i = 0;
    for(auto & atom : atoms)
    {
        atom->setNumber(++i);
    }
    return;
}
