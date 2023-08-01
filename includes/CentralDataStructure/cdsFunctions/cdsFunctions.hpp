#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP_

#include "includes/CentralDataStructure/atom.hpp"
#include <vector>

namespace cds
{

    double getCharge(std::vector<cds::Atom*> atoms);
    void EnsureIntegralCharge(std::vector<cds::Atom*> atoms);

    // Templated functions
    template<typename T> void serializeNumbers(std::vector<T*> elements)
    {
        unsigned int i = 0;
        for (auto& element : elements)
        {
            element->setNumber(++i);
        }
        return;
    }
} // namespace cds
#endif /* INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP_ */
