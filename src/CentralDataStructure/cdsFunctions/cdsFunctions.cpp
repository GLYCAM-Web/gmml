#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CodeUtils/constants.hpp" //dNotSet
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CodeUtils/numbers.hpp" // isNumberIntegral
#include "includes/CodeUtils/logging.hpp"
#include <iomanip> // setprecision

double cds::getCharge(std::vector<cds::Atom*> atoms)
{
    double totalCharge = 0.0;
    for (auto& atom : atoms)
    {
        if (atom->getCharge() != constants::dNotSet) // Better to throw/warn here?
        {
            totalCharge += atom->getCharge();
        }
    }
    return totalCharge;
}

void cds::EnsureIntegralCharge(std::vector<cds::Atom*> atoms)
{
    double charge = cds::getCharge(atoms);
    std::stringstream ss;
    ss << std::fixed;
    ss << "Total charge is: " << std::setprecision(5) << charge << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    if (!codeUtils::isNumberIntegral(charge))
    {
        std::stringstream errorMessage;
        errorMessage << "Non-integral charge (" << charge << "). You cannot run MD with this.\n";
        throw errorMessage.str();
    }
    return;
}
