#ifndef INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP
#define INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP

#include <utility>
#include <map>
#include <string>

namespace atomicBonds
{
    const double maxCutOff = 1.65;
    const double minCutOff = 0.7;

    const std::map<std::string, std::pair<double, double>> bondLengthMap = {
        {"CC", std::make_pair(1.22,  1.67)},
        {"CO", std::make_pair(1.08,  1.68)},
        {"OC", std::make_pair(1.08,  1.68)},
        {"CN", std::make_pair(1.26,  1.55)},
        {"NC", std::make_pair(1.26,  1.55)},
        {"OP", std::make_pair(1.35, 1.776)},
        {"PO", std::make_pair(1.35, 1.776)},
        {"OS", std::make_pair(1.43,  1.78)},
        {"SO", std::make_pair(1.43,  1.78)},
        {"NS", std::make_pair(1.62,  1.77)},
        {"SN", std::make_pair(1.62,  1.77)},
        {"CS", std::make_pair(1.50,  1.91)},
        {"SC", std::make_pair(1.50,  1.91)},
        {"SS", std::make_pair(1.50,  2.10)}
    };
    // FUNCTIONS
    std::pair<double, double> getBondLengthByAtomType(const std::string& atom1Element, const std::string& atom2Element);

    double getMaxBondLengthByAtomType(const std::string& atom1Element, const std::string& atom2Element);

    template<class atomT> inline bool bondAtomsIfClose(atomT* atom1, atomT* atom2)
    {
        // std::mutex mtx;           // mutex for critical section
        double maxLength = atomicBonds::getMaxBondLengthByAtomType(atom1->getElement(), atom2->getElement());
        if (atom1->getCoordinate()->withinDistance(atom2->getCoordinate(), maxLength))
        {
            // std::lock_guard<std::mutex> guard(mtx); // pre C++17 version
            // std::lock_guard guard(mtx);//RAII, the mutex will be unlocked upon guard destruction. Exception safe.
            atom1->addBond(atom2);
            // std::cout << "Bonded " << atom1->getName() << "_" << atom1->getIndex() << " to " << atom2->getName() <<
            // "_" << atom2->getIndex() << std::endl; std::stringstream ss; std::cout << "Bonded " << atom1->getName()
            // << "_" << atom1->getIndex() << " to " << atom2->getName() << "_" << atom2->getIndex() << "\n";
            // gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
            return true;
        }
        return false;
    }

} // namespace atomicBonds

#endif /* INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP_ */
